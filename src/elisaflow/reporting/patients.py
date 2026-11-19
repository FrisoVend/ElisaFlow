"""Render the patients section from stored analysis results."""

from collections.abc import Callable
from dataclasses import dataclass
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_RIGHT
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import BaseDocTemplate, CondPageBreak, Flowable, Frame, FrameSplitter, PageTemplate, Paragraph, Table

from elisaflow.formatting import format_extinction_scale, format_number, format_result, format_text, is_finite_number, metric_decimals, metric_label
from elisaflow.models import Sample, SampleType, ScaledStatistics, Statistics, ValidationMetric
from elisaflow.reporting.common import ReportArea, ReportContext
from elisaflow.reporting.tables import TableColumn, paragraph_width, table_headers
from elisaflow.sample_selection import select_samples

_COLUMN_COUNT = 3
_COLUMN_GAP = 2.5 * mm
_PAGE_TEMPLATE = "Patients"


@dataclass(frozen=True)
class PatientPages:
    """Matching patient page template and section content."""

    template: PageTemplate
    flowables: tuple[Flowable, ...]


@dataclass(frozen=True)
class _PatientColumn:
    """Keep a patient column's heading and ordered cells together."""

    header: TableColumn
    cells: tuple[Paragraph, ...]


def patient_pages(
    context: ReportContext, area: ReportArea, *, on_page: Callable[[Canvas, BaseDocTemplate], None],
) -> PatientPages:
    """Build matching frames and tables for the fixed three-column patient section."""
    table_width = (area.width - (_COLUMN_COUNT - 1) * _COLUMN_GAP) / _COLUMN_COUNT
    if table_width <= 0:
        raise ValueError("Het rapportoppervlak is te smal voor de patiëntkolommen.")
    frame_ids = tuple(f"patient-{index}" for index in range(_COLUMN_COUNT))
    frames = [
        Frame(
            area.left + index * (table_width + _COLUMN_GAP), area.bottom, table_width, area.height,
            id=frame_id, leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0,
        )
        for index, frame_id in enumerate(frame_ids)
    ]
    template = PageTemplate(id=_PAGE_TEMPLATE, frames=frames, onPage=on_page)
    flowables = []
    patients = tuple(
        (sample, context.result.patient_extinction(sample)) for sample in select_samples(context.result.samples, SampleType.UNKNOWN)
    )
    if not patients:
        return PatientPages(template, ())

    table = _patient_table(context, patients, table_width)
    section_style = ParagraphStyle("PatientSection", parent=context.styles["Heading3"], keepWithNext=False)
    flowables.append(CondPageBreak(60))
    flowables.append(Paragraph("Patiëntresultaten", section_style))
    # Fill each column to the page bottom before moving to the next column.
    flowables.append(FrameSplitter(_PAGE_TEMPLATE, nextFrames=list(frame_ids), gap=0, required=30))
    flowables.append(table)

    return PatientPages(template, tuple(flowables))


def _patient_table(context: ReportContext, patients: tuple[tuple[Sample, ScaledStatistics | None], ...], table_width: float) -> Table:
    """Choose, render and size each column from the same patient data."""
    cell_style = ParagraphStyle("PatientCell", parent=context.styles["Normal"], fontSize=7, leading=8.5)
    name_style = ParagraphStyle("PatientName", parent=cell_style, fontName="Helvetica-Bold")
    number_style = ParagraphStyle("PatientNumber", parent=cell_style, fontSize=6.5, leading=8, alignment=TA_RIGHT)
    cv_warning_style = ParagraphStyle("PatientCVWarning", parent=number_style, textColor=colors.firebrick)
    header_style = ParagraphStyle("PatientHeader", parent=cell_style, fontName="Helvetica-Bold", fontSize=6.5, leading=8)
    extinction_labels = {format_extinction_scale(extinction.scale) for _, extinction in patients if extinction is not None}
    extinction_label = next(iter(extinction_labels)) if len(extinction_labels) == 1 else "Extinctie"
    columns = [
        _PatientColumn(
            TableColumn("sample_name", "Patiënt", border_after=True),
            tuple(Paragraph(escape(format_text(sample.name, missing="")), name_style) for sample, _ in patients),
        ),
        _PatientColumn(
            TableColumn("used_wells", "Coörd.", border_after=True),
            tuple(Paragraph(escape(", ".join(sample.used_wells)), cell_style) for sample, _ in patients),
        ),
    ]
    has_dilution = any(sample.dilution != 1 for sample, _ in patients)
    if context.multiple_measurements:
        columns.append(_PatientColumn(
            TableColumn("selected_wavelength_nm", "nm", border_after=not has_dilution),
            tuple(Paragraph(format_number(sample.selected_wavelength_nm, missing=""), cell_style) for sample, _ in patients),
        ))
    if has_dilution:
        columns.append(_PatientColumn(
            TableColumn("sample_dilution", "Verd.", border_after=True),
            tuple(Paragraph(format_number(sample.dilution, missing=""), cell_style) for sample, _ in patients),
        ))
    for metric, summaries in (
        (ValidationMetric.EXTINCTION, tuple(extinction.statistics if extinction is not None else Statistics() for _, extinction in patients)),
        (ValidationMetric.RESULT, tuple(sample.result for sample, _ in patients)),
    ):
        if metric == ValidationMetric.RESULT:
            columns.append(_PatientColumn(
                TableColumn("result_gap", "", border_after=True),
                tuple(Paragraph("", cell_style) for _ in patients),
            ))
        label = (
            extinction_label if metric == ValidationMetric.EXTINCTION
            else metric_label(metric, context.result.settings.result.label, context.result.settings.result.unit)
        )
        has_cv = (metric == ValidationMetric.EXTINCTION or context.report.show_result_cv) and any(is_finite_number(value.cv_percent) for value in summaries)
        decimals = metric_decimals(metric, context.report)
        mean_cells = []
        for (sample, extinction), statistics in zip(patients, summaries, strict=True):
            if metric == ValidationMetric.RESULT:
                text = format_result(statistics.mean, sample.result_status.value, sample.result_limit, decimals)
            else:
                text = format_number(statistics.mean, decimals)
                if len(extinction_labels) > 1 and extinction is not None:
                    text = f"{format_extinction_scale(extinction.scale)} {text}"
            mean_cells.append(Paragraph(escape(text), number_style))
        columns.append(_PatientColumn(
            TableColumn(
                f"{metric.value}_mean",
                "Gem." if has_cv else label,
                group=label if has_cv else None,
                border_after=not has_cv,
                header_alignment=TA_CENTER,
                numeric=True,
            ),
            tuple(mean_cells),
        ))
        if has_cv:
            columns.append(_PatientColumn(
                TableColumn(
                    f"{metric.value}_cv_percent", "CV%", group=label,
                    border_after=True, numeric=True,
                ),
                tuple(
                    Paragraph(
                        format_number(value.cv_percent, context.report.percentage_decimals, missing=""),
                        cv_warning_style if is_finite_number(value.cv_percent) and value.cv_percent > context.report.cv_warning_percent else number_style,
                    )
                    for value in summaries
                ),
            ))

    headers = table_headers([column.header for column in columns], header_style, rows=2, border=context.border)
    group_headers, subheaders = headers.rows
    cell_padding = 1.5
    result_gap_index = next(index for index, column in enumerate(columns) if column.header.key == "result_gap")
    column_widths = []
    for index, column in enumerate(columns):
        if index == result_gap_index:
            column_widths.append(0)
            continue
        needed = max(
            paragraph_width(group_headers[index]) if column.header.group is None else 0,
            paragraph_width(subheaders[index]),
            max(paragraph_width(cell) for cell in column.cells),
        ) + 2 * cell_padding + 1
        # Metadata uses text sizing; right-aligned statistics get more room.
        limit = table_width * (0.30 if column.header.numeric else 0.22)
        column_widths.append(needed if index == 0 else min(needed, limit))
    # Reuse the exact group ranges that determine the header spans.
    for start, end in headers.groups:
        needed = paragraph_width(group_headers[start]) + 2 * cell_padding + 1
        extra = max(0, needed - sum(column_widths[start:end])) / (end - start)
        for index in range(start, end):
            column_widths[index] += extra
    # Short names leave more space for complete numeric values and their labels.
    minimum_result_gap = 6
    other_width = sum(column_widths[1:])
    other_scale = min(1, (table_width - min(column_widths[0], table_width * 0.20) - minimum_result_gap) / other_width)
    column_widths[1:] = [width * other_scale for width in column_widths[1:]]
    column_widths[0] = min(column_widths[0], table_width - sum(column_widths[1:]) - minimum_result_gap)
    column_widths[result_gap_index] = table_width - sum(column_widths)

    patient_rows = [list(row) for row in zip(*(column.cells for column in columns), strict=True)]
    return context.table(
        [*headers.rows, *patient_rows],
        column_widths,
        header_rows=2,
        padding=(cell_padding, 1.5),
        commands=[
            ("LINEBELOW", (0, 2), (-1, -1), 0.25, context.border),
            ("ROWBACKGROUNDS", (0, 2), (-1, -1), [colors.white, context.stripe]),
            *headers.commands,
        ],
    )
