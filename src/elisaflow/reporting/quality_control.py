"""Render the quality control section from stored analysis results."""

from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.enums import TA_RIGHT
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.platypus import Paragraph, Table

from elisaflow.configuration.assay_settings import ResultSettings
from elisaflow.formatting import format_formula, format_number, format_result, is_finite_number, metric_decimals
from elisaflow.models import CalculationMethod, ResponseScale, SampleType, ValidationMetric, ValidationResult, ValidationStatus
from elisaflow.quality_control.summary import QualityControlGroup, quality_control_groups
from elisaflow.reporting.common import ReportContext
from elisaflow.reporting.tables import TableColumn, paragraph_width, table_headers
from elisaflow.sample_selection import select_samples


def _result_settings(context: ReportContext, result: ValidationResult) -> ResultSettings:
    settings = context.result.settings
    if result.sample_type == SampleType.CONTROL and settings.calculation.method == CalculationMethod.FORMULA:
        return settings.calculation.control_results.get(result.sample_name, settings.result)
    return settings.result


def _validation_criterion(context: ReportContext, result: ValidationResult) -> str:
    """Build escaped criterion text with wavelength subscript markup."""
    label = (
        context.extinction_label(raw=result.response_scale == ResponseScale.RAW)
        if result.metric == ValidationMetric.EXTINCTION else _result_settings(context, result).label
    )
    criterion = escape(label)
    wavelengths = [result.wavelength_nm] if result.wavelength_nm is not None else []
    if not wavelengths and context.multiple_measurements:
        if result.sample_type == SampleType.STANDARD:
            samples = select_samples(context.result.samples, result.sample_type, concentration=result.nominal_concentration)
        else:
            samples = select_samples(
                context.result.samples, result.sample_type, name=result.sample_name if result.sample_type == SampleType.CONTROL else None
            )
        wavelengths = sorted({sample.selected_wavelength_nm for sample in samples if sample.selected_wavelength_nm is not None})
    if wavelengths:
        criterion += context.wavelength_subscript(*(int(value) for value in wavelengths))
    if (
        context.result.settings.calculation.method == CalculationMethod.FORMULA
        and result.sample_type == SampleType.CONTROL
        and result.metric == ValidationMetric.RESULT
    ):
        override = context.result.settings.calculation.overrides.get(result.sample_name)
        if override is not None:
            criterion += f" ({escape(format_formula(override, result.sample_name))})"
    decimals = metric_decimals(result.metric, context.report)
    if result.minimum is not None:
        criterion = f"{format_number(result.minimum, decimals)} ≤ {criterion}"
    if result.maximum is not None:
        criterion = f"{criterion} ≤ {format_number(result.maximum, decimals)}"
    if result.minimum is None and result.maximum is None:
        criterion += ": geen grenzen"
    return criterion


def _validation_cells(
    context: ReportContext,
    groups: tuple[QualityControlGroup, ...],
) -> tuple[list[dict[str, Paragraph]], list[tuple[int, int]], bool]:
    status_labels = {
        ValidationStatus.PASSED: ("OK", colors.darkgreen),
        ValidationStatus.FAILED: ("Afgekeurd", colors.firebrick),
        ValidationStatus.NOT_EVALUATED: ("Niet bepaald", colors.dimgrey),
    }
    cell_style = ParagraphStyle("ValidationCell", parent=context.styles["Normal"], fontSize=7.5, leading=9.5)
    sample_style = ParagraphStyle("ValidationSample", parent=cell_style, fontName="Helvetica-Bold")
    number_style = ParagraphStyle("ValidationNumber", parent=cell_style, alignment=TA_RIGHT)
    cv_warning_style = ParagraphStyle("ValidationCVWarning", parent=number_style, textColor=colors.firebrick)
    show_cv = False
    data = []
    group_rows = []
    for group in groups:
        results = group.checks
        first_row = len(data)
        # Separate rows keep values and statuses aligned when a criterion wraps.
        for index, result in enumerate(results or [None]):
            criterion, value_text, unit, status = "", "", "", "-"
            value_cv = None
            if result is not None:
                criterion = _validation_criterion(context, result)
                decimals = metric_decimals(result.metric, context.report)
                value_text = format_number(result.value, decimals)
                if result.metric == ValidationMetric.RESULT:
                    unit = _result_settings(context, result).unit or ""
                    result_status = result.result_status.value if result.result_status is not None else None
                    value_text = format_result(result.value, result_status, result.result_limit, decimals)
                if is_finite_number(result.value) and (result.metric == ValidationMetric.EXTINCTION or context.report.show_result_cv):
                    value_cv = result.cv_percent
                show_cv |= is_finite_number(value_cv)
                label, color = status_labels[result.status]
                status = f'<font color="{color.hexval()}">{label}</font>'
            data.append(
                {
                    "sample": Paragraph(escape(group.name) if index == 0 else "", sample_style),
                    "coordinates": Paragraph(escape(", ".join(group.coordinates)) if index == 0 else "", cell_style),
                    "criterion": Paragraph(criterion, cell_style),
                    "value_gap": Paragraph("", cell_style),
                    "value": Paragraph(escape(value_text), number_style),
                    "value_cv": Paragraph(
                        format_number(value_cv, context.report.percentage_decimals, missing=""),
                        cv_warning_style if is_finite_number(value_cv) and value_cv > context.report.cv_warning_percent else number_style,
                    ),
                    "unit": Paragraph(escape(unit), cell_style),
                    "status": Paragraph(status, cell_style),
                }
            )
        group_rows.append((first_row, len(data) - 1))

    return data, group_rows, show_cv


def _validation_widths(data: list[list[Paragraph]], columns: list[TableColumn], table_width: float) -> list[float]:
    total_weight = sum(column.weight for column in columns)
    widths = [table_width * column.weight / total_weight for column in columns]
    numeric = [index for index, column in enumerate(columns) if column.numeric]
    for index in numeric:
        column = columns[index]
        sizing_rows = data[1:] if column.minimum_width_label else data
        needed = max(paragraph_width(row[index]) for row in sizing_rows) + 7
        if column.minimum_width_label:
            needed = max(needed, stringWidth(column.minimum_width_label, "Helvetica-Bold", 7.5) + 7)
        widths[index] = needed
    fixed_text = [index for index, column in enumerate(columns) if not column.numeric and column.minimum_width_label]
    for index in fixed_text:
        widths[index] = (
            max(stringWidth(word, row[index].style.fontName, row[index].style.fontSize) for row in data for word in row[index].getPlainText().split())
            + 7
        )
        if columns[index].key in ("sample", "coordinates"):
            limit = table_width * (0.22 if columns[index].key == "sample" else 0.16)
            widths[index] = min(max(paragraph_width(row[index]) for row in data) + 7, limit)
    numeric_width = sum(widths[index] for index in numeric)
    scale = min(1, table_width * 0.50 / numeric_width)
    for index in numeric:
        widths[index] *= scale
    remaining = table_width - numeric_width * scale - sum(widths[index] for index in fixed_text)
    text_weight = sum(column.weight for index, column in enumerate(columns) if not column.numeric and index not in fixed_text)
    for index, column in enumerate(columns):
        if not column.numeric and index not in fixed_text:
            widths[index] = remaining * column.weight / text_weight
    return widths


def validation_table(context: ReportContext, width: float) -> Table | None:
    """Render all stored QC together, with blanks before the plate-ordered groups."""
    groups = (
        *quality_control_groups(context.result, blank_only=True),
        *quality_control_groups(context.result),
    )
    if not groups:
        return None

    cells, group_rows, show_cv = _validation_cells(context, groups)
    cell_style = cells[0]["sample"].style
    header_style = ParagraphStyle("ValidationHeader", parent=cell_style, fontName="Helvetica-Bold")
    header_rows = 2 if show_cv else 1
    columns = [
        TableColumn("sample", "Monster", shared=True, border_after=True, minimum_width_label="Monster"),
        TableColumn("coordinates", "Coörd.", shared=True, border_after=True, minimum_width_label="Coörd."),
        TableColumn("criterion", "Criterium", weight=0.90),
        TableColumn("value_gap", "", weight=0.10, shared=True, border_after=True),
        TableColumn(
            "value",
            "Gem.",
            weight=0.14,
            group="Waarde",
            numeric=True,
            border_after=not show_cv,
            minimum_width_label="Gem." if show_cv else "Waarde",
        ),
        *([TableColumn("value_cv", "CV%", weight=0.10, group="Waarde", numeric=True, border_after=True)] if show_cv else []),
        TableColumn("unit", "Eenheid", weight=0.20, minimum_width_label="Eenheid", border_after=True),
        TableColumn("status", "Status", weight=0.18, minimum_width_label="Status", border_after=True),
    ]
    headers = table_headers(columns, header_style, rows=header_rows, border=context.border)
    data = [*headers.rows, *[[row[column.key] for column in columns] for row in cells]]
    column_widths = _validation_widths(data, columns, width)

    group_commands = []
    for group_index, (start, end) in enumerate(group_rows):
        start += header_rows
        end += header_rows
        group_commands.extend(
            [
                (
                    "BACKGROUND",
                    (0, start),
                    (-1, end),
                    colors.white if group_index % 2 == 0 else context.stripe,
                ),
                ("LINEBELOW", (0, end), (-1, end), 0.20, context.border),
            ]
        )
        if start != end:
            for index, column in enumerate(columns):
                if column.shared:
                    group_commands.append(("SPAN", (index, start), (index, end)))
    return context.table(
        data,
        column_widths,
        header_rows=header_rows,
        padding=(3, 2),
        commands=[*headers.commands, *group_commands],
    )
