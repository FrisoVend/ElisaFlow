from datetime import datetime
from math import isnan
from pathlib import Path
from tempfile import NamedTemporaryFile
from xml.sax.saxutils import escape
from zoneinfo import ZoneInfo

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import (
    BaseDocTemplate,
    Flowable,
    Frame,
    HRFlowable,
    Image,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)

from elisaflow.branding import ACCENT, APP_NAME, BORDER, BRAND, HEADER, INK, LOGO_PATH, MUTED
from elisaflow.formatting import (
    format_formula,
    format_number,
    format_protocol_name,
    format_text,
)
from elisaflow.models import (
    AnalysisResult,
    MeasurementResult,
    ValidationStatus,
)
from elisaflow.reporting.calibration import calibration_section
from elisaflow.reporting.common import ReportArea, ReportContext
from elisaflow.reporting.common import cell_commands as _table_cell_commands
from elisaflow.reporting.patients import patient_pages
from elisaflow.reporting.quality_control import validation_table
from elisaflow.runs.analyst import validate_analyst_name


def _format_session_time(value: str) -> str:
    """Display a measurement's timestamp in Dutch local time without a zone label."""
    timestamp = datetime.fromisoformat(value)
    if timestamp.tzinfo is None:
        raise ValueError("Het meettijdstip moet een UTC-offset bevatten om de Nederlandse tijd te bepalen.")
    return timestamp.astimezone(ZoneInfo("Europe/Amsterdam")).strftime("%d-%m-%Y %H:%M:%S")


class PdfReportWriter:
    _INK = colors.HexColor(INK)
    _MUTED = colors.HexColor(MUTED)
    _BRAND = colors.HexColor(BRAND)
    _ACCENT = colors.HexColor(ACCENT)
    _BORDER = colors.HexColor(BORDER)
    _HEADER = colors.HexColor(HEADER)
    _LOGO_PATH = LOGO_PATH
    _FRAME_PADDING = 6

    def __init__(self, result: AnalysisResult, *, analyst: str):
        self.result = result
        self._analyst = validate_analyst_name(analyst)
        self.report = result.settings.report
        self.measurements = tuple(result.measurement_results.values())
        self._multiple_measurements = len(self.measurements) > 1

    def write(self, output_path: Path) -> Path:
        """Build a report in a sibling temporary file, then publish it atomically."""
        file_path = Path(output_path)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        temporary_path = None
        try:
            with NamedTemporaryFile(dir=file_path.parent, prefix=f".{file_path.name}.", suffix=".tmp", delete=False) as temporary:
                temporary_path = Path(temporary.name)

                self._initialize_report(temporary_path)
                self._write_title()
                self._write_metadata()

                for measurement in self.measurements:
                    self._write_extinctions(measurement)
                for measurement in self.measurements:
                    self._write_replicate_exclusions(measurement)
                for measurement in self.measurements:
                    self.story.extend(calibration_section(self.context, measurement))
                patients = patient_pages(self.context, self.report_area, on_page=self._draw_footer)
                self.document.addPageTemplates(patients.template)
                self.story.extend(patients.flowables)

                # Manage the stream ourselves so failed builds close it before cleanup.
                self.document.build(self.story, filename=temporary)
            temporary_path.replace(file_path)
        finally:
            if temporary_path is not None:
                temporary_path.unlink(missing_ok=True)
        return file_path

    def _initialize_report(self, output_path: Path) -> None:
        self.document = BaseDocTemplate(
            str(output_path),
            pagesize=A4,
            leftMargin=12 * mm,
            rightMargin=12 * mm,
            topMargin=14 * mm,
            bottomMargin=14 * mm,
            title=f"{APP_NAME} report | {self.result.metadata.assay_name}",
            author=self._analyst,
        )
        # Keep the standard six-point frame padding around full-width sections.
        self.content_width = self.document.width - 2 * self._FRAME_PADDING
        self.report_area = ReportArea(
            self.document.leftMargin + self._FRAME_PADDING,
            self.document.bottomMargin + self._FRAME_PADDING,
            self.content_width,
            self.document.height - 2 * self._FRAME_PADDING,
        )
        full_width_frame = Frame(
            self.document.leftMargin,
            self.document.bottomMargin,
            self.document.width,
            self.document.height,
            id="report",
        )
        self.document.addPageTemplates(PageTemplate(id="Report", frames=[full_width_frame], onPage=self._draw_footer))
        self.styles = getSampleStyleSheet()
        self.styles["Normal"].fontSize = 8
        self.styles["Normal"].leading = 10
        self.styles["Normal"].textColor = self._INK
        self.styles["Title"].fontName = "Helvetica-Bold"
        self.styles["Title"].fontSize = 17
        self.styles["Title"].leading = 20
        self.styles["Title"].alignment = TA_LEFT
        self.styles["Title"].textColor = self._INK
        self.styles["Title"].spaceAfter = 0
        self.styles["Heading3"].fontName = "Helvetica-Bold"
        self.styles["Heading3"].fontSize = 9
        self.styles["Heading3"].leading = 11
        self.styles["Heading3"].textColor = self._BRAND
        self.styles["Heading3"].spaceBefore = 3 * mm
        self.styles["Heading3"].spaceAfter = 1.5 * mm
        self.styles["Heading3"].keepWithNext = True
        self.styles.add(
            ParagraphStyle(
                "ReportLabel",
                parent=self.styles["Normal"],
                fontName="Helvetica-Bold",
                fontSize=8,
                leading=11,
                textColor=self._BRAND,
                spaceAfter=2 * mm,
            )
        )
        self.styles.add(
            ParagraphStyle(
                "BlockHeading",
                parent=self.styles["Heading3"],
                spaceBefore=0,
            )
        )
        self.styles.add(
            ParagraphStyle(
                "Metadata",
                parent=self.styles["Normal"],
                spaceAfter=1,
            )
        )
        self.context = ReportContext(self.result, self.styles, self.content_width)
        self.story: list[Flowable] = []

    def _draw_footer(self, canvas: Canvas, document: BaseDocTemplate) -> None:
        left = document.leftMargin + self._FRAME_PADDING
        right = A4[0] - document.rightMargin - self._FRAME_PADDING
        canvas.saveState()
        canvas.setStrokeColor(self._BORDER)
        canvas.setLineWidth(0.35)
        canvas.line(left, 10 * mm, right, 10 * mm)
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(self._MUTED)
        label = f"Reinier Diagnostiek | {format_protocol_name(self.result.metadata.assay_name)}"
        available_width = right - left - 30 * mm
        while stringWidth(label, "Helvetica", 7) > available_width:
            label = label[:-2].rstrip() + "\u2026"
        canvas.drawString(left, 7 * mm, label)
        canvas.drawRightString(right, 7 * mm, f"Pagina {document.page}")
        if document.page > 1:
            canvas.setFont("Helvetica-Bold", 7)
            canvas.setFillColor(self._BRAND)
            canvas.drawString(left, A4[1] - 9 * mm, f"{APP_NAME} analyserapport | vervolg")
            canvas.setStrokeColor(self._BORDER)
            canvas.line(left, A4[1] - 11 * mm, right, A4[1] - 11 * mm)
        status_label, status_color = {
            ValidationStatus.PASSED: ("Geldig", colors.darkgreen),
            ValidationStatus.FAILED: ("Ongeldig", colors.firebrick),
            ValidationStatus.NOT_EVALUATED: ("Niet bepaald", colors.dimgrey),
        }[self.result.validation_status]
        canvas.setFont("Helvetica-Bold", 9)
        canvas.setFillColor(status_color)
        canvas.drawRightString(right, A4[1] - 9 * mm, f"Teststatus: {status_label}")
        canvas.restoreState()

    def _write_title(self) -> None:
        """Pair the assay title with the bundled company logo."""
        title = [
            Paragraph(f"{APP_NAME} analyserapport", self.styles["ReportLabel"]),
            Paragraph(escape(format_protocol_name(self.result.metadata.assay_name)), self.styles["Title"]),
        ]
        logo = Image(str(self._LOGO_PATH))
        logo_aspect_ratio = logo.imageWidth / logo.imageHeight
        logo.drawHeight = 20 * mm
        logo.drawWidth = logo.drawHeight * logo_aspect_ratio
        logo.hAlign = "RIGHT"
        logo_width = 22 * mm
        header = Table([[title, logo]], colWidths=[self.content_width - logo_width, logo_width])
        header.setStyle(
            TableStyle(
                [
                    *_table_cell_commands(valign="MIDDLE"),
                    ("ALIGN", (1, 0), (1, 0), "RIGHT"),
                ]
            )
        )
        self.story.extend(
            [
                header,
                Spacer(1, 2.5 * mm),
                HRFlowable(width="100%", thickness=0.8, color=self._ACCENT),
                Spacer(1, 3 * mm),
            ]
        )

    def _write_metadata(self) -> None:
        """Group run and assay details beside the quality-control table."""
        metadata = self.result.metadata
        calculation_details = self.result.calculation_details
        settings = self.result.settings

        base_fields = [
            ("Datum/tijd", _format_session_time(metadata.session_time)),
            ("Analist", self._analyst),
            ("Instrument", format_text(metadata.instrument_name)),
        ]

        def field_paragraph(label: str, value: object, *, markup: bool = False) -> Paragraph:
            contents = str(value) if markup else escape(str(value))
            return Paragraph(
                f'<font color="{self._MUTED.hexval()}">{escape(label)}:</font> {contents}',
                self.styles["Metadata"],
            )

        metadata_block = [Paragraph("Rungegevens", self.styles["BlockHeading"])]
        metadata_block.extend(field_paragraph(label, value) for label, value in base_fields)

        assay_fields = []
        if settings.lot_number is not None:
            assay_fields.append(("Lotnummer kit", settings.lot_number))
        assay_fields.extend(
            (f"Lotnummer {name}", control.lot_number) for name, control in settings.controls.items() if control.lot_number is not None
        )
        assay_fields.append(("Golflengte", " / ".join(measurement.label for measurement in self.measurements)))
        if self._multiple_measurements:
            preferred = settings.measurement.preferred_wavelength_nm
            assay_fields.append(("Voorkeursgolflengte", f"{preferred} nm"))
            cutoff = settings.measurement.cutoff_concentration
            unit = settings.result.unit or ""
            assay_fields.append(("Omschakelconcentratie", f"{cutoff:g} {unit}"))
            for measurement in self.measurements:
                if not measurement.usable_for_results:
                    assay_fields.append((f"Curve {measurement.wavelength_nm} nm", "Niet bruikbaar: kalibratie niet beschikbaar"))
        if calculation_details.formula_settings is not None:
            formula = format_formula(calculation_details.formula_settings.formula, "Monster")
            assay_fields.append(("Berekening", f"{settings.result.label} ({formula})"))

        metadata_block.extend([Spacer(1, 2 * mm), Paragraph("Assaygegevens", self.styles["BlockHeading"])])
        metadata_block.extend(field_paragraph(label, value) for label, value in assay_fields)
        for wavelength, rule in settings.measurement.saturation.items():
            if rule.mode == "calculate":
                target = settings.measurement.wavelength_label(wavelength)
                target_subscript = self.context.wavelength_subscript(wavelength)
                source_subscript = self.context.wavelength_subscript(rule.source_wavelength_nm)
                equation = (
                    f"Ext{target_subscript} = {rule.intercept:.3f} + {rule.slope:.3f} "
                    f"× Ext{source_subscript}"
                )
                metadata_block.append(field_paragraph(f"Verzadiging {target}", equation, markup=True))

        validation_width = self.content_width * (0.65 if self._multiple_measurements else 0.60)
        table = validation_table(self.context, validation_width)
        if table is None:
            self.story.extend(metadata_block)
            return
        gap = 3 * mm
        layout = Table(
            [
                [
                    metadata_block,
                    "",
                    [
                        Paragraph("Kwaliteitscontrole", self.styles["BlockHeading"]),
                        table,
                    ],
                ]
            ],
            colWidths=[self.content_width - validation_width - gap, gap, validation_width],
            hAlign="LEFT",
        )
        layout.setStyle(TableStyle(_table_cell_commands()))
        self.story.append(layout)

    def _write_extinctions(self, measurement: MeasurementResult) -> None:
        """Display measured and reconstructed extinctions in plate coordinates."""
        values = {}
        has_converted = False
        for sample in measurement.samples:
            for well in sample.measurements:
                text = format_number(well.source_extinction, self.report.extinction_decimals, missing="")
                for flag, marker in ((well.saturated, "SAT"), (well.disabled, "UIT"), (well.reference_invalid, "REF")):
                    if flag:
                        text = marker
                if well.converted_extinction is not None and not isnan(well.converted_extinction) and well.included:
                    text = format_number(well.converted_extinction, self.report.extinction_decimals) + "*"
                    has_converted = True
                values[well.well] = text
        columns = range(1, 13)

        self.story.append(Paragraph(
            'Extincties<br/>'
            f'<font name="Helvetica" size="8" color="{MUTED}">{escape(measurement.label)}</font>',
            self.styles["Heading3"],
        ))
        data = [["", *[str(column) for column in columns]]]
        for row in "ABCDEFGH":
            cells = [row]
            for column in columns:
                cells.append(values.get(f"{row}{column}", ""))
            data.append(cells)

        label_width = 8 * mm
        table = Table(
            data,
            colWidths=[label_width] + [(self.content_width - label_width) / 12] * 12,
            rowHeights=4.5 * mm,
            hAlign="LEFT",
            splitByRow=0,
        )
        table.setStyle(
            TableStyle(
                [
                    ("GRID", (1, 1), (-1, -1), 0.35, self._BORDER),
                    ("BACKGROUND", (0, 0), (-1, 0), self._HEADER),
                    ("BACKGROUND", (0, 1), (0, -1), self._HEADER),
                    ("TEXTCOLOR", (0, 0), (-1, -1), self._INK),
                    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                    ("FONTNAME", (0, 1), (0, -1), "Helvetica-Bold"),
                    ("FONTSIZE", (0, 0), (-1, -1), 7.5),
                    ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ]
            )
        )
        self.story.append(table)
        if has_converted:
            self.story.append(Paragraph("* Berekend uit de ingestelde brongolflengte; de oorspronkelijke meting was verzadigd.", self.styles["Normal"]))

    def _write_replicate_exclusions(self, measurement: MeasurementResult) -> None:
        """Identify retained raw readings omitted from calculations and QC."""
        excluded = sorted(
            ((sample, well) for sample in measurement.samples for well in sample.measurements if not well.included),
            key=lambda reading: reading[1].order,
        )
        if not excluded:
            return
        data = [["Well", "Monster", self.context.extinction_label(), "Afwijking (%)", "Reden"]]
        for sample, well in excluded:
            data.append(
                [
                    Paragraph(escape(well.well), self.styles["Normal"]),
                    Paragraph(escape(sample.name), self.styles["Normal"]),
                    format_number(well.extinction, self.report.extinction_decimals),
                    format_number(well.replicate_deviation_percent, self.report.percentage_decimals),
                    Paragraph(escape(well.exclusion_reason), self.styles["Normal"]),
                ]
            )
        table = Table(data, colWidths=[self.content_width * weight for weight in (0.07, 0.17, 0.12, 0.15, 0.49)], repeatRows=1)
        table.setStyle(
            TableStyle(
                [
                    *_table_cell_commands(2, 3),
                    ("BACKGROUND", (0, 0), (-1, 0), self._HEADER),
                    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                    ("FONTSIZE", (0, 0), (-1, -1), 8),
                    ("LINEBELOW", (0, 0), (-1, -1), 0.25, self._BORDER),
                ]
            )
        )
        self.story.extend([Paragraph(f"Uitgesloten replicaten{self.context.wavelength_suffix(measurement)}", self.styles["Heading3"]), table])
