"""Render the calibration section from stored analysis results."""

from io import BytesIO
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.enums import TA_RIGHT
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import Flowable, Image, KeepTogether, Paragraph, Table, TableStyle

from elisaflow.calibration.plot import plot_calibration
from elisaflow.formatting import format_extinction_scale, format_number, is_finite_number
from elisaflow.models import ExtinctionScale, MeasurementResult, Sample, SampleType
from elisaflow.reporting.common import ReportContext, cell_commands
from elisaflow.sample_selection import select_samples


def calibration_section(context: ReportContext, measurement: MeasurementResult) -> list[Flowable]:
    """Render fitted and fixed curves, with recovery diagnostics for fitted ones."""
    flowables = []
    details = measurement.calculation_details
    if details.calibration_error:
        flowables.extend(
            [
                Paragraph(f"Kalibratiecurve{context.wavelength_suffix(measurement)}", context.styles["Heading3"]),
                Paragraph(f"Niet beschikbaar: {escape(details.calibration_error)}", context.styles["Normal"]),
            ]
        )
        return flowables
    curve = details.curve
    if curve is None:
        return flowables
    settings = details.calibration_settings
    scale = measurement.calculation_input_scale
    if settings is None or scale is None:
        raise ValueError("Een kalibratierapport vereist de gebruikte curve-instellingen en invoerschaal.")
    plot_samples = measurement.samples
    if curve.was_fitted:
        plot_samples = tuple(sorted(select_samples(plot_samples, SampleType.STANDARD), key=lambda sample: sample.concentration))

    gap = 3 * mm
    plot_width = (context.content_width - gap) * 0.55
    table_width = context.content_width - plot_width - gap
    plot_height = 62 * mm
    figure = plot_calibration(
        plot_samples,
        curve.model,
        curve.parameters,
        model_name=settings.fit_model,
        r_squared=details.curve_r_squared,
        weighting_exponent=settings.weighting_exponent,
        extinction_label=format_extinction_scale(scale),
        concentration_unit=context.result.settings.result.unit,
        figsize=(plot_width / 72, plot_height / 72),
        input_statistics=tuple(scale.statistics(sample) for sample in select_samples(plot_samples, SampleType.STANDARD)),
        preprocessing=settings.preprocessing,
        curve_was_fitted=curve.was_fitted,
        concentration_range=curve.concentration_range if not curve.was_fitted else None,
    )
    buffer = BytesIO()
    figure.savefig(buffer, format="png", dpi=250)
    buffer.seek(0)
    image = Image(buffer, width=plot_width, height=plot_height)
    if not curve.was_fitted:
        return [KeepTogether([
            Paragraph(f"Vaste kalibratiecurve{context.wavelength_suffix(measurement)}", context.styles["Heading3"]),
            image,
        ])]
    table = calibration_table(context, plot_samples, table_width, scale=scale)
    layout = Table([[image, "", table]], colWidths=[plot_width, gap, table_width], hAlign="LEFT")
    layout.setStyle(TableStyle(cell_commands()))
    flowables.append(
        KeepTogether(
            [
                Paragraph(f"Kalibratiecurve{context.wavelength_suffix(measurement)}", context.styles["Heading3"]),
                layout,
            ]
        )
    )

    return flowables


def calibration_table(context: ReportContext, standards: tuple[Sample, ...], width: float, *, scale: ExtinctionScale) -> Table:
    """Display the already averaged fit inputs and diagnostic recoveries."""
    cell_style = ParagraphStyle("CalibrationCell", parent=context.styles["Normal"], fontSize=7, leading=9)
    number_style = ParagraphStyle("CalibrationNumber", parent=cell_style, alignment=TA_RIGHT)
    recovery_warning_style = ParagraphStyle("CalibrationRecoveryWarning", parent=number_style, textColor=colors.firebrick)
    header_style = ParagraphStyle("CalibrationHeader", parent=number_style, fontName="Helvetica-Bold")
    columns = [
        (context.concentration_label("Standaard Conc"), context.report.calibration_decimals),
        (format_extinction_scale(scale), context.report.extinction_decimals),
        (context.concentration_label("Berekende Conc"), context.report.calibration_decimals),
        ("Terugvinding%", context.report.percentage_decimals),
    ]
    data = [[Paragraph(escape(label), header_style) for label, _ in columns]]
    for sample in standards:
        values = (sample.concentration, scale.statistics(sample).mean, sample.back_calculated_concentration.mean, sample.recovery_percent)
        row = []
        for index, (value, (_, decimals)) in enumerate(zip(values, columns, strict=True)):
            text = format_number(value, decimals)
            style = (
                recovery_warning_style
                if index == 3 and is_finite_number(value) and abs(value - 100) > context.report.recovery_warning_percent
                else number_style
            )
            row.append(Paragraph(text, style))
        data.append(row)
    return context.table(
        data,
        [width * weight for weight in (0.23, 0.19, 0.28, 0.30)],
        commands=[
            ("LINEBELOW", (0, 1), (-1, -1), 0.25, context.border),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, context.stripe]),
        ],
    )
