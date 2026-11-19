"""Shared report formatting and layout settings."""

from dataclasses import dataclass
from math import isfinite
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.styles import StyleSheet1
from reportlab.platypus import Table, TableStyle

from elisaflow.branding import BORDER, HEADER, STRIPE
from elisaflow.configuration.assay_settings import ReportSettings
from elisaflow.models import AnalysisResult, MeasurementResult


def cell_commands(horizontal_padding: float = 0, vertical_padding: float = 0, *, valign: str = "TOP") -> list[tuple]:
    return [
        ("VALIGN", (0, 0), (-1, -1), valign),
        ("LEFTPADDING", (0, 0), (-1, -1), horizontal_padding),
        ("RIGHTPADDING", (0, 0), (-1, -1), horizontal_padding),
        ("TOPPADDING", (0, 0), (-1, -1), vertical_padding),
        ("BOTTOMPADDING", (0, 0), (-1, -1), vertical_padding),
    ]


@dataclass(frozen=True)
class ReportArea:
    """Usable report rectangle, in ReportLab points."""

    left: float
    bottom: float
    width: float
    height: float

    def __post_init__(self) -> None:
        if not all(isfinite(value) for value in (self.left, self.bottom, self.width, self.height)) or self.width <= 0 or self.height <= 0:
            raise ValueError("Het rapportoppervlak moet eindige coördinaten en positieve afmetingen hebben.")


@dataclass(frozen=True)
class ReportContext:
    result: AnalysisResult
    styles: StyleSheet1
    content_width: float

    border = colors.HexColor(BORDER)
    stripe = colors.HexColor(STRIPE)
    header = colors.HexColor(HEADER)

    @property
    def report(self) -> ReportSettings:
        return self.result.settings.report

    @property
    def multiple_measurements(self) -> bool:
        return len(self.result.measurement_results) > 1

    def wavelength_suffix(self, measurement: MeasurementResult) -> str:
        return f" ({measurement.label})" if self.multiple_measurements or self.result.settings.measurement.reference_wavelength_nm is not None else ""

    def wavelength_subscript(self, *wavelengths: int) -> str:
        """Format measured or reference-corrected wavelengths as PDF subscripts."""
        labels = (self.result.settings.measurement.wavelength_label(value).removesuffix(" nm") for value in wavelengths)
        return f"<sub>{escape(' / '.join(labels))}</sub>"

    def concentration_label(self, label: str) -> str:
        unit = self.result.settings.result.unit
        return f"{label} [{unit}]" if unit else label

    def extinction_label(self, *, raw: bool = False) -> str:
        """Label the measured scale, including a correction by a zero blank."""
        return "ΔExt" if self.result.blank_value is not None and not raw else "Ext"

    def table(self, data, widths, *, header_rows: int = 1, padding=(2, 3), commands=()) -> Table:
        """Apply shared header, border and padding styles to result tables."""
        table = Table(data, colWidths=widths, repeatRows=header_rows, hAlign="LEFT")
        table.setStyle(
            TableStyle(
                [
                    ("LINEBELOW", (0, header_rows - 1), (-1, header_rows - 1), 0.5, self.border),
                    ("BACKGROUND", (0, 0), (-1, header_rows - 1), self.header),
                    *cell_commands(*padding),
                    *commands,
                ]
            )
        )
        return table
