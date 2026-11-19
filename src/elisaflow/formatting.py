"""Display helpers shared by desktop and report interfaces."""
from dataclasses import fields
from math import isfinite
from typing import TYPE_CHECKING

from elisaflow.formulas.names import FormulaNames

if TYPE_CHECKING:
    from elisaflow.calibration.curves import CalibrationParameters
    from elisaflow.configuration.assay_settings import ReportSettings
    from elisaflow.formulas.expressions import Formula
    from elisaflow.models import ExtinctionScale, ValidationMetric


def format_formula(formula: Formula, sample_name: str) -> str:
    """Use the same formula notation as the equation editor."""
    return FormulaNames(sample_name, formula.references).format(formula.expression)


def format_curve_parameters(parameters: CalibrationParameters, *, lines: int = 1) -> str:
    """Keep coefficient precision consistent in metadata and calibration graphs."""
    values = [f"{parameter.name}={getattr(parameters, parameter.name):.3g}" for parameter in fields(parameters)]
    per_line = (len(values) + lines - 1) // lines
    return "\n".join("; ".join(values[index:index + per_line]) for index in range(0, len(values), per_line))


def format_extinction_scale(scale: ExtinctionScale) -> str:
    """Label blank correction; preprocessing is shown with the calibration plot."""
    return "ΔExt" if scale.blank_subtracted is not None else "Ext"


def format_protocol_name(name: str) -> str:
    return name.rsplit(".", 1)[0].replace("_", " ")

def is_finite_number(value: float | None) -> bool:
    """Accept finite scalar numbers, treating None and pandas NA as missing."""
    try:
        return isfinite(value)
    except (TypeError, ValueError):
        return False


def format_number(value: float | None, decimals: int | None = None, missing: str = "-") -> str:
    """Format a finite report value; each table chooses its missing-value marker."""
    if not is_finite_number(value):
        return missing
    return f"{value:g}" if decimals is None else f"{value:.{decimals}f}"

def format_text(value: str | float | None, missing: str = "-") -> str:
    return str(value) if isinstance(value, str) or is_finite_number(value) else missing


def metric_decimals(metric: ValidationMetric, settings: ReportSettings) -> int:
    return settings.extinction_decimals if metric.value == "extinction" else settings.result_decimals


def metric_label(metric: ValidationMetric, label: str = "Resultaat", unit: str | None = None) -> str:
    if metric.value == "extinction":
        return "Ext"
    return f"{label} [{unit}]" if unit else label


def format_result(value: float | None, status: str | None, limit: float | None, decimals: int) -> str:
    """Display a reporting bound without substituting it into calculations."""
    from elisaflow.models import ResultStatus

    qualifier = {
        ResultStatus.BELOW_RANGE.value: "<",
        ResultStatus.ABOVE_RANGE.value: ">",
    }.get(status)
    if qualifier is not None and is_finite_number(limit):
        return f"{qualifier}{format_number(limit, decimals)}"
    return format_number(value, decimals)
