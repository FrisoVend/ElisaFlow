"""Draw fixed and fitted calibration curves with available standard measurements."""

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from elisaflow.calibration.fitting import CurveModel
from elisaflow.calibration.logistic import FivePLParameters, FourPLParameters
from elisaflow.formatting import format_curve_parameters, format_formula, format_number
from elisaflow.formulas.expressions import Formula
from elisaflow.models import Sample, SampleType, Statistics
from elisaflow.sample_selection import select_samples


def plot_calibration[P](
    samples: tuple[Sample, ...],
    model: CurveModel[P],
    parameters: P,
    *,
    input_statistics: tuple[Statistics, ...],
    model_name: str | None = None,
    r_squared: float | None = None,
    weighting_exponent: float | None = None,
    extinction_label: str = "Extinctie",
    concentration_unit: str | None = None,
    figsize: tuple[float, float] = (4.2, 2.6),
    preprocessing: Formula | None = None,
    curve_was_fitted: bool = True,
    concentration_range: tuple[float, float] | None = None,
) -> Figure:
    """Return a figure with sample means, available extinction SDs, and the fit.

    Fit-input summaries account for curve preprocessing. Missing SDs leave the
    corresponding mean visible without an error bar. The logarithmic
    concentration axis displays positive standards; zero standards remain in
    the analysis and calibration table. The fit legend renders the R² supplied
    by analysis, or n/a when it is unavailable. A supplied model name adds the
    model and coefficients to the legend, inside the existing plot area.
    When supplied, the weighting exponent appears below the coefficients.
    Preprocessing uses the editor's formula notation inside the legend box.
    Fixed curves need no standards and show no fit statistic. Their
    display domain uses the configured range, available concentrations, or a
    default display interval; it never changes analytical reporting limits.
    Input statistics are required and used as supplied, including unavailable values.
    """
    standards = select_samples(samples, SampleType.STANDARD)
    if len(input_statistics) != len(standards):
        raise ValueError("Elke standaard vereist precies één invoerstatistiek voor de kalibratiegrafiek.")
    nominal = np.array([sample.concentration for sample in standards], dtype=float)
    extinctions = np.array([summary.mean for summary in input_statistics])
    valid = np.isfinite(nominal) & (nominal > 0) & np.isfinite(extinctions)
    if curve_was_fitted and not valid.any():
        raise ValueError("Een logaritmische kalibratiegrafiek vereist eindige standaardgemiddelden bij positieve concentraties.")

    r_squared_label = format_number(r_squared, 4, missing="n/a")

    figure = Figure(figsize=figsize, layout="constrained")
    FigureCanvasAgg(figure)
    axes = figure.add_subplot()
    axes.set_xscale("log")
    lower, upper = _display_range(samples, nominal[valid], parameters, concentration_range)
    grid = np.geomspace(lower, upper, 600)
    axes.plot(
        grid,
        model.calculate_extinctions(grid, parameters),
        color="#087C8E",
        linewidth=1.3,
        label=f"Fit (R² = {r_squared_label})" if curve_was_fitted else "Vaste curve",
    )
    if valid.any():
        axes.scatter(
            nominal[valid],
            extinctions[valid],
            s=15,
            color="#263442",
            zorder=3,
            label="Standaardgemiddelden",
        )
    if input_statistics:
        sd = np.array([summary.sd for summary in input_statistics])
        with_sd = valid & np.isfinite(sd) & (sd >= 0)
        if with_sd.any():
            axes.errorbar(
                nominal[with_sd],
                extinctions[with_sd],
                yerr=sd[with_sd],
                fmt="none",
                ecolor="#12A59F",
                elinewidth=0.8,
                capsize=2,
                label="±1 SD",
            )
    concentration_label = "Standaard concentratie" if curve_was_fitted else "Concentratie"
    if concentration_unit:
        concentration_label += f" [{concentration_unit}]"
    axes.set_xlabel(concentration_label, fontsize=8)
    axes.set_ylabel(extinction_label, fontsize=8)
    axes.tick_params(labelsize=7)
    axes.grid(color="#CBD2D9", linewidth=0.4, alpha=0.6)
    axes.set_axisbelow(True)
    legend_title = f"{model_name.upper()}\n{format_curve_parameters(parameters, lines=2)}" if model_name else None
    if curve_was_fitted and weighting_exponent is not None:
        weighting_label = f"Wegingsexponent: {format_number(weighting_exponent)}"
        legend_title = f"{legend_title}\n{weighting_label}" if legend_title else weighting_label
    if preprocessing is not None:
        preprocessing_label = f"Voorbewerking: {format_formula(preprocessing, 'Monster')}"
        legend_title = f"{legend_title}\n{preprocessing_label}" if legend_title else preprocessing_label
    legend = axes.legend(
        title=legend_title,
        title_fontsize=6.5,
        fontsize=6.5,
        frameon=bool(legend_title),
        facecolor="white",
        edgecolor="#CBD2D9",
        framealpha=0.92,
        alignment="left",
        labelspacing=0.25,
        loc="best",
    )
    legend.get_title().set_parse_math(False)
    for spine in axes.spines.values():
        spine.set_color("#CBD2D9")
    return figure


def _display_range(samples: tuple[Sample, ...], nominal: np.ndarray, parameters: object, configured: tuple[float, float] | None) -> tuple[float, float]:
    """Choose a positive plotting interval without modifying the calculation."""
    if configured is not None:
        lower, upper = configured
        if not np.isfinite([lower, upper]).all() or lower < 0 or lower >= upper:
            raise ValueError("Het concentratiebereik van de grafiek moet eindig, niet-negatief en oplopend zijn.")
    else:
        values = nominal if nominal.size else np.array([
            well.result / sample.dilution if sample.sample_type == SampleType.UNKNOWN else well.result
            for sample in samples for well in sample.result_measurements if np.isfinite(well.result) and well.result > 0
        ])
        if values.size:
            lower, upper = float(values.min()), float(values.max())
        elif isinstance(parameters, (FourPLParameters, FivePLParameters)):
            lower, upper = parameters.C / 100, parameters.C * 100
        else:
            lower, upper = 0.01, 100.0
        if lower == upper:
            lower, upper = lower / 10, upper * 10
    return (lower if lower > 0 else upper / 1000), upper
