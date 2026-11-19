"""Fit through explicit model adapters; share validation, weighting and diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
from numpy.typing import ArrayLike, NDArray

from elisaflow.calibration.errors import CalibrationDataError, CalibrationFitError
from elisaflow.configuration.validation import finite_number

NumericResult = float | NDArray[np.float64]


class CurveModel[P](Protocol):
    """Prediction and inversion, independently of how a curve was fitted."""

    def calculate_extinctions(self, x: ArrayLike, params: P) -> NumericResult: ...

    def calculate_concentrations(self, y: ArrayLike, params: P) -> NumericResult: ...


@dataclass(frozen=True)
class ParameterFit[P]:
    """A solver outcome and the model domain used for its diagnostics."""

    model: CurveModel[P]
    parameters: P
    success: bool
    message: str


class FittingModel[P](CurveModel[P], Protocol):
    minimum_levels: int

    def fit(self, x: NDArray[np.float64], y: NDArray[np.float64], weights: NDArray[np.float64]) -> ParameterFit[P]: ...


class CalibrationModel[P](FittingModel[P], Protocol):
    def for_curve(self, params: P, concentration_range: tuple[float, float] | None) -> CalibrationModel[P]:
        """Select and validate the inverse domain without mutating this model."""
        ...


class NonlinearModel[P](CurveModel[P], Protocol):
    """Initialization, bounded parameters and scaling for nonlinear fitting."""

    def starting_points(self, x: ArrayLike, y: ArrayLike) -> list[P]: ...

    def pack(self, params: P) -> NDArray[np.float64]: ...

    def unpack(self, values: ArrayLike) -> P: ...

    def scale_parameters(self, params: P, response_scale: float) -> P: ...

    def bounds_for_start(self, x: ArrayLike, y: ArrayLike, start: P) -> tuple[NDArray[np.float64], NDArray[np.float64]]: ...


def response_weights(responses: ArrayLike, weighting_exponent: float = 0.0) -> NDArray[np.float64]:
    """Compute observed-response weights y**(-p); p=0 is unweighted."""
    exponent = finite_number(weighting_exponent, "weighting_exponent")
    y = np.asarray(responses, dtype=float)
    if not np.isfinite(y).all():
        raise CalibrationDataError("Weging vereist eindige responswaarden.")
    if exponent == 0:
        return np.ones_like(y)
    if np.any(y <= 0):
        raise CalibrationDataError("Een wegingsexponent ongelijk aan nul vereist positieve responswaarden.")
    with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
        weights = np.power(y, -exponent)
    if not np.isfinite(weights).all() or np.any(weights <= 0):
        raise CalibrationDataError("Gewichten moeten eindig en positief zijn; de wegingsexponent levert gewichten buiten dit bereik op.")
    return weights


@dataclass(frozen=True)
class FitDiagnostics:
    """Own read-only copies of the arrays supplied by the fitter."""

    residuals: NDArray[np.float64]
    predictions: NDArray[np.float64]
    back_calculated: NDArray[np.float64]
    relative_errors: NDArray[np.float64]
    weights: NDArray[np.float64]

    def __post_init__(self) -> None:
        for name in ("residuals", "predictions", "back_calculated", "relative_errors", "weights"):
            values = np.array(getattr(self, name), dtype=float, copy=True)
            values.setflags(write=False)
            object.__setattr__(self, name, values)


@dataclass(frozen=True)
class FitResult[P]:
    parameters: P
    diagnostics: FitDiagnostics
    success: bool
    cost: float
    message: str


def fit_curve[P](
    concentrations: ArrayLike,
    extinctions: ArrayLike,
    model: FittingModel[P],
    weighting_exponent: float = 0.0,
) -> FitResult[P]:
    """Validate observations, fit with the model's adapter, and record diagnostics.

    Every retained replicate contributes independently. Diagnostics use the
    fitting domain, before applying any reporting range. Nonlinear adapters use
    multistart; polynomial adapters solve their monotonic constrained problem.
    A finite unconverged attempt returns success=False; no usable attempt raises
    CalibrationFitError. Unexpected implementation errors propagate unchanged.
    """
    x, y = np.asarray(concentrations, dtype=float), np.asarray(extinctions, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or x.shape != y.shape:
        raise CalibrationDataError("Standaarden moeten overeenkomende eendimensionale reeksen zijn.")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise CalibrationDataError("Standaardconcentraties en responswaarden moeten eindig zijn.")
    if np.any(x < 0):
        raise CalibrationDataError("Standaardconcentraties mogen niet negatief zijn.")
    if np.unique(x).size < model.minimum_levels:
        raise CalibrationDataError(f"Curvepassing vereist ten minste {model.minimum_levels} verschillende concentraties.")
    if np.ptp(y) == 0:
        raise CalibrationDataError("Curvepassing vereist verschillende responswaarden.")
    weights = response_weights(y, weighting_exponent)

    fitted = model.fit(x, y, weights)
    predictions = np.asarray(fitted.model.calculate_extinctions(x, fitted.parameters))
    if not np.isfinite(predictions).all():
        raise CalibrationFitError("De kalibratie leverde niet-eindige voorspellingen in de oorspronkelijke responseenheid op.")
    back_calculated = np.asarray(fitted.model.calculate_concentrations(y, fitted.parameters))
    relative_errors = np.full_like(x, np.nan)
    np.divide(back_calculated - x, x, out=relative_errors, where=x > 0)
    diagnostics = FitDiagnostics(y - predictions, predictions, back_calculated, relative_errors, weights)
    return FitResult(fitted.parameters, diagnostics, fitted.success, float(0.5 * np.sum(weights * diagnostics.residuals**2)), fitted.message)

def fit_nonlinear[P](
    model: NonlinearModel[P], x: NDArray[np.float64], y: NDArray[np.float64], weights: NDArray[np.float64],
) -> ParameterFit[P]:
    """Normalize, try model-provided starts, and restore response units."""
    from scipy.optimize import least_squares

    response_scale = float(np.max(np.abs(y)))
    normalized_y = y / response_scale
    sqrt_weights = np.sqrt(weights / np.max(weights))
    starts = model.starting_points(x, normalized_y)
    if not starts:
        raise CalibrationDataError("Curvepassing vereist ten minste één beginpunt.")

    def residuals(values: NDArray[np.float64]) -> NDArray[np.float64]:
        predictions = model.calculate_extinctions(x, model.unpack(values))
        return sqrt_weights * (predictions - normalized_y)

    result = None
    last_failure = "Geen eindig optimalisatieresultaat."
    for start in starts:
        try:
            candidate = least_squares(
                residuals,
                model.pack(start),
                bounds=model.bounds_for_start(x, normalized_y, start),
                x_scale="jac",
                ftol=1e-12,
                xtol=1e-12,
                gtol=1e-12,
                max_nfev=5000,
            )
            if not np.isfinite(candidate.cost) or not np.isfinite(candidate.x).all() or not np.isfinite(candidate.fun).all():
                last_failure = "De optimalisatie leverde niet-eindige parameters, residuen of kosten op."
                continue
            candidate_params = model.scale_parameters(model.unpack(candidate.x), response_scale)
            if not np.isfinite(model.calculate_extinctions(x, candidate_params)).all():
                last_failure = "De voorspellingen na curvepassing zijn niet-eindig in de oorspronkelijke responseenheid."
                continue
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError) as exc:
            last_failure = str(exc)
            continue
        # A converged fit takes priority over a lower-cost unfinished attempt.
        if result is None or (not candidate.success, candidate.cost) < (not result.success, result.cost):
            result = candidate
            params = candidate_params

    if result is None:
        raise CalibrationFitError(f"De kalibratie convergeerde vanuit geen enkel beginpunt: {last_failure}")
    return ParameterFit(model, params, bool(result.success), str(result.message))
