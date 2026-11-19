"""Constructed calibration curves, reporting ranges and model registration."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import ArrayLike, NDArray

from elisaflow.calibration.errors import CalibrationDataError, CalibrationFitError
from elisaflow.calibration.fitting import CalibrationModel, FitDiagnostics, fit_curve
from elisaflow.calibration.logistic import FivePLModel, FivePLParameters, FourPLModel, FourPLParameters, LogisticParameters
from elisaflow.calibration.polynomial import CubicModel, CubicParameters, QuadraticModel, QuadraticParameters
from elisaflow.calibration.power_law import PowerLawModel, PowerLawParameters
from elisaflow.models import ResultStatus

if TYPE_CHECKING:
    from elisaflow.configuration.assay_settings import CalibrationSettings

CalibrationParameters = LogisticParameters | QuadraticParameters | CubicParameters | PowerLawParameters
NumericResult = float | NDArray[np.float64]


def _numeric_result(values: NDArray[np.float64]) -> NumericResult:
    return float(values) if values.ndim == 0 else values


@dataclass(frozen=True)
class CalibratedValue:
    """Concentration and reporting qualification before patient dilution."""

    concentration: float
    status: ResultStatus
    limit: float = float("nan")


@dataclass(frozen=True)
class CalibrationCurve[P]:
    model: CalibrationModel[P]
    parameters: P
    concentration_range: tuple[float, float] | None = None
    diagnostics: FitDiagnostics | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "model", self.model.for_curve(self.parameters, self.concentration_range))

    @property
    def was_fitted(self) -> bool:
        return self.diagnostics is not None

    def calculate_extinction(self, concentrations: ArrayLike) -> NumericResult:
        return self.model.calculate_extinctions(concentrations, self.parameters)

    def calculate_concentration(self, extinctions: ArrayLike) -> NumericResult:
        """Invert responses, limiting numeric results to the reporting range."""
        responses = np.asarray(extinctions, dtype=float)
        concentrations = np.asarray(self.model.calculate_concentrations(responses, self.parameters))

        if self.concentration_range is not None:
            lower, upper = self.concentration_range
            low_response, high_response = self.calculate_extinction([lower, upper])
            # Exact endpoint responses have known concentrations, independent of inversion round-off.
            if np.isfinite([low_response, high_response]).all() and low_response != high_response:
                concentrations = np.where(responses == low_response, lower, concentrations)
                concentrations = np.where(responses == high_response, upper, concentrations)
            # Include endpoints despite small floating-point inversion errors.
            below = (concentrations < lower) & ~np.isclose(concentrations, lower, rtol=1e-12, atol=0)
            above = (concentrations > upper) & ~np.isclose(concentrations, upper, rtol=1e-12, atol=0)
            concentrations = np.where(below | above, np.nan, concentrations)
        return _numeric_result(concentrations)

    def calculate_results(self, extinctions: ArrayLike) -> tuple[CalibratedValue, ...]:
        """Invert and qualify responses on this curve's reporting range."""
        responses = np.asarray(extinctions, dtype=float)
        if responses.ndim != 1:
            raise ValueError("Kalibratie-invoer moet een eendimensionale reeks zijn.")
        concentrations = self.calculate_concentration(responses)
        bounds = None
        if self.concentration_range is not None:
            lower, upper = self.concentration_range
            low_response, high_response = self.calculate_extinction([lower, upper])
            bounds = lower, upper, low_response, high_response
        results = []
        for response, concentration in zip(responses, concentrations, strict=True):
            status = ResultStatus.IN_RANGE if np.isfinite(concentration) else ResultStatus.UNAVAILABLE
            limit = float("nan")
            if bounds is not None and not np.isfinite(concentration) and np.isfinite(response):
                lower, upper, low_response, high_response = bounds
                increasing = low_response < high_response
                below = response < low_response if increasing else response > low_response
                above = response > high_response if increasing else response < high_response
                if lower > 0 and below:
                    status, limit = ResultStatus.BELOW_RANGE, lower
                if above:
                    status, limit = ResultStatus.ABOVE_RANGE, upper
            results.append(CalibratedValue(float(concentration), status, limit))
        return tuple(results)

    def calculate_bound(self, response: float, *, above: bool) -> CalibratedValue:
        """Qualify a strict response bound without treating it as a measurement."""
        unavailable = CalibratedValue(float("nan"), ResultStatus.UNAVAILABLE)
        if self.concentration_range is None or not np.isfinite(response):
            return unavailable
        lower, upper = self.concentration_range
        low_response, high_response = self.calculate_extinction([lower, upper])
        if not np.isfinite([low_response, high_response]).all() or low_response == high_response:
            return unavailable
        above_concentration = above == (low_response < high_response)
        endpoint = high_response if above_concentration else low_response
        beyond = response > endpoint if above else response < endpoint
        # Allow numerical fitting error when the bound is the endpoint standard.
        if not beyond and not np.isclose(response, endpoint, rtol=1e-8, atol=0):
            return unavailable
        if not above_concentration and lower <= 0:
            return unavailable
        return CalibratedValue(
            float("nan"),
            ResultStatus.ABOVE_RANGE if above_concentration else ResultStatus.BELOW_RANGE,
            upper if above_concentration else lower,
        )


def create_calibration_curve(
    settings: CalibrationSettings,
    *,
    concentrations: ArrayLike | None = None,
    extinction_values: ArrayLike | None = None,
) -> CalibrationCurve[CalibrationParameters]:
    """Construct a fixed curve or fit all standards with a bounded reporting range."""
    model = CURVE_MODELS[settings.fit_model].create_model()
    if settings.parameters is not None:
        return CalibrationCurve(model, settings.parameters, settings.concentration_range)
    if concentrations is None or extinction_values is None:
        raise CalibrationDataError("Curvepassing vereist standaardconcentraties en extinctiewaarden.")

    concentrations = np.asarray(concentrations, dtype=float)
    fit = fit_curve(concentrations, extinction_values, model, settings.weighting_exponent)
    if not fit.success:
        raise CalibrationFitError(f"De kalibratie convergeerde niet: {fit.message}")
    positive = concentrations[concentrations > 0]
    if not positive.size:
        raise CalibrationDataError("Een rapportagebereik vereist behouden positieve standaardconcentraties.")
    concentration_range = (float(np.min(positive)), float(np.max(positive)))
    if settings.concentration_range is not None:
        concentration_range = (
            max(concentration_range[0], settings.concentration_range[0]),
            min(concentration_range[1], settings.concentration_range[1]),
        )
    if concentration_range[0] >= concentration_range[1]:
        raise CalibrationDataError("Het ingestelde concentratiebereik overlapt niet met het bereik van de overgebleven standaarden.")
    return CalibrationCurve(model, fit.parameters, concentration_range, fit.diagnostics)


@dataclass(frozen=True)
class CurveModelDefinition[P]:
    create_model: Callable[[], CalibrationModel[P]]
    parse_parameters: Callable[[Mapping[str, object]], P]
    parameter_type: type[P]


LogisticModelDefinition = CurveModelDefinition[FourPLParameters] | CurveModelDefinition[FivePLParameters]
CalibrationModelDefinition = (
    LogisticModelDefinition | CurveModelDefinition[QuadraticParameters] | CurveModelDefinition[CubicParameters] | CurveModelDefinition[PowerLawParameters]
)


CURVE_MODELS: Mapping[str, CalibrationModelDefinition] = MappingProxyType(
    {
        "4pl": CurveModelDefinition(FourPLModel, FourPLParameters.from_mapping, FourPLParameters),
        "5pl": CurveModelDefinition(FivePLModel, FivePLParameters.from_mapping, FivePLParameters),
        "quadratic": CurveModelDefinition(QuadraticModel, QuadraticParameters.from_mapping, QuadraticParameters),
        "cubic": CurveModelDefinition(CubicModel, CubicParameters.from_mapping, CubicParameters),
        "power_law": CurveModelDefinition(PowerLawModel, PowerLawParameters.from_mapping, PowerLawParameters),
    }
)
