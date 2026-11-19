"""Offset power-law calibration on the original weighted response scale."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import astuple, dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from elisaflow.calibration.errors import CalibrationDataError
from elisaflow.calibration.fitting import ParameterFit, fit_nonlinear
from elisaflow.configuration.validation import finite_number

NumericResult = float | NDArray[np.float64]


@dataclass(frozen=True)
class PowerLawParameters:
    """Coefficients of y = C + A*x**B in the original units."""

    A: float
    B: float
    C: float

    def __post_init__(self) -> None:
        for name in self.__dataclass_fields__:
            object.__setattr__(self, name, finite_number(getattr(self, name), name))
        if self.A == 0 or self.B == 0:
            raise ValueError("A en B van een machtsfunctie mogen niet nul zijn.")

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> PowerLawParameters:
        if not isinstance(values, Mapping) or set(values) != {"A", "B", "C"}:
            raise ValueError("De parameters van een machtsfunctie moeten A, B en C zijn.")
        return cls(**values)


class PowerLawModel:
    minimum_levels = 3

    def for_curve(self, params: PowerLawParameters, concentration_range: tuple[float, float] | None) -> PowerLawModel:
        if params.B < 0 and concentration_range is not None and concentration_range[0] <= 0:
            raise CalibrationDataError("Een machtsfunctie met negatieve exponent vereist positieve concentraties.")
        return self

    def calculate_extinctions(self, concentrations: ArrayLike, params: PowerLawParameters) -> NumericResult:
        x = np.asarray(concentrations, dtype=float)
        if np.any(x < 0):
            raise ValueError("Concentraties mogen niet negatief zijn.")
        if params.B < 0 and np.any(x == 0):
            raise ValueError("Een machtsfunctie met negatieve exponent vereist positieve concentraties.")
        with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
            # Combine units in log space without first forming x**B.
            values = params.C + np.sign(params.A) * np.exp(np.log(abs(params.A)) + params.B * np.log(x))
        return float(values) if values.ndim == 0 else values

    def calculate_concentrations(self, extinctions: ArrayLike, params: PowerLawParameters) -> NumericResult:
        y = np.asarray(extinctions, dtype=float)
        with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
            delta = y - params.C
            # A changes with concentration units; only response units set this tolerance.
            tolerance = 32 * np.finfo(float).eps * np.maximum(np.abs(y), abs(params.C))
            outside = np.sign(params.A) * delta < 0
            delta = np.where(outside & (np.abs(delta) <= tolerance), 0, delta)
            valid = ((np.sign(params.A) * delta > 0) | ((params.B > 0) & (delta == 0))) & np.isfinite(y) & np.isfinite(delta)
            values = np.exp((np.log(np.abs(delta)) - np.log(abs(params.A))) / params.B)
        values = np.where(valid & np.isfinite(values), values, np.nan)
        return float(values) if values.ndim == 0 else values

    def fit(
        self, x: NDArray[np.float64], y: NDArray[np.float64], weights: NDArray[np.float64],
    ) -> ParameterFit[PowerLawParameters]:
        x_scale = float(x.max())
        fitted = fit_nonlinear(self, x / x_scale, y, weights)
        params = fitted.parameters
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            amplitude = np.sign(params.A) * np.exp(np.log(abs(params.A)) - params.B * np.log(x_scale))
        try:
            restored = PowerLawParameters(amplitude, params.B, params.C)
        except ValueError as error:
            raise CalibrationDataError("Machtsfunctieparameters moeten eindig en niet nul zijn in de oorspronkelijke concentratie-eenheid.") from error
        return ParameterFit(self, restored, fitted.success, fitted.message)

    def starting_points(self, concentrations: ArrayLike, extinctions: ArrayLike) -> list[PowerLawParameters]:
        """Estimate A and C linearly for several signed exponent seeds."""
        x, y = np.asarray(concentrations, dtype=float), np.asarray(extinctions, dtype=float)
        limit = self._exponent_limit(x)
        magnitudes = dict.fromkeys(min(value, limit / 2) for value in (0.25, 0.5, 1, 2, 4))
        signs = (1, -1) if np.all(x > 0) else (1,)
        candidates = []
        for sign in signs:
            for magnitude in magnitudes:
                b = sign * magnitude
                with np.errstate(over="ignore", under="ignore", invalid="ignore"):
                    design = np.column_stack((x ** b, np.ones_like(x)))
                if not np.isfinite(design).all():
                    continue
                try:
                    (a, c), _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
                    if rank != 2:
                        continue
                    params = PowerLawParameters(a, b, c)
                except (ValueError, np.linalg.LinAlgError):
                    continue
                cost = float(np.sum((self.calculate_extinctions(x, params) - y) ** 2))
                if np.isfinite(cost):
                    candidates.append((cost, params))
        if not candidates:
            raise CalibrationDataError("Er konden geen eindige beginwaarden voor de machtsfunctie worden bepaald.")
        return [params for _, params in sorted(candidates, key=lambda candidate: candidate[0])]

    @staticmethod
    def _exponent_limit(x: NDArray[np.float64]) -> float:
        log_span = float(np.max(np.abs(np.log(x[x > 0]))))
        # Keep every trial's powered concentrations finite after scaling x.
        return min(100.0, 500.0 / log_span) if log_span else 100.0

    def pack(self, params: PowerLawParameters) -> NDArray[np.float64]:
        return np.asarray(astuple(params))

    def unpack(self, values: ArrayLike) -> PowerLawParameters:
        return PowerLawParameters(*np.asarray(values, dtype=float))

    def scale_parameters(self, params: PowerLawParameters, response_scale: float) -> PowerLawParameters:
        return PowerLawParameters(params.A * response_scale, params.B, params.C * response_scale)

    def bounds_for_start(
        self, concentrations: ArrayLike, extinctions: ArrayLike, start: PowerLawParameters,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        limit = self._exponent_limit(np.asarray(concentrations, dtype=float))
        a_bounds = (0, np.inf) if start.A > 0 else (-np.inf, 0)
        b_bounds = (1e-4, limit) if start.B > 0 else (-limit, -1e-4)
        return np.array([a_bounds[0], b_bounds[0], -np.inf]), np.array([a_bounds[1], b_bounds[1], np.inf])
