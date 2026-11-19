"""Weighted polynomial calibration and inversion on a monotonic dose interval."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import astuple, dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from elisaflow.calibration.errors import CalibrationDataError, CalibrationFitError
from elisaflow.calibration.fitting import ParameterFit
from elisaflow.configuration.validation import finite_number


@dataclass(frozen=True)
class QuadraticParameters:
    """Coefficients of y = A + B*x + C*x**2 in the original units."""

    A: float
    B: float
    C: float

    def __post_init__(self) -> None:
        for name in self.__dataclass_fields__:
            object.__setattr__(self, name, finite_number(getattr(self, name), name))
        if not any(astuple(self)[1:]):
            raise ValueError("Een polynomiale kalibratie mag niet constant zijn.")

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> QuadraticParameters:
        if not isinstance(values, Mapping) or set(values) != {"A", "B", "C"}:
            raise ValueError("De kwadratische parameters moeten A, B en C zijn.")
        return cls(**values)


@dataclass(frozen=True)
class CubicParameters:
    """Coefficients of y = A + B*x + C*x**2 + D*x**3."""

    A: float
    B: float
    C: float
    D: float

    def __post_init__(self) -> None:
        for name in self.__dataclass_fields__:
            object.__setattr__(self, name, finite_number(getattr(self, name), name))
        if not any(astuple(self)[1:]):
            raise ValueError("Een polynomiale kalibratie mag niet constant zijn.")

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> CubicParameters:
        if not isinstance(values, Mapping) or set(values) != {"A", "B", "C", "D"}:
            raise ValueError("De kubische parameters moeten A, B, C en D zijn.")
        return cls(**values)


NumericResult = float | NDArray[np.float64]


class _PolynomialModel[P: (QuadraticParameters, CubicParameters)]:
    parameter_type: type[P]
    minimum_levels: int

    def __init__(self, concentration_range: tuple[float, float] | None = None):
        if concentration_range is not None:
            lower, upper = concentration_range
            if not np.isfinite([lower, upper]).all() or lower < 0 or upper <= lower:
                raise ValueError("Het polynomiale concentratiebereik moet eindig, niet-negatief en oplopend zijn.")
        self.concentration_range = concentration_range

    def bind_concentration_range(self, concentration_range: tuple[float, float]) -> _PolynomialModel[P]:
        """Return a separate model so fitting never changes an existing curve."""
        return type(self)(concentration_range)

    def for_curve(self, params: P, concentration_range: tuple[float, float] | None) -> _PolynomialModel[P]:
        if concentration_range is None:
            raise CalibrationDataError("Polynomiale kalibratie vereist een concentratiebereik om een monotone tak te selecteren.")
        model = self.bind_concentration_range(concentration_range)
        model.validate_parameters(params)
        return model

    def fit(self, x: NDArray[np.float64], y: NDArray[np.float64], weights: NDArray[np.float64]) -> ParameterFit[P]:
        """Polynomial adapter: fit and diagnose on the retained standard domain."""
        lower, upper = float(x.min()), float(x.max())
        x_scale = upper
        response_scale = float(np.max(np.abs(y)))
        normalized_x, normalized_y = x / x_scale, y / response_scale
        weights = np.sqrt(weights)
        weights = weights / weights.max()
        coefficients = self._solve_monotonic(normalized_x, normalized_y, weights, lower / x_scale, 1.0)
        coefficients *= response_scale
        # Repeated division avoids forming x_scale**degree, which can overflow.
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            for index in range(1, len(coefficients)):
                for _ in range(index):
                    coefficients[index] /= x_scale
        if not np.isfinite(coefficients).all():
            raise CalibrationDataError("Polynomiale coëfficiënten moeten eindig zijn in de oorspronkelijke concentratie-eenheid.")
        return ParameterFit(
            self.bind_concentration_range((lower, upper)), self.unpack(coefficients), True,
            "Monotonic weighted polynomial least-squares fit.",
        )

    def pack(self, params: P) -> NDArray[np.float64]:
        return np.array(astuple(params))

    def unpack(self, values: ArrayLike) -> P:
        return self.parameter_type(*np.asarray(values, dtype=float))

    @staticmethod
    def _minimum_derivative(coefficients: NDArray[np.float64], lower: float, upper: float) -> tuple[float, float]:
        """Find the exact minimum slope, including a cubic derivative's vertex."""
        derivative = np.polynomial.polynomial.polyder(coefficients)
        points = [lower, upper]
        if len(derivative) == 3 and derivative[2] > 0:
            vertex = -derivative[1] / (2 * derivative[2])
            if lower < vertex < upper:
                points.append(float(vertex))
        return min((float(np.polynomial.polynomial.polyval(point, derivative)), point) for point in points)

    def _solve_monotonic(
        self, x: NDArray[np.float64], y: NDArray[np.float64], sqrt_weights: NDArray[np.float64],
        lower: float, upper: float,
    ) -> NDArray[np.float64]:
        """Always solve the same constrained problem in both monotonic directions."""
        from scipy.optimize import minimize

        objective_scale = float(np.linalg.norm(y * sqrt_weights))
        design = np.polynomial.polynomial.polyvander(x, self.minimum_levels - 1) * sqrt_weights[:, None] / objective_scale
        column_norms = np.linalg.norm(design, axis=0)
        try:
            if np.any(column_norms == 0) or np.linalg.matrix_rank(design / column_norms) < self.minimum_levels:
                raise CalibrationDataError("De standaarden ondersteunen geen numeriek stabiele polynomiale curvepassing; gebruik een lagere graad of grotere concentratieverschillen.")
            # QR conditions the optimization coordinates; every fit still runs
            # through the constrained optimizer, including exact monotonic data.
            orthogonal, triangular = np.linalg.qr(design, mode="reduced")
            transform = np.linalg.solve(triangular, np.eye(self.minimum_levels))
        except np.linalg.LinAlgError as error:
            raise CalibrationFitError("De monotone polynomiale kalibratie kon het gewogen kleinste-kwadratenprobleem niet opstellen.") from error

        def slope(values):
            return self._minimum_derivative(transform @ values, lower, upper)[0]

        def slope_gradient(values):
            _, point = self._minimum_derivative(transform @ values, lower, upper)
            derivative = np.array([0.0, *(power * point ** (power - 1) for power in range(1, len(values)))])
            return derivative @ transform

        candidates = []
        for direction in (1, -1):
            target = orthogonal.T @ (direction * y * sqrt_weights / objective_scale)
            start = np.zeros(self.minimum_levels)
            start[0] = np.average(direction * y, weights=sqrt_weights**2)

            def objective(values):
                residuals = values - target
                return 0.5 * float(residuals @ residuals)

            def gradient(values):
                return values - target

            result = minimize(
                objective, triangular @ start, jac=gradient, method="SLSQP",
                constraints={"type": "ineq", "fun": slope, "jac": slope_gradient},
                options={"ftol": 1e-12, "maxiter": 1000},
            )
            fitted = transform @ result.x
            minimum, _ = self._minimum_derivative(fitted, lower, upper)
            tolerance = 1e-9 * max(1.0, float(np.max(np.abs(fitted[1:]))))
            if not result.success or not np.isfinite(result.x).all() or minimum < -tolerance:
                raise CalibrationFitError("De monotone polynomiale kalibratie convergeerde niet naar een geldige oplossing.")
            # Remove only numerical constraint violations, including double roots
            # that would otherwise appear as two artificial turning points.
            roundoff = 64 * np.finfo(float).eps * max(1.0, float(np.max(np.abs(fitted[1:]))))
            if minimum <= roundoff:
                fitted[1] += max(0.0, -minimum) + roundoff
            candidates.append((objective(triangular @ fitted), direction * fitted))
        return min(candidates, key=lambda candidate: candidate[0])[1]

    def calculate_extinctions(self, concentrations: ArrayLike, params: P) -> NumericResult:
        x = np.asarray(concentrations, dtype=float)
        if np.any(x < 0):
            raise ValueError("Concentraties mogen niet negatief zijn.")
        with np.errstate(over="ignore", invalid="ignore"):
            values = np.asarray(np.polynomial.polynomial.polyval(x, np.trim_zeros(self.pack(params), "b")))
        return float(values) if values.ndim == 0 else values

    def _scaled_coefficients(self, params: P, x_scale: float) -> NDArray[np.float64]:
        coefficients = self.pack(params)
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            for index in range(1, len(coefficients)):
                for _ in range(index):
                    coefficients[index] *= x_scale
        if not np.isfinite(coefficients).all():
            raise CalibrationDataError("Polynomiale parameters veroorzaken numerieke overloop binnen het concentratiebereik.")
        return coefficients

    @staticmethod
    def _turning_points(coefficients: NDArray[np.float64]) -> list[float]:
        derivative = np.polynomial.polynomial.polyder(coefficients)
        derivative /= np.max(np.abs(derivative))
        if len(derivative) == 2 or derivative[-1] == 0:
            roots = [-derivative[0] / derivative[1]] if derivative[1] != 0 else []
        else:
            c, b, a = derivative
            discriminant = b * b - 4 * a * c
            if discriminant <= 0:
                # A double derivative root is a stationary inflection, not a turn.
                return []
            q = -0.5 * (b + np.copysign(np.sqrt(discriminant), b))
            roots = [q / a, c / q]
        return sorted(root for root in roots if np.isfinite(root) and root > 0)

    def validate_parameters(self, params: P) -> None:
        if self.concentration_range is None:
            return
        lower, upper = self.concentration_range
        coefficients = self._scaled_coefficients(params, upper)
        tolerance = 16 * np.finfo(float).eps
        if any(lower / upper + tolerance < root < 1 - tolerance for root in self._turning_points(coefficients)):
            raise CalibrationDataError("De polynomiale kalibratie moet monotoon zijn over het concentratiebereik; verlaag de graad of beperk het bereik.")
        endpoints = self.calculate_extinctions([lower, upper], params)
        if not np.isfinite(endpoints).all() or endpoints[0] == endpoints[1]:
            raise CalibrationDataError("Polynomiale kalibratie vereist verschillende eindige responswaarden op de grenspunten van het bereik.")

    def calculate_concentrations(self, extinctions: ArrayLike, params: P) -> NumericResult:
        y = np.asarray(extinctions, dtype=float)
        values = np.full(y.shape, np.nan)
        if self.concentration_range is None:
            raise CalibrationDataError("Polynomiale inversie vereist een concentratiebereik om een monotone tak te selecteren.")
        from scipy.optimize import brentq

        self.validate_parameters(params)
        lower, upper = self.concentration_range
        coefficients = self._scaled_coefficients(params, upper)
        response_scale = float(np.max(np.abs(coefficients)))
        coefficients /= response_scale
        # Extend only along the same monotonic branch, stopping at a turn.
        turns = self._turning_points(coefficients)
        tolerance = 16 * np.finfo(float).eps
        left = max((point for point in turns if point <= lower / upper + tolerance), default=0.0)
        right = min((point for point in turns if point >= 1 - tolerance), default=np.inf)

        def evaluate(z: float) -> float:
            with np.errstate(over="ignore", invalid="ignore"):
                return float(np.polynomial.polynomial.polyval(z, coefficients))

        y_left = evaluate(left)
        increasing = self.calculate_extinctions(upper, params) > self.calculate_extinctions(lower, params)
        for index in np.ndindex(y.shape):
            if not np.isfinite(y[index]):
                continue
            with np.errstate(over="ignore", invalid="ignore"):
                target = float(y[index] / response_scale)
            if not np.isfinite(target):
                continue
            top = right if np.isfinite(right) else max(1.0, left + 1)
            y_top = evaluate(top)
            while not np.isfinite(right) and (y_top < target if increasing else y_top > target) and top < np.finfo(float).max / 2:
                top *= 2
                y_top = evaluate(top)
            endpoint_tolerance = 4 * np.finfo(float).eps * max(abs(y_left), abs(y_top) if np.isfinite(y_top) else 0, abs(target))
            # Clamp only just-outside responses; preserve small interior roots.
            below = target < min(y_left, y_top)
            above = target > max(y_left, y_top)
            if below or above:
                endpoint = min(y_left, y_top) if below else max(y_left, y_top)
                if abs(target - endpoint) > endpoint_tolerance:
                    continue
                target = endpoint
            if target == y_left:
                values[index] = left * upper
            elif target == y_top:
                values[index] = top * upper
            else:
                values[index] = brentq(
                    lambda z: evaluate(z) - target, left, top,
                    xtol=np.nextafter(0.0, 1.0), rtol=8 * np.finfo(float).eps, maxiter=1000,
                ) * upper
        return float(values) if values.ndim == 0 else values


class QuadraticModel(_PolynomialModel[QuadraticParameters]):
    minimum_levels = 3
    parameter_type = QuadraticParameters


class CubicModel(_PolynomialModel[CubicParameters]):
    minimum_levels = 4
    parameter_type = CubicParameters
