"""Four- and five-parameter logistic calibration models and their parameters."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from elisaflow.calibration.errors import CalibrationDataError
from elisaflow.calibration.fitting import ParameterFit, fit_nonlinear
from elisaflow.configuration.validation import finite_number

NumericResult = float | NDArray[np.float64]


def _numeric_result(values: NDArray[np.float64]) -> NumericResult:
    return float(values) if values.ndim == 0 else values


@dataclass(frozen=True)
class FourPLParameters:
    """Parameters for D + (A-D) / (1 + (x/C)**B)."""
    A: float
    B: float
    C: float
    D: float

    def __post_init__(self) -> None:
        for name in ("A", "B", "C", "D"):
            object.__setattr__(self, name, finite_number(getattr(self, name), name))
        if self.B == 0:
            raise ValueError("B mag niet nul zijn.")
        if self.C <= 0:
            raise ValueError("C moet positief zijn.")
        if self.A == self.D:
            raise ValueError("A en D moeten verschillen.")

    def __str__(self) -> str:
        return f"A={self.A}, B={self.B}, C={self.C}, D={self.D}"

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> FourPLParameters:
        if not isinstance(values, Mapping) or set(values) != {"A", "B", "C", "D"}:
            raise ValueError("De 4PL-parameters moeten A, B, C en D zijn.")
        return cls(**values)


@dataclass(frozen=True)
class FivePLParameters:
    """Parameters for D + (A-D) / (1 + (x/C)**B)**E."""

    A: float
    B: float
    C: float
    D: float
    E: float

    def __post_init__(self) -> None:
        base = FourPLParameters(self.A, self.B, self.C, self.D)
        for name in ("A", "B", "C", "D"):
            object.__setattr__(self, name, getattr(base, name))
        object.__setattr__(self, "E", finite_number(self.E, "E"))
        if self.E <= 0:
            raise ValueError("E moet positief zijn.")

    def __str__(self) -> str:
        return f"A={self.A}, B={self.B}, C={self.C}, D={self.D}, E={self.E}"

    @classmethod
    def from_mapping(cls, values: Mapping[str, object]) -> FivePLParameters:
        if not isinstance(values, Mapping) or set(values) != {"A", "B", "C", "D", "E"}:
            raise ValueError("De 5PL-parameters moeten A, B, C, D en E zijn.")
        return cls(**values)


LogisticParameters = FourPLParameters | FivePLParameters


class FourPLModel:
    minimum_levels = 4

    def fit(self, x: NDArray[np.float64], y: NDArray[np.float64], weights: NDArray[np.float64]) -> ParameterFit[FourPLParameters]:
        return fit_nonlinear(self, x, y, weights)

    def for_curve(self, params: FourPLParameters, concentration_range: tuple[float, float] | None) -> FourPLModel:
        return self

    def calculate_extinctions(self, concentrations: ArrayLike, params: FourPLParameters) -> NumericResult:
        from scipy.special import expit

        x = np.asarray(concentrations, dtype=float)
        if np.any(x < 0):
            raise ValueError("Concentraties mogen niet negatief zijn.")
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            z = params.B * (np.log(x) - np.log(params.C))
            # Computing each tail separately avoids cancellation near an asymptote.
            values = params.A * expit(-z) + params.D * expit(z)
        return _numeric_result(np.asarray(values))

    def calculate_concentrations(self, extinctions: ArrayLike, params: FourPLParameters) -> NumericResult:
        y = np.asarray(extinctions, dtype=float)
        lower, upper = min(params.A, params.D), max(params.A, params.D)
        tolerance = 4 * np.finfo(float).eps * max(abs(params.A), abs(params.D))
        with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
            # A fitted asymptote may miss an exact endpoint by round-off.
            # Clamp only values just outside, preserving tiny interior tails.
            y = np.where((y < lower) & np.isclose(y, lower, rtol=0, atol=tolerance), lower, y)
            y = np.where((y > upper) & np.isclose(y, upper, rtol=0, atol=tolerance), upper, y)
            valid = (y >= lower) & (y <= upper)
            log_ratio = np.log(np.abs(y - params.A)) - np.log(np.abs(params.D - y))
            values = np.exp(np.log(params.C) + log_ratio / params.B)
        return _numeric_result(np.where(valid, values, np.nan))

    def _initial_guesses(self, concentrations: ArrayLike, extinctions: ArrayLike) -> list[FourPLParameters]:
        """Rank Hill, Mead-style, and fallback guesses in response space.

        See sections 3.2 and 5.2 of https://journal.r-project.org/articles/RJ-2019-003/.
        Mead's two-column regression is evaluated over a bounded Cauchy slope
        grid. The best candidate per asymptote padding is retained, limiting
        the number of subsequent nonlinear multistart fits.
        """
        x = np.asarray(concentrations, dtype=float)
        y = np.asarray(extinctions, dtype=float)
        if x.ndim != 1 or y.ndim != 1 or x.shape != y.shape or not x.size:
            raise ValueError("Standaarden moeten overeenkomende, niet-lege eendimensionale reeksen zijn.")
        if not np.isfinite(x).all() or not np.isfinite(y).all():
            raise ValueError("Standaardconcentraties en responswaarden moeten eindig zijn.")
        if np.any(x < 0):
            raise ValueError("Standaardconcentraties mogen niet negatief zijn.")
        response_scale = float(np.max(np.abs(y)))
        if response_scale == 0 or np.ptp(y / response_scale) == 0:
            raise ValueError("Curvepassing vereist verschillende responswaarden.")
        normalized_y = y / response_scale

        levels = np.unique(x)
        responses = np.array([normalized_y[x == level].mean() for level in levels])
        positive = levels > 0
        if not np.any(positive):
            raise ValueError("Curvepassing vereist positieve standaardconcentraties.")
        log_x = np.log(levels[positive])
        positive_responses = responses[positive]
        centered_log_x = log_x - log_x.mean()
        log_x_variance = float(centered_log_x @ centered_log_x)
        lower_bounds, upper_bounds = self.bounds(x, normalized_y)
        log_c_min = max(lower_bounds[2], np.log(np.nextafter(0.0, 1.0)))
        log_c_max = min(upper_bounds[2], np.log(np.finfo(float).max))

        span = float(np.ptp(responses))
        response_min, response_max = float(responses.min()), float(responses.max())
        if span == 0:
            # Replicates can vary even when every concentration has the same mean.
            response_min, response_max = float(normalized_y.min()), float(normalized_y.max())
            span = response_max - response_min

        def make_guess(a: float, b: float, log_c: float, d: float) -> FourPLParameters | None:
            if not np.isfinite([a, b, log_c, d]).all() or a <= d or b == 0:
                return None
            magnitude = float(np.exp(np.clip(np.log(abs(b)), lower_bounds[1], upper_bounds[1])))
            c = float(np.exp(np.clip(log_c, log_c_min, log_c_max)))
            return FourPLParameters(a, float(np.copysign(magnitude, b)), c, d)

        def error(params: FourPLParameters) -> float:
            residuals = np.asarray(self.calculate_extinctions(x, params)) - normalized_y
            with np.errstate(over="ignore", invalid="ignore"):
                return float(residuals @ residuals)

        # Safe fallback uses concentration order and a trend across all positive
        # levels; response sorting can invent crossings in noisy data.
        trend = float(centered_log_x @ (positive_responses - positive_responses.mean()))
        if trend == 0:
            trend = float(responses[-1] - responses[0])
        lower = response_min - 0.1 * span
        upper = response_max + 0.1 * span
        midpoint = (lower + upper) / 2
        nearest = int(np.argmin(np.abs(positive_responses - midpoint)))
        fallback = make_guess(upper, -1.0 if trend >= 0 else 1.0, float(log_x[nearest]), lower)
        guesses = [fallback] if fallback is not None else []

        # Cauchy quantiles cover both signs and a broad range of slopes without
        # performing nonlinear optimization during initialization.
        slopes = np.tan(np.pi * (np.arange(1, 25) / 25 - 0.5))
        magnitudes = np.exp([lower_bounds[1], upper_bounds[1]])
        slopes = np.unique(np.concatenate((slopes, magnitudes, -magnitudes)))
        slopes = slopes[(np.abs(slopes) >= magnitudes[0]) & (np.abs(slopes) <= magnitudes[1])]

        for padding in (0.001, 0.1, 0.5):
            a = response_max + padding * span
            d = response_min - padding * span
            # In the reported A > D convention, log((y-D)/(A-y)) has slope -B.
            with np.errstate(divide="ignore", invalid="ignore"):
                logit_y = np.log(positive_responses - d) - np.log(a - positive_responses)
            fraction = (positive_responses - d) / (a - d)
            best_hill = None
            best_error = np.inf
            # Logits amplify errors near plateaus. Also try just the transition
            # region, retaining the better guess by error on all observations.
            for mask in (np.isfinite(logit_y), np.isfinite(logit_y) & (fraction >= 0.05) & (fraction <= 0.95)):
                if np.count_nonzero(mask) < 2:
                    continue
                hill_x, hill_y = log_x[mask], logit_y[mask]
                centered_hill_x = hill_x - hill_x.mean()
                variance = float(centered_hill_x @ centered_hill_x)
                if variance == 0:
                    continue
                slope = float(centered_hill_x @ (hill_y - hill_y.mean()) / variance)
                if not np.isfinite(slope) or abs(slope) <= np.finfo(float).eps:
                    continue
                hill = make_guess(a, -slope, float(hill_x.mean() - hill_y.mean() / slope), d)
                if hill is not None and (cost := error(hill)) < best_error:
                    best_hill, best_error = hill, cost
            if best_hill is not None:
                guesses.append(best_hill)

            shifted_y = positive_responses - d
            if log_x_variance == 0 or np.any(shifted_y <= 0):
                continue
            best_mead = None
            best_error = np.inf
            for b in slopes:
                # Mead: (y-D) * (u + v * exp(B*(log(x)-reference))) ~= 1.
                # Normalize the exponential column to prevent overflow; undo
                # that shift analytically when recovering log(C).
                exponent = b * centered_log_x
                shift = float(exponent.max())
                with np.errstate(under="ignore"):
                    design = np.column_stack((shifted_y, shifted_y * np.exp(exponent - shift)))
                column_scales = np.max(np.abs(design), axis=0)
                if np.any(column_scales == 0):
                    continue
                try:
                    coefficients, _, rank, _ = np.linalg.lstsq(design / column_scales, np.ones_like(shifted_y), rcond=None)
                except np.linalg.LinAlgError:
                    continue
                with np.errstate(over="ignore", invalid="ignore"):
                    coefficients /= column_scales
                if rank < 2 or not np.isfinite(coefficients).all() or np.any(coefficients <= 0):
                    continue
                u, v = coefficients
                with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
                    mead = make_guess(float(d + 1 / u), float(b), float(log_x.mean() + (np.log(u) - np.log(v) + shift) / b), d)
                if mead is not None and (cost := error(mead)) < best_error:
                    best_mead, best_error = mead, cost
            if best_mead is not None:
                guesses.append(best_mead)

        ranked = sorted((guess for guess in dict.fromkeys(guesses) if np.isfinite(error(guess))), key=error)
        restored = []
        for guess in ranked:
            try:
                restored.append(self.scale_parameters(guess, response_scale))
            except (ValueError, OverflowError):
                continue
        if not restored:
            raise CalibrationDataError("Er konden geen eindige beginwaarden voor de 4PL-parameters worden bepaald.")
        return restored

    def pack(self, params: FourPLParameters) -> NDArray[np.float64]:
        # A signed slope has an equivalent positive slope with swapped asymptotes.
        a, b, d = (params.A, params.B, params.D) if params.B > 0 else (params.D, -params.B, params.A)
        return np.array([a, np.log(b), np.log(params.C), d])

    def unpack(self, values: ArrayLike) -> FourPLParameters:
        a, log_b, log_c, d = np.asarray(values, dtype=float)
        b = np.exp(log_b)
        # Fitted parameters use a consistent asymptote order without changing the curve.
        if a < d:
            a, b, d = d, -b, a
        return FourPLParameters(a, b, np.exp(log_c), d)

    def scale_parameters(self, params: FourPLParameters, response_scale: float) -> FourPLParameters:
        return FourPLParameters(params.A * response_scale, params.B, params.C, params.D * response_scale)

    def bounds(self, concentrations: ArrayLike, extinctions: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        x = np.asarray(concentrations, dtype=float)
        positive = x[x > 0]
        if positive.size == 0:
            raise ValueError("Curvepassing vereist positieve standaardconcentraties.")
        log_c_min = np.log(np.min(positive)) - np.log(100.0)
        log_c_max = np.log(np.max(positive)) + np.log(100.0)
        return (
            np.array([-np.inf, np.log(0.05), log_c_min, -np.inf]),
            np.array([np.inf, np.log(20.0), log_c_max, np.inf]),
        )

    def bounds_for_start(
        self, concentrations: ArrayLike, extinctions: ArrayLike, start: FourPLParameters,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        return self.bounds(concentrations, extinctions)

    def starting_points(self, concentrations: ArrayLike, extinctions: ArrayLike) -> list[FourPLParameters]:
        """Combine ranked initial estimates and a small log-dose grid.

        Hill/logit initialization follows the approach discussed in
        https://journal.r-project.org/articles/RJ-2019-003/ and used by drc.
        The padding and grid below are deterministic coverage heuristics.
        Replicates are averaged only to construct guesses, not for fitting.
        """
        x = np.asarray(concentrations, dtype=float)
        y = np.asarray(extinctions, dtype=float)
        starts = self._initial_guesses(x, y)
        initial = starts[0]
        log_x = np.log(np.unique(x[x > 0]))

        # Vary slope and midpoint independently to cover shallow/steep curves
        # and basins away from the observed response midpoint.
        log_midpoints = np.unique(np.append(np.linspace(log_x.min(), log_x.max(), 5)[1:-1], np.log(initial.C)))
        for b in (0.5, 1.0, 4.0):
            for log_c in log_midpoints:
                starts.append(FourPLParameters(initial.A, float(np.copysign(b, initial.B)), float(np.exp(log_c)), initial.D))
        return list(dict.fromkeys(starts))


class FivePLModel:
    """Asymmetric logistic curve with stable tails and two signed-slope branches.

    Unlike 4PL, swapping A/D and negating B changes a 5PL curve when E != 1.
    Parameters therefore retain their asymptote order throughout optimization.
    """

    minimum_levels = 5

    def fit(self, x: NDArray[np.float64], y: NDArray[np.float64], weights: NDArray[np.float64]) -> ParameterFit[FivePLParameters]:
        return fit_nonlinear(self, x, y, weights)

    def for_curve(self, params: FivePLParameters, concentration_range: tuple[float, float] | None) -> FivePLModel:
        return self

    def calculate_extinctions(self, concentrations: ArrayLike, params: FivePLParameters) -> NumericResult:
        x = np.asarray(concentrations, dtype=float)
        if np.any(x < 0):
            raise ValueError("Concentraties mogen niet negatief zijn.")
        with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
            z = params.B * (np.log(x) - np.log(params.C))
            log_fraction = -params.E * np.logaddexp(0.0, z)
            # Compute both tails separately, preserving small responses near A.
            values = params.A * np.exp(log_fraction) + params.D * -np.expm1(log_fraction)
        return _numeric_result(np.asarray(values))

    def calculate_concentrations(self, extinctions: ArrayLike, params: FivePLParameters) -> NumericResult:
        y = np.asarray(extinctions, dtype=float)
        lower, upper = min(params.A, params.D), max(params.A, params.D)
        tolerance = 4 * np.finfo(float).eps * max(abs(params.A), abs(params.D))
        with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
            y = np.where((y < lower) & np.isclose(y, lower, rtol=0, atol=tolerance), lower, y)
            y = np.where((y > upper) & np.isclose(y, upper, rtol=0, atol=tolerance), upper, y)
            valid = (y >= lower) & (y <= upper)
            span = params.A - params.D
            # log1p retains the small distance from A; log differences retain
            # tiny fractions near D without first dividing them to zero.
            log_fraction = np.where(
                np.abs(y - params.A) <= np.abs(y - params.D),
                np.log1p((y - params.A) / span),
                np.log(np.abs(y - params.D)) - np.log(abs(span)),
            )
            power = -log_fraction / params.E
            # log(expm1(power)), avoiding overflow even for very large power.
            log_ratio = power + np.log(-np.expm1(-power))
            values = np.exp(np.log(params.C) + log_ratio / params.B)
        return _numeric_result(np.where(valid, values, np.nan))

    def pack(self, params: FivePLParameters) -> NDArray[np.float64]:
        return np.array([params.A, params.B, np.log(params.C), params.D, np.log(params.E)])

    def unpack(self, values: ArrayLike) -> FivePLParameters:
        a, b, log_c, d, log_e = np.asarray(values, dtype=float)
        return FivePLParameters(a, b, np.exp(log_c), d, np.exp(log_e))

    def scale_parameters(self, params: FivePLParameters, response_scale: float) -> FivePLParameters:
        return FivePLParameters(params.A * response_scale, params.B, params.C, params.D * response_scale, params.E)

    def bounds(self, concentrations: ArrayLike, extinctions: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        lower, upper = FourPLModel().bounds(concentrations, extinctions)
        lower[1], upper[1] = -20.0, 20.0
        return np.append(lower, np.log(0.05)), np.append(upper, np.log(20.0))

    def bounds_for_start(
        self, concentrations: ArrayLike, extinctions: ArrayLike, start: FivePLParameters
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Keep each attempt on its nonzero slope branch; try both via starts."""
        lower, upper = self.bounds(concentrations, extinctions)
        if start.B > 0:
            lower[1] = 0.05
        else:
            upper[1] = -0.05
        return lower, upper

    def starting_points(self, concentrations: ArrayLike, extinctions: ArrayLike) -> list[FivePLParameters]:
        """Extend 4PL guesses with both orientations and midpoint-matched E seeds.

        At E=1 the reflected orientation is identical. Once E is fitted, the
        two branches describe different asymmetric shapes and must both be tried.
        Replicates are retained by the shared fitter.
        """
        x, y = np.asarray(concentrations, dtype=float), np.asarray(extinctions, dtype=float)
        bases = FourPLModel().starting_points(x, y)
        lower, upper = self.bounds(x, y)
        log_c_min = max(lower[2], np.log(np.nextafter(0.0, 1.0)))
        log_c_max = min(upper[2], np.log(np.finfo(float).max))
        starts = []
        for index, base in enumerate(bases):
            for a, b, d in ((base.A, base.B, base.D), (base.D, -base.B, base.A)):
                for e in ((1.0, 0.3, 3.0) if index < 3 else (1.0,)):
                    power = np.log(2.0) / e
                    log_ratio = power + np.log(-np.expm1(-power))
                    log_c = np.clip(np.log(base.C) - log_ratio / b, log_c_min, log_c_max)
                    starts.append(FivePLParameters(a, b, np.exp(log_c), d, e))
        response_scale = float(np.max(np.abs(y)))

        def error(params: FivePLParameters) -> float:
            residuals = np.asarray(self.calculate_extinctions(x, params)) / response_scale - y / response_scale
            return float(residuals @ residuals)

        return sorted(dict.fromkeys(starts), key=error)
