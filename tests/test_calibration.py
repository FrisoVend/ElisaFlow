import unittest
from dataclasses import astuple
from itertools import product
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
from fixtures import assay_metadata, samples_from_records
from scipy.optimize import least_squares

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.calibration.curves import create_calibration_curve
from elisaflow.calibration.fitting import fit_curve, fit_nonlinear
from elisaflow.calibration.logistic import FourPLModel, FourPLParameters
from elisaflow.configuration.assay_settings import AssaySettings, CalibrationSettings
from elisaflow.models import AssayData, SampleType


class CalibrationModelTests(unittest.TestCase):
    def test_models_round_trip_ascending_and_descending_curves(self) -> None:
        concentrations = np.array([0.0, 0.2, 1.0, 5.0, 20.0])
        cases = [
            (FourPLModel(), FourPLParameters(A=2.0, B=-1.4, C=2.5, D=0.1)),
            (FourPLModel(), FourPLParameters(A=2.0, B=1.4, C=2.5, D=0.1)),
            (FourPLModel(), FourPLParameters(A=0.1, B=-1.4, C=2.5, D=2.0)),
            (FourPLModel(), FourPLParameters(A=0.1, B=1.4, C=2.5, D=2.0)),
        ]

        for model, params in cases:
            with self.subTest(model=type(model).__name__, B=params.B):
                responses = model.calculate_extinctions(concentrations, params)
                calculated = model.calculate_concentrations(responses, params)
                np.testing.assert_allclose(calculated, concentrations)

    def test_negative_concentrations_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, 'niet.?negatief'):
            FourPLModel().calculate_extinctions(-1.0, FourPLParameters(A=2.0, B=1.0, C=2.0, D=0.1))

    def test_extreme_tails_round_trip(self) -> None:
        for slope, concentration in ((1.0, 1e300), (-1.0, 1e-300)):
            with self.subTest(B=slope):
                params = FourPLParameters(A=1.0, B=slope, C=1.0, D=0.0)
                response = FourPLModel().calculate_extinctions(concentration, params)
                np.testing.assert_allclose(response, 1e-300, rtol=1e-12, atol=0)
                np.testing.assert_allclose(
                    FourPLModel().calculate_concentrations(response, params),
                    concentration,
                    rtol=1e-12,
                    atol=0,
                )

    def test_small_response_near_a_round_trips_without_cancellation(self) -> None:
        params = FourPLParameters(A=0.0, B=1.0, C=1.0, D=1.0)
        # y = x / (1 + x) is approximately x at tiny x.
        response = FourPLModel().calculate_extinctions(1e-20, params)
        np.testing.assert_allclose(response, 1e-20, rtol=1e-12, atol=0)
        np.testing.assert_allclose(FourPLModel().calculate_concentrations(response, params), 1e-20, rtol=1e-12, atol=0)

    def test_asymptotes_and_invalid_responses_for_all_four_parameter_forms(self) -> None:
        for (a, d), b in product(((2.0, 0.1), (0.1, 2.0)), (-1.4, 1.4)):
            with self.subTest(A=a, B=b, D=d):
                params = FourPLParameters(a, b, 2.5, d)
                model = FourPLModel()
                endpoints = [a, d] if b > 0 else [d, a]
                np.testing.assert_array_equal(model.calculate_extinctions([0, np.inf], params), endpoints)
                np.testing.assert_array_equal(model.calculate_concentrations(endpoints, params), [0, np.inf])
                self.assertTrue(np.isnan(model.calculate_concentrations([0, 3, np.nan], params)).all())

    def test_roundoff_just_outside_asymptotes_maps_to_endpoints(self) -> None:
        for (a, d), b in product(((2.0, 0.1), (0.1, 2.0), (0.0, 2.0)), (-1.4, 1.4)):
            with self.subTest(A=a, B=b, D=d):
                params = FourPLParameters(a, b, 2.5, d)
                endpoints = [a, d] if b > 0 else [d, a]
                outside = [np.nextafter(value, -np.inf if value == min(a, d) else np.inf) for value in endpoints]
                model = FourPLModel()
                np.testing.assert_array_equal(model.calculate_concentrations(outside, params), [0, np.inf])
                self.assertTrue(np.isnan(model.calculate_concentrations([min(a, d) - 1e-10, max(a, d) + 1e-10], params)).all())

    def test_starting_points_are_distinct_finite_and_within_packed_bounds(self) -> None:
        model = FourPLModel()
        x = np.repeat(np.concatenate(([0.0], np.logspace(-5, 8, 12))), 2)[::-1]
        for b in (-3.4, 3.4):
            with self.subTest(B=b):
                y = model.calculate_extinctions(x, FourPLParameters(3.2, b, 4_000.0, 0.04))
                starts = model.starting_points(x, y)
                self.assertGreater(len(starts), 1)
                self.assertEqual(len(starts), len(set(starts)))
                self.assertGreater(len({abs(start.B) for start in starts}), 1)
                self.assertGreater(len({start.C for start in starts}), 1)
                self.assertGreater(len({(start.A, start.D) for start in starts}), 1)
                lower, upper = model.bounds(x, y)
                for start in starts:
                    packed = model.pack(start)
                    self.assertTrue(np.isfinite(packed).all())
                    self.assertTrue(np.all(packed >= lower) and np.all(packed <= upper))
                    self.assertTrue(np.isfinite(model.calculate_extinctions(x, start)).all())

    def test_starting_points_do_not_depend_on_order_or_identical_replicates(self) -> None:
        model = FourPLModel()
        x = np.array([0.0, 0.25, 0.75, 2.0, 6.0, 20.0])
        y = model.calculate_extinctions(x, FourPLParameters(2.2, -1.3, 2.0, 0.2))
        expected = model.starting_points(x, y)
        self.assertEqual(model.starting_points(x[::-1], y[::-1]), expected)
        self.assertEqual(model.starting_points(np.repeat(x, 2), np.repeat(y, 2)), expected)

    def test_first_start_estimate_signed_slope_and_midpoint(self) -> None:
        model = FourPLModel()
        x = np.logspace(-3, 3, 12)
        for b in (-4.0, -1.3, -0.3, 0.3, 1.3, 4.0):
            with self.subTest(B=b):
                expected = FourPLParameters(2.5, b, 1.0, 0.05)
                y = model.calculate_extinctions(x, expected)
                initial = model.starting_points(x, y)[0]
                self.assertGreater(initial.A, initial.D)
                np.testing.assert_allclose(initial.B, b, rtol=0.2)
                np.testing.assert_allclose(initial.C, expected.C, rtol=0.1)
                rmse = np.sqrt(np.mean((model.calculate_extinctions(x, initial) - y) ** 2))
                self.assertLess(rmse, 0.02 * (expected.A - expected.D))

    def test_mead_initialization_improves_partial_curves_over_hill_only(self) -> None:
        model = FourPLModel()
        for b, limits in product((-1.3, 1.3), ((-2, -0.5), (0.5, 2))):
            with self.subTest(B=b, limits=limits):
                x = np.logspace(*limits, 12)
                y = model.calculate_extinctions(x, FourPLParameters(2.5, b, 1.0, 0.05))
                initial = model.starting_points(x, y)[0]
                # Disable the auxiliary Mead solves to compare with the guarded
                # Hill/fallback initializer on the same incomplete curve.
                with patch("elisaflow.calibration.logistic.np.linalg.lstsq", side_effect=np.linalg.LinAlgError("Singular")):
                    hill_only = model.starting_points(x, y)[0]
                initial_error = np.sum((model.calculate_extinctions(x, initial) - y) ** 2)
                hill_error = np.sum((model.calculate_extinctions(x, hill_only) - y) ** 2)
                self.assertLess(initial_error, 0.3 * hill_error)

    def test_first_start_survive_failed_auxiliary_regressions(self) -> None:
        model = FourPLModel()
        x = np.logspace(-2, 2, 12)
        y = model.calculate_extinctions(x, FourPLParameters(2.5, -1.3, 1.0, 0.05))
        with patch("elisaflow.calibration.logistic.np.linalg.lstsq", side_effect=np.linalg.LinAlgError("Singular")):
            initial = model.starting_points(x, y)[0]
        self.assertTrue(np.isfinite(astuple(initial)).all())
        self.assertLess(initial.B, 0)
        self.assertTrue(np.isfinite(model.calculate_extinctions(x, initial)).all())

    def test_first_start_handle_flat_replicate_means(self) -> None:
        model = FourPLModel()
        x = np.repeat([0.0, 1.0, 2.0, 3.0], 2)
        y = np.tile([0.5, 1.5], 4)
        initial = model.starting_points(x, y)[0]
        self.assertTrue(np.isfinite(astuple(initial)).all())
        self.assertGreater(initial.A, initial.D)
        self.assertTrue(np.isfinite(model.calculate_extinctions(x, initial)).all())

    def test_first_start_handle_one_positive_level_without_hill_regression(self) -> None:
        initial = FourPLModel().starting_points([0, 1], [0.2, 1.0])[0]
        self.assertTrue(np.isfinite(astuple(initial)).all())
        self.assertLess(initial.B, 0)

    def test_first_start_preserve_response_and_concentration_units(self) -> None:
        model = FourPLModel()
        x = np.concatenate(([0.0], np.logspace(-2, 2, 12)))
        y = np.asarray(model.calculate_extinctions(x, FourPLParameters(2.5, -1.3, 1.0, 0.05)))
        reference = model.starting_points(x, y)[0]
        for factor in (1e-100, 1e100):
            with self.subTest(response_scale=factor):
                initial = model.starting_points(x, y * factor)[0]
                np.testing.assert_allclose(astuple(initial), np.asarray(astuple(reference)) * [factor, 1, 1, factor], rtol=1e-10)
            with self.subTest(concentration_scale=factor):
                initial = model.starting_points(x * factor, y)[0]
                np.testing.assert_allclose(astuple(initial), np.asarray(astuple(reference)) * [1, 1, factor, 1], rtol=1e-10)

    def test_first_start_reject_invalid_input_directly(self) -> None:
        cases = (
            ([], [], "niet-lege"),
            ([0, 1], [1], "overeenkomende"),
            ([[0, 1]], [[1, 2]], "eendimensionale"),
            ([0, np.inf], [1, 2], 'eindig'),
            ([0, 1], [1, np.nan], 'eindig'),
            ([-1, 1], [1, 2], 'niet.?negatief'),
            ([0, 0], [1, 2], "positieve standaard"),
            ([0, 1], [1, 1], "verschillende responswaarden"),
        )
        for x, y, message in cases:
            with self.subTest(x=x, y=y), self.assertRaisesRegex(ValueError, message):
                FourPLModel().starting_points(x, y)[0]


class CalibrationSettingsTests(unittest.TestCase):
    def test_configured_4pl_uses_lettered_parameters_and_range(self) -> None:
        settings = CalibrationSettings.from_dict(
            {"response_scale": "blank_corrected", "method": "calibration", "fit_model": "4pl", "parameters": {"A": 2.0, "B": -1.4, "C": 2.5, "D": 0.1}, "concentration_range": [0.2, 20.0]}
        )

        curve = create_calibration_curve(settings)

        self.assertIsInstance(curve.parameters, FourPLParameters)
        self.assertEqual(curve.concentration_range, (0.2, 20.0))
        below_range_response = curve.calculate_extinction(0.1)
        self.assertTrue(np.isnan(curve.calculate_concentration(below_range_response)))

    def test_old_parameter_names_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, 'A, B, C en D'):
            CalibrationSettings.from_dict(
                {"response_scale": "blank_corrected", "method": "calibration", "fit_model": "4pl", "parameters": {"bottom": 0.1, "top": 2.0, "scale": 2.5, "hill_slope": 1.4}}
            )

    def test_range_includes_endpoints_for_both_slope_directions(self) -> None:
        for slope in (-1.4, 1.4):
            with self.subTest(B=slope):
                curve = create_calibration_curve(
                    CalibrationSettings(
                        "4pl",
                        FourPLParameters(A=2.0, B=slope, C=2.5, D=0.1),
                        concentration_range=(0.2, 20.0), response_scale="blank_corrected",
                    )
                )
                x = np.array([0.1, 0.2, 2.0, 20.0, 30.0])
                y = curve.calculate_extinction(x)
                np.testing.assert_allclose(curve.calculate_concentration(y), [np.nan, 0.2, 2.0, 20.0, np.nan])
                self.assertIsInstance(curve.calculate_extinction(0.2), float)
                self.assertIsInstance(curve.calculate_concentration(y[1]), float)
                self.assertAlmostEqual(curve.calculate_concentration(y[1]), 0.2)

    def test_fitting_requires_both_standard_arrays(self) -> None:
        for arguments in (
            {},
            {"concentrations": [1, 2, 3, 4]},
            {"extinction_values": [1, 2, 3, 4]},
        ):
            with self.subTest(arguments=arguments):
                with self.assertRaisesRegex(ValueError, 'vereist standaard'):
                    create_calibration_curve(CalibrationSettings("4pl", response_scale="raw"), **arguments)

    def test_configured_curves_allow_either_slope_direction(self) -> None:
        for slope in (-1.4, 1.4):
            with self.subTest(B=slope):
                settings = CalibrationSettings(
                    "4pl",
                    FourPLParameters(A=2.0, B=slope, C=2.5, D=0.1), response_scale="blank_corrected",
                )
                self.assertEqual(settings.parameters.B, slope)

    def test_weighting_setting_is_validated(self) -> None:
        with self.assertRaisesRegex(ValueError, "weighting_exponent"):
            CalibrationSettings("4pl", weighting_exponent="automatic", response_scale="raw")

    def test_unsupported_model_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, 'Niet-ondersteund kalibratiemodel'):
            CalibrationSettings("6pl", response_scale="raw")


class CurveFittingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.parameters = FourPLParameters(A=2.2, B=-1.3, C=2.0, D=0.2)
        self.concentrations = np.array([0.0, 0.25, 0.75, 2.0, 6.0, 20.0])
        self.responses = np.asarray(FourPLModel().calculate_extinctions(self.concentrations, self.parameters), dtype=float)

    @staticmethod
    def _optimizer_at_start(residuals, initial, **kwargs):
        errors = residuals(initial)
        return SimpleNamespace(x=initial, fun=errors, cost=float(0.5 * (errors @ errors)), success=True, message="Converged.")

    def test_multistart_selects_weighted_cost_and_uses_the_winners_diagnostics(self) -> None:
        model = FourPLModel()
        p = self.parameters
        bottom_error = FourPLParameters(p.A, p.B, p.C, p.D + 0.08)
        top_error = FourPLParameters(p.A + 0.3, p.B, p.C, p.D)
        scale = float(self.responses.max())
        starts = [model.scale_parameters(params, 1 / scale) for params in (bottom_error, top_error)]

        for exponent, expected in ((0.0, bottom_error), (2.0, top_error)):
            with self.subTest(weighting_exponent=exponent):
                with (
                    patch.object(model, "starting_points", return_value=starts) as generate,
                    patch("scipy.optimize.least_squares", side_effect=self._optimizer_at_start) as optimize,
                ):
                    fit = fit_curve(self.concentrations, self.responses, model, exponent)
                self.assertEqual(optimize.call_count, len(starts))
                np.testing.assert_array_equal(generate.call_args.args[1], self.responses / scale)
                np.testing.assert_allclose(astuple(fit.parameters), astuple(expected))
                residuals = self.responses - model.calculate_extinctions(self.concentrations, expected)
                np.testing.assert_allclose(fit.diagnostics.residuals, residuals)
                weights = self.responses ** -exponent
                np.testing.assert_allclose(fit.diagnostics.weights, weights)
                self.assertAlmostEqual(fit.cost, 0.5 * np.sum(weights * residuals**2))

    def test_multistart_continues_after_numerical_errors_and_nonfinite_results(self) -> None:
        model = FourPLModel()
        start = model.starting_points(self.concentrations, self.responses / self.responses.max())[0]
        valid = self._optimizer_at_start(
            lambda values: model.calculate_extinctions(self.concentrations, model.unpack(values)) - self.responses / self.responses.max(),
            model.pack(start),
        )
        invalid_cost = SimpleNamespace(**(vars(valid) | {"cost": np.nan}))
        invalid_parameters = SimpleNamespace(**(vars(valid) | {"x": np.full(4, np.inf)}))
        invalid_residuals = SimpleNamespace(**(vars(valid) | {"fun": np.full_like(self.responses, np.nan)}))
        attempts = [
            ValueError("Bad start"),
            FloatingPointError("Numerical failure"),
            OverflowError("Overflow"),
            np.linalg.LinAlgError("Singular Jacobian"),
            invalid_cost,
            invalid_parameters,
            valid,
            invalid_residuals,
        ]
        with (
            patch.object(model, "starting_points", return_value=[start] * len(attempts)),
            patch("scipy.optimize.least_squares", side_effect=attempts) as optimize,
        ):
            fit = fit_curve(self.concentrations, self.responses, model)
        self.assertTrue(fit.success)
        self.assertTrue(np.isfinite(fit.cost))
        self.assertEqual(optimize.call_count, len(attempts))

    def test_multistart_prefers_convergence_over_a_lower_unfinished_cost(self) -> None:
        model = FourPLModel()
        scale = self.responses.max()
        exact = model.scale_parameters(self.parameters, 1 / scale)
        approximate = model.starting_points(self.concentrations, self.responses / scale)[0]

        def solve(residuals, initial, **kwargs):
            result = self._optimizer_at_start(residuals, initial)
            if np.array_equal(initial, model.pack(exact)):
                result.success = False
                result.message = "Evaluation limit."
            return result

        with (
            patch.object(model, "starting_points", return_value=[approximate, exact]),
            patch("scipy.optimize.least_squares", side_effect=solve),
        ):
            fit = fit_curve(self.concentrations, self.responses, model)
        self.assertTrue(fit.success)
        self.assertEqual(fit.message, "Converged.")
        np.testing.assert_allclose(
            fit.diagnostics.predictions,
            model.calculate_extinctions(self.concentrations, model.scale_parameters(approximate, scale)),
        )

    def test_all_unconverged_attempts_return_the_best_finite_failed_fit(self) -> None:
        model = FourPLModel()
        normalized_y = self.responses / self.responses.max()
        exact = model.scale_parameters(self.parameters, 1 / self.responses.max())
        approximate = model.starting_points(self.concentrations, normalized_y)[0]

        def stop(residuals, initial, **kwargs):
            result = self._optimizer_at_start(residuals, initial)
            result.success = False
            result.message = "Evaluation limit."
            return result

        with (
            patch.object(model, "starting_points", return_value=[approximate, exact]),
            patch("scipy.optimize.least_squares", side_effect=stop),
        ):
            fit = fit_curve(self.concentrations, self.responses, model)
        self.assertFalse(fit.success)
        np.testing.assert_allclose(astuple(fit.parameters), astuple(self.parameters))
        self.assertLess(fit.cost, 1e-28)

    def test_all_unusable_starts_raise_a_clear_convergence_error(self) -> None:
        model = FourPLModel()
        start = model.starting_points(self.concentrations, self.responses / self.responses.max())[0]
        with (
            patch.object(model, "starting_points", return_value=[start, start]),
            patch("scipy.optimize.least_squares", side_effect=ValueError("Bad start")) as optimize,
            self.assertRaisesRegex(RuntimeError, 'convergeerde vanuit geen enkel beginpunt: Bad start'),
        ):
            fit_curve(self.concentrations, self.responses, model)
        self.assertEqual(optimize.call_count, 2)

    def test_empty_starting_points_are_rejected(self) -> None:
        model = FourPLModel()
        with patch.object(model, "starting_points", return_value=[]), self.assertRaisesRegex(ValueError, 'ten minste één beginpunt'):
            fit_curve(self.concentrations, self.responses, model)

    def test_multistart_loop_supports_a_different_model_and_parameter_type(self) -> None:
        class LinearModel:
            minimum_levels = 2

            def fit(self, x, y, weights):
                return fit_nonlinear(self, x, y, weights)

            def starting_points(self, x, y):
                return [0.1, 0.5, 2.0]

            def calculate_extinctions(self, x, params):
                return np.asarray(x) * params

            def calculate_concentrations(self, y, params):
                return np.asarray(y) / params

            def pack(self, params):
                return np.array([params])

            def unpack(self, values):
                return float(values[0])

            def scale_parameters(self, params, response_scale):
                return params * response_scale

            def bounds(self, x, y):
                return np.array([0.0]), np.array([np.inf])

            def bounds_for_start(self, x, y, start):
                return self.bounds(x, y)

        fit = fit_curve([1, 2, 3], [1.7, 3.4, 5.1], LinearModel())
        self.assertTrue(fit.success, fit.message)
        self.assertAlmostEqual(fit.parameters, 1.7)
        np.testing.assert_allclose(fit.diagnostics.back_calculated, [1, 2, 3])

    def test_multistart_recovers_a_low_weighted_cost_for_a_noisy_curve(self) -> None:
        # A steep response with a noisy low plateau can attract the original
        # starting guess to a worse minimum under inverse-response weighting.
        x = np.logspace(-2, 2, 8)
        y = np.array(
            [
                0.07071158722623438,
                0.10177810130104592,
                0.09006972690113677,
                0.012346457865549831,
                0.005,
                0.018075901585293187,
                0.12736821982579524,
                2.2370070230675982,
            ]
        )
        multi = fit_curve(x, y, FourPLModel(), 2.0)
        self.assertTrue(multi.success, multi.message)
        # A single ranked start gets stuck at weighted cost 1.4940584;
        # the production multistart fit reaches approximately 1.2293426.
        self.assertLess(multi.cost, 1.3)

    def test_zero_standard_back_calculation_survives_fitted_endpoint_roundoff(self) -> None:
        x = np.repeat([0, 1, 5, 25, 100], 2)
        y = FourPLModel().calculate_extinctions(x, self.parameters)
        fit = fit_curve(x, y, FourPLModel())
        self.assertTrue(fit.success, fit.message)
        np.testing.assert_allclose(fit.diagnostics.back_calculated, x, rtol=1e-6, atol=1e-9)

    def test_invalid_standards_are_rejected(self) -> None:
        cases = (
            ([], [], "4 verschillende"),
            ([1, 2, 3, 4], [1, 2, 3], "overeenkomende eendimensionale"),
            ([[1, 2, 3, 4]], [[1, 2, 3, 4]], "overeenkomende eendimensionale"),
            ([0, 1, 2, np.nan], [1, 2, 3, 4], 'eindig'),
            ([0, 1, 2, 3], [1, 2, 3, np.inf], 'eindig'),
            ([-1, 0, 1, 2], [1, 2, 3, 4], 'niet.?negatief'),
            ([0, 1, 2, 2], [1, 2, 3, 4], "4 verschillende"),
            ([0, 1, 2, 3], [1, 1, 1, 1], "verschillende responswaarden"),
        )
        for x, y, message in cases:
            with self.subTest(x=x, y=y):
                with self.assertRaisesRegex(ValueError, message):
                    fit_curve(x, y, FourPLModel(), weighting_exponent=0.0)

    def test_legacy_weighting_keyword_and_custom_functions_are_rejected(self) -> None:
        with self.assertRaisesRegex(TypeError, "weighting"):
            fit_curve(self.concentrations, self.responses, FourPLModel(), weighting="1/y")
        with self.assertRaisesRegex(ValueError, "weighting_exponent"):
            fit_curve(self.concentrations, self.responses, FourPLModel(), lambda y: np.ones_like(y))

    def test_failed_optimizer_convergence_is_rejected(self) -> None:
        def stop_early(*args, **kwargs):
            return least_squares(*args, **(kwargs | {"max_nfev": 1}))

        with patch("scipy.optimize.least_squares", side_effect=stop_early):
            with self.assertRaisesRegex(RuntimeError, 'convergeerde'):
                create_calibration_curve(
                    CalibrationSettings("4pl", response_scale="raw"),
                    concentrations=self.concentrations,
                    extinction_values=self.responses,
                )

    def test_configured_range_is_applied_to_fitted_curves(self) -> None:
        curve = create_calibration_curve(
            CalibrationSettings("4pl", concentration_range=(0.75, 6.0), response_scale="raw"),
            concentrations=self.concentrations.astype(str),
            extinction_values=self.responses,
        )
        self.assertEqual(curve.concentration_range, (0.75, 6.0))
        y = curve.calculate_extinction([0.25, 0.75, 2.0, 6.0, 20.0])
        np.testing.assert_allclose(curve.calculate_concentration(y), [np.nan, 0.75, 2, 6, np.nan])

    def test_fit_recovers_equivalent_parameters_for_both_asymptote_orders(self) -> None:
        for (a, d), b, include_zero in product(((2.6, 0.15), (0.15, 2.6)), (-1.4, 1.4), (False, True)):
            with self.subTest(A=a, B=b, D=d, include_zero=include_zero):
                x = np.logspace(-2, 2, 15)
                if include_zero:
                    x = np.append(x, 0.0)
                x = x[::-1]  # Fitting must not depend on sorted input or a zero standard.
                expected = FourPLParameters(max(a, d), b if a > d else -b, 2.0, min(a, d))
                # Generate independently from the defining equation.
                with np.errstate(divide="ignore"):
                    y = d + (a - d) / (1 + (x / 2.0) ** b)
                fit = fit_curve(x, y, FourPLModel(), weighting_exponent=0.0)

                self.assertTrue(fit.success, fit.message)
                np.testing.assert_allclose(astuple(fit.parameters), astuple(expected), rtol=1e-5)
                np.testing.assert_allclose(fit.diagnostics.residuals, 0, atol=1e-8)
                np.testing.assert_allclose(fit.diagnostics.back_calculated[x > 0], x[x > 0], rtol=1e-5)

    def test_noisy_four_parameter_fit_preserves_diagnostics_for_all_forms(self) -> None:
        x = np.repeat(np.logspace(-2, 2, 15), 2)
        model = FourPLModel()
        for (a, d), b in product(((2.6, 0.15), (0.15, 2.6)), (-1.4, 1.4)):
            with self.subTest(A=a, B=b, D=d):
                params = FourPLParameters(a, b, 2.0, d)
                expected = np.asarray(model.calculate_extinctions(x, params))
                y = expected * np.tile([0.995, 1.005], len(x) // 2)
                fit = fit_curve(x, y, FourPLModel(), weighting_exponent=2.0)
                self.assertTrue(fit.success, fit.message)
                self.assertGreater(fit.parameters.A, fit.parameters.D)
                self.assertEqual(np.sign(fit.parameters.B), np.sign(b if a > d else -b))
                predicted = model.calculate_extinctions(x, fit.parameters)
                np.testing.assert_allclose(predicted, expected, rtol=1e-3)
                diagnostics = fit.diagnostics
                np.testing.assert_allclose(diagnostics.residuals, y - predicted)
                np.testing.assert_allclose(diagnostics.weights, 1 / y**2)
                np.testing.assert_allclose(diagnostics.back_calculated, model.calculate_concentrations(y, fit.parameters))
                np.testing.assert_allclose(diagnostics.relative_errors, (diagnostics.back_calculated - x) / x)

    def test_response_units_do_not_change_fit(self) -> None:
        x = np.repeat(np.concatenate(([0.0], np.logspace(-2, 2, 12))), 2)
        y = np.asarray(FourPLModel().calculate_extinctions(x, self.parameters))
        y += np.tile([-0.002, 0.002], len(x) // 2)
        reference = fit_curve(x, y, FourPLModel(), weighting_exponent=0.0)
        for factor in (1e-9, 1e9):
            with self.subTest(factor=factor):
                fit = fit_curve(x, y * factor, FourPLModel(), weighting_exponent=0.0)
                conversion = np.array([factor, 1, 1, factor])
                np.testing.assert_allclose(
                    astuple(fit.parameters),
                    np.asarray(astuple(reference.parameters)) * conversion,
                    rtol=1e-5,
                )

    def test_four_standards_fit_and_zero_has_undefined_relative_error(self) -> None:
        x = np.array([0.0, 0.5, 2.0, 8.0])
        y = FourPLModel().calculate_extinctions(x, self.parameters)
        fit = fit_curve(x, y, FourPLModel(), weighting_exponent=0.0)
        self.assertTrue(fit.success, fit.message)
        self.assertTrue(np.isnan(fit.diagnostics.relative_errors[0]))

    def test_fitter_rejects_unknown_weighting(self) -> None:
        with self.assertRaisesRegex(ValueError, "weighting_exponent"):
            fit_curve(self.concentrations, self.responses, FourPLModel(), weighting_exponent="bad")

    def test_fit_detects_ascending_curve(self) -> None:
        curve = create_calibration_curve(
            CalibrationSettings("4pl", response_scale="raw"),
            concentrations=self.concentrations,
            extinction_values=self.responses,
        )

        self.assertLess(curve.parameters.B, 0)
        self.assertIsNotNone(curve.diagnostics)
        self.assertEqual(curve.concentration_range, (0.25, 20.0))
        np.testing.assert_allclose(
            curve.diagnostics.back_calculated[1:],
            self.concentrations[1:],
            rtol=1e-5,
            atol=1e-8,
        )

    def test_fit_detects_descending_curve(self) -> None:
        params = FourPLParameters(A=2.2, B=1.3, C=2.0, D=0.2)
        responses = FourPLModel().calculate_extinctions(self.concentrations, params)

        curve = create_calibration_curve(
            CalibrationSettings("4pl", response_scale="raw"),
            concentrations=self.concentrations,
            extinction_values=responses,
        )

        self.assertGreater(curve.parameters.B, 0)
        np.testing.assert_allclose(
            curve.calculate_extinction(self.concentrations),
            responses,
            rtol=1e-6,
            atol=1e-8,
        )

    def test_four_parameter_fit_is_stable_across_wide_concentration_scales(self) -> None:
        concentrations = np.concatenate(([0.0], np.logspace(-5, 8, 18)))
        params = FourPLParameters(A=3.2, B=-3.4, C=4_000.0, D=0.04)
        responses = FourPLModel().calculate_extinctions(concentrations, params)

        curve = create_calibration_curve(
            CalibrationSettings("4pl", response_scale="raw"),
            concentrations=concentrations,
            extinction_values=responses,
        )

        np.testing.assert_allclose(
            curve.calculate_extinction(concentrations),
            responses,
            rtol=1e-8,
            atol=1e-10,
        )

    def test_weighted_fit_is_deterministic(self) -> None:
        concentrations = np.repeat(np.concatenate(([0.0], np.logspace(-3, 3, 9))), 2)
        params = FourPLParameters(A=2.4, B=-3.2, C=17.0, D=0.08)
        responses = np.asarray(FourPLModel().calculate_extinctions(concentrations, params), dtype=float)
        responses += np.tile([-0.015, 0.015], len(responses) // 2)

        first = fit_curve(concentrations, responses, FourPLModel(), weighting_exponent=2.0)
        second = fit_curve(concentrations, responses, FourPLModel(), weighting_exponent=2.0)

        self.assertEqual(first.parameters, second.parameters)
        np.testing.assert_array_equal(first.diagnostics.residuals, second.diagnostics.residuals)

    def test_four_distinct_levels_are_sufficient_with_replicates(self) -> None:
        concentrations = np.repeat([0.0, 0.5, 2.0, 8.0], 2)
        responses = FourPLModel().calculate_extinctions(concentrations, self.parameters)

        curve = create_calibration_curve(
            CalibrationSettings("4pl", response_scale="raw"),
            concentrations=concentrations,
            extinction_values=responses,
        )

        self.assertIsNotNone(curve.diagnostics)

    def test_configured_weighting_is_applied(self) -> None:
        for exponent, expected in (
            (0.0, np.ones_like(self.responses)),
            (1.0, 1 / self.responses),
            (1.5, 1 / self.responses**1.5),
            (2.0, 1 / self.responses**2),
        ):
            with self.subTest(weighting_exponent=exponent):
                curve = create_calibration_curve(
                    CalibrationSettings("4pl", weighting_exponent=exponent, response_scale="raw"),
                    concentrations=self.concentrations,
                    extinction_values=self.responses,
                )
                np.testing.assert_allclose(curve.diagnostics.weights, expected)

    def test_inverse_response_weighting_requires_positive_responses(self) -> None:
        responses = self.responses.copy()
        responses[0] = 0

        with self.assertRaisesRegex(ValueError, 'vereist positieve responswaarden'):
            create_calibration_curve(
                CalibrationSettings("4pl", weighting_exponent=1.0, response_scale="raw"),
                concentrations=self.concentrations,
                extinction_values=responses,
            )

    def test_relative_errors_are_diagnostics_without_a_fit_acceptance_threshold(self) -> None:
        responses = self.responses.copy()
        responses[2] += 0.15

        curve = create_calibration_curve(
            CalibrationSettings("4pl", response_scale="raw"), concentrations=self.concentrations, extinction_values=responses
        )
        self.assertTrue(np.any(np.abs(curve.diagnostics.relative_errors[1:]) > 0.01))

    def test_analysis_returns_curve_diagnostics(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": [f"standard-{x}" for x in self.concentrations],
                "sample_name": [f"Standard {x}" for x in self.concentrations],
                "sample_type": [SampleType.STANDARD.value] * len(self.concentrations),
                "sample_concentration": self.concentrations,
                "source_extinction": self.responses,
            }
        )
        analyzer = AssayAnalyzer(
            AssayData(samples_from_records(frame.to_dict("records")), assay_metadata()),
            AssaySettings("test.skax", "lot", calculation=CalibrationSettings("4pl", response_scale="raw")),
        )

        details = analyzer.analyze().calculation_details

        self.assertIsNotNone(details.curve.diagnostics)
        self.assertEqual(details.curve.concentration_range, (0.25, 20.0))


if __name__ == "__main__":
    unittest.main()
