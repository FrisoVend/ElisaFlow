import unittest
from dataclasses import asdict, astuple, replace
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
from fixtures import assay_metadata, samples_from_records
from table_helpers import sample_table, well_table

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.calibration.curves import create_calibration_curve
from elisaflow.calibration.errors import CalibrationDataError, CalibrationFitError
from elisaflow.calibration.fitting import fit_curve
from elisaflow.calibration.plot import plot_calibration
from elisaflow.calibration.polynomial import CubicModel, CubicParameters, QuadraticModel, QuadraticParameters
from elisaflow.configuration.assay_settings import AssaySettings, CalibrationSettings
from elisaflow.models import AssayData, SampleType
from elisaflow.sample_selection import select_samples


class PolynomialTests(unittest.TestCase):
    cases = (
        ("quadratic", QuadraticModel, QuadraticParameters(0.1, 0.3, 0.02)),
        ("cubic", CubicModel, CubicParameters(0.1, 0.3, 0.02, 0.003)),
    )

    def test_forward_inverse_scalars_and_arrays_for_both_directions(self) -> None:
        x = np.array([[0, 0.01, 0.5], [1, 5, 10]])
        for name, model_type, parameters in self.cases:
            for direction in (1, -1):
                with self.subTest(model=name, direction=direction):
                    coefficients = np.asarray(astuple(parameters)) * direction
                    model = model_type((0, 10))
                    params = model.unpack(coefficients)
                    y = sum(coefficient * x**index for index, coefficient in enumerate(coefficients))
                    np.testing.assert_allclose(model.calculate_extinctions(x, params), y)
                    np.testing.assert_allclose(model.calculate_concentrations(y, params), x, atol=1e-13)
                    self.assertIsInstance(model.calculate_extinctions(1, params), float)
                    self.assertIsInstance(model.calculate_concentrations(float(y[1, 0]), params), float)
                    self.assertTrue(np.isnan(model.calculate_concentrations([np.nan, np.inf], params)).all())

    def test_exact_weighted_fits_recover_coefficients_and_diagnostics(self) -> None:
        x = np.repeat([0, 0.1, 0.5, 1, 3, 5, 10], 2)[::-1]
        for name, model_type, params in self.cases:
            for exponent in (0.0, 1.0, 1.5, 2.0):
                with self.subTest(model=name, weighting_exponent=exponent):
                    y = model_type().calculate_extinctions(x, params)
                    with patch("scipy.optimize.least_squares", side_effect=AssertionError("Polynomial fitting uses its constrained solver")):
                        fit = fit_curve(x, y, model_type(), exponent)
                    self.assertTrue(fit.success)
                    np.testing.assert_allclose(astuple(fit.parameters), astuple(params), rtol=1e-11, atol=1e-13)
                    np.testing.assert_allclose(fit.diagnostics.back_calculated, x, atol=1e-12)
                    np.testing.assert_allclose(fit.diagnostics.predictions, y, atol=1e-12)
                    np.testing.assert_allclose(fit.diagnostics.weights, y ** -exponent)
                    self.assertTrue(np.isnan(fit.diagnostics.relative_errors[x == 0]).all())
                    self.assertFalse(fit.diagnostics.predictions.flags.writeable)

    def test_noisy_weighted_fit_matches_independent_linear_algebra(self) -> None:
        x = np.repeat(np.linspace(0, 10, 9), 2)
        noise = np.array([0.003, -0.002, 0.001, 0.002, -0.001, 0.004, 0.002, -0.002, 0.001])
        for name, model_type, params in self.cases:
            y = model_type().calculate_extinctions(x, params) + np.repeat(noise, 2) + np.tile([-0.001, 0.001], 9)
            for exponent in (0.0, 1.0, 1.5, 2.0):
                with self.subTest(model=name, weighting_exponent=exponent):
                    weights = y ** -exponent
                    design = np.column_stack([x**index for index in range(model_type.minimum_levels)])
                    expected = np.linalg.lstsq(design * np.sqrt(weights)[:, None], y * np.sqrt(weights), rcond=None)[0]
                    fit = fit_curve(x, y, model_type(), exponent)
                    np.testing.assert_allclose(astuple(fit.parameters), expected, rtol=1e-10, atol=1e-12)
                    np.testing.assert_allclose(fit.diagnostics.residuals, y - design @ expected, atol=1e-12)
                    self.assertAlmostEqual(fit.cost, 0.5 * np.sum(weights * fit.diagnostics.residuals**2))

    def test_concentration_and_response_units_do_not_change_the_fit(self) -> None:
        x = np.array([0, 0.1, 0.5, 1, 3, 5, 10])
        for name, model_type, params in self.cases:
            y = model_type().calculate_extinctions(x, params)
            for x_scale, y_scale in ((1e-80, 1), (1e80, 1), (1, 1e-100), (1, 1e100)):
                with self.subTest(model=name, x_scale=x_scale, y_scale=y_scale):
                    fit = fit_curve(x * x_scale, y * y_scale, model_type(), 2.0)
                    np.testing.assert_allclose(fit.diagnostics.predictions / y_scale, y, rtol=1e-11)
                    np.testing.assert_allclose(fit.diagnostics.back_calculated / x_scale, x, atol=1e-12)
                    converted = np.asarray(astuple(fit.parameters)).copy() / y_scale
                    for index in range(1, len(converted)):
                        for _ in range(index):
                            converted[index] *= x_scale
                    np.testing.assert_allclose(converted, astuple(params), rtol=1e-10, atol=1e-12)

    def test_minimum_distinct_levels_and_degenerate_linear_polynomials(self) -> None:
        for name, model_type, params in self.cases:
            count = model_type.minimum_levels
            with self.subTest(model=name):
                x = np.arange(count, dtype=float)
                fit = fit_curve(x, model_type().calculate_extinctions(x, params), model_type())
                self.assertTrue(fit.success)
                with self.assertRaisesRegex(ValueError, f"{count} verschillende"):
                    fit_curve(np.repeat(x[:-1], 3), np.tile([1, 2, 3], count - 1), model_type())
                linear = model_type().unpack([0.1, 0.3] + [0] * (count - 2))
                bound = model_type((0, 10))
                np.testing.assert_allclose(bound.calculate_concentrations(bound.calculate_extinctions(x, linear), linear), x)

    def test_fixed_settings_serialization_and_matching_parameter_validation(self) -> None:
        for name, model_type, params in self.cases:
            config = {
                "method": "calibration",
                "fit_model": name,
                "parameters": asdict(params),
                "concentration_range": [0, 10],
                "response_scale": "raw",
            }
            settings = CalibrationSettings.from_dict(config)
            self.assertEqual(CalibrationSettings.from_dict(settings.to_dict()), settings)
            curve = create_calibration_curve(settings)
            self.assertIsInstance(curve.model, model_type)
            self.assertEqual(curve.parameters, params)
            with self.assertRaisesRegex(ValueError, 'vereist een concentratiebereik'):
                CalibrationSettings(name, params, response_scale="blank_corrected")
            with self.assertRaisesRegex(ValueError, 'parameters moeten'):
                CalibrationSettings("4pl", params, response_scale="blank_corrected")
            with self.assertRaises(ValueError):
                CalibrationSettings.from_dict(config | {"parameters": {"A": 1}})
            for changes in ({"A": np.nan}, {"B": True}, {"B": 0, "C": 0, **({"D": 0} if name == "cubic" else {})}):
                with self.subTest(model=name, changes=changes), self.assertRaises(ValueError):
                    replace(params, **changes)

    def test_range_selects_a_unique_branch_and_limits_results(self) -> None:
        params = QuadraticParameters(0, -2, 1)  # Minimum at x=1.
        left = create_calibration_curve(CalibrationSettings("quadratic", params, concentration_range=(0, 0.5), response_scale="blank_corrected"))
        right = create_calibration_curve(CalibrationSettings("quadratic", params, concentration_range=(2, 3), response_scale="blank_corrected"))
        self.assertEqual(left.calculate_concentration(0), 0)
        self.assertAlmostEqual(right.calculate_concentration(0), 2)
        self.assertTrue(np.isnan(left.calculate_concentration(1)))
        self.assertTrue(np.isnan(right.calculate_concentration(-0.5)))

    def test_inversion_requires_a_selected_monotone_branch(self) -> None:
        for name, model_type, params in self.cases:
            with self.subTest(model=name), self.assertRaisesRegex(CalibrationDataError, "monotone tak"):
                model_type().calculate_concentrations([0.1, 0.5], params)

    def test_turning_data_is_fitted_monotonically_and_fixed_turning_curves_are_rejected(self) -> None:
        x = np.linspace(0, 4, 9)
        cases = (
            ("quadratic", QuadraticModel, QuadraticParameters(5, -4, 1)),
            ("cubic", CubicModel, CubicParameters(5, -3, 0, 1)),
        )
        for name, model_type, params in cases:
            with self.subTest(model=name):
                fit = fit_curve(x, model_type().calculate_extinctions(x, params), model_type())
                model = model_type((0, 4))
                model.validate_parameters(fit.parameters)
                values = model.calculate_extinctions(np.linspace(0, 4, 1001), fit.parameters)
                differences = np.diff(values)
                self.assertTrue(np.all(differences >= -1e-12) or np.all(differences <= 1e-12))
                self.assertGreater(fit.cost, 0)
                with self.assertRaisesRegex(ValueError, 'monotoon'):
                    create_calibration_curve(CalibrationSettings(name, params, concentration_range=(0, 4), response_scale="blank_corrected"))

    def test_wide_configured_range_does_not_extend_the_fit_or_reporting_domain(self) -> None:
        x = np.linspace(0, 0.5, 6)
        for name in ("quadratic", "cubic"):
            with self.subTest(model=name):
                curve = create_calibration_curve(
                    CalibrationSettings(name, concentration_range=(0, 3), response_scale="raw"), concentrations=x, extinction_values=1 + x - x * x
                )
                self.assertEqual(curve.concentration_range, (0.1, 0.5))
                values = curve.calculate_extinction(np.linspace(0, 0.5, 1001))
                self.assertTrue(np.all(np.diff(values) >= -1e-12))
                np.testing.assert_allclose(curve.calculate_extinction(x), 1 + x - x * x, atol=1e-8)
                self.assertTrue(np.isnan(curve.calculate_concentration(curve.calculate_extinction(0))))

    def test_constrained_quadratic_matches_an_independent_bounded_least_squares_solution(self) -> None:
        from scipy.optimize import lsq_linear

        x = np.array([10.0, 20.0, 75.0, 500.0])
        z = (x - x.min()) / np.ptp(x)
        # These coefficients are the intercept and the slopes at both interval
        # endpoints. Non-negative endpoint slopes characterize every increasing quadratic.
        design = np.column_stack([np.ones_like(z), z - z**2 / 2, z**2 / 2])
        for direction in (1, -1):
            y = np.array([0.0778, 0.1576, 0.6106, 1.7979])
            y = y if direction == 1 else 3 - y
            for exponent in (0.0, 1.5, 2.0):
                with self.subTest(direction=direction, exponent=exponent):
                    sqrt_weights = y ** (-exponent / 2)
                    expected = lsq_linear(
                        design * sqrt_weights[:, None], direction * y * sqrt_weights,
                        bounds=([-np.inf, 0, 0], np.inf), tol=1e-13,
                    )
                    fit = fit_curve(x, y, QuadraticModel(), exponent)
                    np.testing.assert_allclose(fit.diagnostics.predictions, direction * (design @ expected.x), rtol=1e-6, atol=1e-8)
                    self.assertGreater(fit.cost, 0)
                    QuadraticModel((10, 500)).validate_parameters(fit.parameters)

    def test_cubic_constraint_finds_a_turn_between_monotonic_standard_points(self) -> None:
        x = np.array([0.0, 0.2, 0.8, 1.0])
        # Derivative (x - 0.5)**2 - 0.01 is positive at every standard, but negative between them.
        y = 2 + 0.24 * x - 0.5 * x**2 + x**3 / 3
        self.assertTrue(np.all(np.diff(y) > 0))
        for direction in (1, -1):
            with self.subTest(direction=direction):
                fit = fit_curve(x, direction * y, CubicModel())
                CubicModel((0, 1)).validate_parameters(fit.parameters)
                derivative = np.polynomial.polynomial.polyder(direction * np.array(astuple(fit.parameters)))
                vertex = -derivative[1] / (2 * derivative[2])
                self.assertGreater(vertex, 0)
                self.assertLess(vertex, 1)
                self.assertGreaterEqual(np.polynomial.polynomial.polyval(vertex, derivative), -1e-12)
                # Symmetry makes the exact constrained solution an intercept plus (x - 0.5)**3.
                design = np.column_stack([np.ones_like(x), (x - 0.5)**3])
                expected = design @ np.linalg.lstsq(design, direction * y, rcond=None)[0]
                np.testing.assert_allclose(fit.diagnostics.predictions, expected, atol=1e-7)
                self.assertGreater(fit.cost, 0)

    def test_constrained_fits_are_invariant_to_concentration_and_response_units(self) -> None:
        x = np.array([10.0, 20.0, 75.0, 500.0])
        y = np.array([0.0778, 0.1576, 0.6106, 1.7979])
        for model in (QuadraticModel, CubicModel):
            expected = fit_curve(x, y, model(), 2.0)
            for x_scale, y_scale in ((1e-80, 1), (1e80, 1), (1, 1e-100), (1, 1e100)):
                with self.subTest(model=model.__name__, x_scale=x_scale, y_scale=y_scale):
                    fit = fit_curve(x * x_scale, y * y_scale, model(), 2.0)
                    np.testing.assert_allclose(fit.diagnostics.predictions / y_scale, expected.diagnostics.predictions, rtol=1e-7)
                    np.testing.assert_allclose(fit.diagnostics.back_calculated / x_scale, expected.diagnostics.back_calculated, rtol=1e-6)

    def test_noisy_cubic_converges_at_the_monotonicity_boundary(self) -> None:
        x = np.linspace(0, 1, 8)
        y = np.array([
            1.7240179310917423, 1.922787187790152, 2.6817843697212145, 2.2153351761885105,
            1.7709364252285016, 2.3504633336320366, 2.204396696077302, 2.0615011878867247,
        ])
        fit = fit_curve(x, y, CubicModel(), 2.0)
        self.assertTrue(fit.success)
        CubicModel((0, 1)).validate_parameters(fit.parameters)
        self.assertTrue(np.isfinite(fit.diagnostics.predictions).all())

    def test_constrained_solver_failure_is_a_calibration_failure(self) -> None:
        x = np.array([10.0, 20.0, 75.0, 500.0])
        for _, model_type, parameters in self.cases:
            for y in (np.array([0.0778, 0.1576, 0.6106, 1.7979]), model_type().calculate_extinctions(x, parameters)):
                with self.subTest(model=model_type.__name__, responses=y.tolist()):
                    failure = SimpleNamespace(success=False, x=np.zeros(model_type.minimum_levels))
                    with patch("scipy.optimize.minimize", return_value=failure):
                        with self.assertRaises(CalibrationFitError):
                            fit_curve(x, y, model_type(), 2.0)

    def test_rank_deficient_standard_design_is_rejected(self) -> None:
        for _, model_type, _ in self.cases:
            with self.subTest(model=model_type.__name__):
                x = 1 + np.arange(model_type.minimum_levels) * 1e-8
                y = np.arange(model_type.minimum_levels) + 1.0
                with self.assertRaisesRegex(ValueError, "numeriek stabiele"):
                    fit_curve(x, y, model_type())

    def test_stationary_inflection_and_boundary_turns_remain_invertible(self) -> None:
        params = CubicParameters(-1, 3, -3, 1)  # (x-1)**3, monotonic.
        model = CubicModel((0, 2))
        x = np.array([0, 0.5, 1, 1.5, 2])
        np.testing.assert_allclose(model.calculate_concentrations(model.calculate_extinctions(x, params), params), x, atol=1e-12)
        model = QuadraticModel((0, 1))
        params = QuadraticParameters(0, -2, 1)
        np.testing.assert_allclose(model.calculate_concentrations([0, -1], params), [0, 1])

    def test_range_endpoints_and_missing_values(self) -> None:
        for name, model_type, params in self.cases:
            curve = create_calibration_curve(CalibrationSettings(name, params, concentration_range=(0.1, 10), response_scale="blank_corrected"))
            x = [0, 0.1, 1, 10, 11]
            y = curve.calculate_extinction(x)
            np.testing.assert_allclose(curve.calculate_concentration(y), [np.nan, 0.1, 1, 10, np.nan])
            with self.assertRaisesRegex(ValueError, 'niet.?negatief'):
                model_type().calculate_extinctions(-1, params)

    def test_raw_analysis_recovery_dilution_and_plot_use_polynomial_model(self) -> None:
        x = np.repeat([0, 0.1, 0.5, 1, 3, 5, 10], 2)
        for name, model_type, params in self.cases:
            with self.subTest(model=name):
                frame = pd.DataFrame(
                    {
                        "sample_id": [f"standard-{v}" for v in x] + ["patient", "blank"],
                        "sample_name": [f"Standard {v}" for v in x] + ["Patient", "Blank"],
                        "sample_type": ["Standard"] * len(x) + ["Unknown", "Blank"],
                        "sample_concentration": np.r_[x, np.nan, np.nan],
                        "sample_dilution": [1] * len(x) + [3, 1],
                        "source_extinction": np.r_[model_type().calculate_extinctions(np.r_[x, 2], params), 0.1],
                    }
                )
                settings = AssaySettings(
                    "test.skax", "lot", calculation=CalibrationSettings(name, response_scale="raw", weighting_exponent=2.0)
                )
                result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), assay_metadata()), settings).analyze()
                self.assertEqual(result.calculation_details.calibration_settings.fit_model, name)
                np.testing.assert_allclose(well_table(result.samples).loc[well_table(result.samples)["sample_type"] == "Unknown", "result"], 6)
                samples = sample_table(result.samples)
                standards = samples[(samples["sample_type"] == "Standard") & (samples["sample_concentration"] > 0)]
                np.testing.assert_allclose(standards["recovery_percent"], 100, rtol=1e-6)
                np.testing.assert_allclose(result.calculation_details.curve.diagnostics.weights, 1 / frame.iloc[: len(x)]["source_extinction"] ** 2)
                figure = plot_calibration(
                    result.samples, result.calculation_details.curve.model, result.calculation_details.curve.parameters,
                    input_statistics=tuple(sample.curve_extinction for sample in select_samples(result.samples, SampleType.STANDARD)),
                )
                line = figure.axes[0].lines[0]
                np.testing.assert_allclose(line.get_ydata(), model_type().calculate_extinctions(line.get_xdata(), params), rtol=1e-10)

    def test_model_is_not_mutated_when_fitted_on_another_concentration_range(self) -> None:
        model = QuadraticModel()
        x = np.arange(4, dtype=float)
        fit_curve(x, 1 + x + x * x, model)
        self.assertIsNone(model.concentration_range)
        fit_curve(x * 10, 1 + x + x * x, model)
        self.assertIsNone(model.concentration_range)


if __name__ == "__main__":
    unittest.main()
