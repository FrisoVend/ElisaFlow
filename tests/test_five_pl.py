import unittest
from dataclasses import asdict, astuple
from itertools import product

import numpy as np
import pandas as pd
from fixtures import assay_metadata, samples_from_records
from table_helpers import well_table

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.calibration.curves import create_calibration_curve
from elisaflow.calibration.fitting import fit_curve
from elisaflow.calibration.logistic import FivePLModel, FivePLParameters, FourPLModel, FourPLParameters
from elisaflow.calibration.plot import plot_calibration
from elisaflow.configuration.assay_settings import AssaySettings, CalibrationSettings
from elisaflow.models import AssayData, SampleType
from elisaflow.sample_selection import select_samples


class FivePLModelTests(unittest.TestCase):
    def test_forward_inverse_and_packing_for_both_asymptote_orders_and_slope_signs(self) -> None:
        x = np.array([0, 0.2, 1, 5, 20, np.inf])
        model = FivePLModel()
        for (a, d), b, e in product(((2.5, 0.1), (0.1, 2.5)), (-1.3, 1.3), (0.3, 1, 3)):
            with self.subTest(A=a, B=b, E=e):
                params = FivePLParameters(a, b, 2, d, e)
                with np.errstate(divide="ignore", invalid="ignore"):
                    expected = d + (a - d) / (1 + (x / 2) ** b) ** e
                responses = model.calculate_extinctions(x, params)
                np.testing.assert_allclose(responses, expected)
                np.testing.assert_allclose(model.calculate_concentrations(responses, params), x, atol=1e-14)
                np.testing.assert_allclose(astuple(model.unpack(model.pack(params))), astuple(params))
                endpoints = [a, d] if b > 0 else [d, a]
                np.testing.assert_array_equal(model.calculate_concentrations(endpoints, params), [0, np.inf])
                self.assertTrue(np.isnan(model.calculate_concentrations([-1, 4, np.nan, np.inf], params)).all())
                self.assertIsInstance(model.calculate_extinctions(1, params), float)
                self.assertIsInstance(model.calculate_concentrations(1, params), float)

    def test_e_one_matches_four_pl(self) -> None:
        x = np.r_[0, np.logspace(-100, 100, 50), np.inf]
        for (a, d), b in product(((2.5, 0.1), (0.1, 2.5)), (-1.3, 1.3)):
            with self.subTest(A=a, B=b):
                four = FourPLParameters(a, b, 2, d)
                five = FivePLParameters(a, b, 2, d, 1)
                np.testing.assert_allclose(FivePLModel().calculate_extinctions(x, five), FourPLModel().calculate_extinctions(x, four))
                y = np.linspace(min(a, d), max(a, d), 100)
                np.testing.assert_allclose(
                    FivePLModel().calculate_concentrations(y, five), FourPLModel().calculate_concentrations(y, four), rtol=1e-12
                )

    def test_extreme_tails_preserve_representable_responses_and_inverse(self) -> None:
        model = FivePLModel()
        cases = (
            (FivePLParameters(1, 1, 1, 0, 2), 1e150, 1e-300),
            (FivePLParameters(1, -1, 1, 0, 2), 1e-150, 1e-300),
            (FivePLParameters(0, 1, 1, 1, 2), 5e-301, 1e-300),
            (FivePLParameters(1, 1, 1, 0, 0.05), 1e300, 1e-15),
        )
        for params, x, expected in cases:
            with self.subTest(params=params):
                with np.errstate(all="raise"):
                    response = model.calculate_extinctions(x, params)
                    calculated = model.calculate_concentrations(response, params)
                np.testing.assert_allclose(response, expected, rtol=1e-12, atol=0)
                np.testing.assert_allclose(calculated, x, rtol=1e-12, atol=0)

    def test_roundoff_outside_asymptotes_is_clamped_but_interior_tails_are_preserved(self) -> None:
        model = FivePLModel()
        for (a, d), b in product(((2, 0.1), (0.1, 2), (0, 2)), (-1.3, 1.3)):
            params = FivePLParameters(a, b, 2, d, 0.4)
            endpoints = [a, d] if b > 0 else [d, a]
            outside = [np.nextafter(value, -np.inf if value == min(a, d) else np.inf) for value in endpoints]
            np.testing.assert_array_equal(model.calculate_concentrations(outside, params), [0, np.inf])
            self.assertTrue(np.isnan(model.calculate_concentrations([min(a, d) - 1e-10, max(a, d) + 1e-10], params)).all())

    def test_parameter_and_concentration_validation(self) -> None:
        valid = dict(A=2, B=-1.3, C=2, D=0.1, E=0.4)
        for changes in (dict(E=0), dict(E=-1), dict(E=np.inf), dict(E=np.nan), dict(E=True), dict(B=0), dict(C=0), dict(A=0.1)):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                FivePLParameters(**(valid | changes))
        with self.assertRaisesRegex(ValueError, 'niet.?negatief'):
            FivePLModel().calculate_extinctions(-1, FivePLParameters(**valid))

    def test_configured_parameters_and_settings_round_trip(self) -> None:
        params = FivePLParameters(2.5, -1.3, 2, 0.1, 0.4)
        settings = CalibrationSettings.from_dict(
            {"response_scale": "blank_corrected", "method": "calibration", "fit_model": "5pl", "parameters": asdict(params), "concentration_range": [0.2, 20]}
        )
        self.assertEqual(settings.parameters, params)
        self.assertEqual(CalibrationSettings.from_dict(settings.to_dict()), settings)
        self.assertEqual(
            CalibrationSettings.from_dict({"response_scale": "blank_corrected", "method": "calibration", "fit_model": "5pl", "parameters": asdict(params)}).parameters, params
        )
        curve = create_calibration_curve(settings)
        x = np.array([0.1, 0.2, 2, 20, 30])
        y = curve.calculate_extinction(x)
        np.testing.assert_allclose(curve.calculate_concentration(y), [np.nan, 0.2, 2, 20, np.nan])
        with self.assertRaisesRegex(ValueError, '5PL-parameters'):
            CalibrationSettings.from_dict({"response_scale": "blank_corrected", "method": "calibration", "fit_model": "5pl", "parameters": {"A": 2, "B": 1, "C": 1, "D": 0}})
        for model_name, wrong in (("4pl", params), ("5pl", FourPLParameters(2, 1, 1, 0))):
            with self.subTest(model=model_name), self.assertRaisesRegex(ValueError, 'parameters moeten'):
                CalibrationSettings(model_name, wrong, response_scale="blank_corrected")


class FivePLFittingTests(unittest.TestCase):
    def test_exact_asymmetric_fits_recover_both_asymptote_orders_and_slope_signs(self) -> None:
        x = np.r_[0, np.logspace(-3, 3, 15)]
        for (a, d), b, e in product(((2.5, 0.1), (0.1, 2.5)), (-1.3, 1.3), (0.4, 2.5)):
            with self.subTest(A=a, B=b, E=e):
                with np.errstate(divide="ignore"):
                    y = d + (a - d) / (1 + (x / 2) ** b) ** e
                fit = fit_curve(x, y, FivePLModel())
                self.assertTrue(fit.success, fit.message)
                np.testing.assert_allclose(astuple(fit.parameters), [a, b, 2, d, e], rtol=1e-5)
                np.testing.assert_allclose(fit.diagnostics.predictions, y, atol=1e-8)
                np.testing.assert_allclose(fit.diagnostics.back_calculated[1:], x[1:], rtol=1e-5)
                self.assertTrue(np.isnan(fit.diagnostics.relative_errors[0]))

    def test_noisy_replicates_weighting_diagnostics_and_determinism(self) -> None:
        model = FivePLModel()
        x = np.repeat(np.r_[0, np.logspace(-2, 2, 10)], 2)
        expected = model.calculate_extinctions(x, FivePLParameters(2.5, -1.3, 2, 0.1, 0.4))
        y = expected * np.tile([0.99, 1.01], len(x) // 2)
        for exponent in (0.0, 1.0, 1.5, 2.0):
            with self.subTest(weighting_exponent=exponent):
                curve = create_calibration_curve(
                    CalibrationSettings("5pl", weighting_exponent=exponent, response_scale="raw"), concentrations=x, extinction_values=y
                )
                diagnostics = curve.diagnostics
                self.assertEqual(curve.concentration_range, (0.01, 100))
                np.testing.assert_allclose(diagnostics.predictions, expected, rtol=0.001)
                np.testing.assert_allclose(diagnostics.residuals, y - diagnostics.predictions)
                np.testing.assert_allclose(diagnostics.weights, y**-exponent)
                np.testing.assert_allclose(diagnostics.back_calculated, model.calculate_concentrations(y, curve.parameters))
                np.testing.assert_allclose(diagnostics.relative_errors[2:], (diagnostics.back_calculated[2:] - x[2:]) / x[2:])
                self.assertEqual(len(diagnostics.predictions), len(x))
                self.assertFalse(diagnostics.predictions.flags.writeable)
                repeated = create_calibration_curve(
                    CalibrationSettings("5pl", weighting_exponent=exponent, response_scale="raw"), concentrations=x, extinction_values=y
                )
                self.assertEqual(repeated.parameters, curve.parameters)

    def test_units_and_wide_concentration_ranges(self) -> None:
        model = FivePLModel()
        x = np.r_[0, np.logspace(-5, 8, 18)]
        params = FivePLParameters(3.2, -1.4, 4000, 0.04, 0.6)
        y = model.calculate_extinctions(x, params)
        for concentration_scale, response_scale in ((1, 1), (1e-100, 1e-100), (1e100, 1e100)):
            with self.subTest(concentration_scale=concentration_scale, response_scale=response_scale):
                fit = fit_curve(x * concentration_scale, y * response_scale, model)
                self.assertTrue(fit.success, fit.message)
                np.testing.assert_allclose(fit.diagnostics.predictions / response_scale, y, rtol=1e-7)
                np.testing.assert_allclose(
                    astuple(fit.parameters), np.asarray(astuple(params)) * [response_scale, 1, concentration_scale, response_scale, 1], rtol=1e-5
                )

    def test_five_distinct_levels_required_and_sufficient(self) -> None:
        model = FivePLModel()
        params = FivePLParameters(2.5, -1.3, 2, 0.1, 0.4)
        x = np.repeat([0, 0.2, 1, 5], 3)
        with self.assertRaisesRegex(ValueError, '5 verschillende'):
            fit_curve(x, model.calculate_extinctions(x, params), model)
        x = np.array([0, 0.2, 1, 5, 20])
        y = model.calculate_extinctions(x, params)
        curve = create_calibration_curve(CalibrationSettings("5pl", response_scale="raw"), concentrations=x, extinction_values=y)
        np.testing.assert_allclose(curve.diagnostics.back_calculated[1:], x[1:], rtol=1e-5)

    def test_starts_cover_asymmetry_and_both_slope_branches_with_valid_bounds(self) -> None:
        model = FivePLModel()
        x = np.r_[0, np.logspace(-3, 3, 10)]
        y = model.calculate_extinctions(x, FivePLParameters(2.5, -1.3, 2, 0.1, 0.4))
        starts = model.starting_points(x, y)
        self.assertEqual(len(starts), len(set(starts)))
        self.assertEqual({np.sign(start.B) for start in starts}, {-1, 1})
        self.assertEqual({start.E for start in starts}, {0.3, 1, 3})
        for start in starts:
            lower, upper = model.bounds_for_start(x, y, start)
            packed = model.pack(start)
            self.assertTrue(np.isfinite(packed).all())
            self.assertTrue(np.all(lower <= packed) and np.all(packed <= upper))
            self.assertTrue(lower[1] > 0 or upper[1] < 0)

    def test_analysis_and_plot_use_five_pl_parameters(self) -> None:
        x = np.r_[0, np.logspace(-2, 2, 8)]
        model = FivePLModel()
        params = FivePLParameters(2.5, -1.3, 2, 0.1, 0.4)
        concentrations = np.r_[x, 3.0]
        frame = pd.DataFrame(
            {
                "sample_id": [f"standard-{value}" for value in x] + ["patient"],
                "sample_name": [f"Standard {value}" for value in x] + ["Patient"],
                "sample_type": [SampleType.STANDARD.value] * len(x) + [SampleType.UNKNOWN.value],
                "sample_concentration": np.r_[x, np.nan],
                "source_extinction": model.calculate_extinctions(concentrations, params),
            }
        )
        assay = AssayData(samples_from_records(frame.to_dict("records")), assay_metadata())
        for fixed in (None, params):
            with self.subTest(fixed=fixed is not None):
                settings = AssaySettings(
                    "test.skax", "lot", calculation=CalibrationSettings("5pl", fixed, response_scale="raw" if fixed is None else "blank_corrected"),
                )
                result = AssayAnalyzer(assay, settings).analyze()
                details = result.calculation_details
                self.assertIsInstance(details.curve.parameters, FivePLParameters)
                self.assertEqual(details.calibration_settings.fit_model, "5pl")
                self.assertEqual(details.curve.was_fitted, fixed is None)
                patient = well_table(result.samples)[well_table(result.samples)["sample_id"] == "patient"]
                np.testing.assert_allclose(patient["result"], 3, rtol=1e-5)
                figure = plot_calibration(
                    result.samples, model, details.curve.parameters,
                    input_statistics=tuple(sample.curve_extinction for sample in select_samples(result.samples, SampleType.STANDARD)),
                )
                line = figure.axes[0].lines[0]
                np.testing.assert_allclose(line.get_ydata(), model.calculate_extinctions(line.get_xdata(), params), rtol=1e-6)


if __name__ == "__main__":
    unittest.main()
