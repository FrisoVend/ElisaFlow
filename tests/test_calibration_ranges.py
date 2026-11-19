"""Fitted reporting ranges follow retained standards for every curve model."""

import unittest

import numpy as np

from elisaflow.calibration.curves import CalibrationCurve, create_calibration_curve
from elisaflow.calibration.errors import CalibrationDataError
from elisaflow.calibration.logistic import FivePLModel, FivePLParameters, FourPLModel, FourPLParameters
from elisaflow.calibration.polynomial import CubicModel, CubicParameters, QuadraticModel, QuadraticParameters
from elisaflow.calibration.power_law import PowerLawModel, PowerLawParameters
from elisaflow.configuration.assay_settings import CalibrationSettings
from elisaflow.models import ResultStatus


class FittedReportingRangeTests(unittest.TestCase):
    cases = (
        ("4pl", FourPLModel(), FourPLParameters(4, -1.2, 4, 0.1)),
        ("5pl", FivePLModel(), FivePLParameters(4, -1.2, 4, 0.1, 1.4)),
        ("quadratic", QuadraticModel(), QuadraticParameters(0.1, 0.3, 0.02)),
        ("cubic", CubicModel(), CubicParameters(0.1, 0.3, 0.02, 0.001)),
        ("power_law", PowerLawModel(), PowerLawParameters(0.3, 1.7, 0.1)),
    )

    def test_wider_and_narrower_settings_bound_results_without_dropping_fitting_observations(self):
        x = np.repeat([0, 1, 2, 4, 8, 16], 2)
        for name, model, parameters in self.cases:
            y = model.calculate_extinctions(x, parameters)
            for configured, expected in (((0.25, 64), (1, 16)), ((2, 8), (2, 8))):
                with self.subTest(model=name, configured=configured):
                    curve = create_calibration_curve(
                        CalibrationSettings(name, concentration_range=configured, response_scale="raw"), concentrations=x, extinction_values=y,
                    )
                    self.assertEqual(curve.concentration_range, expected)
                    self.assertEqual(len(curve.diagnostics.predictions), len(x))
                    np.testing.assert_allclose(curve.diagnostics.predictions, y, rtol=1e-6, atol=1e-8)
                    np.testing.assert_allclose(curve.diagnostics.back_calculated, x, rtol=1e-6, atol=1e-8)
                    lower, upper = expected
                    responses = curve.calculate_extinction([lower / 2, lower, upper, upper * 2])
                    calculated = curve.calculate_results(responses)
                    self.assertEqual([value.status for value in calculated], [
                        ResultStatus.BELOW_RANGE, ResultStatus.IN_RANGE, ResultStatus.IN_RANGE, ResultStatus.ABOVE_RANGE,
                    ])
                    self.assertEqual(calculated[0].limit, lower)
                    self.assertEqual(calculated[-1].limit, upper)
                    np.testing.assert_allclose([calculated[1].concentration, calculated[2].concentration], expected, rtol=1e-8)
                    self.assertTrue(np.isnan(calculated[0].concentration))
                    self.assertTrue(np.isnan(calculated[-1].concentration))

    def test_empty_or_single_point_overlap_is_a_recognized_calibration_failure(self):
        for configured in ((20, 40), (16, 20), (0.1, 1)):
            with self.subTest(configured=configured), self.assertRaisesRegex(CalibrationDataError, 'overlapt niet'):
                create_calibration_curve(
                    CalibrationSettings("quadratic", concentration_range=configured, response_scale="raw"),
                    concentrations=[0, 1, 2, 4, 8, 16], extinction_values=[1, 2, 3, 5, 9, 17],
                )


class CalibrationEndpointTests(unittest.TestCase):
    def test_endpoint_roundoff_preserves_exact_bounds_without_extrapolating(self):
        for direction in (1, -1):
            with self.subTest(direction=direction):
                curve = CalibrationCurve(
                    QuadraticModel(), QuadraticParameters(0.2, direction * 0.0002, -direction * 0.001), (0.01, 0.1),
                )
                low_response, middle_response, high_response = curve.calculate_extinction([0.01, 0.05, 0.1])
                responses = [low_response - direction * 1e-9, low_response, middle_response, high_response, high_response + direction * 1e-9]
                results = curve.calculate_results(responses)
                self.assertEqual([value.status for value in results], [
                    ResultStatus.BELOW_RANGE, ResultStatus.IN_RANGE, ResultStatus.IN_RANGE, ResultStatus.IN_RANGE, ResultStatus.ABOVE_RANGE,
                ])
                self.assertEqual(results[1].concentration, 0.01)
                self.assertEqual(results[3].concentration, 0.1)
                self.assertAlmostEqual(results[2].concentration, 0.05)
                self.assertEqual(curve.calculate_concentration(low_response), 0.01)
                self.assertEqual(results[0].limit, 0.01)
                self.assertEqual(results[-1].limit, 0.1)
                self.assertTrue(np.isnan(results[0].concentration))
                self.assertTrue(np.isnan(results[-1].concentration))


if __name__ == "__main__":
    unittest.main()
