import unittest

import numpy as np

from elisaflow.calibration.fitting import response_weights
from elisaflow.configuration.assay_settings import CalibrationSettings, WavelengthCalibrationSettings


class WeightingExponentTests(unittest.TestCase):
    def test_integer_fractional_and_negative_exponents(self):
        responses = [0.25, 1.0, 4.0]
        for exponent, expected in (
            (0, [1, 1, 1]),
            (0.5, [2, 1, 0.5]),
            (1, [4, 1, 0.25]),
            (1.5, [8, 1, 0.125]),
            (2, [16, 1, 0.0625]),
            (-0.5, [0.5, 1, 2]),
        ):
            with self.subTest(exponent=exponent):
                np.testing.assert_allclose(response_weights(responses, exponent), expected)

    def test_zero_exponent_accepts_zero_and_negative_responses(self):
        np.testing.assert_array_equal(response_weights([-2.0, 0.0, 3.0], 0), [1, 1, 1])
        for responses in ([0, 1], [-1, 1]):
            for exponent in (-0.5, 0.5, 1, 2):
                with self.subTest(responses=responses, exponent=exponent), self.assertRaisesRegex(ValueError, 'positieve responswaarden'):
                    response_weights(responses, exponent)

    def test_invalid_exponents_are_rejected_by_settings_and_weights(self):
        for value in (None, True, False, "1", "1/y", "1/y^2", [], {}, np.nan, np.inf, -np.inf, 10**1000):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "weighting_exponent"):
                    CalibrationSettings("4pl", weighting_exponent=value, response_scale="raw")
                with self.assertRaisesRegex(ValueError, "weighting_exponent"):
                    CalibrationSettings.from_dict({"response_scale": "raw", "method": "calibration", "fit_model": "4pl", "weighting_exponent": value})
                with self.assertRaisesRegex(ValueError, "weighting_exponent"):
                    response_weights([1.0, 2.0], value)

    def test_nonfinite_responses_and_unrepresentable_weights_are_rejected(self):
        for response in (np.nan, np.inf, -np.inf):
            with self.subTest(response=response), self.assertRaisesRegex(ValueError, 'eindige responswaarden'):
                response_weights([response, 1.0], 0)
        for response, exponent in ((1e-200, 2), (1e200, 2), (1e-200, -2), (1e200, -2)):
            with self.subTest(response=response, exponent=exponent), self.assertRaisesRegex(ValueError, 'Gewichten moeten eindig en positief zijn'):
                response_weights([response, 1.0], exponent)

    def test_default_and_numeric_exponents_round_trip_as_numbers(self):
        self.assertEqual(CalibrationSettings("4pl", response_scale="raw").weighting_exponent, 0.0)
        for exponent in (0, 1, 2, 1.375, -0.5):
            with self.subTest(exponent=exponent):
                settings = CalibrationSettings.from_dict({"response_scale": "raw", "method": "calibration", "fit_model": "4pl", "weighting_exponent": exponent})
                self.assertIsInstance(settings.weighting_exponent, float)
                self.assertEqual(settings.to_dict()["weighting_exponent"], exponent)
                self.assertEqual(CalibrationSettings.from_dict(settings.to_dict()), settings)

    def test_old_weighting_field_is_rejected_in_single_and_wavelength_configs(self):
        for value in ("1", "1/y", "1/y^2", 1.5):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, 'Onbekende velden.*weighting'):
                    CalibrationSettings.from_dict({"response_scale": "raw", "method": "calibration", "fit_model": "4pl", "weighting": value})
                with self.assertRaisesRegex(ValueError, 'Onbekende velden.*weighting'):
                    WavelengthCalibrationSettings.from_dict({
                        "method": "calibration",
                        "wavelengths": {
                            "405": {"response_scale": "raw", "fit_model": "quadratic", "weighting_exponent": 1.5},
                            "450": {"response_scale": "raw", "fit_model": "4pl", "weighting": value},
                        },
                    })


if __name__ == "__main__":
    unittest.main()
