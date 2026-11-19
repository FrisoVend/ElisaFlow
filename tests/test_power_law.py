"""Offset power-law calibration shares the assay and reporting pipeline."""

import unittest
from dataclasses import asdict, astuple, replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
from fixtures import assay_metadata, samples_from_records
from scipy.optimize import least_squares
from table_helpers import sample_table, well_table
from test_multi_wavelength import dual_data, dual_settings, saturate

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.calibration.curves import create_calibration_curve
from elisaflow.calibration.errors import CalibrationDataError
from elisaflow.calibration.fitting import fit_curve
from elisaflow.calibration.power_law import PowerLawModel, PowerLawParameters
from elisaflow.configuration.assay_settings import AssaySettings, CalibrationSettings
from elisaflow.models import AssayData, ResultStatus, SampleType
from elisaflow.reporting.pdf_report import PdfReportWriter


class PowerLawTests(unittest.TestCase):
    cases = (
        PowerLawParameters(0.3, 1.7, 0.1),
        PowerLawParameters(-0.3, 0.7, 4),
        PowerLawParameters(0.3, -0.7, 0.1),
        PowerLawParameters(-0.3, -1.2, 4),
    )

    def test_forward_and_inverse_support_signed_parameters_and_missing_values(self):
        model = PowerLawModel()
        x = np.array([[0.1, 0.5, 1], [2, 5, 10]])
        for parameters in self.cases:
            with self.subTest(parameters=parameters):
                response = parameters.C + parameters.A * x ** parameters.B
                np.testing.assert_allclose(model.calculate_extinctions(x, parameters), response)
                np.testing.assert_allclose(model.calculate_concentrations(response, parameters), x, rtol=1e-11)
                self.assertIsInstance(model.calculate_extinctions(1, parameters), float)
                self.assertIsInstance(model.calculate_concentrations(response[0, 0], parameters), float)
                self.assertTrue(np.isnan(model.calculate_concentrations(np.nan, parameters)))
                self.assertTrue(np.isnan(model.calculate_concentrations(parameters.C - parameters.A, parameters)))
                with self.assertRaisesRegex(ValueError, "negatief"):
                    model.calculate_extinctions(-1, parameters)

    def test_zero_concentration_and_offset_have_explicit_domains(self):
        model = PowerLawModel()
        increasing = self.cases[0]
        self.assertEqual(model.calculate_extinctions(0, increasing), increasing.C)
        self.assertEqual(model.calculate_concentrations(increasing.C, increasing), 0)
        decreasing = self.cases[2]
        with self.assertRaisesRegex(ValueError, "positieve concentraties"):
            model.calculate_extinctions(0, decreasing)
        self.assertTrue(np.isnan(model.calculate_concentrations(decreasing.C, decreasing)))
        with self.assertRaisesRegex(CalibrationDataError, "positieve concentraties"):
            create_calibration_curve(CalibrationSettings("power_law", decreasing, concentration_range=(0, 10), response_scale="blank_corrected"))

    def test_inverse_domain_tolerance_is_independent_of_amplitude_units(self):
        model = PowerLawModel()
        parameters = PowerLawParameters(1e100, 1.7, 0.1)
        self.assertTrue(np.isnan(model.calculate_concentrations(-1, parameters)))
        for response in (float("inf"), float("-inf")):
            self.assertTrue(np.isnan(model.calculate_concentrations(response, parameters)))
        expected = np.exp((np.log(0.3) - np.log(parameters.A)) / parameters.B)
        self.assertAlmostEqual(model.calculate_concentrations(0.4, parameters) / expected, 1)

    def test_invalid_parameters_and_wrong_model_parameter_types_are_rejected(self):
        for changes in ({"A": 0}, {"B": 0}, {"A": True}, {"B": np.inf}, {"C": "0.1"}, {"C": np.nan}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                replace(self.cases[0], **changes)
        for mapping in ({"A": 1, "B": 2}, {"A": 1, "B": 2, "C": 3, "D": 4}, []):
            with self.subTest(mapping=mapping), self.assertRaises(ValueError):
                PowerLawParameters.from_mapping(mapping)
        with self.assertRaisesRegex(ValueError, "parameters moeten"):
            CalibrationSettings("4pl", self.cases[0], response_scale="blank_corrected")

    def test_exact_fits_recover_both_curvature_families_and_all_replicates(self):
        model = PowerLawModel()
        for parameters in self.cases:
            x = np.repeat(np.geomspace(0.1, 10, 7), 2)
            y = model.calculate_extinctions(x, parameters)
            for weighting in (0, 1.5, 2):
                if weighting and np.any(y <= 0):
                    continue
                with self.subTest(parameters=parameters, weighting=weighting):
                    fitted = fit_curve(x, y, model, weighting)
                    self.assertTrue(fitted.success, fitted.message)
                    np.testing.assert_allclose(astuple(fitted.parameters), astuple(parameters), rtol=1e-6, atol=1e-9)
                    np.testing.assert_allclose(fitted.diagnostics.back_calculated, x, rtol=1e-6)
                    np.testing.assert_allclose(fitted.diagnostics.predictions, y, rtol=1e-8, atol=1e-9)
                    self.assertEqual(len(fitted.diagnostics.weights), len(x))

    def test_zero_standard_is_retained_in_fitting_and_minimum_levels_are_enforced(self):
        model, parameters = PowerLawModel(), self.cases[0]
        x = np.repeat([0, 1, 2, 4, 8], 2)
        fitted = fit_curve(x, model.calculate_extinctions(x, parameters), model)
        self.assertTrue(fitted.success)
        np.testing.assert_allclose(astuple(fitted.parameters), astuple(parameters), rtol=1e-7)
        self.assertEqual(len(fitted.diagnostics.predictions), len(x))
        x = np.array([0, 1, 4])
        self.assertTrue(fit_curve(x, model.calculate_extinctions(x, parameters), model).success)
        with self.assertRaisesRegex(CalibrationDataError, "3 verschillende"):
            fit_curve([1, 1, 2, 2], [1, 1, 2, 2], model)

    def test_weighted_noisy_fit_matches_independent_optimization_in_response_space(self):
        x = np.repeat(np.geomspace(0.2, 12, 8), 2)
        parameters = self.cases[0]
        y = parameters.C + parameters.A * x ** parameters.B
        y += np.tile([0.008, -0.005], 8) * np.repeat(np.linspace(0.5, 2, 8), 2)
        for exponent in (0, 0.75, 2, -1.25):
            with self.subTest(exponent=exponent):
                weights = y ** -exponent
                expected = least_squares(
                    lambda p: np.sqrt(weights) * (p[2] + p[0] * x ** p[1] - y),
                    astuple(parameters), ftol=1e-13, xtol=1e-13, gtol=1e-13,
                )
                fitted = fit_curve(x, y, PowerLawModel(), exponent)
                self.assertTrue(fitted.success, fitted.message)
                np.testing.assert_allclose(astuple(fitted.parameters), expected.x, rtol=1e-5, atol=1e-7)
                np.testing.assert_allclose(fitted.diagnostics.weights, weights)
                self.assertAlmostEqual(fitted.cost, expected.cost, places=10)

    def test_fitting_is_independent_of_concentration_and_response_units(self):
        x = np.repeat([0, 0.1, 0.5, 1, 3, 5, 10], 2)
        parameters = self.cases[0]
        y = PowerLawModel().calculate_extinctions(x, parameters)
        for x_scale, y_scale in ((1e-80, 1), (1e80, 1), (1, 1e-100), (1, 1e100)):
            with self.subTest(x_scale=x_scale, y_scale=y_scale):
                fitted = fit_curve(x * x_scale, y * y_scale, PowerLawModel(), 2)
                self.assertTrue(fitted.success, fitted.message)
                np.testing.assert_allclose(fitted.diagnostics.predictions / y_scale, y, rtol=1e-7)
                np.testing.assert_allclose(fitted.diagnostics.back_calculated[x > 0] / x_scale, x[x > 0], rtol=1e-6)
                self.assertAlmostEqual(fitted.parameters.B, parameters.B, places=6)

    def test_settings_round_trip_fixed_curves_and_out_of_range_results(self):
        for parameters in self.cases:
            with self.subTest(parameters=parameters):
                settings = CalibrationSettings.from_dict({
                    "method": "calibration", "fit_model": "power_law", "parameters": asdict(parameters),
                    "concentration_range": [0.2, 8], "response_scale": "raw",
                })
                self.assertEqual(CalibrationSettings.from_dict(settings.to_dict()), settings)
                curve = create_calibration_curve(settings)
                responses = curve.calculate_extinction([0.1, 0.2, 1, 8, 10])
                calculated = curve.calculate_results(responses)
                self.assertEqual([item.status for item in calculated], [
                    ResultStatus.BELOW_RANGE, ResultStatus.IN_RANGE, ResultStatus.IN_RANGE, ResultStatus.IN_RANGE, ResultStatus.ABOVE_RANGE,
                ])
                self.assertEqual(calculated[0].limit, 0.2)
                self.assertEqual(calculated[-1].limit, 8)
                np.testing.assert_allclose(curve.calculate_concentration(responses), [np.nan, 0.2, 1, 8, np.nan])

    def test_raw_fitting_preserves_background_and_dilution_and_pdf_supports_fitted_and_fixed_curves(self):
        import pymupdf

        model, parameters = PowerLawModel(), PowerLawParameters(0.3, 1.7, 0.12)
        x = np.repeat([0, 0.1, 0.5, 1, 3, 5, 10], 2)
        frame = pd.DataFrame({
            "sample_id": [f"standard-{v}" for v in x] + ["patient", "blank"],
            "sample_name": [f"Standaard {v}" for v in x] + ["Patient", "Blanco"],
            "sample_type": ["Standard"] * len(x) + ["Unknown", "Blank"],
            "sample_concentration": np.r_[x, np.nan, np.nan],
            "sample_dilution": [1] * len(x) + [3, 1],
            "source_extinction": np.r_[model.calculate_extinctions(np.r_[x, 2], parameters), 0.1],
            "well": [f"{chr(65 + i % 8)}{i // 8 + 1}" for i in range(len(x) + 2)],
        })
        data = AssayData(samples_from_records(frame.to_dict("records")), assay_metadata())
        for calculation in (
            CalibrationSettings("power_law", weighting_exponent=2, response_scale="raw"),
            CalibrationSettings("power_law", replace(parameters, C=0.02), concentration_range=(0.1, 10), response_scale="blank_corrected"),
        ):
            with self.subTest(fixed=calculation.parameters is not None), TemporaryDirectory() as directory:
                result = AssayAnalyzer(data, AssaySettings("test.skax", "lot", calculation=calculation)).analyze()
                patient = next(sample for sample in result.samples if sample.sample_id == "patient")
                self.assertAlmostEqual(patient.result.mean, 6, places=5)
                self.assertAlmostEqual(result.calculation_details.curve.parameters.C, 0.02 if calculation.parameters else 0.12, places=6)
                path = PdfReportWriter(result, analyst="operator").write(Path(directory) / "power-law.pdf")
                with pymupdf.open(path) as document:
                    text = " ".join(page.get_text() for page in document)
                if calculation.parameters is None:
                    self.assertIn("Kalibratiecurve", text)
                    self.assertIn("100.0", text)
                else:
                    self.assertIn("Vaste kalibratiecurve", text)
                    self.assertNotIn("Curvemodel:", text)
                    self.assertNotIn("Parameters:", text)

    def test_power_law_can_be_an_independent_wavelength_with_saturation_fallback(self):
        data = dual_data()
        parameters = PowerLawParameters(0.003, 1.2, 0.03)
        data = replace(data, samples=tuple(replace(sample, measurements=tuple(
            replace(well, source_extinction=float(PowerLawModel().calculate_extinctions(sample.concentration, parameters)))
            if well.wavelength_nm == 405 and sample.sample_type != SampleType.BLANK else well
            for well in sample.measurements
        )) for sample in data.samples))
        config = dual_settings().to_dict()
        config["calculation"]["wavelengths"]["405"] = {"response_scale": "raw", "fit_model": "power_law", "weighting_exponent": 1.5}
        config["validation"] = {}
        settings = AssaySettings.from_dict("test.skax", config)
        self.assertEqual(AssaySettings.from_dict("test.skax", settings.to_dict()), settings)
        data = saturate(data, 450, "p1", all_replicates=True)
        result = AssayAnalyzer(data, settings).analyze()
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 405)
        self.assertAlmostEqual(patient["result_mean"], 20, places=5)
        standards = well_table(result.measurement_results[405].samples)
        standards = standards.loc[standards["sample_type"].eq("Standard")]
        np.testing.assert_allclose(standards["back_calculated_concentration"], standards["sample_concentration"], rtol=1e-6)


if __name__ == "__main__":
    unittest.main()
