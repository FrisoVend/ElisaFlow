import unittest
from dataclasses import asdict, replace
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from fixtures import assay_metadata, processed_samples, samples_from_records, with_validation
from table_helpers import sample_table, well_table

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.analysis.calculation import calculate_samples
from elisaflow.assays.registry import AssayRegistry
from elisaflow.calibration.curves import CalibrationCurve
from elisaflow.calibration.fitting import FitDiagnostics
from elisaflow.calibration.logistic import FivePLModel, FivePLParameters, FourPLModel, FourPLParameters
from elisaflow.calibration.polynomial import QuadraticParameters
from elisaflow.configuration.assay_settings import AssaySettings, CalibrationSettings, ReplicateFilterSettings, SampleSettings
from elisaflow.formulas.expressions import Formula
from elisaflow.models import AssayData, SampleType, ValidationStatus


class CalibrationQualityTests(unittest.TestCase):
    @staticmethod
    def calculate(records):
        """Assess stored R² through analysis with a known curve, independently of fitting."""
        samples = processed_samples(records)
        count = sum(len(sample.retained) for sample in samples if sample.sample_type == SampleType.STANDARD)
        curve = CalibrationCurve(
            FourPLModel(), FourPLParameters(2, -1, 1, 0),
            diagnostics=FitDiagnostics(*(np.zeros(count) for _ in range(5))),
        )
        with patch("elisaflow.analysis.calculation.create_calibration_curve", return_value=curve):
            return calculate_samples(
                samples, CalibrationSettings("4pl", response_scale="raw"),
                wavelength_nm=450, reference_wavelength_nm=None, blank_value=0.2, allow_fit_failure=False,
            )

    def test_r_squared_uses_retained_standard_means_on_the_recorded_curve_input_scale(self):
        records = [
            {
                "sample_id": f"standard-{nominal}", "sample_type": "Standard", "sample_concentration": nominal,
                "source_extinction": response, "extinction": 100,
            }
            for nominal, responses in ((0, [0]), (1, [0.8, 1.2]), (3, [1.8, 2, 2.2]), (5, [np.nan]))
            for response in responses
        ]
        records += [
            {
                "sample_id": "standard-1", "sample_type": "Standard", "sample_concentration": 1,
                "source_extinction": 500, "included_in_analysis": False,
            },
            {"sample_type": "Unknown", "source_extinction": 999},
        ]
        outcome = self.calculate(records)
        # Standard means [0, 1, 2] and predictions [0, 1, 1.5]: SSE=0.25, SST=2.
        self.assertAlmostEqual(outcome.details.curve_r_squared, 0.875)

    def test_r_squared_for_perfect_and_poor_fits(self):
        for responses, expected in (([0, 1, 1.5], 1), ([1, 2, 3], -1.125)):
            with self.subTest(responses=responses):
                outcome = self.calculate([
                    {"sample_type": "Standard", "sample_concentration": nominal, "source_extinction": response}
                    for nominal, response in zip([0, 1, 3], responses, strict=True)
                ])
                self.assertAlmostEqual(outcome.details.curve_r_squared, expected)

    def test_r_squared_is_unavailable_without_finite_varying_standard_means(self):
        for responses in ([], [np.nan, np.inf], [1], [1, 1, 1], [1, np.nan]):
            with self.subTest(responses=responses):
                outcome = self.calculate([
                    {"sample_type": "Standard", "sample_concentration": nominal, "source_extinction": response}
                    for nominal, response in enumerate(responses)
                ])
                self.assertIsNone(outcome.details.curve_r_squared)


class CalibrationResponseTests(unittest.TestCase):
    @staticmethod
    def data(model, parameters, *, blank=0.1):
        nominal = np.repeat([0, 0.1, 0.4, 1, 4, 10, 50], 3)
        concentrations = np.r_[nominal, 2.5, 1, 1]
        frame = pd.DataFrame(
            {
                "sample_id": [f"standard-{x}" for x in nominal] + ["patient", "control", "blank"],
                "sample_name": [f"Standard {x}" for x in nominal] + ["Patient", "Control", "Blank"],
                "sample_type": ["Standard"] * len(nominal) + ["Unknown", "Control", "Blank"],
                "sample_concentration": np.r_[nominal, np.nan, np.nan, np.nan],
                "sample_dilution": [1] * len(nominal) + [3, 1, 1],
                "source_extinction": model.calculate_extinctions(concentrations, parameters),
                "well": [f"A{index + 1}" for index in range(len(concentrations))],
            }
        )
        frame.loc[frame["sample_type"] == "Blank", "source_extinction"] = blank
        return AssayData(samples_from_records(frame.to_dict('records')), assay_metadata())

    @staticmethod
    def settings(model_name, *, response_scale, **options):
        return AssaySettings("test.skax", "lot", calculation=CalibrationSettings(model_name, response_scale=response_scale, **options))

    def test_fitted_and_fixed_curves_save_their_explicit_response_scale(self) -> None:
        for model_name in ("4pl", "5pl"):
            parameters = FourPLParameters(2.5, -1.3, 2, 0.1) if model_name == "4pl" else FivePLParameters(2.5, -1.3, 2, 0.1, 0.4)
            for scale in ("raw", "blank_corrected"):
                settings = CalibrationSettings.from_dict(
                    {"method": "calibration", "fit_model": model_name, "response_scale": scale, "parameters": asdict(parameters)}
                )
                self.assertEqual(settings.response_scale, scale)
                self.assertEqual(CalibrationSettings.from_dict(settings.to_dict()), settings)
                fitted = CalibrationSettings.from_dict({"method": "calibration", "fit_model": model_name, "response_scale": scale})
                self.assertEqual(fitted.response_scale, scale)
                self.assertEqual(CalibrationSettings.from_dict(fitted.to_dict()), fitted)
            fitted = CalibrationSettings.from_dict({"response_scale": "raw", "method": "calibration", "fit_model": model_name, "weighting_exponent": 2.0})
            self.assertEqual(fitted.response_scale, "raw")
            self.assertEqual(fitted.to_dict()["response_scale"], "raw")
            self.assertEqual(CalibrationSettings.from_dict(fitted.to_dict()), fitted)
        for value in (None, True, [], {}, "", "corrected"):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "response_scale"):
                CalibrationSettings.from_dict({"method": "calibration", "fit_model": "4pl", "response_scale": value})

    def test_raw_fitting_weights_inversion_and_qc_for_four_and_five_pl(self) -> None:
        cases = (
            ("4pl", FourPLModel(), FourPLParameters(2.5, -1.3, 2, 0.1)),
            ("4pl", FourPLModel(), FourPLParameters(2.5, 1.3, 2, 0.1)),
            ("5pl", FivePLModel(), FivePLParameters(2.5, -1.3, 2, 0.1, 0.4)),
            ("5pl", FivePLModel(), FivePLParameters(2.5, 1.3, 2, 0.1, 0.4)),
        )
        for name, model, parameters in cases:
            with self.subTest(model=name, slope=parameters.B):
                data = self.data(model, parameters)
                original = well_table(data.samples).copy(deep=True)
                settings = self.settings(name, weighting_exponent=2.0, response_scale="raw")
                result = AssayAnalyzer(data, settings).analyze()
                wells = well_table(result.samples)
                standards = wells[wells["sample_type"] == "Standard"]
                diagnostics = result.calculation_details.curve.diagnostics
                np.testing.assert_allclose(diagnostics.weights, 1 / standards["source_extinction"] ** 2)
                np.testing.assert_allclose(diagnostics.predictions, standards["source_extinction"], atol=1e-8)
                np.testing.assert_allclose(standards["curve_extinction"], standards["source_extinction"])
                np.testing.assert_allclose(wells["extinction"], wells["source_extinction"] - 0.1)
                np.testing.assert_allclose(wells.loc[wells["sample_type"] == "Unknown", "result"], 7.5, rtol=1e-6)
                np.testing.assert_allclose(wells.loc[wells["sample_type"] == "Control", "result"], 1, rtol=1e-6)
                samples = sample_table(result.samples)
                standards_summary = samples[samples["sample_type"] == "Standard"]
                np.testing.assert_allclose(standards_summary["curve_extinction_mean"], standards_summary["extinction_mean"] + 0.1)
                positive = standards_summary["sample_concentration"] > 0
                np.testing.assert_allclose(standards_summary.loc[positive, "recovery_percent"], 100, rtol=1e-6)
                self.assertAlmostEqual(result.calculation_details.curve_r_squared, 1)
                self.assertEqual(result.blank_value, 0.1)
                pd.testing.assert_frame_equal(well_table(data.samples), original)
                # A different measured blank changes QC extinctions, but not
                # the raw fit, its weights, or the patient calculation.
                changed = self.data(model, parameters, blank=0.08)
                other = AssayAnalyzer(changed, settings).analyze()
                self.assertEqual(other.calculation_details.curve.parameters, result.calculation_details.curve.parameters)
                np.testing.assert_allclose(well_table(other.samples)["result"], wells["result"], equal_nan=True)

    def test_fixed_curves_use_the_explicit_response_scale(self) -> None:
        model = FourPLModel()
        raw_parameters = FourPLParameters(2.5, -1.3, 2, 0.1)
        data = self.data(model, raw_parameters)
        for scale, parameters in (
            ("raw", raw_parameters),
            ("blank_corrected", replace(raw_parameters, A=2.4, D=0)),
        ):
            with self.subTest(scale=scale):
                settings = AssaySettings("test.skax", "lot", calculation=CalibrationSettings("4pl", parameters, response_scale=scale))
                result = AssayAnalyzer(data, settings).analyze()
                patient = well_table(result.samples)[well_table(result.samples)["sample_type"] == "Unknown"]
                np.testing.assert_allclose(patient["result"], 7.5)
                expected_input = patient["source_extinction"] - (result.blank_value if scale == "blank_corrected" else 0)
                np.testing.assert_allclose(patient["curve_extinction"], expected_input)
                self.assertIsNone(result.calculation_details.curve.diagnostics)

    def test_blank_corrected_four_and_five_pl_fits_use_corrected_weights_and_inversion(self) -> None:
        for name, model, parameters in (
            ("4pl", FourPLModel(), FourPLParameters(2.5, -1.3, 2, 0.1)),
            ("5pl", FivePLModel(), FivePLParameters(2.5, -1.3, 2, 0.1, 0.4)),
        ):
            with self.subTest(model=name):
                result = AssayAnalyzer(
                    self.data(model, parameters, blank=0.08),
                    self.settings(name, response_scale="blank_corrected", weighting_exponent=2),
                ).analyze()
                wells = well_table(result.samples)
                standards = wells[wells["sample_type"] == "Standard"]
                diagnostics = result.calculation_details.curve.diagnostics
                np.testing.assert_allclose(diagnostics.weights, standards["extinction"] ** -2)
                np.testing.assert_allclose(diagnostics.predictions, standards["extinction"], atol=1e-8)
                np.testing.assert_allclose(wells.loc[wells["sample_type"] == "Unknown", "result"], 7.5, rtol=1e-6)
                np.testing.assert_allclose(wells.loc[wells["sample_type"] == "Control", "result"], 1, rtol=1e-6)
                patient = next(sample for sample in result.samples if sample.sample_type == SampleType.UNKNOWN)
                self.assertAlmostEqual(result.patient_extinction(patient).scale.blank_subtracted, 0.08)

    def test_raw_curve_preprocessing_uses_original_raw_reference_means(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": ["blank", "reference", "reference", "control", "patient"],
                "sample_name": ["Blank", "Reference", "Reference", "Control", "Patient"],
                "sample_type": ["Blank", "Control", "Control", "Control", "Unknown"],
                "source_extinction": [0.1, 0.8, 1, 1.8, 1.1],
            }
        )
        controls = {
            "Reference": SampleSettings(),
            "Control": SampleSettings(),
        }
        parameters = FourPLParameters(2.5, -1.3, 2, 0.1)
        settings = self.settings(
            "4pl", parameters=parameters, response_scale="raw", preprocessing=Formula("(sample * reference)", {"reference": "Reference"})
        )
        for ordered in (controls, dict(reversed(list(controls.items())))):
            result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), assay_metadata()), replace(settings, controls=ordered)).analyze()
            inputs = frame["source_extinction"].to_numpy() * 0.9
            inputs[frame["sample_type"].eq("Blank")] = np.nan
            expected = FourPLModel().calculate_concentrations(inputs, parameters)
            np.testing.assert_allclose(well_table(result.samples)["result"], expected, equal_nan=True)
            np.testing.assert_allclose(well_table(result.samples)["curve_extinction"], inputs, equal_nan=True)
            np.testing.assert_allclose(well_table(result.samples)["extinction"], [0, 0.7, 0.9, 1.7, 1])

    def test_reciprocal_preprocessing_skips_blanks_for_fixed_curves(self) -> None:
        frame = pd.DataFrame({
            "sample_id": ["blank", "s1", "s2", "s4", "control", "patient"],
            "sample_name": ["Blank", "Standard 1", "Standard 2", "Standard 4", "Control", "Patient"],
            "sample_type": ["Blank", "Standard", "Standard", "Standard", "Control", "Unknown"],
            "sample_concentration": [np.nan, 1, 2, 4, np.nan, np.nan],
            "source_extinction": [0.1, 1.1, 0.6, 0.35, 0.35, 0.6],
        })
        data = AssayData(samples_from_records(frame.to_dict("records")), assay_metadata())
        settings = with_validation(
            self.settings(
                "quadratic", parameters=QuadraticParameters(0, 1, 0), concentration_range=(1, 4),
                preprocessing=Formula("1 / sample"), response_scale="blank_corrected",
            ),
            blank={"extinction": {"response_scale": "raw", "max": 0.15}},
            controls={"Control": {"result": {"min": 3.9, "max": 4.1}}},
        )
        result = AssayAnalyzer(data, settings).analyze()
        self.assertAlmostEqual(result.blank_value, 0.1)
        self.assertEqual(result.validation_status, ValidationStatus.PASSED)
        samples = sample_table(result.samples).set_index("sample_id")
        np.testing.assert_allclose(samples.loc[["s1", "s2", "s4", "control", "patient"], "result_mean"], [1, 2, 4, 4, 2])
        self.assertEqual(samples.loc["blank", "extinction_mean"], 0)
        self.assertTrue(np.isnan(samples.loc["blank", "result_mean"]))
        blank = next(sample for sample in result.samples if sample.sample_type == SampleType.BLANK)
        self.assertTrue(all(well.included and well.curve_extinction is None for well in blank.measurements))

    def test_replicate_filter_uses_final_extinction_while_fitting_uses_raw(self) -> None:
        model = FourPLModel()
        parameters = FourPLParameters(2.5, -1.3, 2, 0.1)
        data = self.data(model, parameters, blank=0.08)
        _readings = well_table(data.samples)
        _readings.loc[2, "source_extinction"] = 0.108
        data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
        settings = replace(
            self.settings("4pl", weighting_exponent=2.0, response_scale="raw"),
            standards={
                0: SampleSettings(replicate_filter=ReplicateFilterSettings(15, 1, response_scale="blank_corrected")),
            },
        )
        result = AssayAnalyzer(data, settings).analyze()
        excluded = well_table(result.samples)[~well_table(result.samples)["included_in_analysis"]]
        self.assertEqual(excluded.index.tolist(), [2])
        final_mean = np.mean([0.02, 0.02, 0.028])
        self.assertAlmostEqual(excluded.iloc[0]["replicate_mean"], final_mean)
        self.assertAlmostEqual(excluded.iloc[0]["replicate_deviation_percent"], 100 * (0.028 - final_mean) / final_mean)
        # On the raw scale the outlier is only about 5% from the mean.
        self.assertLess(100 * (0.108 - np.mean([0.1, 0.1, 0.108])) / np.mean([0.1, 0.1, 0.108]), 15)
        self.assertEqual(len(result.calculation_details.curve.diagnostics.weights), 20)
        included = well_table(result.samples)[well_table(result.samples)["included_in_analysis"] & well_table(result.samples)["sample_type"].eq("Standard")]
        np.testing.assert_allclose(result.calculation_details.curve.diagnostics.weights, 1 / included["source_extinction"] ** 2)
        zero = sample_table(result.samples)[sample_table(result.samples)["sample_concentration"] == 0].iloc[0]
        self.assertEqual(zero["count"], 2)
        self.assertAlmostEqual(zero["curve_extinction_mean"], 0.1)
        self.assertAlmostEqual(zero["extinction_mean"], 0.02)
        for scale in ("raw", "blank_corrected"):
            fixed_parameters = parameters if scale == "raw" else replace(parameters, A=2.42, D=0.02)
            fixed_settings = replace(settings, calculation=CalibrationSettings("4pl", fixed_parameters, response_scale=scale))
            fixed = AssayAnalyzer(data, fixed_settings).analyze()
            np.testing.assert_array_equal(well_table(fixed.samples)["included_in_analysis"], well_table(result.samples)["included_in_analysis"])

    def test_raw_calibration_preserves_raw_blank_validation(self) -> None:
        data = self.data(FourPLModel(), FourPLParameters(2.5, -1.3, 2, 0.1))
        settings = with_validation(
            replace(self.settings("4pl", response_scale="raw"), blank=SampleSettings()), blank={"extinction": {"response_scale": "raw", "max": 0.09}}
        )
        result = AssayAnalyzer(data, settings).analyze()
        self.assertEqual(result.validation_status, ValidationStatus.FAILED)
        self.assertAlmostEqual(result.validation_results[0].value, 0.1)

    def test_numeric_source_strings_and_instrument_exclusions_are_supported(self) -> None:
        data = self.data(FourPLModel(), FourPLParameters(2.5, -1.3, 2, 0.1))
        _readings = well_table(data.samples)
        _readings["disabled"] = False
        data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
        _readings = well_table(data.samples)
        _readings.loc[2, "disabled"] = True
        data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
        _readings = well_table(data.samples)
        _readings["source_extinction"] = _readings["source_extinction"].astype(str)
        data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
        result = AssayAnalyzer(data, self.settings("4pl", weighting_exponent=2.0, response_scale="raw")).analyze()
        self.assertFalse(well_table(result.samples).loc[2, "included_in_analysis"])
        self.assertTrue(np.isnan(well_table(result.samples).loc[2, "result"]))
        self.assertEqual(len(result.calculation_details.curve.diagnostics.weights), 20)

    def test_bundled_fitted_assays_explicitly_use_raw_response(self) -> None:
        project = Path(__file__).resolve().parents[1]
        registry = AssayRegistry(project / "configs")
        for name in ("IA-2_as", "GAD_as", "ZnT8_as", "Hib_IgG", "Hib_IgG_duplo"):
            with self.subTest(assay=name):
                self.assertEqual(registry.get(project / "protocols" / f"{name}.skax").calculation.response_scale, "raw")
        fixed = registry.get(project / "protocols" / "Qkoorts_IgG.skax")
        self.assertEqual(fixed.calculation.response_scale, "blank_corrected")


if __name__ == "__main__":
    unittest.main()
