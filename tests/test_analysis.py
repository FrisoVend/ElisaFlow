import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
from fixtures import assay_metadata, borrelia_settings, reference_parameters, samples_from_records
from table_helpers import sample_table, well_table

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.assays.registry import AssayRegistry
from elisaflow.calibration.logistic import FourPLModel, FourPLParameters
from elisaflow.configuration.assay_settings import (
    AssaySettings,
    CalibrationSettings,
    FormulaOverrideSettings,
    FormulaSettings,
    MeasurementSettings,
    SampleSettings,
    ValidationPlan,
)
from elisaflow.formulas.expressions import Formula
from elisaflow.models import (
    AnalysisResult,
    AssayData,
    ResultStatus,
    Sample,
    SampleType,
    ValidationMetric,
    ValidationStatus,
    WellMeasurement,
)


class AnalysisTests(unittest.TestCase):
    def setUp(self) -> None:
        self.metadata = assay_metadata(name="test", session_time="2026-10-05T10:00:00+02:00")
        self.parameters = reference_parameters()

    def test_missing_control_replicate_cannot_pass_extinction_qc(self) -> None:
        settings = AssayRegistry(Path(__file__).resolve().parents[1] / "configs").get("Qkoorts_IgM.skax")
        for flag in (None, "disabled", "saturated"):
            with self.subTest(exclusion=flag):
                config = settings.to_dict()
                assay = AssaySettings.from_dict(settings.protocol_name, config)
                frame = pd.DataFrame({
                    "sample_id": ["blank", "positive", "positive", "cutoff", "negative", "patient"],
                    "sample_name": ["Blank", "Positief", "Positief", "Grens", "Negatief", "Patient"],
                    "sample_type": ["Blank", "Control", "Control", "Control", "Control", "Unknown"],
                    "well": ["A1", "A2", "A3", "A4", "A5", "A6"],
                    "source_extinction": [0.15, 0.85, np.nan, 0.65, 0.35, 0.75],
                })
                if flag:
                    frame[flag] = frame["well"].eq("A3")
                reference = frame.assign(wavelength_nm=620, source_extinction=0.05)
                wells = pd.concat([frame.assign(wavelength_nm=450), reference], ignore_index=True)
                result = AssayAnalyzer(AssayData(samples_from_records(wells.to_dict('records')), self.metadata), assay).analyze()
                check = next(check for check in result.validation_results if check.sample_name == "Positief")
                if flag:
                    self.assertEqual(result.validation_status, ValidationStatus.PASSED)
                    self.assertAlmostEqual(check.value, 0.7)
                else:
                    self.assertEqual(result.validation_status, ValidationStatus.NOT_EVALUATED)
                    self.assertEqual(check.status, ValidationStatus.NOT_EVALUATED)
                    self.assertIsNone(check.value)
                    positive = sample_table(result.samples).loc[sample_table(result.samples)["sample_name"].eq("Positief")].iloc[0]
                    for column in ("extinction_mean", "extinction_sd", "extinction_cv_percent"):
                        self.assertTrue(np.isnan(positive[column]), column)

    def test_only_empty_unknown_wells_are_excluded_before_calculations(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": ["blank", "blank", "control", "control", "patient", "patient", "empty", "standard"],
                "sample_name": ["Blank", "Blank", "Grens", "Grens", "Patient", "Patient", "Empty", "Standard"],
                "sample_type": ["Blank", "Blank", "Control", "Control", "Unknown", "Unknown", "Unknown", "Standard"],
                "source_extinction": [0.008, 0.1, 0.008, 0.6, 0.008, 0.009, -0.1, 0],
                "well": [f"A{index + 1}" for index in range(8)],
            }
        )
        original = frame.copy(deep=True)
        settings = AssaySettings(
            "test.skax", "lot", calculation=FormulaSettings(Formula("sample / grens", {"grens": "Grens"}), {}, response_scale="blank_corrected")
        )
        analyzer = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings)
        result = analyzer.analyze()
        self.assertEqual(well_table(result.samples)["well"].tolist(), ["A1", "A2", "A3", "A4", "A6", "A8"])
        self.assertAlmostEqual(result.blank_value, 0.054)

        patient = sample_table(result.samples).loc[sample_table(result.samples)["sample_id"].eq("patient")].iloc[0]
        self.assertEqual(patient["count"], 1)
        self.assertEqual(patient["used_wells"], "A6")
        self.assertAlmostEqual(patient["extinction_mean"], -0.045)
        self.assertAlmostEqual(patient["result_mean"], -0.18)
        self.assertTrue(np.isnan(patient["extinction_sd"]))
        self.assertTrue(np.isnan(patient["result_sd"]))
        pd.testing.assert_frame_equal(frame, original)
        pd.testing.assert_frame_equal(well_table(analyzer.analyze().samples), well_table(result.samples))

    def test_empty_well_filter_preserves_missing_values_and_infinite_value_errors(self) -> None:
        settings = AssaySettings(
            "test.skax", "lot", calculation=FormulaSettings(Formula("sample / grens", {"grens": "Grens"}), {}, response_scale="blank_corrected")
        )
        for value in (np.nan, np.inf, -np.inf):
            with self.subTest(value=value):
                frame = pd.DataFrame(
                    {
                        "sample_id": ["control", "patient"],
                        "sample_name": ["Grens", "Patient"],
                        "sample_type": ["Control", "Unknown"],
                        "source_extinction": [0.5, value],
                    }
                )
                analyzer = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings)
                if np.isinf(value):
                    with self.assertRaisesRegex(ValueError, 'Bronmetingen mogen niet oneindig zijn'):
                        analyzer.analyze()
                else:
                    result = analyzer.analyze()
                    self.assertEqual(len(well_table(result.samples)), 2)
                    patient = sample_table(result.samples).iloc[1]
                    self.assertTrue(np.isnan(patient["result_mean"]))

    def test_json_acceptance_rules_are_attached_to_analysis_results(self) -> None:
        registry = AssayRegistry(Path(__file__).resolve().parents[1] / "configs")
        settings = replace(registry.get("Qkoorts_IgM.skax"), measurement=MeasurementSettings())
        frame = pd.DataFrame(
            {
                "sample_id": ["blank", "positive", "cutoff", "negative", "patient"],
                "sample_name": ["Blank", "Positief", "Grens", "Negatief", "Patient"],
                "sample_type": ["Blank", "Control", "Control", "Control", "Unknown"],
                "source_extinction": [0.3, 0.9, 0.7, 0.4, 0.8],
            }
        )
        analyzer = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings)
        result = analyzer.analyze()
        self.assertEqual(len(result.validation_results), 4)
        blank = next(check for check in result.validation_results if check.sample_type == SampleType.BLANK)
        self.assertAlmostEqual(blank.value, 0.3)
        self.assertEqual(blank.status, ValidationStatus.FAILED)
        self.assertEqual(result.validation_status, ValidationStatus.FAILED)
        positive = next(check for check in result.validation_results if check.sample_name == "Positief")
        self.assertAlmostEqual(positive.value, 0.6)
        self.assertEqual(positive.status, ValidationStatus.PASSED)
        self.assertEqual(analyzer.analyze().validation_results, result.validation_results)

        missing_negative = frame.loc[frame["sample_name"].ne("Negatief")]
        failed_and_incomplete = AssayAnalyzer(AssayData(samples_from_records(missing_negative.to_dict('records')), self.metadata), settings).analyze()
        self.assertEqual(failed_and_incomplete.validation_status, ValidationStatus.FAILED)
        self.assertIn(
            ValidationStatus.NOT_EVALUATED,
            [check.status for check in failed_and_incomplete.validation_results],
        )

        # Reducing the raw blank leaves all configured control criteria passing.
        frame.loc[0, "source_extinction"] = 0.1
        passed = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
        self.assertEqual(passed.validation_status, ValidationStatus.PASSED)

        missing_negative = frame.loc[frame["sample_name"].ne("Negatief")]
        incomplete = AssayAnalyzer(AssayData(samples_from_records(missing_negative.to_dict('records')), self.metadata), settings).analyze()
        self.assertEqual(incomplete.validation_status, ValidationStatus.NOT_EVALUATED)

        no_rules = replace(settings, validation=ValidationPlan())
        unchecked = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), no_rules).analyze()
        self.assertEqual(unchecked.validation_results, ())
        self.assertEqual(unchecked.validation_status, ValidationStatus.NOT_EVALUATED)

    def test_gad_rules_validate_fitted_concentrations_and_nominal_standard(self) -> None:
        registry = AssayRegistry(Path(__file__).resolve().parents[1] / "configs")
        settings = registry.get("test_GAD_as.skax")
        concentrations = np.repeat([0, 5, 25, 50, 100, 250], 2)
        frame = pd.DataFrame(
            {
                "sample_id": [f"standard-{x}" for x in concentrations],
                "sample_name": [f"Calibrator {x}" for x in concentrations],
                "sample_type": "Standard",
                "sample_concentration": concentrations,
                "source_extinction": FourPLModel().calculate_extinctions(concentrations, self.parameters) + 0.05,
            }
        )
        frame = pd.concat(
            [
                frame,
                pd.DataFrame(
                    {
                        "sample_id": ["blank", "positive"],
                        "sample_name": ["Blank", "Positief"],
                        "sample_type": ["Blank", "Control"],
                        "sample_concentration": [np.nan, np.nan],
                        "source_extinction": [
                            0.05,
                            FourPLModel().calculate_extinctions(51, self.parameters) + 0.05,
                        ],
                    }
                ),
            ],
            ignore_index=True,
        )
        frame["well"] = [f"{'ABCDEFGH'[index % 8]}{index // 8 + 1}" for index in range(len(frame))]
        data = pd.concat(
            [
                frame.assign(wavelength_nm=450),
                frame.assign(wavelength_nm=405, source_extinction=frame["source_extinction"] * 0.5),
            ],
            ignore_index=True,
        )
        result = AssayAnalyzer(AssayData(samples_from_records(data.to_dict('records')), self.metadata), settings).analyze()
        self.assertEqual(set(result.measurement_results), {405, 450})
        self.assertTrue(all(measurement.calculation_details.curve.was_fitted for measurement in result.measurement_results.values()))
        self.assertEqual(len(result.validation_results), 3)
        self.assertEqual(result.validation_status, ValidationStatus.PASSED)
        control = next(check for check in result.validation_results if check.sample_type == SampleType.CONTROL)
        self.assertAlmostEqual(control.value, 51, places=5)
        standard = next(check for check in result.validation_results if check.sample_type == SampleType.STANDARD)
        self.assertEqual(standard.sample_name, "Standaard 250")
        self.assertEqual(standard.wavelength_nm, 450)
        self.assertGreaterEqual(standard.value, 2)

    def test_concentration_validation_uses_precomputed_sample_means(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": ["control-a", "control-b", "control-b"],
                "sample_name": ["Positive"] * 3,
                "sample_type": ["Control"] * 3,
                "source_extinction": FourPLModel().calculate_extinctions([1, 4, 4], self.parameters),
            }
        )
        settings = AssaySettings.from_dict(
            "test.skax",
            {
                "lot_number": "lot",
                "controls": {"Positive": {}},
                "validation": {"result": {"controls": {"Positive": {"result": {"min": 2.4, "max": 2.6}}}}},
                "calculation": {"response_scale": "blank_corrected", "method": "calibration", "fit_model": "4pl", "parameters": {"A": 2.2, "B": -1.3, "C": 2, "D": 0.2}},
                "result": {"label": "Conc"},
            },
        )
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
        (check,) = result.validation_results
        self.assertEqual(check.metric, ValidationMetric.RESULT)
        self.assertAlmostEqual(check.value, 2.5)
        np.testing.assert_allclose(sample_table(result.samples)["result_mean"], [1, 4])
        self.assertAlmostEqual(well_table(result.samples)["result"].mean(), 3)
        self.assertEqual(result.validation_status, ValidationStatus.PASSED)

        frame.loc[1, "source_extinction"] = self.parameters.A + 1
        incomplete = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
        self.assertEqual(incomplete.validation_results[0].status, ValidationStatus.NOT_EVALUATED)
        self.assertEqual(incomplete.validation_status, ValidationStatus.NOT_EVALUATED)

    def test_fitted_analysis_returns_results_without_mutating_input(self) -> None:
        concentrations = np.repeat([0.0, 0.25, 0.75, 2.0, 6.0, 20.0], 2)
        responses = FourPLModel().calculate_extinctions(concentrations, self.parameters)
        frame = pd.DataFrame(
            {
                "sample_id": [f"standard-{x}" for x in concentrations],
                "sample_name": [f"Standard {x}" for x in concentrations],
                "sample_type": SampleType.STANDARD.value,
                "sample_concentration": concentrations,
                "source_extinction": responses + 0.05,
            }
        )
        frame = pd.concat(
            [
                frame,
                pd.DataFrame(
                    {
                        "sample_id": ["blank", "blank", "patient-a", "patient-b"],
                        "sample_name": ["Blank", "Blank", "Same name", "Same name"],
                        "sample_type": ["Blank", "Blank", "Unknown", "Unknown"],
                        "sample_concentration": [np.nan] * 4,
                        "source_extinction": [
                            0.04,
                            0.06,
                            *(FourPLModel().calculate_extinctions([1.0, 4.0], self.parameters) + 0.05),
                        ],
                    }
                ),
            ],
            ignore_index=True,
        )
        original = frame.copy(deep=True)
        settings = AssaySettings(
            "test.skax", "lot", calculation=CalibrationSettings("4pl", response_scale="raw"), standards={250.0: SampleSettings()}
        )
        analyzer = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings)
        result = analyzer.analyze()
        self.assertIsInstance(result, AnalysisResult)
        pd.testing.assert_frame_equal(frame, original)
        self.assertAlmostEqual(result.blank_value, 0.05)
        details = result.calculation_details
        self.assertTrue(details.curve.was_fitted)
        self.assertIsNotNone(details.curve.diagnostics)
        self.assertEqual(details.curve.concentration_range, (0.25, 20.0))
        patients = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")]
        np.testing.assert_allclose(patients["result_mean"], [1.0, 4.0], rtol=1e-6)
        self.assertEqual(len(patients), 2)
        second = analyzer.analyze()
        view = well_table(result.samples)
        view.loc[0, "extinction"] = 999
        self.assertNotEqual(well_table(second.samples).loc[0, "extinction"], 999)

    def test_configured_analysis_needs_no_standards_and_requires_complete_replicates(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": ["patient"] * 3,
                "sample_name": ["Patient"] * 3,
                "sample_type": ["Unknown"] * 3,
                "source_extinction": [
                    FourPLModel().calculate_extinctions(1.0, self.parameters),
                    FourPLModel().calculate_extinctions(4.0, self.parameters),
                    np.nan,
                ],
            }
        )
        settings = AssaySettings("test.skax", "lot", calculation=CalibrationSettings("4pl", self.parameters, response_scale="blank_corrected"))
        for invalid_response, count in ((np.nan, 2), (self.parameters.A, 3)):
            with self.subTest(invalid_response=invalid_response):
                frame.loc[2, "source_extinction"] = invalid_response
                result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
                self.assertFalse(result.calculation_details.curve.was_fitted)
                self.assertIs(result.calculation_details.curve.parameters, self.parameters)
                sample = sample_table(result.samples).iloc[0]
                self.assertEqual(sample["count"], count)
                for column in (
                    "result_mean",
                    "result_sd",
                    "result_cv_percent",
                ):
                    self.assertTrue(np.isnan(sample[column]), column)

    def test_calibration_keeps_results_with_their_wells_after_reordering(self) -> None:
        frame = pd.DataFrame({
            "sample_id": ["paired", "other", "paired", "disabled", "blank"],
            "sample_name": ["Patient", "Patient", "Patient", "Disabled", "Blank"],
            "sample_type": ["Unknown"] * 4 + ["Blank"],
            "sample_dilution": [3, 10, 3, 1, 1],
            "well": ["H12", "A3", "C7", "D2", "F6"],
            "source_extinction": [*FourPLModel().calculate_extinctions([2, 1, 5, 7], self.parameters) + 0.05, 0.05],
            "disabled": [False, False, False, True, False],
        })
        data = AssayData(samples_from_records(frame.to_dict("records")), self.metadata)
        settings = AssaySettings("test.skax", "lot", calculation=CalibrationSettings("4pl", self.parameters, response_scale="blank_corrected"))
        shuffled = replace(data, samples=tuple(replace(sample, measurements=sample.measurements[::-1]) for sample in data.samples[::-1]))
        for inputs in (data, shuffled):
            with self.subTest(reordered=inputs is shuffled):
                result = AssayAnalyzer(inputs, settings).analyze()
                wells = well_table(result.samples).set_index("well")
                np.testing.assert_allclose(wells.loc[["H12", "A3", "C7"], "result"], [6, 10, 15], rtol=1e-12)
                self.assertTrue(np.isnan(wells.loc["D2", "result"]))
                self.assertFalse(wells.loc["D2", "included_in_analysis"])
                np.testing.assert_allclose(wells.loc[frame["well"], "source_extinction"], frame["source_extinction"])
        self.assertTrue(all(np.isnan(well.result) for sample in data.samples for well in sample.measurements))

    def test_equal_readings_in_distinct_wells_keep_independent_results(self) -> None:
        first, second = (
            WellMeasurement(source_extinction=float(value), well=f"A{index}", wavelength_nm=450)
            for index, value in enumerate(FourPLModel().calculate_extinctions([1, 2], self.parameters), 1)
        )
        data = AssayData((
            Sample("paired", "Patient", SampleType.UNKNOWN, (first, second), dilution=3),
            Sample("single", "Patient", SampleType.UNKNOWN, (replace(first, well="B1"),), dilution=10),
        ), self.metadata)
        settings = AssaySettings.from_dict("test.skax", {
            "lot_number": "lot",
            "calculation": CalibrationSettings("4pl", self.parameters, response_scale="blank_corrected").to_dict(),
            "patients": {"replicate_filter": {"response_scale": "blank_corrected", "max_deviation_percent": 100}},
        })
        result = AssayAnalyzer(data, settings).analyze()
        np.testing.assert_allclose([well.result for sample in result.samples for well in sample.measurements], [3, 6, 10], rtol=1e-12)

    def test_single_wavelength_rejects_duplicate_wells_across_or_within_samples(self) -> None:
        well = WellMeasurement(source_extinction=0.5, well="A1", wavelength_nm=450)
        first = Sample("p1", "Patient", SampleType.UNKNOWN, (well,))
        settings = AssaySettings("test.skax", "lot", calculation=FormulaSettings(Formula("sample"), response_scale="blank_corrected"))
        for samples in (
            (replace(first, measurements=(well, replace(well, order=1))),),
            (first, replace(first, sample_id="p2", measurements=(replace(well, order=1),))),
        ):
            with self.subTest(samples=len(samples)), self.assertRaisesRegex(ValueError, "Dubbele wellmetingen"):
                AssayAnalyzer(AssayData(samples, self.metadata), settings).analyze()

    def test_concentration_statistics_use_calculated_replicates(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": ["patient", "patient", "single"],
                "sample_name": ["Patient", "Patient", "Single"],
                "sample_type": ["Unknown"] * 3,
                "source_extinction": FourPLModel().calculate_extinctions([1.0, 4.0, 2.0], self.parameters),
            }
        )
        settings = AssaySettings("test.skax", "lot", calculation=CalibrationSettings("4pl", self.parameters, response_scale="blank_corrected"))
        samples = sample_table(AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze().samples)
        self.assertNotIn("concentration_count", samples.columns)
        self.assertNotIn("result", samples.columns)
        patient = samples.loc[samples["sample_id"].eq("patient")].iloc[0]
        self.assertEqual(patient["count"], 2)
        self.assertAlmostEqual(patient["result_mean"], 2.5)
        self.assertAlmostEqual(patient["result_sd"], 3 / np.sqrt(2))
        self.assertAlmostEqual(patient["result_cv_percent"], 100 * 3 / np.sqrt(2) / 2.5)
        single = samples.loc[samples["sample_id"].eq("single")].iloc[0]
        self.assertEqual(single["count"], 1)
        self.assertAlmostEqual(single["result_mean"], 2.0)
        self.assertTrue(np.isnan(single["result_sd"]))
        self.assertTrue(np.isnan(single["result_cv_percent"]))

    def test_ratio_statistics_use_blank_corrected_replicates(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": ["blank", "blank", "control", "patient", "patient", "missing"],
                "sample_name": ["Blank", "Blank", "Grens", "Patient", "Patient", "Missing"],
                "sample_type": ["Blank", "Blank", "Control", "Unknown", "Unknown", "Unknown"],
                "source_extinction": [0.1, 0.1, 0.6, 1.1, 2.1, np.nan],
            }
        )
        settings = AssaySettings(
            "test.skax", "lot", calculation=FormulaSettings(Formula("sample / grens", {"grens": "Grens"}), {}, response_scale="blank_corrected")
        )
        samples = sample_table(AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze().samples)
        self.assertNotIn("ratio_count", samples.columns)
        self.assertNotIn("result", samples.columns)
        patient = samples.loc[samples["sample_id"].eq("patient")].iloc[0]
        self.assertEqual(patient["count"], 2)
        self.assertAlmostEqual(patient["extinction_mean"], 1.5)
        self.assertAlmostEqual(patient["extinction_sd"], 1 / np.sqrt(2))
        self.assertAlmostEqual(patient["extinction_cv_percent"], 100 / np.sqrt(2) / 1.5)
        self.assertAlmostEqual(patient["result_mean"], 3.0)
        self.assertAlmostEqual(patient["result_sd"], np.sqrt(2))
        self.assertAlmostEqual(patient["result_cv_percent"], 100 * np.sqrt(2) / 3)
        blank = samples.loc[samples["sample_id"].eq("blank")].iloc[0]
        self.assertEqual(blank["count"], 2)
        self.assertAlmostEqual(blank["result_mean"], 0.0)
        self.assertAlmostEqual(blank["result_sd"], 0.0)
        self.assertTrue(np.isnan(blank["result_cv_percent"]))
        missing = samples.loc[samples["sample_id"].eq("missing")].iloc[0]
        self.assertEqual(missing["count"], 0)
        for column in ("result_mean", "result_sd", "result_cv_percent"):
            self.assertTrue(np.isnan(missing[column]), column)

    def test_ratio_cutoff_uses_control_wells_and_returns_calculation_details(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": ["blank", "control", "control", "patient"],
                "sample_name": ["Blank", "Grens", "Grens", "Grens"],
                "sample_type": ["Blank", "Control", "Control", "Unknown"],
                "source_extinction": [0.1, 0.3, 0.5, 1.1],
            }
        )
        settings = AssaySettings(
            "test.skax", "lot", calculation=FormulaSettings(Formula("sample / grens", {"grens": "Grens"}), {}, response_scale="blank_corrected")
        )
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()

        self.assertAlmostEqual(well_table(result.samples).loc[3, "result"], 1 / 0.3)

    def test_missing_standard_measurements_are_rejected(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": ["s0", "s1", "s2", "s3"],
                "sample_name": ["S"] * 4,
                "sample_type": ["Standard"] * 4,
                "sample_concentration": [0.0, 1.0, 2.0, 3.0],
                "source_extinction": [0.1, 0.7, np.nan, 1.6],
            }
        )
        settings = AssaySettings("test.skax", "lot", calculation=CalibrationSettings("4pl", response_scale="raw"))
        with self.assertRaisesRegex(ValueError, 'eindig'):
            AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()

    def test_borrelia_formula_corrects_cutoff_after_blank_subtraction(self) -> None:
        settings = borrelia_settings()
        frame = pd.DataFrame(
            {
                "sample_id": ["blank", "control", "control", "control", "patient"],
                "sample_name": ["Blank", "Kalibratie", "Kalibratie", "Kalibratie", "Kalibratie"],
                "sample_type": ["Blank", "Control", "Control", "Control", "Unknown"],
                "source_extinction": [0.1, 0.39, 0.41, 0.4, 1.1],
            }
        )
        original = frame.copy(deep=True)
        analyzer = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings)
        result = analyzer.analyze()

        self.assertAlmostEqual(well_table(result.samples).loc[4, "result"], 1 / 0.129)
        np.testing.assert_allclose(well_table(result.samples)["extinction"], [0, 0.29, 0.31, 0.3, 1])
        control = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Control")]
        self.assertAlmostEqual(control.iloc[0]["extinction_mean"], 0.3)
        pd.testing.assert_frame_equal(frame, original)
        pd.testing.assert_frame_equal(well_table(analyzer.analyze().samples), well_table(result.samples))

        zero_settings = replace(
            settings, calculation=replace(settings.calculation, formula=Formula("sample / (kalibratie * 0)", {"kalibratie": "Kalibratie"}))
        )
        with self.assertRaisesRegex(ValueError, 'door nul delen'):
            AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), zero_settings).analyze()

    def test_borrelia_control_ratios_use_calibrator_mean_independently_of_cutoff_factor(self) -> None:
        settings = borrelia_settings()
        frame = pd.DataFrame(
            {
                "sample_id": ["blank", "cal", "cal", "cal", "positive", "negative", "patient"],
                "sample_name": ["Blank", "Kalibratie", "Kalibratie", "Kalibratie", "Positief", "Negatief", "Positief"],
                "sample_type": ["Blank", "Control", "Control", "Control", "Control", "Control", "Unknown"],
                "source_extinction": [0.05, 0.44, 0.46, 0.45, 0.65, 0.25, 0.45],
            }
        )
        original = frame.copy(deep=True)
        analyzer = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings)
        result = analyzer.analyze()

        np.testing.assert_allclose(well_table(result.samples)["result"].iloc[4:], [1.5, 0.5, 1 / 0.43])
        np.testing.assert_allclose(well_table(result.samples)["extinction"], [0, 0.39, 0.41, 0.4, 0.6, 0.2, 0.4])
        self.assertEqual(result.validation_status, ValidationStatus.PASSED)
        checks = {check.sample_name: check for check in result.validation_results if check.metric == ValidationMetric.RESULT}
        self.assertAlmostEqual(checks["Positief"].value, 1.5)
        self.assertAlmostEqual(checks["Negatief"].value, 0.5)
        changed_factor = replace(
            settings, calculation=replace(settings.calculation, formula=Formula("sample / (kalibratie * 0.5)", {"kalibratie": "Kalibratie"}))
        )
        changed = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), changed_factor).analyze()

        np.testing.assert_allclose(well_table(changed.samples)["result"].iloc[4:], [1.5, 0.5, 2])
        pd.testing.assert_frame_equal(frame, original)
        pd.testing.assert_frame_equal(well_table(analyzer.analyze().samples), well_table(result.samples))

        # QC fails when a control exceeds its calibrator-relative limit.
        frame.loc[5, "source_extinction"] = 0.45
        failed = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
        negative = next(check for check in failed.validation_results if check.sample_name == "Negatief" and check.metric == ValidationMetric.RESULT)
        self.assertAlmostEqual(negative.value, 1)
        self.assertEqual(negative.status, ValidationStatus.FAILED)

        # This positive control passes the extinction minimum and would pass
        # the cutoff ratio minimum, but fails the calibrator ratio minimum.
        frame.loc[5, "source_extinction"] = 0.25
        frame.loc[[1, 2, 3, 4], "source_extinction"] = [0.64, 0.66, 0.65, 0.57]
        failed = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
        positive = next(check for check in failed.validation_results if check.sample_name == "Positief" and check.metric == ValidationMetric.RESULT)
        self.assertAlmostEqual(positive.value, 0.52 / 0.6)
        self.assertEqual(positive.status, ValidationStatus.FAILED)

    def test_control_overrides_include_their_own_corrections(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": ["cal", "positive", "patient"],
                "sample_name": ["Cal", "Positive", "Patient"],
                "sample_type": ["Control", "Control", "Unknown"],
                "source_extinction": [0.4, 0.2, 0.4],
            }
        )
        settings = AssaySettings(
            "test.skax",
            "lot",
            calculation=FormulaSettings(
                Formula("sample / (cal * 0.5)", {"cal": "Cal"}),
                {
                    "Cal": FormulaOverrideSettings("sample * 2 / (cal * 0.5)", response_scale="blank_corrected"),
                    "Positive": FormulaOverrideSettings("sample * 3 / cal", response_scale="blank_corrected"),
                },
                response_scale="blank_corrected",
            ),
        )
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()

        np.testing.assert_allclose(well_table(result.samples)["extinction"], [0.4, 0.2, 0.4])
        np.testing.assert_allclose(well_table(result.samples)["result"], [4, 1.5, 2])

    def test_ratio_references_reject_missing_nonfinite_zero_and_noncontrol_targets(self) -> None:
        settings = AssaySettings(
            "test.skax",
            "lot",
            calculation=FormulaSettings(
                Formula("sample / cutoff", {"cutoff": "Cutoff", "reference": "Reference"}),
                {"Positive": FormulaOverrideSettings("sample / reference", response_scale="blank_corrected")},
                response_scale="blank_corrected",
            ),
            controls={"Positive": SampleSettings()},
        )
        for reference_type, reference_values, message in (
            ("Unknown", [0.4, 0.4], 'Geen metingen'),
            ("Control", [0.4, np.nan], 'moeten eindig zijn'),
            ("Control", [0, 0], 'door nul delen'),
        ):
            frame = pd.DataFrame(
                {
                    "sample_id": ["cutoff", "positive", "ref", "ref"],
                    "sample_name": ["Cutoff", "Positive", "Reference", "Reference"],
                    "sample_type": ["Control", "Control", reference_type, reference_type],
                    "source_extinction": [0.5, 0.8, *reference_values],
                }
            )
            with self.subTest(reference=reference_values), self.assertRaisesRegex(ValueError, message):
                AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()

    def test_control_overrides_leave_reference_means_unchanged(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": ["cal", "reference", "patient"],
                "sample_name": ["Cal", "Reference", "Patient"],
                "sample_type": ["Control", "Control", "Unknown"],
                "source_extinction": [0.4, 0.8, 0.6],
            }
        )
        settings = AssaySettings(
            "test.skax",
            "lot",
            calculation=FormulaSettings(
                Formula("sample / ((cal * 0.5) / reference)", {"cal": "Cal", "reference": "Reference"}),
                {"Reference": FormulaOverrideSettings("sample * 10 / ((cal * 0.5) / reference)", response_scale="blank_corrected")},
                response_scale="blank_corrected",
            ),
        )
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()

        self.assertAlmostEqual(well_table(result.samples).iloc[-1]["result"], 2.4)
        self.assertAlmostEqual(well_table(result.samples).iloc[1]["result"], 32)

    def test_igg_preprocessing_uses_original_control_mean_for_configured_curve(self) -> None:
        registry = AssayRegistry(Path(__file__).resolve().parents[1] / "configs")
        settings = replace(registry.get("Qkoorts_IgG.skax"), measurement=MeasurementSettings())
        frame = pd.DataFrame(
            {
                "sample_id": ["blank", "control", "control", "patient", "same-name"],
                "sample_name": ["Blank", "Standaard", "Standaard", "Patient", "Standaard"],
                "sample_type": ["Blank", "Control", "Control", "Unknown", "Unknown"],
                "source_extinction": [0.1, 0.7, 1.1, 0.5, 5.1],
            }
        )
        original = frame.copy(deep=True)
        analyzer = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings)
        result = analyzer.analyze()
        parameters = settings.calculation.parameters
        # The patient's corrected curve input is (0.5 - 0.1) * 0.86 / 0.8 = 0.43.
        expected = parameters.C * ((0.43 - parameters.A) / (parameters.D - 0.43)) ** (1 / parameters.B)
        self.assertAlmostEqual(well_table(result.samples).loc[3, "result"], expected)
        np.testing.assert_allclose(well_table(result.samples)["extinction"], [0, 0.6, 1, 0.4, 5])
        control = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Control")]
        self.assertAlmostEqual(control.iloc[0]["extinction_mean"], 0.8)
        pd.testing.assert_frame_equal(frame, original)
        pd.testing.assert_frame_equal(well_table(analyzer.analyze().samples), well_table(result.samples))

    def test_control_overrides_use_original_reference_means_independently_of_order(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": ["blank", "cutoff", "positive", "positive", "same-name"],
                "sample_name": ["Blank", "Grens", "Positive", "Positive", "Positive"],
                "sample_type": ["Blank", "Control", "Control", "Control", "Unknown"],
                "source_extinction": [0.1, 0.6, 0.3, 0.5, 0.8],
            }
        )
        original = frame.copy(deep=True)
        settings = AssaySettings.from_dict(
            "test.skax",
            {
                "lot_number": "lot",
                "controls": {
                    "Grens": {},
                    "Positive": {},
                },
                "validation": {"result": {"controls": {"Positive": {"result": {"min": 24, "max": 24.1}}}}},
                "calculation": {
                    "response_scale": "blank_corrected",
                    "method": "formula", "expression": "sample / (grens * 0.5)", "references": {"grens": "Grens"},
                    "control_overrides": {"Grens": {"expression": "sample / grens", "response_scale": "blank_corrected"}, "Positive": {"expression": "sample * 10 / (grens * grens * 0.5)", "response_scale": "blank_corrected"}},
                },
                "result": {"label": "Ratio"},
            },
        )
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()

        self.assertNotIn("preprocessed_extinction", well_table(result.samples))
        np.testing.assert_allclose(well_table(result.samples)["result"], [0, 1, 16, 32, 2.8])
        positive = sample_table(result.samples).loc[sample_table(result.samples)["sample_id"].eq("positive")].iloc[0]
        self.assertAlmostEqual(positive["extinction_mean"], 0.3)
        self.assertAlmostEqual(positive["result_mean"], 24)
        self.assertEqual(result.validation_status, ValidationStatus.PASSED)
        reversed_settings = replace(settings, calculation=replace(
            settings.calculation, control_overrides=dict(reversed(list(settings.calculation.control_overrides.items()))),
        ))
        reordered = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), reversed_settings).analyze()
        pd.testing.assert_frame_equal(well_table(reordered.samples), well_table(result.samples))
        pd.testing.assert_frame_equal(frame, original)

    def test_shared_curve_preprocessing_applies_once_to_controls_and_patients(self) -> None:
        responses = FourPLModel().calculate_extinctions([1, 4], self.parameters)
        frame = pd.DataFrame(
            {
                "sample_id": ["blank", "positive", "negative", "patient"],
                "sample_name": ["Blank", "Positive", "Negative", "Patient"],
                "sample_type": ["Blank", "Control", "Control", "Unknown"],
                "source_extinction": [0.05, responses[0] / 2 + 0.05, responses[1] / 2 + 0.05, responses[0] / 2 + 0.05],
            }
        )
        settings = AssaySettings(
            "test.skax",
            "lot",
            calculation=CalibrationSettings("4pl", self.parameters, preprocessing=Formula("(sample * 2)", {}), response_scale="blank_corrected"),
        )
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
        np.testing.assert_allclose(well_table(result.samples)["result"].iloc[1:], [1, 4, 1])

    def test_patient_concentrations_are_dilution_corrected_and_grouped_by_factor(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": ["patient"] * 4 + ["control"],
                "sample_name": ["Patient"] * 4 + ["Positive"],
                "sample_type": ["Unknown"] * 4 + ["Control"],
                "sample_dilution": [2, 2, 10, 10, 0],
                "source_extinction": FourPLModel().calculate_extinctions([1, 4, 1, 4, 2], self.parameters),
                "well": ["A1", "A2", "A3", "A4", "A5"],
            }
        )
        original = frame.copy(deep=True)
        settings = AssaySettings(
            "test.skax", "lot", calculation=CalibrationSettings("4pl", self.parameters, concentration_range=(1, 4), response_scale="blank_corrected")
        )
        analyzer = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings)
        result = analyzer.analyze()
        np.testing.assert_allclose(well_table(result.samples)["result"], [2, 8, 10, 40, 2])
        patients = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")]
        self.assertEqual(patients["sample_dilution"].tolist(), [2, 10])
        self.assertEqual(patients["count"].tolist(), [2, 2])
        self.assertEqual(patients["used_wells"].tolist(), ["A1, A2", "A3, A4"])
        np.testing.assert_allclose(patients["result_mean"], [5, 25])
        np.testing.assert_allclose(patients["result_sd"], [6 / np.sqrt(2), 30 / np.sqrt(2)])
        pd.testing.assert_frame_equal(frame, original)
        pd.testing.assert_frame_equal(well_table(analyzer.analyze().samples), well_table(result.samples))

        # Dilution cannot rescue a reading outside the curve's reporting range.
        frame.loc[0, "source_extinction"] = FourPLModel().calculate_extinctions(0.5, self.parameters)
        incomplete = sample_table(AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze().samples).iloc[0]
        self.assertTrue(np.isnan(incomplete["result_mean"]))
        self.assertTrue(np.isnan(incomplete["result_sd"]))

    def test_concentration_reporting_bounds_support_both_curve_directions_and_asymptote_orders(self) -> None:
        for a, d in ((2.2, 0.2), (0.2, 2.2)):
            for b in (-1.3, 1.3):
                for sample_type in ("Unknown", "Control"):
                    with self.subTest(a=a, d=d, b=b, sample_type=sample_type):
                        parameters = FourPLParameters(a, b, 2, d)
                        model = FourPLModel()
                        responses = list(model.calculate_extinctions([0.5, 1, 2, 4, 8], parameters))
                        increasing = responses[0] < responses[-1]
                        extremes = [0.1, 2.3] if increasing else [2.3, 0.1]
                        frame = pd.DataFrame(
                            {
                                "sample_id": list(range(8)),
                                "sample_name": [f"Sample {index}" for index in range(8)],
                                "sample_type": sample_type,
                                "source_extinction": [*responses, *extremes, np.nan],
                            }
                        )
                        settings = AssaySettings(
                            "test.skax",
                            "lot",
                            calculation=CalibrationSettings("4pl", parameters, concentration_range=(1, 4), response_scale="blank_corrected"),
                        )
                        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
                        expected = ["below_range", "in_range", "in_range", "in_range", "above_range", "below_range", "above_range", "unavailable"]
                        self.assertEqual(well_table(result.samples)["result_status"].tolist(), expected)
                        self.assertEqual(sample_table(result.samples)["result_status"].tolist(), expected)
                        np.testing.assert_allclose(well_table(result.samples)["result"], [np.nan, 1, 2, 4, np.nan, np.nan, np.nan, np.nan])
                        np.testing.assert_allclose(sample_table(result.samples)["result_limit"], [1, np.nan, np.nan, np.nan, 4, 1, 4, np.nan])

    def test_bounded_replicates_require_matching_statuses_and_receive_patient_dilution(self) -> None:
        groups = {
            "Below": [0.5, 0.75],
            "Above": [5, 8],
            "Inside": [1, 4],
            "Mixed": [0.5, 2],
            "Missing": [0.5, np.nan],
            "Opposite": [0.5, 8],
        }
        names = [name for name in groups for _ in groups[name]]
        concentrations = [value for values in groups.values() for value in values]
        frame = pd.DataFrame(
            {
                "sample_id": names,
                "sample_name": names,
                "sample_type": "Unknown",
                "sample_dilution": 2.5,
                "source_extinction": FourPLModel().calculate_extinctions(concentrations, self.parameters),
            }
        )
        settings = AssaySettings(
            "test.skax", "lot", calculation=CalibrationSettings("4pl", self.parameters, concentration_range=(1, 4), response_scale="blank_corrected")
        )
        analyzer = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings)
        result = analyzer.analyze()
        samples = sample_table(result.samples).set_index("sample_name")
        self.assertEqual(samples["result_status"].tolist(), ["below_range", "above_range", "in_range", "unavailable", "unavailable", "unavailable"])
        np.testing.assert_allclose(samples["result_limit"], [2.5, 10, np.nan, np.nan, np.nan, np.nan])
        self.assertAlmostEqual(samples.loc["Inside", "result_mean"], 6.25)
        for name in ("Below", "Above", "Mixed", "Missing", "Opposite"):
            for metric in ("mean", "sd", "cv_percent"):
                self.assertTrue(np.isnan(samples.loc[name, f"result_{metric}"]))
        pd.testing.assert_frame_equal(sample_table(analyzer.analyze().samples), sample_table(result.samples))

    def test_zero_lower_reporting_limit_does_not_report_negative_concentration(self) -> None:
        frame = pd.DataFrame({"sample_id": [1], "sample_name": ["Control"], "sample_type": "Control", "source_extinction": [0.1]})
        settings = AssaySettings(
            "test.skax", "lot", calculation=CalibrationSettings("4pl", self.parameters, concentration_range=(0, 4), response_scale="blank_corrected")
        )
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
        self.assertEqual(sample_table(result.samples).iloc[0]["result_status"], ResultStatus.UNAVAILABLE.value)
        self.assertTrue(np.isnan(sample_table(result.samples).iloc[0]["result_limit"]))

    def test_dilution_corrected_reporting_bound_overflow_is_rejected(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": [1],
                "sample_name": ["Patient"],
                "sample_type": "Unknown",
                "sample_dilution": [1e308],
                "source_extinction": FourPLModel().calculate_extinctions([8], self.parameters),
            }
        )
        settings = AssaySettings(
            "test.skax", "lot", calculation=CalibrationSettings("4pl", self.parameters, concentration_range=(1, 4), response_scale="blank_corrected")
        )
        with self.assertRaisesRegex(ValueError, 'Voor verdunning gecorrigeerde concentratiegrenzen moeten eindig zijn'):
            AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()

    def test_dilution_factors_are_validated_before_filtering_empty_wells(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": ["cutoff", "patient"],
                "sample_name": ["Grens", "Patient"],
                "sample_type": ["Control", "Unknown"],
                "source_extinction": [0.5, 0],
                "sample_dilution": [0, 1],
            }
        )
        concentration = AssaySettings("test.skax", "lot", calculation=CalibrationSettings("4pl", self.parameters, response_scale="blank_corrected"))
        ratio = AssaySettings(
            "test.skax", "lot", calculation=FormulaSettings(Formula("sample / grens", {"grens": "Grens"}), {}, response_scale="blank_corrected")
        )
        for factor in (0, -1, 0.5, np.nan, np.inf, -np.inf, "bad", True, None):
            frame["sample_dilution"] = pd.Series([0, factor], dtype=object)
            with self.subTest(factor=factor), self.assertRaises(ValueError):
                AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), concentration).analyze()
        frame["sample_dilution"] = [0, 2]
        with self.assertRaisesRegex(ValueError, 'assayspecifieke verdunningscorrectie'):
            AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), ratio).analyze()
        frame["sample_dilution"] = [2, 1]
        with self.assertRaisesRegex(ValueError, 'alleen ondersteund voor patiëntmonsters'):
            AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), concentration).analyze()

    def test_dilution_correction_rejects_numeric_overflow(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": ["patient"],
                "sample_name": ["Patient"],
                "sample_type": ["Unknown"],
                "sample_dilution": [1e308],
                "source_extinction": [FourPLModel().calculate_extinctions(4, self.parameters)],
            }
        )
        settings = AssaySettings("test.skax", "lot", calculation=CalibrationSettings("4pl", self.parameters, response_scale="blank_corrected"))
        with self.assertRaisesRegex(ValueError, 'Voor verdunning gecorrigeerde concentraties moeten eindig zijn'):
            AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()

    def test_fitted_controls_and_patient_dilution_leave_calibration_intact(self) -> None:
        nominal = np.repeat([0, 1, 5, 25, 100], 2)
        responses = FourPLModel().calculate_extinctions(nominal, self.parameters)
        control_response = FourPLModel().calculate_extinctions(2, self.parameters)
        patient_response = FourPLModel().calculate_extinctions(3, self.parameters)
        frame = pd.DataFrame(
            {
                "sample_id": [f"standard-{x}" for x in nominal] + ["blank", "positive", "patient"],
                "sample_name": [f"Standard {x}" for x in nominal] + ["Blank", "Positive", "Patient"],
                "sample_type": ["Standard"] * len(nominal) + ["Blank", "Control", "Unknown"],
                "sample_concentration": [*nominal, np.nan, np.nan, np.nan],
                "sample_dilution": [0] * len(nominal) + [0, 0, 10],
                "source_extinction": [*(responses + 0.05), 0.05, control_response + 0.05, patient_response + 0.05],
            }
        )
        settings = AssaySettings.from_dict(
            "test.skax",
            {
                "lot_number": "lot",
                "controls": {"Positive": {}},
                "validation": {"result": {"controls": {"Positive": {"result": {"min": 1.9, "max": 2.1}}}}},
                "calculation": {"response_scale": "raw", "method": "calibration", "fit_model": "4pl"},
                "result": {"label": "Conc"},
            },
        )
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
        self.assertEqual(result.validation_status, ValidationStatus.PASSED)
        np.testing.assert_allclose(well_table(result.samples)["result"].iloc[-2:], [2, 30], rtol=1e-6)
        diagnostics = result.calculation_details.curve.diagnostics
        np.testing.assert_allclose(diagnostics.predictions, responses + 0.05, rtol=1e-6)
        np.testing.assert_allclose(diagnostics.back_calculated, nominal, rtol=1e-6, atol=1e-9)
        standards = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Standard")]
        np.testing.assert_allclose(standards.iloc[1:]["recovery_percent"], 100, rtol=1e-6)

    def test_flagged_low_readings_are_excluded_instead_of_treated_as_empty_wells(self) -> None:
        settings = AssaySettings(
            "test.skax", "lot", calculation=FormulaSettings(Formula("sample / grens", {"grens": "Grens"}), {}, response_scale="blank_corrected")
        )
        for name in ("saturated", "disabled"):
            frame = pd.DataFrame(
                {
                    "sample_id": ["cutoff", "patient"],
                    "sample_name": ["Grens", "Patient"],
                    "sample_type": ["Control", "Unknown"],
                    "source_extinction": [0.5, 0],
                    "well": ["A1", "A2"],
                    name: [False, True],
                }
            )
            with self.subTest(flag=name):
                result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
                self.assertEqual(well_table(result.samples)["well"].tolist(), ["A1", "A2"])
                self.assertFalse(well_table(result.samples).iloc[1]["included_in_analysis"])
                self.assertTrue(np.isnan(well_table(result.samples).iloc[1]["extinction"]))
                patients = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")]
                if name == "disabled":
                    self.assertTrue(patients.empty)
                else:
                    self.assertEqual(patients.iloc[0]["count"], 0)
                    self.assertEqual(patients.iloc[0]["used_wells"], "")
                    self.assertTrue(np.isnan(patients.iloc[0]["result_mean"]))
            frame[name] = [False, "false"]
            with self.subTest(malformed_flag=name), self.assertRaisesRegex(ValueError, 'waar of onwaar'):
                AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()

    def test_instrument_exclusions_use_healthy_replicates_for_blank_cutoff_qc_and_patients(self) -> None:
        settings = AssaySettings.from_dict(
            "test.skax",
            {
                "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / grens", "references": {"grens": "Grens"}},
                "lot_number": "lot",
                "blank": {},
                "controls": {
                    name: {"replicate_filter": {"response_scale": "blank_corrected", "max_deviation_percent": 15, "max_exclusions": 0}}
                    for name, minimum, maximum in (("Grens", 0.4, 0.6), ("Positive", 0.6, 0.8))
                },
                "validation": {
                    "result": {
                        "controls": {
                            name: {"extinction": {"response_scale": "blank_corrected", "min": minimum, "max": maximum}}
                            for name, minimum, maximum in (("Grens", 0.4, 0.6), ("Positive", 0.6, 0.8))
                        },
                        "blank": {"extinction": {"response_scale": "raw", "max": 0.2}},
                    }
                },
            },
        )
        for flag in ("disabled", "saturated"):
            for separate_ids in (False, True):
                with self.subTest(flag=flag, separate_ids=separate_ids):
                    frame = pd.DataFrame(
                        {
                            "sample_id": ["blank", "blank", "cutoff", "cutoff", "positive", "positive", "patient", "patient", "flagged", "other"],
                            "sample_name": ["Blank", "Blank", "Grens", "Grens", "Positive", "Positive", "Patient", "Patient", "Flagged", "Other"],
                            "sample_type": ["Blank"] * 2 + ["Control"] * 4 + ["Unknown"] * 4,
                            "source_extinction": [0.1, 6, 0.6, 6, 0.8, 6, 1.1, 6, 6, 1.6],
                            "well": [f"A{index + 1}" for index in range(10)],
                            flag: [False, True, False, True, False, True, False, True, True, False],
                        }
                    )
                    if separate_ids:
                        frame.loc[[1, 3, 5], "sample_id"] += "-flagged"
                    original = frame.copy(deep=True)
                    analyzer = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings)
                    result = analyzer.analyze()
                    self.assertAlmostEqual(result.blank_value, 0.1)

                    self.assertEqual(result.validation_status, ValidationStatus.PASSED)
                    samples = sample_table(result.samples).set_index("sample_id")
                    self.assertEqual(samples.loc["patient", "count"], 1)
                    self.assertEqual(samples.loc["patient", "used_wells"], "A7")
                    self.assertAlmostEqual(samples.loc["patient", "result_mean"], 2)
                    self.assertAlmostEqual(samples.loc["other", "result_mean"], 3)
                    positive = next(
                        check for check in result.validation_results
                        if check.sample_key == (SampleType.CONTROL, "Positive") and check.metric == ValidationMetric.EXTINCTION
                    )
                    self.assertAlmostEqual(positive.value, 0.7)
                    self.assertIsNone(positive.cv_percent)
                    excluded = well_table(result.samples).loc[frame[flag]]
                    self.assertFalse(excluded["included_in_analysis"].any())
                    self.assertTrue(excluded[["extinction", "result"]].isna().all().all())
                    reason = "uitgeschakeld" if flag == "disabled" else "verzadigd"
                    self.assertTrue(excluded["exclusion_reason"].str.contains(reason).all())
                    if flag == "saturated":
                        self.assertEqual(samples.loc["flagged", "count"], 0)
                        self.assertTrue(np.isnan(samples.loc["flagged", "result_mean"]))
                    else:
                        self.assertNotIn("flagged", samples.index)
                    pd.testing.assert_frame_equal(frame, original)
                    pd.testing.assert_frame_equal(well_table(analyzer.analyze().samples), well_table(result.samples))

    def test_flagged_standards_are_excluded_before_fitting(self) -> None:
        nominal = np.repeat([0.25, 0.75, 2, 6, 20], 2)
        responses = FourPLModel().calculate_extinctions(nominal, self.parameters)
        settings = AssaySettings("test.skax", "lot", calculation=CalibrationSettings("4pl", response_scale="raw"))
        for flag in ("disabled", "saturated"):
            with self.subTest(flag=flag):
                frame = pd.DataFrame(
                    {
                        "sample_id": [f"standard-{value}" for value in nominal],
                        "sample_name": [f"Standard {value}" for value in nominal],
                        "sample_type": "Standard",
                        "sample_concentration": nominal,
                        "source_extinction": responses,
                        flag: [True] + [False] * 9,
                    }
                )
                frame.loc[0, "source_extinction"] = 6
                result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
                self.assertEqual(len(result.calculation_details.curve.diagnostics.predictions), 9)
                np.testing.assert_allclose(result.calculation_details.curve.diagnostics.predictions, responses[1:], rtol=1e-6)
                self.assertTrue(np.isnan(well_table(result.samples).loc[0, "result"]))
                self.assertTrue(np.isnan(well_table(result.samples).loc[0, "back_calculated_concentration"]))
                self.assertEqual(sample_table(result.samples).iloc[0]["count"], 1)

    def test_unbounded_curve_cannot_report_a_surviving_duplicate_when_another_is_saturated(self) -> None:
        settings = AssaySettings("test.skax", "lot", calculation=CalibrationSettings("4pl", self.parameters, response_scale="blank_corrected"))
        frame = pd.DataFrame(
            {
                "sample_id": ["paired", "paired", "saturated", "disabled", "both"],
                "sample_name": ["Paired", "Paired", "Saturated", "Disabled", "Both"],
                "sample_type": "Unknown",
                "sample_dilution": 10,
                "source_extinction": [FourPLModel().calculate_extinctions(3, self.parameters), 6, 6, 6, 6],
                "saturated": [False, True, True, False, True],
                "disabled": [False, False, False, True, True],
            }
        )
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
        samples = sample_table(result.samples).set_index("sample_id")
        self.assertEqual(list(samples.index), ["paired", "saturated"])
        self.assertEqual(samples.loc["paired", "count"], 1)
        self.assertTrue(np.isnan(samples.loc["paired", "result_mean"]))
        self.assertEqual(samples.loc["paired", "result_status"], ResultStatus.UNAVAILABLE.value)
        self.assertEqual(samples.loc["saturated", "count"], 0)
        self.assertTrue(np.isnan(samples.loc["saturated", "result_mean"]))
        self.assertEqual(samples.loc["saturated", "result_status"], ResultStatus.UNAVAILABLE.value)
        self.assertEqual(samples.loc["saturated", "sample_dilution"], 10)

    def test_fully_flagged_optional_controls_are_unavailable_and_required_inputs_still_fail(self) -> None:
        settings = AssaySettings.from_dict(
            "test.skax",
            {
                "lot_number": "lot",
                "controls": {"Optional": {"replicate_filter": {"response_scale": "blank_corrected", "max_deviation_percent": 15}}},
                "validation": {"result": {"controls": {"Optional": {"extinction": {"response_scale": "blank_corrected", "min": 0.1}}}}},
                "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / grens", "references": {"grens": "Grens"}},
                "result": {"label": "Ratio"},
            },
        )
        for flag in ("disabled", "saturated"):
            with self.subTest(flag=flag):
                frame = pd.DataFrame(
                    {
                        "sample_id": ["blank", "cutoff", "optional", "patient"],
                        "sample_name": ["Blank", "Grens", "Optional", "Patient"],
                        "sample_type": ["Blank", "Control", "Control", "Unknown"],
                        "source_extinction": [0.1, 0.6, 6, 1.1],
                        flag: [False, False, True, False],
                    }
                )
                result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
                self.assertEqual(result.validation_status, ValidationStatus.NOT_EVALUATED)
                optional = next(check for check in result.validation_results if check.sample_key == (SampleType.CONTROL, "Optional"))
                self.assertIsNone(optional.value)
                self.assertEqual(optional.status, ValidationStatus.NOT_EVALUATED)
                self.assertAlmostEqual(sample_table(result.samples).iloc[-1]["result_mean"], 2)
                for index, message in ((0, 'geen bruikbare blancometingen'), (1, "Geen metingen gevonden voor controle 'Grens'")):
                    with self.subTest(required=index):
                        frame.loc[index, flag] = True
                        with self.assertRaisesRegex(ValueError, message):
                            AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
                        frame.loc[index, flag] = False

    def test_fixed_curve_preprocessing_corrects_inputs_before_concentration(self) -> None:
        concentrations = np.repeat([0, 0.25, 0.75, 2, 6, 20], 2)
        responses = FourPLModel().calculate_extinctions(concentrations, self.parameters)
        frame = pd.DataFrame(
            {
                "sample_id": [f"standard-{x}" for x in concentrations],
                "sample_name": [f"Standard {x}" for x in concentrations],
                "sample_type": "Standard",
                "sample_concentration": concentrations,
                "source_extinction": responses * 0.8 / 0.86 + 0.05,
            }
        )
        frame = pd.concat(
            [
                frame,
                pd.DataFrame(
                    {
                        "sample_id": ["blank", "control", "patient-a", "patient-b"],
                        "sample_name": ["Blank", "Standaard", "Patient A", "Patient B"],
                        "sample_type": ["Blank", "Control", "Unknown", "Unknown"],
                        "sample_concentration": [np.nan] * 4,
                        "source_extinction": [
                            0.05,
                            0.85,
                            *(FourPLModel().calculate_extinctions([1, 4], self.parameters) * 0.8 / 0.86 + 0.05),
                        ],
                    }
                ),
            ],
            ignore_index=True,
        )
        settings = AssaySettings(
            "test.skax",
            "lot",
            controls={"Standaard": SampleSettings()},
            calculation=CalibrationSettings(
                "4pl", self.parameters,
                preprocessing=Formula("((sample * 0.86) / standaard)", {"standaard": "Standaard"}), response_scale="blank_corrected",
            ),
        )
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
        patients = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")]
        np.testing.assert_allclose(patients["result_mean"], [1, 4], rtol=1e-6)
        standards = well_table(result.samples).loc[well_table(result.samples)["sample_type"].eq("Standard")]
        np.testing.assert_allclose(standards["curve_extinction"], responses, rtol=1e-6)
        self.assertIsNone(result.calculation_details.curve.diagnostics)

    def test_preprocessing_rejects_missing_incomplete_and_zero_control_measurements(self) -> None:
        registry = AssayRegistry(Path(__file__).resolve().parents[1] / "configs")
        settings = replace(registry.get("Qkoorts_IgG.skax"), measurement=MeasurementSettings())
        for control_types, control_values, error in (
            (["Unknown", "Unknown"], [0.8, 0.8], 'Geen metingen'),
            (["Control", "Control"], [0.8, np.nan], 'moeten eindig zijn'),
            (["Control", "Control"], [0, 0], 'door nul delen'),
        ):
            frame = pd.DataFrame(
                {
                    "sample_id": ["blank", "control-a", "control-b", "patient"],
                    "sample_name": ["Blank", "Standaard", "Standaard", "Patient"],
                    "sample_type": ["Blank", *control_types, "Unknown"],
                    "source_extinction": [0.05, *(value + 0.05 for value in control_values), 0.45],
                }
            )
            with self.subTest(values=control_values), self.assertRaisesRegex(ValueError, error):
                AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()

    def test_fitted_sample_summaries_preserve_input_sds_and_diagnostic_recovery(self) -> None:
        concentrations = np.repeat([0, 1, 5, 25, 100], 2)
        responses = FourPLModel().calculate_extinctions(concentrations, self.parameters)
        responses += np.tile([-0.002, 0.002], 5)
        frame = pd.DataFrame(
            {
                "sample_id": [f"standard-{x}" for x in concentrations] + ["blank"],
                "sample_name": [f"Standard {x}" for x in concentrations] + ["Blank"],
                "sample_type": ["Standard"] * len(concentrations) + ["Blank"],
                "sample_concentration": [*concentrations, np.nan],
                "source_extinction": [*(responses + 0.05), 0.05],
            }
        )
        original = frame.copy(deep=True)
        settings = AssaySettings(
            "test.skax", "lot", calculation=CalibrationSettings("4pl", response_scale="raw")
        )
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.metadata), settings).analyze()
        standards = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Standard")]
        np.testing.assert_allclose(standards["curve_extinction_mean"], responses.reshape(-1, 2).mean(axis=1) + 0.05)
        np.testing.assert_allclose(standards["curve_extinction_sd"], np.sqrt(2) * 0.002)
        np.testing.assert_allclose(standards["extinction_sd"], np.sqrt(2) * 0.002)
        diagnostic_values = result.calculation_details.curve.diagnostics.back_calculated.reshape(-1, 2)
        for index, (_, sample) in enumerate(standards.iterrows()):
            values = diagnostic_values[index]
            if np.isfinite(values).all():
                self.assertAlmostEqual(sample["back_calculated_concentration_mean"], values.mean())
                self.assertAlmostEqual(sample["back_calculated_concentration_sd"], values.std(ddof=1))
            else:
                self.assertTrue(np.isnan(sample["back_calculated_concentration_mean"]))
        self.assertTrue(np.isnan(standards.iloc[0]["recovery_percent"]))
        np.testing.assert_allclose(
            standards.iloc[1:]["recovery_percent"],
            100 * diagnostic_values[1:].mean(axis=1) / [1, 5, 25, 100],
        )
        # A low-end replicate is outside the reporting range, but its fit
        # diagnostic remains available for the standard's calibration recovery.
        self.assertTrue(np.isnan(standards.iloc[1]["result_mean"]))
        self.assertTrue(np.isfinite(standards.iloc[1]["back_calculated_concentration_mean"]))
        pd.testing.assert_frame_equal(frame, original)


if __name__ == "__main__":
    unittest.main()
