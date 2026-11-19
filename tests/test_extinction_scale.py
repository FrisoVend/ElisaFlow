"""Calculation input semantics are recorded by analysis, independently of display."""

import unittest
from dataclasses import replace

import numpy as np
import test_calibration_response as response
from fixtures import assay_metadata
from test_multi_wavelength import dual_data, dual_settings, saturate

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.calibration.logistic import FourPLModel, FourPLParameters
from elisaflow.configuration.assay_settings import AssaySettings, CalibrationSettings, FormulaOverrideSettings, FormulaSettings
from elisaflow.formatting import format_extinction_scale
from elisaflow.formulas.expressions import Formula
from elisaflow.models import AssayData, ExtinctionScale, ResponseScale, Sample, SampleType, WellMeasurement


class ExtinctionScaleTests(unittest.TestCase):
    def test_formula_records_the_applied_blank_including_zero_and_absence(self):
        patient = Sample("p", "Patient", SampleType.UNKNOWN, tuple(
            WellMeasurement(source_extinction=value, wavelength_nm=450, well=f"A{index}")
            for index, value in enumerate((1.2, 1.8), 1)
        ))
        settings = AssaySettings("test.skax", "lot", calculation=FormulaSettings(Formula("sample"), response_scale="blank_corrected"))
        for blank in (None, 0.0, 0.25):
            samples = (patient,) if blank is None else (
                Sample("b", "Blank", SampleType.BLANK, (WellMeasurement(source_extinction=blank, wavelength_nm=450, well="B1"),)),
                patient,
            )
            with self.subTest(blank=blank):
                result = AssayAnalyzer(AssayData(samples, assay_metadata()), settings).analyze()
                selected = next(sample for sample in result.samples if sample.sample_type == SampleType.UNKNOWN)
                value = result.patient_extinction(selected)
                self.assertEqual(value.scale.stage, "formula_input")
                self.assertEqual(value.scale.wavelength_nm, 450)
                self.assertEqual(value.scale.blank_subtracted, blank)
                self.assertIsNone(value.scale.preprocessing_expression)
                self.assertIsNone(selected.selected_wavelength_nm)
                self.assertAlmostEqual(value.statistics.mean, 1.5 - (blank or 0))
                self.assertEqual(format_extinction_scale(value.scale), "Ext" if blank is None else "ΔExt")
                self.assertTrue(all(well.curve_extinction is None for well in selected.measurements))

    def test_fixed_and_fitted_curves_record_actual_inputs(self):
        parameters = FourPLParameters(2.5, -1.3, 2, 0.1)
        data = response.CalibrationResponseTests.data(FourPLModel(), parameters)
        for fixed, scale in ((False, "raw"), (True, "raw"), (True, "blank_corrected")):
            calculation = CalibrationSettings(
                "4pl", parameters if fixed else None, response_scale=scale, preprocessing=Formula("sample * 2") if fixed else None,
            )
            settings = AssaySettings("test.skax", "lot", calculation=calculation)
            with self.subTest(fixed=fixed, scale=scale):
                result = AssayAnalyzer(data, settings).analyze()
                patient = next(sample for sample in result.samples if sample.sample_type == SampleType.UNKNOWN)
                value = result.patient_extinction(patient)
                blank = None if scale == "raw" else 0.1
                expected = [(2 if fixed else 1) * (well.source_extinction - (blank or 0)) for well in patient.retained]
                self.assertEqual(value.scale.stage, "curve_input")
                self.assertEqual(value.scale.blank_subtracted, blank)
                self.assertEqual(value.scale.preprocessing_expression, calculation.preprocessing.expression if fixed else None)
                self.assertAlmostEqual(value.statistics.mean, np.mean(expected))
                # QC still sees the independently blank-corrected measurements.
                self.assertAlmostEqual(
                    ResponseScale.BLANK_CORRECTED.statistics(patient).mean,
                    np.mean([well.source_extinction for well in patient.retained]) - 0.1,
                )

    def test_selected_scale_and_statistics_follow_the_selected_current_sample(self):
        config = dual_settings().to_dict()
        config["calculation"]["wavelengths"]["405"] = {
            "fit_model": "4pl", "response_scale": "blank_corrected",
            "parameters": {"A": 1.395, "B": -1.2, "C": 30, "D": 0.020},
        }
        data = saturate(dual_data(), 450, "p1", all_replicates=True)
        result = AssayAnalyzer(data, AssaySettings.from_dict("test.skax", config)).analyze()
        patient = next(sample for sample in result.samples if sample.sample_id == "p1")
        value = result.patient_extinction(patient)
        self.assertEqual(value.scale.wavelength_nm, 405)
        self.assertAlmostEqual(value.scale.blank_subtracted, 0.005)
        changed = replace(patient, name="Renamed patient", measurements=tuple(
            replace(well, curve_extinction=2 * well.curve_extinction) for well in patient.measurements
        ))
        self.assertAlmostEqual(result.patient_extinction(changed).statistics.mean, 2 * value.statistics.mean)
        # A result query describes executed analysis, even if a caller replaces settings.
        config["calculation"]["wavelengths"]["405"]["response_scale"] = "raw"
        changed_result = replace(result, settings=AssaySettings.from_dict("test.skax", config))
        self.assertEqual(changed_result.patient_extinction(patient), value)

    def test_missing_actual_inputs_remain_unavailable_despite_finite_extinction(self):
        result = AssayAnalyzer(dual_data(), dual_settings()).analyze()
        patient = next(sample for sample in result.samples if sample.sample_id == "p1")
        missing = replace(patient, measurements=tuple(replace(well, curve_extinction=None) for well in patient.measurements))
        self.assertTrue(np.isfinite(ResponseScale.BLANK_CORRECTED.statistics(missing).mean))
        self.assertIsNone(result.patient_extinction(missing))
        self.assertIsNone(result.patient_extinction(replace(patient, selected_wavelength_nm=None)))

    def test_failed_curve_has_no_calculation_input_scale(self):
        data = dual_data()
        for identity in ("std-30", "std-100", "std-300"):
            data = saturate(data, 450, identity, all_replicates=True)
        result = AssayAnalyzer(data, dual_settings()).analyze()
        failed = result.measurement_results[450]
        self.assertIsNotNone(failed.calculation_details.calibration_error)
        self.assertIsNone(failed.calculation_input_scale)
        patient = next(sample for sample in failed.samples if sample.sample_id == "p1")
        self.assertTrue(np.isfinite(ResponseScale.BLANK_CORRECTED.statistics(patient).mean))
        self.assertIsNone(failed.calculation_extinction(patient))

    def test_control_keeps_executed_scale_when_renamed_or_formula_details_are_replaced(self):
        settings = AssaySettings("test.skax", "lot", calculation=FormulaSettings(
            Formula("sample"), response_scale="blank_corrected",
            control_overrides={"QC": FormulaOverrideSettings(response_scale="raw")},
        ))
        data = AssayData((
            Sample("b", "Blank", SampleType.BLANK, (WellMeasurement(well="C1", wavelength_nm=450, source_extinction=0.1),)),
            Sample("c", "QC", SampleType.CONTROL, (
                WellMeasurement(well="D1", wavelength_nm=450, source_extinction=0.5),
                WellMeasurement(well="E1", wavelength_nm=450, source_extinction=0.7),
            )),
        ), assay_metadata())
        result = AssayAnalyzer(data, settings).analyze()
        measurement = next(iter(result.measurement_results.values()))
        control = next(sample for sample in measurement.samples if sample.sample_type == SampleType.CONTROL)
        actual = measurement.calculation_extinction(control)
        self.assertEqual(control.calculation_input_scale.response_scale, ResponseScale.RAW)
        self.assertAlmostEqual(actual.statistics.mean, 0.6)
        self.assertIsNone(actual.scale.blank_subtracted)
        changed = replace(measurement, calculation_details=replace(
            measurement.calculation_details, formula_settings=FormulaSettings(Formula("sample * 2"), response_scale="blank_corrected"),
        ))
        self.assertEqual(changed.calculation_extinction(replace(control, name="New name")), actual)

    def test_partial_and_nan_curve_inputs_never_fall_back_to_other_extinctions(self):
        result = AssayAnalyzer(dual_data(), dual_settings()).analyze()
        patient = next(sample for sample in result.samples if sample.sample_id == "p1")
        for missing in (None, float("nan")):
            with self.subTest(missing=missing):
                changed = replace(patient, measurements=(replace(patient.measurements[0], curve_extinction=missing), *patient.measurements[1:]))
                actual = result.patient_extinction(changed)
                self.assertIsNotNone(actual)
                self.assertEqual(actual.scale, patient.calculation_input_scale)
                self.assertTrue(np.isnan(actual.statistics.mean))
        other_channel = replace(patient, selected_wavelength_nm=405)
        self.assertIsNone(result.patient_extinction(other_channel))

    def test_scale_rejects_contradictory_correction_metadata(self):
        values = dict(stage="formula_input", wavelength_nm=450, reference_wavelength_nm=None, blank_subtracted=None, preprocessing_expression=None)
        for invalid in (
            {"response_scale": "raw", "blank_subtracted": 0},
            {"blank_subtracted": float("nan")},
            {"preprocessing_expression": "sample * 2"},
        ):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                ExtinctionScale(**(values | invalid))


if __name__ == "__main__":
    unittest.main()
