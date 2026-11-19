"""Unresolved patient saturation preserves justified concentration bounds."""

import unittest
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
import test_calibration_response as calibration_response
from table_helpers import well_table
from test_multi_wavelength import dual_data, dual_settings, saturate
from test_reference_wavelength import reference_case

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.calibration.logistic import FivePLModel, FivePLParameters, FourPLModel, FourPLParameters
from elisaflow.calibration.polynomial import CubicModel, CubicParameters, QuadraticModel, QuadraticParameters
from elisaflow.calibration.power_law import PowerLawModel, PowerLawParameters
from elisaflow.configuration.assay_settings import AssaySettings, CalibrationSettings, MeasurementSettings
from elisaflow.formatting import format_result
from elisaflow.formulas.expressions import Formula
from elisaflow.models import ResultStatus, SampleType
from elisaflow.reporting.pdf_report import PdfReportWriter


def single_case():
    data = dual_data()
    data = replace(data, samples=tuple(
        replace(sample, measurements=tuple(well for well in sample.measurements if well.wavelength_nm == 450))
        for sample in data.samples
    ))
    config = dual_settings().to_dict()
    config["measurement"] = {"mode": "single"}
    config["calculation"].update(config["calculation"].pop("wavelengths")["450"])
    config["validation"] = {}
    return data, AssaySettings.from_dict("test.skax", config)


def change_patient(data, **changes):
    return replace(data, samples=tuple(
        replace(sample, measurements=tuple(replace(well, **changes) for well in sample.measurements))
        if sample.sample_type == SampleType.UNKNOWN else sample
        for sample in data.samples
    ))


class SaturationBoundTests(unittest.TestCase):
    def patient(self, data, settings, identity="p1"):
        return next(sample for sample in AssayAnalyzer(data, settings).analyze().samples if sample.sample_id == identity)

    def test_all_fitted_models_preserve_saturation_bounds_and_dilution_without_using_flagged_values(self):
        cases = (
            ("4pl", FourPLModel(), FourPLParameters(2.5, -1.3, 2, 0.1)),
            ("4pl", FourPLModel(), FourPLParameters(2.5, 1.3, 2, 0.1)),
            ("5pl", FivePLModel(), FivePLParameters(2.5, -1.3, 2, 0.1, 0.4)),
            ("5pl", FivePLModel(), FivePLParameters(2.5, 1.3, 2, 0.1, 0.4)),
            ("quadratic", QuadraticModel(), QuadraticParameters(0.1, 0.3, 0.02)),
            ("cubic", CubicModel(), CubicParameters(0.1, 0.3, 0.02, 0.001)),
            ("power_law", PowerLawModel(), PowerLawParameters(0.3, 1.2, 0.1)),
            ("power_law", PowerLawModel(), PowerLawParameters(-0.3, 0.7, 6)),
        )
        for name, model, parameters in cases:
            data = calibration_response.CalibrationResponseTests.data(model, parameters)
            increasing = model.calculate_extinctions(50, parameters) > model.calculate_extinctions(0.1, parameters)
            expected_status = ResultStatus.ABOVE_RANGE if increasing else ResultStatus.BELOW_RANGE
            expected_limit = 150 if increasing else 0.3
            for flagged_value in (0, float("nan"), float("inf")):
                with self.subTest(model=name, increasing=increasing, flagged_value=flagged_value):
                    changed = change_patient(data, saturated=True, source_extinction=flagged_value)
                    original = well_table(changed.samples).copy(deep=True)
                    patient = self.patient(changed, calibration_response.CalibrationResponseTests.settings(name, response_scale="raw"), "patient")
                    self.assertEqual(patient.result_status, expected_status)
                    self.assertAlmostEqual(patient.result_limit, expected_limit)
                    self.assertTrue(np.isnan(patient.result.mean))
                    self.assertEqual(len(patient.retained), 0)
                    self.assertTrue(all(not well.included and np.isnan(well.curve_extinction) for well in patient.measurements))
                    self.assertTrue(all(well.saturated and well.exclusion_reason for well in patient.measurements))
                    pd.testing.assert_frame_equal(well_table(changed.samples), original)

    def test_excluded_top_standard_and_configured_range_limit_the_reported_bound(self):
        data, settings = single_case()
        data = saturate(saturate(data, 450, "p1", all_replicates=True), 450, "std-300", all_replicates=True)
        for configured, expected in ((None, 100), ((3, 80), 80)):
            with self.subTest(configured=configured):
                patient = self.patient(data, replace(settings, calculation=replace(settings.calculation, concentration_range=configured)))
                self.assertEqual(patient.result_status, ResultStatus.ABOVE_RANGE)
                self.assertEqual(patient.result_limit, expected)

    def test_one_saturated_duplicate_does_not_silently_report_the_other_duplicate(self):
        data, settings = single_case()
        patient = self.patient(saturate(data, 450, "p1"), settings)
        self.assertEqual(patient.result_status, ResultStatus.UNAVAILABLE)
        self.assertIsNone(patient.result_limit)
        self.assertTrue(np.isnan(patient.result.mean))
        self.assertEqual(len(patient.retained), 1)

    def test_saturated_and_unsaturated_above_range_duplicates_agree(self):
        data, settings = single_case()
        response = FourPLModel().calculate_extinctions(600, FourPLParameters(4.5, -1.2, 30, 0.08))
        data = change_patient(data, source_extinction=float(response))
        patient = self.patient(saturate(data, 450, "p1"), settings)
        self.assertEqual(patient.result_status, ResultStatus.ABOVE_RANGE)
        self.assertEqual(patient.result_limit, 300)
        self.assertTrue(np.isnan(patient.result.mean))

    def test_both_saturated_channels_select_the_stronger_bound(self):
        data = saturate(saturate(dual_data(), 450, "p1", all_replicates=True), 405, "p1", all_replicates=True)
        data = saturate(data, 450, "std-300", all_replicates=True)
        patient = self.patient(data, dual_settings())
        self.assertEqual(patient.result_status, ResultStatus.ABOVE_RANGE)
        self.assertEqual(patient.result_limit, 300)
        self.assertEqual(patient.selected_wavelength_nm, 405)

    def test_both_decreasing_channels_select_the_stronger_lower_bound(self):
        data = dual_data()
        parameters = {450: FourPLParameters(4.5, 1.2, 30, 0.08), 405: FourPLParameters(1.4, 1.2, 30, 0.025)}
        data = replace(data, samples=tuple(replace(sample, measurements=tuple(
            replace(well, source_extinction=float(FourPLModel().calculate_extinctions(sample.concentration, parameters[well.wavelength_nm])))
            if sample.sample_type != SampleType.BLANK else well for well in sample.measurements
        )) for sample in data.samples))
        data = saturate(saturate(saturate(data, 450, "p1", all_replicates=True), 405, "p1", all_replicates=True),
                        450, "std-1", all_replicates=True)
        patient = self.patient(data, dual_settings())
        self.assertEqual(patient.result_status, ResultStatus.BELOW_RANGE)
        self.assertEqual(patient.result_limit, 1)
        self.assertEqual(patient.selected_wavelength_nm, 405)

    def test_valid_alternate_numeric_result_still_has_priority(self):
        patient = self.patient(saturate(dual_data(), 450, "p1", all_replicates=True), dual_settings())
        self.assertEqual(patient.result_status, ResultStatus.IN_RANGE)
        self.assertEqual(patient.selected_wavelength_nm, 405)
        self.assertAlmostEqual(patient.result.mean, 20, places=5)

    def test_failed_curve_does_not_prevent_the_other_channels_saturation_bound(self):
        data = saturate(dual_data(), 405, "p1", all_replicates=True)
        for identity in ("std-30", "std-100", "std-300"):
            data = saturate(data, 450, identity, all_replicates=True)
        result = AssayAnalyzer(data, dual_settings()).analyze()
        self.assertIsNotNone(result.measurement_results[450].calculation_details.calibration_error)
        patient = next(sample for sample in result.samples if sample.sample_id == "p1")
        self.assertEqual(patient.result_status, ResultStatus.ABOVE_RANGE)
        self.assertEqual(patient.result_limit, 300)
        self.assertEqual(patient.selected_wavelength_nm, 405)

    def test_disabled_saturation_and_saturated_controls_do_not_acquire_bounds(self):
        data, settings = single_case()
        disabled = change_patient(saturate(data, 450, "p1", all_replicates=True), disabled=True)
        patient = self.patient(disabled, settings)
        self.assertEqual(patient.result_status, ResultStatus.UNAVAILABLE)
        result = AssayAnalyzer(saturate(data, 450, "qc", all_replicates=True), settings).analyze()
        control = next(sample for sample in result.samples if sample.sample_id == "qc")
        self.assertEqual(control.result_status, ResultStatus.UNAVAILABLE)
        self.assertTrue(all(np.isnan(well.result_limit) for well in control.measurements))

    def test_reference_correction_must_still_support_the_bound(self):
        data, settings = reference_case(switching=False)
        data = replace(data, samples=tuple(replace(sample, measurements=tuple(
            replace(well, source_extinction=0.05) if sample.sample_id == "p1" and well.wavelength_nm == 620 else well
            for well in sample.measurements
        )) for sample in data.samples))
        data = saturate(data, 450, "p1", all_replicates=True)
        patient = self.patient(data, settings)
        self.assertEqual(patient.result_status, ResultStatus.ABOVE_RANGE)
        for changes in ({"source_extinction": 100}, {"saturated": True}, {"disabled": True}):
            changed = replace(data, samples=tuple(
                replace(sample, measurements=tuple(
                    replace(well, **changes) if sample.sample_id == "p1" and well.wavelength_nm == 620 else well
                    for well in sample.measurements
                )) for sample in data.samples
            ))
            with self.subTest(changes=changes):
                patient = self.patient(changed, settings)
                self.assertEqual(patient.result_status, ResultStatus.UNAVAILABLE)
                self.assertIsNone(patient.result_limit)

    def test_affine_preprocessing_preserves_direction_but_nonlinear_preprocessing_does_not_invent_a_bound(self):
        data, settings = single_case()
        data = saturate(data, 450, "p1", all_replicates=True)
        for expression, parameters, status in (
            ("sample * 2", FourPLParameters(9, -1.2, 30, 0.16), ResultStatus.ABOVE_RANGE),
            ("10 - sample", FourPLParameters(5.5, -1.2, 30, 9.92), ResultStatus.ABOVE_RANGE),
            ("sample * sample", FourPLParameters(20.25, -1.2, 30, 0.0064), ResultStatus.UNAVAILABLE),
        ):
            with self.subTest(expression=expression):
                calculation = replace(settings.calculation, parameters=parameters, weighting_exponent=0, preprocessing=Formula(expression))
                patient = self.patient(data, replace(settings, calculation=calculation))
                self.assertEqual(patient.result_status, status)

        control_mean = float(FourPLModel().calculate_extinctions(40, FourPLParameters(4.5, -1.2, 30, 0.08)))
        calculation = replace(
            settings.calculation, parameters=FourPLParameters(9 / control_mean - 0.001, -1.2, 30, 0.16 / control_mean - 0.001),
            weighting_exponent=0, preprocessing=Formula("(+sample * 2) / qc + -0.001", {"qc": "QC"}),
        )
        patient = self.patient(data, replace(settings, calculation=calculation))
        self.assertEqual(patient.result_status, ResultStatus.ABOVE_RANGE)
        self.assertEqual(patient.result_limit, 300)

    def test_fixed_curves_require_reliable_standard_evidence_and_a_reporting_range(self):
        data, settings = single_case()
        data = saturate(data, 450, "p1", all_replicates=True)
        parameters = FourPLParameters(4.5, -1.2, 30, 0.08)
        for scale in ("raw", "blank_corrected"):
            parameters_for_scale = parameters if scale == "raw" else replace(parameters, A=4.49, D=0.07)
            calculation = CalibrationSettings("4pl", parameters_for_scale, concentration_range=(1, 100), response_scale=scale)
            fixed = replace(settings, calculation=calculation)
            with self.subTest(scale=scale):
                self.assertEqual(self.patient(data, fixed).result_limit, 100)
                no_standards = replace(data, samples=tuple(sample for sample in data.samples if sample.sample_type != SampleType.STANDARD))
                self.assertEqual(self.patient(no_standards, fixed).result_status, ResultStatus.UNAVAILABLE)
                unbounded = replace(fixed, calculation=replace(calculation, concentration_range=None))
                self.assertEqual(self.patient(data, unbounded).result_status, ResultStatus.UNAVAILABLE)
                unsupported = replace(fixed, calculation=replace(calculation, concentration_range=(1, 600)))
                self.assertEqual(self.patient(data, unsupported).result_status, ResultStatus.UNAVAILABLE)

    def test_reconstructed_standard_signals_do_not_establish_a_detector_bound(self):
        from test_saturation_conversion import conversion_settings

        data = saturate(saturate(dual_data(), 450, "p1", all_replicates=True), 405, "p1", all_replicates=True)
        data = saturate(data, 450, "std-300", all_replicates=True)
        result = AssayAnalyzer(data, conversion_settings()).analyze()
        for wavelength, expected in ((450, ResultStatus.UNAVAILABLE), (405, ResultStatus.ABOVE_RANGE)):
            patient = next(sample for sample in result.measurement_results[wavelength].samples if sample.sample_id == "p1")
            self.assertEqual(patient.result_status, expected)
        self.assertEqual(next(sample for sample in result.samples if sample.sample_id == "p1").result_limit, 300)

    def test_decreasing_fixed_curve_does_not_report_a_bound_below_zero(self):
        model, parameters = FourPLModel(), FourPLParameters(2.5, 1.3, 2, 0.1)
        data = change_patient(calibration_response.CalibrationResponseTests.data(model, parameters), saturated=True)
        calculation = CalibrationSettings("4pl", parameters, response_scale="raw", concentration_range=(0, 50))
        patient = self.patient(data, AssaySettings("test.skax", "lot", calculation=calculation), "patient")
        self.assertEqual(patient.result_status, ResultStatus.UNAVAILABLE)
        self.assertIsNone(patient.result_limit)

    def test_pdf_displays_the_diluted_bound_and_preserves_saturation_markers(self):
        import pymupdf

        data, settings = single_case()
        data = change_patient(data, saturated=True)
        data = replace(data, samples=tuple(replace(sample, dilution=10) if sample.sample_id == "p1" else sample for sample in data.samples))
        result = AssayAnalyzer(data, replace(settings, measurement=MeasurementSettings())).analyze()
        patient = next(sample for sample in result.samples if sample.sample_id == "p1")
        self.assertEqual(format_result(patient.result.mean, patient.result_status.value, patient.result_limit, 1), ">3000.0")
        with TemporaryDirectory() as directory:
            path = PdfReportWriter(result, analyst="operator").write(Path(directory) / "report.pdf")
            with pymupdf.open(path) as document:
                text = " ".join(page.get_text() for page in document)
        self.assertIn(">3000.0", text)
        self.assertIn("SAT", text)


if __name__ == "__main__":
    unittest.main()
