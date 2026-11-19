"""Select usable calibration coverage without extrapolation or mixed replicates."""

import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
from fixtures import assay_metadata, samples_from_records
from table_helpers import sample_table, well_table
from test_multi_wavelength import dual_data, dual_settings, saturate

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.assays.registry import AssayRegistry
from elisaflow.calibration.logistic import FourPLModel, FourPLParameters
from elisaflow.configuration.assay_settings import AssaySettings
from elisaflow.models import AssayData, ResultStatus, ValidationStatus
from elisaflow.readers.skax_reader import SkaxReader


def patient_readings(data, concentrations, *, dilution=1):
    values = {
        wavelength: iter(FourPLModel().calculate_extinctions(np.broadcast_to(concentrations, (2,)), FourPLParameters(top, -1.2, 30, baseline)))
        for wavelength, top, baseline in ((450, 4.5, 0.08), (405, 1.4, 0.025))
    }
    return replace(data, samples=tuple(
        replace(sample, dilution=dilution, measurements=tuple(
            replace(well, source_extinction=float(next(values[well.wavelength_nm]))) for well in sample.measurements
        )) if sample.sample_id == 'p1' else sample for sample in data.samples
    ))


def range_settings():
    config = dual_settings().to_dict()
    config["validation"]["wavelengths"]["450"].pop("standards")
    return AssaySettings.from_dict("test.skax", config)


class CurveSelectionTests(unittest.TestCase):
    def test_disabled_or_saturated_top_level_restricts_primary_and_switches_the_whole_sample(self):
        for flags in (("disabled",), ("saturated",), ("disabled", "saturated")):
            with self.subTest(flags=flags):
                data = dual_data()
                mask = well_table(data.samples)["wavelength_nm"].eq(450) & well_table(data.samples)["sample_id"].eq("std-300")
                _readings = well_table(data.samples)
                _readings.loc[mask, list(flags)] = True
                data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
                data = patient_readings(data, [90, 150])
                result = AssayAnalyzer(data, range_settings()).analyze()
                patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
                self.assertEqual(result.measurement_results[450].calculation_details.curve.concentration_range, (1, 100))
                self.assertEqual(patient["selected_wavelength_nm"], 405)
                self.assertAlmostEqual(patient["result_mean"], 120, places=4)
                wells = well_table(result.samples).loc[well_table(result.samples)["sample_id"].eq("p1")]
                self.assertEqual(wells["selected_wavelength_nm"].tolist(), [405, 405])

    def test_both_above_range_report_the_widest_supported_dilution_corrected_bound(self):
        data = dual_data()
        data = saturate(data, 450, "std-300", all_replicates=True)
        data = patient_readings(data, 400, dilution=10)
        result = AssayAnalyzer(data, range_settings()).analyze()
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 405)
        self.assertEqual(patient["result_status"], ResultStatus.ABOVE_RANGE.value)
        self.assertEqual(patient["result_limit"], 3000)
        self.assertTrue(np.isnan(patient["result_mean"]))
        selected = next(sample for sample in result.samples if sample.sample_id == "p1")
        extinction = result.patient_extinction(selected)
        self.assertIsNotNone(extinction)
        self.assertEqual(extinction.scale.wavelength_nm, 405)
        self.assertTrue(np.isfinite(extinction.statistics.mean))

    def test_low_range_switching_and_lowest_supported_bound(self):
        for concentration, expected_status in ((2, ResultStatus.IN_RANGE), (0.5, ResultStatus.BELOW_RANGE)):
            with self.subTest(concentration=concentration):
                data = dual_data()
                for identity in ("std-1", "std-3"):
                    data = saturate(data, 450, identity, all_replicates=True)
                data = patient_readings(data, concentration)
                result = AssayAnalyzer(data, range_settings()).analyze()
                patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
                self.assertEqual(patient["selected_wavelength_nm"], 405)
                self.assertEqual(patient["result_status"], expected_status.value)
                if expected_status == ResultStatus.IN_RANGE:
                    self.assertAlmostEqual(patient["result_mean"], 2, places=5)
                else:
                    self.assertEqual(patient["result_limit"], 1)
                    self.assertTrue(np.isnan(patient["result_mean"]))

    def test_unusable_fallback_does_not_block_a_usable_preferred_curve(self):
        data = dual_data()
        for identity in ("std-30", "std-100", "std-300"):
            data = saturate(data, 405, identity, all_replicates=True)
        result = AssayAnalyzer(data, dual_settings()).analyze()
        self.assertIsNotNone(result.measurement_results[405].calculation_details.calibration_error)
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 450)
        self.assertAlmostEqual(patient["result_mean"], 20, places=5)
        self.assertEqual(result.validation_status, ValidationStatus.PASSED)

    def test_both_failed_curves_keep_patients_unavailable_with_diagnostics(self):
        data = dual_data()
        for wavelength in (450, 405):
            for identity in ("std-30", "std-100", "std-300"):
                data = saturate(data, wavelength, identity, all_replicates=True)
        result = AssayAnalyzer(data, dual_settings()).analyze()
        patients = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")]
        self.assertTrue(patients["result_mean"].isna().all())
        self.assertTrue(patients["selected_wavelength_nm"].isna().all())
        self.assertTrue(patients["result_status"].eq(ResultStatus.UNAVAILABLE.value).all())
        self.assertTrue(all(measurement.calculation_details.calibration_error for measurement in result.measurement_results.values()))
        self.assertEqual(result.validation_status, ValidationStatus.NOT_EVALUATED)

    def test_failed_preferred_curve_qc_does_not_trigger_switching(self):
        config = dual_settings().to_dict()
        config["validation"]["wavelengths"]["450"]["controls"]["QC"]["extinction"] = {"min": 5, "response_scale": "blank_corrected"}
        result = AssayAnalyzer(dual_data(), AssaySettings.from_dict("test.skax", config)).analyze()
        self.assertTrue(result.measurement_results[450].usable_for_results)
        samples = sample_table(result.samples)
        self.assertTrue(samples.loc[samples["sample_type"].ne("Blank"), "selected_wavelength_nm"].eq(450).all())
        self.assertTrue(samples.loc[samples["sample_type"].eq("Blank"), "selected_wavelength_nm"].isna().all())
        patients = samples.set_index("sample_id").loc[["p1", "p2"]]
        self.assertTrue(patients["selection_reason"].eq("preferred").all())
        np.testing.assert_allclose(patients["result_mean"], [20, 2], rtol=1e-6)
        self.assertTrue(any(check.wavelength_nm == 450 and check.status == ValidationStatus.FAILED for check in result.validation_results))
        self.assertEqual(result.validation_status, ValidationStatus.FAILED)

    def test_failed_fallback_qc_does_not_block_switching_for_range(self):
        config = range_settings().to_dict()
        config["validation"]["wavelengths"]["405"]["blank"]["extinction"] = {"max": 0.001, "response_scale": "raw"}
        data = dual_data()
        data = saturate(data, 450, "std-300", all_replicates=True)
        data = patient_readings(data, 150)
        result = AssayAnalyzer(data, AssaySettings.from_dict("test.skax", config)).analyze()
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertTrue(result.measurement_results[405].usable_for_results)
        self.assertEqual(patient["selected_wavelength_nm"], 405)
        self.assertEqual(patient["selection_reason"], "preferred_out_of_range")
        self.assertEqual(patient["result_status"], ResultStatus.IN_RANGE.value)
        self.assertAlmostEqual(patient["result_mean"], 150, places=4)
        self.assertEqual(result.validation_status, ValidationStatus.FAILED)

    def test_shared_control_failure_still_fails_the_run(self):
        config = dual_settings().to_dict()
        config["validation"]["result"]["controls"]["QC"]["result"] = {"max": 1}
        result = AssayAnalyzer(dual_data(), AssaySettings.from_dict("test.skax", config)).analyze()
        self.assertTrue(all(measurement.usable_for_results for measurement in result.measurement_results.values()))
        self.assertEqual(result.validation_status, ValidationStatus.FAILED)

    def test_both_failed_wavelength_qc_leave_preferred_results_unchanged(self):
        config = dual_settings().to_dict()
        for wavelength in ("450", "405"):
            config["validation"]["wavelengths"][wavelength]["blank"]["extinction"] = {"max": 0, "response_scale": "raw"}
        result = AssayAnalyzer(dual_data(), AssaySettings.from_dict("test.skax", config)).analyze()
        patients = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")]
        self.assertTrue(patients["selected_wavelength_nm"].eq(450).all())
        np.testing.assert_allclose(patients.set_index("sample_id").loc[["p1", "p2"], "result_mean"], [20, 2], rtol=1e-6)
        self.assertEqual(result.validation_status, ValidationStatus.FAILED)

    def test_conflicting_upper_and_lower_bounds_remain_unavailable(self):
        data = dual_data()
        data = saturate(data, 450, "std-300", all_replicates=True)
        data = patient_readings(data, 150)
        mask = well_table(data.samples)["wavelength_nm"].eq(405) & well_table(data.samples)["sample_id"].eq("p1")
        _readings = well_table(data.samples)
        _readings.loc[mask, "source_extinction"] = FourPLModel().calculate_extinctions(0.5, FourPLParameters(1.4, -1.2, 30, 0.025))
        data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
        result = AssayAnalyzer(data, range_settings()).analyze()
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["result_status"], ResultStatus.UNAVAILABLE.value)
        self.assertTrue(np.isnan(patient["result_mean"]))
        self.assertTrue(np.isnan(patient["result_limit"]))

    def test_low_preferred_readings_are_not_dropped_when_the_fallback_can_quantify_them(self):
        data = dual_data()
        mask = well_table(data.samples)["wavelength_nm"].eq(450) & well_table(data.samples)["sample_id"].eq("p1")
        _readings = well_table(data.samples)
        _readings.loc[mask, "source_extinction"] = 0.005
        data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
        result = AssayAnalyzer(data, dual_settings()).analyze()
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 405)
        self.assertEqual(patient["count"], 2)
        self.assertAlmostEqual(patient["result_mean"], 20, places=5)


class SingleCurveRangeTests(unittest.TestCase):
    def test_losing_a_standard_limits_the_single_curve_and_dilution_corrected_patient_bound(self):
        data = saturate(patient_readings(dual_data(), 150, dilution=10), 450, "std-300", all_replicates=True)
        data = replace(data, samples=tuple(
            replace(sample, measurements=tuple(well for well in sample.measurements if well.wavelength_nm == 450))
            for sample in data.samples
        ))
        config = dual_settings().to_dict()
        config["measurement"] = {"mode": "single"}
        config["calculation"].update(config["calculation"].pop("wavelengths")["450"])
        config["validation"].pop("wavelengths")
        result = AssayAnalyzer(data, AssaySettings.from_dict("test.skax", config)).analyze()
        patient = next(sample for sample in result.samples if sample.sample_id == "p1")
        self.assertEqual(result.calculation_details.curve.concentration_range, (1, 100))
        self.assertEqual(patient.result_status, ResultStatus.ABOVE_RANGE)
        self.assertEqual(patient.result_limit, 1000)
        self.assertTrue(np.isnan(patient.result.mean))


class IA2CurveSelectionTests(unittest.TestCase):
    def test_actual_ia2_configuration_uses_405_after_one_or_two_top_450_levels_are_disabled(self):
        root = Path(__file__).resolve().parents[1]
        protocol = SkaxReader.from_file(root / "protocols" / "test_IA-2_as.skax")
        settings = AssayRegistry(root / "configs").get("test_IA-2_as.skax")
        for cutoff in (4000, 350):
            with self.subTest(excluded_levels_from=cutoff):
                frames = []
                for wavelength, top, baseline, blank in ((450, 4.5, 0.08, 0.01), (405, 1.4, 0.025, 0.005)):
                    frame = protocol.wells.copy()
                    nominal = frame["sample_concentration"].where(frame["sample_type"].eq("Standard"), 600.0)
                    nominal = nominal.where(~frame["sample_type"].eq("Control"), 125.0)
                    frame["source_extinction"] = FourPLModel().calculate_extinctions(nominal, FourPLParameters(top, -1.2, 30, baseline))
                    frame.loc[frame["sample_type"].eq("Blank"), "source_extinction"] = blank
                    frame.loc[frame["sample_type"].eq("Unknown"), "sample_dilution"] = 1
                    frame["wavelength_nm"] = wavelength
                    excluded = frame["sample_type"].eq("Standard") & frame["sample_concentration"].ge(cutoff) & (wavelength == 450)
                    frame["saturated"] = excluded
                    frame["disabled"] = excluded
                    frames.append(frame)
                result = AssayAnalyzer(AssayData(samples_from_records(pd.concat(frames, ignore_index=True).to_dict('records')), assay_metadata(name=settings.protocol_name)), settings).analyze()
                patients = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")]
                self.assertTrue(patients["selected_wavelength_nm"].eq(405).all())
                np.testing.assert_allclose(patients["result_mean"], 600, rtol=1e-6)
                self.assertEqual(result.measurement_results[405].calculation_details.curve.concentration_range, (7.5, 4000))
                if cutoff == 4000:
                    self.assertEqual(result.validation_status, ValidationStatus.PASSED)
                    self.assertEqual(result.measurement_results[450].calculation_details.curve.concentration_range, (7.5, 350))
                else:
                    self.assertEqual(result.validation_status, ValidationStatus.NOT_EVALUATED)
                    self.assertIn('ten minste 4 verschillende', result.measurement_results[450].calculation_details.calibration_error)


class ConcentrationCutoffTests(unittest.TestCase):
    def settings(self, *, preferred=450, cutoff=100):
        settings = range_settings()
        return replace(settings, measurement=replace(settings.measurement, preferred_wavelength_nm=preferred, cutoff_concentration=cutoff))

    def test_cutoff_switches_without_saturation_including_the_boundary(self):
        for concentration, expected in ((99, 450), (100, 405), (101, 405)):
            with self.subTest(concentration=concentration):
                result = AssayAnalyzer(patient_readings(dual_data(), concentration), self.settings()).analyze()
                patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
                self.assertEqual(patient["selected_wavelength_nm"], expected)
                self.assertEqual(patient["selection_reason"], "concentration_cutoff" if expected == 405 else "preferred")
                self.assertAlmostEqual(patient["result_mean"], concentration, places=4)

    def test_a_duplicate_above_cutoff_switches_the_whole_sample(self):
        result = AssayAnalyzer(patient_readings(dual_data(), [90, 110]), self.settings()).analyze()
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 405)
        self.assertEqual(patient["count"], 2)
        self.assertAlmostEqual(patient["result_mean"], 100, places=4)

    def test_cutoff_is_evaluated_before_dilution_correction(self):
        result = AssayAnalyzer(patient_readings(dual_data(), 20, dilution=10), self.settings()).analyze()
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 450)
        self.assertAlmostEqual(patient["result_mean"], 200, places=4)

    def test_preferred_405_can_switch_to_450_at_cutoff(self):
        result = AssayAnalyzer(patient_readings(dual_data(), 150), self.settings(preferred=405)).analyze()
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 450)
        self.assertAlmostEqual(patient["result_mean"], 150, places=4)

    def test_unavailable_selected_channel_does_not_return_preferred_above_cutoff(self):
        data = saturate(patient_readings(dual_data(), 150), 405, "p1", all_replicates=True)
        data = replace(data, samples=tuple(replace(sample, measurements=tuple(
            replace(well, disabled=True) if sample.sample_id == "p1" and well.wavelength_nm == 405 else well
            for well in sample.measurements
        )) for sample in data.samples))
        result = AssayAnalyzer(data, self.settings()).analyze()
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertTrue(np.isnan(patient["result_mean"]))
        self.assertEqual(patient["result_status"], ResultStatus.UNAVAILABLE.value)

    def test_saturated_selected_channel_preserves_its_bound_above_cutoff(self):
        data = saturate(patient_readings(dual_data(), 150), 405, "p1", all_replicates=True)
        result = AssayAnalyzer(data, self.settings()).analyze()
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertTrue(np.isnan(patient["result_mean"]))
        self.assertEqual(patient["result_status"], ResultStatus.ABOVE_RANGE.value)
        self.assertEqual(patient["result_limit"], 300)
        self.assertEqual(patient["selected_wavelength_nm"], 405)
        self.assertEqual(patient["selection_reason"], "concentration_cutoff")

    def test_selected_channel_reporting_bound_is_preserved(self):
        data = saturate(patient_readings(dual_data(), 150), 405, "std-300", all_replicates=True)
        result = AssayAnalyzer(data, self.settings()).analyze()
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 405)
        self.assertEqual(patient["result_status"], ResultStatus.ABOVE_RANGE.value)
        self.assertEqual(patient["result_limit"], 100)


if __name__ == "__main__":
    unittest.main()
