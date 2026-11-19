"""Configuration and measurement contracts across the complete analysis pathway."""

import copy
import unittest
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
from fixtures import borrelia_settings, samples_from_records
from table_helpers import sample_table, well_table
from test_multi_wavelength import dual_data, dual_settings, mixed_calibration, saturate

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.analysis.wavelengths import select_measurement_results
from elisaflow.configuration.assay_settings import AssaySettings, CalibrationSettings
from elisaflow.models import AnalysisResult, AssayData, MeasurementResult, ResultStatus


class ConfigurationContractTests(unittest.TestCase):
    def test_current_config_round_trips_without_mutating_input(self):
        settings = dual_settings()
        config = settings.to_dict()
        original = copy.deepcopy(config)
        restored = AssaySettings.from_dict("test.skax", config)
        self.assertEqual(restored, settings)
        self.assertEqual(config, original)
        self.assertIs(restored.validation, restored.validation)
        self.assertEqual(restored.validation.to_dict(), config["validation"])

    def test_obsolete_config_fields_are_rejected(self):
        config = dual_settings().to_dict() | {"validation_by_wavelength": {}}
        with self.assertRaisesRegex(ValueError, 'Onbekende'):
            AssaySettings.from_dict("test.skax", config)
        for fields in ({"fit_parameters": {}}, {"preferred_wavelength_nm": 450}):
            with self.subTest(fields=fields), self.assertRaisesRegex(ValueError, 'Onbekende'):
                CalibrationSettings.from_dict({"response_scale": "raw", "method": "calibration", "fit_model": "4pl", **fields})
        config = dual_settings().to_dict()
        config["controls"] = {"QC": {"validation": {"result": {"min": 38}}}}
        with self.assertRaisesRegex(ValueError, 'Onbekende'):
            AssaySettings.from_dict("test.skax", config)

    def test_single_assay_uses_the_same_validation_contract(self):
        settings = borrelia_settings()
        config = settings.to_dict()
        self.assertEqual(config["measurement"], {"mode": "single"})
        self.assertIn("result", config["validation"]["result"]["controls"]["Positief"])
        self.assertEqual(AssaySettings.from_dict(settings.protocol_name, config), settings)

    def test_old_selected_scope_is_rejected(self):
        config = dual_settings().to_dict()
        config["validation"]["selected"] = config["validation"].pop("result")
        with self.assertRaisesRegex(ValueError, 'Onbekende.*selected'):
            AssaySettings.from_dict("test.skax", config)

    def test_result_limits_are_rejected_in_wavelength_scopes(self):
        for target in (
            {"controls": {"QC": {"result": {"min": 1}}}}, {"standards": {"300": {"result": {"max": 400}}}}, {"blank": {"result": {"max": 1}}},
        ):
            config = dual_settings().to_dict()
            config["validation"]["wavelengths"]["450"] = target
            with self.subTest(target=target), self.assertRaisesRegex(ValueError, 'Resultaatcriteria moeten in validation.result'):
                AssaySettings.from_dict("test.skax", config)

    def test_single_measurement_rejects_duplicate_wavelength_scopes(self):
        config = borrelia_settings().to_dict()
        config["validation"]["wavelengths"] = {"450": {"blank": {"extinction": {"response_scale": "raw", "max": 0.1}}}}
        with self.assertRaisesRegex(ValueError, 'Criteria voor een enkele meting horen in validation.result'):
            AssaySettings.from_dict("test.skax", config)


class MeasurementPipelineTests(unittest.TestCase):
    def test_contradictory_preferred_result_cannot_override_a_valid_alternate(self):
        settings = dual_settings()
        result = AssayAnalyzer(dual_data(), settings).analyze()
        measurements = dict(result.measurement_results)
        primary = measurements[450]
        measurements[450] = replace(primary, samples=tuple(
            replace(sample, measurements=tuple(
                replace(well, result=12, result_status=ResultStatus.ABOVE_RANGE, result_limit=10) for well in sample.measurements
            )) if sample.sample_id == "p1" else sample for sample in primary.samples
        ))
        patient = next(sample for sample in select_measurement_results(measurements, settings.measurement) if sample.sample_id == "p1")
        self.assertEqual(patient.selected_wavelength_nm, 405)
        self.assertEqual(patient.result_status, ResultStatus.IN_RANGE)
        self.assertAlmostEqual(patient.result.mean, 20, places=5)

    def test_assay_can_be_analyzed_repeatedly_without_mutating_the_inputs(self):
        data = dual_data()
        settings = dual_settings()

        original = data.samples
        analyzer = AssayAnalyzer(data, settings)
        first = analyzer.analyze()
        second = analyzer.analyze()
        self.assertIs(data.samples, original)
        pd.testing.assert_frame_equal(sample_table(first.samples), sample_table(second.samples))
        for wavelength in first.measurement_results:
            pd.testing.assert_frame_equal(
                well_table(first.measurement_results[wavelength].samples), well_table(second.measurement_results[wavelength].samples),
            )
        self.assertTrue(all(np.isnan(well.extinction) for sample in data.samples for well in sample.measurements))

    def test_unexpected_calibration_errors_are_not_treated_as_unusable_curves(self):
        for error in (ValueError("Unexpected input bug"), RuntimeError("Unexpected implementation bug")):
            with self.subTest(error=error), patch("elisaflow.analysis.calculation.create_calibration_curve", side_effect=error):
                with self.assertRaises(type(error)) as raised:
                    AssayAnalyzer(dual_data(), dual_settings()).analyze()
                self.assertIs(raised.exception, error)

    def test_polynomial_solver_failure_retains_the_other_wavelength(self):
        data, settings = mixed_calibration()
        with patch("scipy.optimize.minimize", return_value=SimpleNamespace(success=False, x=np.zeros(3))):
            result = AssayAnalyzer(data, settings).analyze()
        self.assertFalse(result.measurement_results[405].usable_for_results)
        self.assertIn('monotone polynomiale kalibratie', result.measurement_results[405].calculation_details.calibration_error)
        self.assertTrue(result.measurement_results[450].usable_for_results)
        patients = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")]
        self.assertTrue(patients["selected_wavelength_nm"].eq(450).all())

    def test_single_and_dual_measurements_produce_the_same_flat_result_type(self):
        data = dual_data()
        dual = AssayAnalyzer(data, dual_settings()).analyze()
        config = dual_settings().to_dict()
        config["measurement"] = {"mode": "single"}
        config["calculation"].update(config["calculation"].pop("wavelengths")["450"])
        config["validation"].pop("wavelengths")
        single_data = AssayData(samples_from_records(well_table(data.samples)[well_table(data.samples)['wavelength_nm'].eq(450)].to_dict('records')), data.metadata)
        single = AssayAnalyzer(single_data, AssaySettings.from_dict("test.skax", config)).analyze()
        self.assertEqual(len(single.measurement_results), 1)
        self.assertEqual(len(dual.measurement_results), 2)
        for result in (single, dual):
            for measurement in result.measurement_results.values():
                self.assertIsInstance(measurement, MeasurementResult)
                self.assertNotIsInstance(measurement, AnalysisResult)
                self.assertFalse(hasattr(measurement, "measurement_results"))
        self.assertAlmostEqual(sample_table(single.samples).set_index("sample_id").loc["p1", "result_mean"], 20, places=5)

    def test_mismatched_sample_identity_cannot_switch_into_another_sample(self):
        data = dual_data()
        data = saturate(data, 450, "p1")
        mask = well_table(data.samples)["wavelength_nm"].eq(405) & well_table(data.samples)["sample_id"].eq("p1")
        _readings = well_table(data.samples)
        _readings.loc[mask, "sample_id"] = "other patient"
        data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
        with self.assertRaisesRegex(ValueError, 'Monsteridentiteiten'):
            AssayAnalyzer(data, dual_settings()).analyze()


if __name__ == "__main__":
    unittest.main()
