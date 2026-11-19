"""Fixed optical conversions preserve raw data and restore usable curve inputs."""

import unittest
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
from table_helpers import sample_table, well_table
from test_multi_wavelength import dual_data, dual_settings, saturate
from test_reference_wavelength import reference_case

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.analysis.wavelengths import measurement_inputs
from elisaflow.configuration.assay_settings import AssaySettings, MeasurementSettings, SaturationSettings
from elisaflow.models import ResultStatus
from elisaflow.reporting.pdf_report import PdfReportWriter

SLOPE = (4.5 - 0.08) / (1.4 - 0.025)
INTERCEPT = 0.08 - SLOPE * 0.025


def conversion_settings(target=450, **changes):
    config = dual_settings().to_dict()
    config["measurement"]["saturation"] = {
        str(target): {
            "mode": "calculate", "source_wavelength_nm": 405 if target == 450 else 450,
            "slope": SLOPE if target == 450 else 1 / SLOPE,
            "intercept": INTERCEPT if target == 450 else -INTERCEPT / SLOPE,
        }
    }
    config["measurement"].update(changes)
    return AssaySettings.from_dict("test.skax", config)


class SaturationConfigurationTests(unittest.TestCase):
    def test_defaults_and_round_trip_preserve_explicit_coefficients(self):
        self.assertEqual(SaturationSettings().to_dict(), {"mode": "exclude"})
        settings = conversion_settings()
        self.assertEqual(AssaySettings.from_dict("test.skax", settings.to_dict()), settings)
        self.assertEqual(settings.measurement.saturation[450].slope, SLOPE)
        with self.assertRaises(TypeError):
            settings.measurement.saturation[405] = SaturationSettings()

    def test_switching_requires_a_positive_finite_cutoff_and_rejects_old_mode(self):
        base = {"mode": "switch_on_concentration", "preferred_wavelength_nm": 450}
        for cutoff in (None, True, "100", 0, -1, float("nan"), float("inf")):
            with self.subTest(cutoff=cutoff), self.assertRaisesRegex(ValueError, "cutoff_concentration"):
                MeasurementSettings.from_dict(base | {"cutoff_concentration": cutoff})
        with self.assertRaises(ValueError):
            MeasurementSettings.from_dict({"mode": "switch_on_saturation", "preferred_wavelength_nm": 450})
        with self.assertRaises(ValueError):
            MeasurementSettings(cutoff_concentration=100)

    def test_invalid_conversion_parameters_and_sources_are_rejected(self):
        valid = {"mode": "calculate", "source_wavelength_nm": 405, "slope": 3.4, "intercept": 0}
        invalid = [
            {"mode": "automatic"}, {"mode": "exclude", "slope": 3.4},
            *[valid | {name: value} for name in ("slope", "intercept") for value in (None, True, "3.4", float("nan"), float("inf"))],
            valid | {"slope": 0}, valid | {"slope": -1},
            *[valid | {"source_wavelength_nm": value} for value in (None, True, "405", 405.5, 0, -405)],
        ]
        for rule in invalid:
            with self.subTest(rule=rule), self.assertRaises(ValueError):
                SaturationSettings.from_dict(rule)
        for source in (450, 620):
            config = conversion_settings().measurement.to_dict()
            config["saturation"]["450"]["source_wavelength_nm"] = source
            with self.subTest(source=source), self.assertRaises(ValueError):
                MeasurementSettings.from_dict(config).analysis_wavelengths((405, 450))
        for mapping in (None, [], {True: {}}, {"450": {}, "0450": {}}):
            with self.subTest(mapping=mapping), self.assertRaises(ValueError):
                MeasurementSettings.from_dict({"saturation": mapping})

    def test_reference_cannot_be_a_conversion_source_or_target(self):
        for target, source in ((450, 620), (620, 405)):
            with self.subTest(target=target, source=source), self.assertRaisesRegex(ValueError, "referentiegolflengte"):
                MeasurementSettings.from_dict({
                    "mode": "switch_on_concentration", "preferred_wavelength_nm": 450,
                    "cutoff_concentration": 100, "reference_wavelength_nm": 620,
                    "saturation": {str(target): {"mode": "calculate", "source_wavelength_nm": source, "slope": 2, "intercept": 0}},
                })


class SaturationConversionTests(unittest.TestCase):
    def test_configured_coefficients_are_used_even_when_they_differ_from_the_plate_relationship(self):
        data = saturate(dual_data(), 450, "p1", all_replicates=True)
        measurement = conversion_settings().measurement
        measurement = replace(measurement, saturation={450: SaturationSettings("calculate", 405, 2.1, 0.37)})
        inputs = measurement_inputs(data, measurement)
        donor = {well.key: well for sample in inputs[405].samples for well in sample.measurements}
        patient = next(sample for sample in inputs[450].samples if sample.sample_id == "p1")
        for well in patient.measurements:
            self.assertAlmostEqual(well.converted_extinction, 0.37 + 2.1 * donor[well.key].source_extinction)

    def test_reconstructed_standards_restore_full_4pl_without_changing_raw_readings(self):
        data = dual_data()
        for identity in ("std-100", "std-300"):
            data = saturate(data, 450, identity, all_replicates=True)
        original = well_table(data.samples).copy(deep=True)
        result = AssayAnalyzer(data, conversion_settings()).analyze()
        measurement = result.measurement_results[450]
        self.assertEqual(measurement.calculation_details.curve.concentration_range, (1, 300))
        self.assertEqual(len(measurement.calculation_details.curve.diagnostics.weights), 12)
        wells = well_table(measurement.samples)
        converted = wells.loc[wells["saturated"]]
        np.testing.assert_allclose(converted["converted_extinction"], converted["source_extinction"], rtol=1e-12)
        self.assertTrue(converted["included_in_analysis"].all())
        self.assertTrue(converted["conversion_source_wavelength_nm"].eq(405).all())
        standards = wells.loc[wells["sample_type"].eq("Standard")]
        np.testing.assert_allclose(standards["back_calculated_concentration"], standards["sample_concentration"], rtol=1e-6)
        pd.testing.assert_frame_equal(well_table(data.samples), original)

    def test_converted_patient_remains_at_preferred_wavelength_below_cutoff(self):
        data = saturate(dual_data(), 450, "p1")
        result = AssayAnalyzer(data, conversion_settings()).analyze()
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 450)
        self.assertAlmostEqual(patient["result_mean"], 20, places=5)
        self.assertEqual(patient["count"], 2)

    def test_reconstruction_does_not_override_the_concentration_cutoff(self):
        data = saturate(dual_data(), 450, "p1", all_replicates=True)
        result = AssayAnalyzer(data, conversion_settings(cutoff_concentration=10)).analyze()
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 405)
        self.assertEqual(patient["selection_reason"], "concentration_cutoff")
        self.assertAlmostEqual(patient["result_mean"], 20, places=5)

    def test_reverse_conversion_is_selectable(self):
        data = saturate(dual_data(), 405, "p1", all_replicates=True)
        result = AssayAnalyzer(data, conversion_settings(405, preferred_wavelength_nm=405)).analyze()
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 405)
        self.assertAlmostEqual(patient["result_mean"], 20, places=5)

    def test_invalid_source_and_disabled_target_are_never_reconstructed(self):
        for flag in ("saturated", "disabled", "missing", "infinite", "target_disabled"):
            data = saturate(dual_data(), 450, "p1", all_replicates=True)
            samples = []
            for sample in data.samples:
                wells = []
                for well in sample.measurements:
                    if sample.sample_id == "p1":
                        if well.wavelength_nm == 405:
                            if flag in ("saturated", "disabled"):
                                well = replace(well, **{flag: True})
                            elif flag in ("missing", "infinite"):
                                well = replace(well, source_extinction=float("nan") if flag == "missing" else float("inf"))
                        elif flag == "target_disabled":
                            well = replace(well, disabled=True)
                    wells.append(well)
                samples.append(replace(sample, measurements=tuple(wells)))
            inputs = measurement_inputs(replace(data, samples=tuple(samples)), conversion_settings().measurement)
            patient = next(sample for sample in inputs[450].samples if sample.sample_id == "p1")
            with self.subTest(flag=flag):
                self.assertTrue(all(well.converted_extinction is None for well in patient.measurements))

    def test_both_channels_saturated_cannot_convert_each_other(self):
        for reference in (False, True):
            with self.subTest(reference=reference):
                data = reference_case()[0] if reference else dual_data()
                for wavelength in (405, 450):
                    data = saturate(data, wavelength, "p1", all_replicates=True)
                config = conversion_settings().to_dict()
                if reference:
                    config["measurement"]["reference_wavelength_nm"] = 620
                config["measurement"]["saturation"].update(conversion_settings(405).to_dict()["measurement"]["saturation"])
                result = AssayAnalyzer(data, AssaySettings.from_dict("test.skax", config)).analyze()
                patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
                expected_status = ResultStatus.UNAVAILABLE if reference else ResultStatus.ABOVE_RANGE
                self.assertEqual(patient["result_status"], expected_status.value)
                self.assertTrue(np.isnan(patient["result_mean"]))
                for measurement in result.measurement_results.values():
                    self.assertTrue(all(well.converted_extinction is None for sample in measurement.samples for well in sample.measurements))

    def test_exclusion_is_default_and_never_estimates_coefficients(self):
        data = saturate(dual_data(), 450, "std-300", all_replicates=True)
        result = AssayAnalyzer(data, dual_settings()).analyze()
        self.assertEqual(result.measurement_results[450].calculation_details.curve.concentration_range, (1, 100))
        self.assertTrue(all(well.converted_extinction is None for sample in result.measurement_results[450].samples for well in sample.measurements))

    def test_conversion_follows_reference_subtraction_and_precedes_target_blank_subtraction(self):
        for target, source, blank in ((450, 405, 0.01), (405, 450, 0.005)):
            with self.subTest(target=target):
                data, settings = reference_case()
                for identity in ("p1", "std-100", "std-300"):
                    data = saturate(data, target, identity, all_replicates=True)
                rules = conversion_settings(target).measurement.saturation
                settings = replace(settings, measurement=replace(settings.measurement, preferred_wavelength_nm=target, saturation=rules))
                original = {(well.wavelength_nm, well.key): well for sample in data.samples for well in sample.measurements}
                original_table = well_table(data.samples).copy(deep=True)
                result = AssayAnalyzer(data, settings).analyze()
                measurement = result.measurement_results[target]
                self.assertAlmostEqual(measurement.blank_value, blank)
                self.assertEqual(measurement.calculation_details.curve.concentration_range, (1, 300))
                for sample in measurement.samples:
                    for well in sample.measurements:
                        reference = original[620, well.key].source_extinction
                        self.assertEqual(well.raw_extinction, original[target, well.key].source_extinction)
                        self.assertEqual(well.reference_extinction, reference)
                        if well.saturated:
                            donor = original[source, well.key].source_extinction
                            expected = rules[target].intercept + rules[target].slope * (donor - reference)
                            self.assertAlmostEqual(well.converted_extinction, expected)
                            self.assertAlmostEqual(well.extinction, expected - blank)
                            self.assertEqual(well.conversion_source_wavelength_nm, source)
                            self.assertTrue(well.included)
                patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
                self.assertEqual(patient["selected_wavelength_nm"], target)
                self.assertAlmostEqual(patient["result_mean"], 20, places=5)
                self.assertEqual(patient["count"], 2)
                pd.testing.assert_frame_equal(well_table(data.samples), original_table)

    def test_invalid_reference_prevents_saturation_reconstruction(self):
        for changes in ({"saturated": True}, {"disabled": True}, {"source_extinction": float("nan")}, {"source_extinction": float("inf")}):
            with self.subTest(changes=changes):
                data, settings = reference_case()
                data = saturate(data, 450, "p1", all_replicates=True)
                data = replace(data, samples=tuple(
                    replace(sample, measurements=tuple(
                        replace(well, **changes) if sample.sample_id == "p1" and well.wavelength_nm == 620 else well
                        for well in sample.measurements
                    )) for sample in data.samples
                ))
                settings = replace(settings, measurement=replace(settings.measurement, saturation=conversion_settings().measurement.saturation))
                result = AssayAnalyzer(data, settings).analyze()
                for measurement in result.measurement_results.values():
                    patient = next(sample for sample in measurement.samples if sample.sample_id == "p1")
                    self.assertTrue(all(well.reference_invalid and not well.included for well in patient.measurements))
                    self.assertTrue(all(well.converted_extinction is None for well in patient.measurements))
                patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
                self.assertEqual(patient["result_status"], ResultStatus.UNAVAILABLE.value)
                self.assertEqual(patient["count"], 0)

    def test_pdf_identifies_the_cutoff_and_reconstructed_readings(self):
        try:
            import pymupdf
        except ImportError:
            self.skipTest("Install pymupdf to inspect the report")
        data = saturate(dual_data(), 450, "std-300", all_replicates=True)
        result = AssayAnalyzer(data, conversion_settings(cutoff_concentration=100)).analyze()
        with TemporaryDirectory() as temporary:
            output = PdfReportWriter(result, analyst="operator").write(Path(temporary) / "converted.pdf")
            with pymupdf.open(output) as document:
                text = "\n".join(page.get_text() for page in document)
        self.assertIn("Omschakelconcentratie", text)
        self.assertIn("100 U/mL", text)
        self.assertIn("Verzadiging 450 nm", text)
        self.assertIn("oorspronkelijke meting was verzadigd", text)


if __name__ == "__main__":
    unittest.main()
