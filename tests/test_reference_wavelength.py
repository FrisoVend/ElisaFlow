"""Reference correction shares the single and switching analysis pathway."""
import unittest
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
from fixtures import assay_metadata, samples_from_records
from skax_fixtures import write_skax
from table_helpers import sample_table, well_table
from test_multi_wavelength import dual_data, dual_settings, saturate

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.analysis.wavelengths import measurement_inputs
from elisaflow.assays.registry import AssayRegistry
from elisaflow.configuration.assay_settings import AssaySettings, MeasurementSettings
from elisaflow.models import AssayData, ResultStatus
from elisaflow.readers.assay_loader import AssayLoader
from elisaflow.readers.skax_document import SkaxDocument
from elisaflow.readers.skax_reader import SkaxReader
from elisaflow.reporting.pdf_report import PdfReportWriter


def reference_case(*, switching=True):
    data = dual_data()
    if not switching:
        data = AssayData(samples_from_records(well_table(data.samples).loc[well_table(data.samples)['wavelength_nm'].eq(450)].copy().to_dict('records')), data.metadata)
    reference = well_table(data.samples).loc[well_table(data.samples)["wavelength_nm"].eq(450)].copy()
    reference["wavelength_nm"] = 620
    reference["source_extinction"] = np.linspace(0.1, 0.3, len(reference))
    background = reference.set_index("well")["source_extinction"]
    _readings = well_table(data.samples)
    _readings["source_extinction"] += _readings["well"].map(background)
    data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
    data = AssayData(samples_from_records(pd.concat([well_table(data.samples), reference.iloc[::-1]], ignore_index=True).to_dict('records')), data.metadata)
    config = dual_settings().to_dict()
    config["measurement"] = {
        "mode": "switch_on_concentration" if switching else "single",
        "reference_wavelength_nm": 620,
    }
    if switching:
        config["measurement"]["preferred_wavelength_nm"] = 450
        config["measurement"]["cutoff_concentration"] = 1000.0
    else:
        config["calculation"].update(config["calculation"].pop("wavelengths")["450"])
        scope = config["validation"].pop("wavelengths")["450"]
        scope["controls"]["QC"].update(config["validation"]["result"]["controls"]["QC"])
        config["validation"]["result"] = scope
    return data, AssaySettings.from_dict("test.skax", config)


class ReferenceWavelengthTests(unittest.TestCase):
    def test_settings_round_trip_and_reference_is_not_an_analysis_channel(self):
        _, settings = reference_case()
        measurement = settings.measurement
        self.assertEqual(measurement.analysis_wavelengths((405, 450, 620)), (405, 450))
        self.assertEqual(MeasurementSettings.from_dict(measurement.to_dict()), measurement)
        for reference in (True, "620", 620.0, 0, -1, 630, 450):
            with self.subTest(reference=reference), self.assertRaises(ValueError):
                replace(measurement, reference_wavelength_nm=reference).analysis_wavelengths((405, 450, 620))
        config = settings.to_dict()
        config["validation"]["wavelengths"]["620"] = {}
        with self.assertRaisesRegex(ValueError, 'validatiegolflengten'):
            from elisaflow.assays.protocol_samples import ProtocolSamples

            ProtocolSamples((), {}, (405, 450, 620)).validate_measurements(AssaySettings.from_dict("test.skax", config))

    def test_correction_aligns_wells_preserves_raw_values_and_precedes_blank_and_fitting(self):
        for switching in (False, True):
            with self.subTest(switching=switching):
                data, settings = reference_case(switching=switching)
                original = well_table(data.samples).copy(deep=True)
                inputs = measurement_inputs(data, settings.measurement)
                self.assertEqual(set(inputs), {405, 450} if switching else {450})
                result = AssayAnalyzer(data, settings).analyze()
                for wavelength, measurement in result.measurement_results.items():
                    baseline = well_table(dual_data().samples).query("wavelength_nm == @wavelength").set_index("well")
                    actual = well_table(measurement.samples)
                    np.testing.assert_allclose(actual["source_extinction"], baseline.loc[actual["well"], "source_extinction"])
                    np.testing.assert_allclose(actual["raw_extinction"] - actual["reference_extinction"], actual["source_extinction"])
                    self.assertAlmostEqual(measurement.blank_value, 0.01 if wavelength == 450 else 0.005)
                    self.assertEqual(measurement.label, f"{wavelength} − 620 nm")
                patients = sample_table(result.samples).set_index("sample_id")
                self.assertAlmostEqual(patients.loc["p1", "result_mean"], 20, places=5)
                self.assertAlmostEqual(patients.loc["p2", "result_mean"], 2, places=5)
                pd.testing.assert_frame_equal(well_table(data.samples), original)

    def test_raw_saturation_still_switches_the_whole_sample(self):
        data, settings = reference_case()
        data = saturate(data, 450, "p1")
        result = AssayAnalyzer(data, settings).analyze()
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 405)
        self.assertAlmostEqual(patient["result_mean"], 20, places=5)

    def test_disabled_signal_does_not_make_another_replicates_reference_unusable(self):
        data, settings = reference_case()
        index = well_table(data.samples).index[well_table(data.samples)["wavelength_nm"].eq(620) & well_table(data.samples)["sample_id"].eq("p1")][0]
        well = well_table(data.samples).loc[index, "well"]
        _readings = well_table(data.samples)
        _readings.loc[index, "saturated"] = True
        data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
        _readings = well_table(data.samples)
        _readings.loc[_readings["well"].eq(well), "disabled"] = True
        data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
        patient = sample_table(AssayAnalyzer(data, settings).analyze().samples).set_index("sample_id").loc["p1"]
        self.assertAlmostEqual(patient["result_mean"], 20, places=5)

    def test_all_invalid_reference_replicates_keep_an_unavailable_patient(self):
        for switching in (False, True):
            with self.subTest(switching=switching):
                data, settings = reference_case(switching=switching)
                mask = well_table(data.samples)["wavelength_nm"].eq(620) & well_table(data.samples)["sample_id"].eq("p1")
                _readings = well_table(data.samples)
                _readings.loc[mask, "saturated"] = True
                data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
                result = AssayAnalyzer(data, settings).analyze()
                patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
                self.assertTrue(np.isnan(patient["result_mean"]))
                self.assertEqual(patient["count"], 0)
                self.assertEqual(patient["result_status"], ResultStatus.UNAVAILABLE.value)

    def test_optional_flags_and_dilutions_can_be_omitted(self):
        data, settings = reference_case(switching=False)
        data = AssayData(samples_from_records(well_table(data.samples).drop(columns=['saturated', 'disabled', 'sample_dilution']).to_dict('records')), data.metadata)
        result = AssayAnalyzer(data, settings).analyze()
        self.assertAlmostEqual(sample_table(result.samples).set_index("sample_id").loc["p1", "result_mean"], 20, places=5)

    def test_invalid_reference_invalidates_both_channels_without_hiding_the_patient(self):
        for switching in (False, True):
            for failure in ("missing", "saturated", "disabled", 'oneindig'):
                with self.subTest(switching=switching, failure=failure):
                    data, settings = reference_case(switching=switching)
                    index = well_table(data.samples).index[well_table(data.samples)["wavelength_nm"].eq(620) & well_table(data.samples)["sample_id"].eq("p1")][0]
                    if failure in ("missing", 'oneindig'):
                        _readings = well_table(data.samples)
                        _readings.loc[index, "source_extinction"] = np.nan if failure == "missing" else np.inf
                        data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
                    else:
                        _readings = well_table(data.samples)
                        _readings.loc[index, failure] = True
                        data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
                    result = AssayAnalyzer(data, settings).analyze()
                    patients = sample_table(result.samples).set_index("sample_id")
                    self.assertTrue(np.isnan(patients.loc["p1", "result_mean"]))
                    self.assertEqual(patients.loc["p1", "result_status"], ResultStatus.UNAVAILABLE.value)
                    self.assertAlmostEqual(patients.loc["p2", "result_mean"], 2, places=5)
                    for measurement in result.measurement_results.values():
                        self.assertTrue(well_table(measurement.samples)["reference_invalid"].any())

    def test_missing_duplicate_and_mismatched_reference_wells_are_rejected(self):
        for failure in ("missing_channel", "missing_well", "dubbele", "identity"):
            with self.subTest(failure=failure):
                data, settings = reference_case()
                mask = well_table(data.samples)["wavelength_nm"].eq(620)
                index = well_table(data.samples).index[mask][0]
                if failure == "missing_channel":
                    readings = well_table(data.samples).loc[~mask]
                    data = replace(data, samples=samples_from_records(readings.to_dict("records")))
                elif failure == "missing_well":
                    readings = well_table(data.samples).drop(index=index)
                    data = replace(data, samples=samples_from_records(readings.to_dict("records")))
                elif failure == "dubbele":
                    data = AssayData(samples_from_records(pd.concat([well_table(data.samples), well_table(data.samples).loc[[index]]]).to_dict('records')), data.metadata)
                else:
                    _readings = well_table(data.samples)
                    _readings.loc[index, "sample_name"] = "Different patient"
                    data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
                with self.assertRaises(ValueError):
                    measurement_inputs(data, settings.measurement)

    def test_skax_loading_and_pdf_label_both_corrected_curves(self):
        data, settings = reference_case()
        with TemporaryDirectory() as directory:
            path = write_skax(Path(directory) / "measurement.skax", well_table(data.samples).to_dict("records"))
            loaded = AssayLoader.load(path, measurement_settings=settings.measurement)
            self.assertEqual(set(well_table(loaded.samples)["wavelength_nm"]), {405, 450, 620})
            result = AssayAnalyzer(loaded, settings).analyze()
            report = PdfReportWriter(result, analyst="operator").write(Path(directory) / "report.pdf")
            self.assertTrue(report.read_bytes().startswith(b"%PDF"))
            import pymupdf
            with pymupdf.open(report) as document:
                text = " ".join(page.get_text() for page in document)
            self.assertIn("450", text)
            self.assertIn("405", text)
            self.assertIn("620", text)

    def test_reference_protocols_define_raw_signals_with_valid_hashes(self):
        project = Path(__file__).resolve().parents[1]
        registry = AssayRegistry(project / "configs")
        for name, wavelength in (("Dengue_NS1_as", 450), ("Qkoorts_IgG", 405), ("Qkoorts_IgM", 405)):
            with self.subTest(assay=name):
                document = SkaxDocument.from_file(project / "protocols" / f"{name}.skax")
                document.validate_hashes()
                settings = registry.get(f"{name}.skax").measurement
                self.assertEqual(settings.reference_wavelength_nm, 620)
                self.assertEqual(settings.analysis_wavelengths(document.wavelengths_nm), (wavelength,))

    def test_migrated_assays_match_their_previous_precalculated_input(self):
        project = Path(__file__).resolve().parents[1]
        registry = AssayRegistry(project / "configs")
        for name in ("Dengue_NS1_as", "Qkoorts_IgG", "Qkoorts_IgM"):
            with self.subTest(assay=name):
                path = project / "protocols" / f"{name}.skax"
                settings = registry.get(path)
                corrected = SkaxReader.from_file(path).wells.copy()
                corrected["source_extinction"] = corrected["sample_name"].map({
                    "Positief": 0.6, "Grens": 0.2, "Negatief": 0.1, "Standaard": 0.8,
                }).fillna(0.4)
                corrected.loc[corrected["sample_type"].eq("Blank"), "source_extinction"] = 0.05
                metadata = assay_metadata()
                expected = AssayAnalyzer(AssayData(samples_from_records(corrected.to_dict('records')), metadata), replace(settings, measurement=MeasurementSettings())).analyze()
                reference = corrected.copy()
                reference["wavelength_nm"] = 620
                reference["source_extinction"] = np.linspace(0.1, 0.3, len(reference))
                signal = corrected.copy()
                signal["wavelength_nm"] = settings.measurement.analysis_wavelengths(SkaxReader.from_file(path).wavelengths_nm)[0]
                signal["source_extinction"] += reference["source_extinction"]
                raw = AssayData(samples_from_records(pd.concat([signal, reference], ignore_index=True).to_dict('records')), metadata)
                actual = AssayAnalyzer(raw, settings).analyze()
                np.testing.assert_allclose(sample_table(actual.samples)["result_mean"], sample_table(expected.samples)["result_mean"], rtol=1e-10, equal_nan=True)
                np.testing.assert_allclose(sample_table(actual.samples)["extinction_mean"], sample_table(expected.samples)["extinction_mean"], atol=1e-14)
                self.assertEqual(actual.validation_status, expected.validation_status)
