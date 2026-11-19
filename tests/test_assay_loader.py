import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from xml.etree import ElementTree as ET
from zipfile import ZipFile

import numpy as np
from skax_fixtures import basic_records, replace_member, write_skax

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.configuration.assay_settings import AssaySettings, MeasurementSettings
from elisaflow.readers.assay_loader import AssayLoader


class AssayLoaderTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.path = Path(temporary.name) / "run.skax"

    def load(self, records=None, **options):
        write_skax(self.path, records, **options)
        return AssayLoader.load(self.path, measurement_settings=MeasurementSettings())

    @staticmethod
    def qc_records(wavelengths=(450,)):
        rows = basic_records(wavelengths)
        rows.extend(
            dict(well="A4", sample_id="standard", sample_name="Standard", sample_type="Standard", sample_concentration=1,
                 source_extinction=0.8, wavelength_nm=wavelength)
            for wavelength in wavelengths
        )
        return rows + [dict(row, well=row["well"].replace("A", "B")) for row in rows if row["sample_type"] != "Unknown"]

    def test_each_planned_qc_replicate_requires_all_protocol_wavelengths_including_reference(self):
        layout = self.qc_records((450, 405, 620))
        settings = MeasurementSettings("switch_on_concentration", 450, cutoff_concentration=10, reference_wavelength_nm=620)
        for well in ("A1", "A2", "A4"):
            for wavelength in (None, 450, 405, 620):
                rows = [
                    row for row in layout
                    if not (row["well"] == well and (wavelength is None or row["wavelength_nm"] == wavelength))
                ]
                write_skax(self.path, rows, layout=layout)
                with self.subTest(well=well, wavelength=wavelength), self.assertRaisesRegex(
                    ValueError, rf"Meetrecord ontbreekt: {well} bij {wavelength or 450} nm",
                ):
                    AssayLoader.load(self.path, measurement_settings=settings)

    def test_disabled_saturated_and_nonfinite_qc_records_count_as_present(self):
        rows = [
            dict(row, disabled=True, saturated=True, source_extinction=float("nan")) if row["sample_type"] != "Unknown" else row
            for row in self.qc_records((450, 620))
        ]
        write_skax(self.path, rows)
        data = AssayLoader.load(self.path, measurement_settings=MeasurementSettings(reference_wavelength_nm=620))
        qc = [well for sample in data.samples if sample.sample_type.value != "Unknown" for well in sample.measurements]
        self.assertEqual(len(qc), 12)
        self.assertTrue(all(well.disabled and well.saturated and np.isnan(well.source_extinction) for well in qc))

    def test_patient_records_are_not_required_at_every_wavelength(self):
        layout = self.qc_records((450, 620))
        rows = [row for row in layout if row["sample_type"] != "Unknown" or row["wavelength_nm"] == 450]
        write_skax(self.path, rows, layout=layout)
        data = AssayLoader.load(self.path, measurement_settings=MeasurementSettings(reference_wavelength_nm=620))
        patient = next(sample for sample in data.samples if sample.name == "Patient")
        self.assertEqual([well.wavelength_nm for well in patient.measurements], [450])

    def test_qc_roles_are_only_required_when_present_in_the_layout(self):
        for kind in (None, "Blank", "Control", "Standard"):
            rows = [row for row in self.qc_records() if row["sample_type"] in (kind, "Unknown")]
            with self.subTest(kind=kind):
                data = self.load(rows)
                self.assertEqual({sample.sample_type.value for sample in data.samples}, {"Unknown"} | ({kind} if kind else set()))

    def test_complete_session_supplies_samples_readings_and_metadata(self):
        data = self.load(name="My assay")
        self.assertEqual([sample.name for sample in data.samples], ["Blank", "Cutoff", "Patient"])
        self.assertEqual([well.source_extinction for sample in data.samples for well in sample.measurements], [0.1, 0.6, 1.1])
        self.assertEqual(data.metadata.assay_name, "My assay")
        self.assertEqual(data.metadata.measurement_wavelengths, ("450",))
        self.assertEqual(data.metadata.session_time, "2026-10-06T10:00:00+02:00")
        self.assertIsNone(data.metadata.instrument_name)
        settings = AssaySettings.from_dict("run.skax", {
            "lot_number": "lot", "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / cutoff", "references": {"cutoff": "Cutoff"}},
        })
        result = AssayAnalyzer(data, settings).analyze()
        self.assertAlmostEqual(result.samples[-1].result.mean, 2)
        self.assertEqual(list(self.path.parent.iterdir()), [self.path])

    def test_only_recorded_wells_are_loaded_without_missing_value_placeholders(self):
        layout = basic_records()
        data = self.load(layout[:2], layout=layout, wavelengths=(450,))
        self.assertEqual([sample.name for sample in data.samples], ["Blank", "Cutoff"])
        self.assertEqual(sum(len(sample.measurements) for sample in data.samples), 2)

    def test_raw_records_group_replicates_by_identity_and_dilution_in_protocol_order(self):
        rows = [
            dict(well=well, sample_id=identity, sample_name=name, sample_type=kind,
                 sample_dilution=dilution, sample_concentration=concentration,
                 source_extinction=value, wavelength_nm=wavelength)
            for wavelength in (405, 450)
            for well, identity, name, kind, dilution, concentration, value in (
                ("B2", "standard", "Standard", "Standard", 1, 3.2, 0.8),
                ("A2", "patient", "Patient", "Unknown", 5, None, 0.4),
                ("B1", "patient", "Patient", "Unknown", 2.5, None, 0.3),
                ("A1", "patient", "Patient", "Unknown", 2.5, None, 0.2),
            )
        ]
        write_skax(self.path, rows, wavelengths=(450, 405))
        settings = MeasurementSettings("switch_on_concentration", 450, cutoff_concentration=10)
        data = AssayLoader.load(self.path, measurement_settings=settings)
        self.assertEqual([(sample.sample_id, sample.dilution) for sample in data.samples],
                         [("patient", 2.5), ("patient", 5), ("standard", 1)])
        self.assertEqual([sample.concentration for sample in data.samples], [None, None, 3.2])
        self.assertEqual([(well.well, well.wavelength_nm, well.order) for well in data.samples[0].measurements],
                         [("A1", 450, 0), ("B1", 450, 1), ("A1", 405, 4), ("B1", 405, 5)])
        self.assertEqual([well.source_extinction for well in data.samples[0].measurements], [0.2, 0.3, 0.2, 0.3])

    def test_plate_identifiers_are_mapped_without_assuming_zero_based_indices(self):
        data = self.load([dict(row, plate_index=7) for row in basic_records()])
        self.assertTrue(all(well.plate_index == 7 for sample in data.samples for well in sample.measurements))

    def test_missing_plate_names_remain_empty_strings(self):
        write_skax(self.path)
        with ZipFile(self.path) as archive:
            session = ET.fromstring(archive.read("Session.xml"))
        plate = session.find("PlateLayout/Plates/Plate")
        plate.remove(plate.find("Name"))
        replace_member(self.path, "Session.xml", ET.tostring(session))
        data = AssayLoader.load(self.path, measurement_settings=MeasurementSettings())
        self.assertEqual({well.plate_name for sample in data.samples for well in sample.measurements}, {""})

    def test_saturation_and_disabled_indices_stay_with_the_original_measurement_records(self):
        rows = basic_records((450, 405))[::-1]
        rows[0]["saturated"] = True
        rows[2]["disabled"] = True
        write_skax(self.path, rows)
        settings = MeasurementSettings("switch_on_concentration", 450, cutoff_concentration=10)
        data = AssayLoader.load(self.path, measurement_settings=settings)
        flags = {(well.well, well.wavelength_nm): (well.saturated, well.disabled) for sample in data.samples for well in sample.measurements}
        self.assertEqual(flags[(rows[0]["well"], rows[0]["wavelength_nm"])], (True, False))
        self.assertEqual(flags[(rows[2]["well"], rows[2]["wavelength_nm"])], (False, True))
        self.assertEqual(sum(disabled for _, disabled in flags.values()), 1)

    def test_optional_disabled_file_defaults_to_no_disabled_records(self):
        data = self.load()
        with ZipFile(self.path) as archive:
            self.assertNotIn("Disableds.data", archive.namelist())
        self.assertFalse(any(well.disabled for sample in data.samples for well in sample.measurements))

    def test_disabled_indices_must_refer_to_existing_measurements(self):
        for indices in ((-1,), (3,), (999,)):
            with self.subTest(indices=indices), self.assertRaisesRegex(ValueError, "Disableds.data"):
                self.load(disabled_indices=indices)

    def test_archive_member_hashes_are_verified(self):
        for member in ("Session.xml", "Measurements.data", "Disableds.data"):
            with self.subTest(member=member):
                write_skax(self.path, disabled_indices=(1,))
                with ZipFile(self.path) as archive:
                    contents = archive.read(member)
                replace_member(self.path, member, contents + b" ", update_hashes=False)
                with self.assertRaisesRegex(ValueError, "SKAX-.*Hash"):
                    AssayLoader.load(self.path, measurement_settings=MeasurementSettings())

    def test_raw_records_must_have_sixteen_numeric_fields_and_valid_flags(self):
        for value in (b"1,2,3\n", b"1,0,0,0,0,1,0,2,0.5,0,450,0,0,0,0,0\n", b"1,0,0,0,0,1,0,0,bad,0,450,0,0,0,0,0\n"):
            with self.subTest(value=value):
                write_skax(self.path)
                replace_member(self.path, "Measurements.data", value)
                with self.assertRaises(ValueError):
                    AssayLoader.load(self.path, measurement_settings=MeasurementSettings())

    def test_non_endpoint_records_duplicate_measurements_and_unknown_wells_are_rejected(self):
        for change in ({"iteration": 1}, {"well_point_x": 1}, {"well_point_y": 1}, {"well": "H12"}):
            rows = basic_records()
            rows[0] = rows[0] | change
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.load(rows, layout=basic_records())
        rows = basic_records()
        with self.assertRaisesRegex(ValueError, "Dubbele wellmetingen"):
            self.load(rows + [rows[0]])

    def test_wavelength_and_step_mappings_must_match_the_records(self):
        rows = basic_records((450, 405))
        for index in (0, 99):
            rows[0]["step_index"] = index
            with self.subTest(step=index), self.assertRaisesRegex(ValueError, "meetgolflengten"):
                write_skax(self.path, rows)
                AssayLoader.load(self.path, measurement_settings=MeasurementSettings("switch_on_concentration", 450, cutoff_concentration=10))
        rows[0]["step_index"] = 2
        with self.assertRaisesRegex(ValueError, "meetgolflengten"):
            write_skax(self.path, rows)
            AssayLoader.load(self.path, measurement_settings=MeasurementSettings("switch_on_concentration", 450, cutoff_concentration=10))

    def test_every_configured_wavelength_must_have_recorded_measurements(self):
        write_skax(self.path, basic_records(), layout=basic_records((450, 405)), wavelengths=(450, 405))
        with self.assertRaisesRegex(ValueError, "alle ruwe protocolgolflengten"):
            AssayLoader.load(self.path, measurement_settings=MeasurementSettings("switch_on_concentration", 450, cutoff_concentration=10))

    def test_templates_and_unfinished_sessions_are_not_analyzed(self):
        with self.assertRaisesRegex(ValueError, "voltooide SKAX"):
            self.load(completed=False)
        for end in ("", "0001-01-01T00:00:00", "2020-01-01T00:00:00+02:00", "invalid"):
            write_skax(self.path)
            with ZipFile(self.path) as archive:
                session = ET.fromstring(archive.read("Session.xml"))
            session.find("ExecutionEnd").text = end
            replace_member(self.path, "Session.xml", ET.tostring(session))
            with self.subTest(end=end), self.assertRaisesRegex(ValueError, "voltooide SKAX"):
                AssayLoader.load(self.path, measurement_settings=MeasurementSettings())

    def test_multiple_measured_plates_are_rejected(self):
        rows = basic_records()
        rows[0]["plate_index"] = 2
        with self.assertRaisesRegex(ValueError, "precies een gemeten plaat"):
            self.load(rows)

    def test_fitted_curves_reference_correction_and_flags_match_typed_inputs(self):
        from table_helpers import sample_table, well_table
        from test_reference_wavelength import reference_case
        for switching in (False, True):
            original, settings = reference_case(switching=switching)
            write_skax(self.path, well_table(original.samples).to_dict("records"))
            loaded = AssayLoader.load(self.path, measurement_settings=settings.measurement)
            expected = sample_table(AssayAnalyzer(original, settings).analyze().samples).set_index("sample_id").sort_index()
            actual = sample_table(AssayAnalyzer(loaded, settings).analyze().samples).set_index("sample_id").sort_index()
            np.testing.assert_allclose(actual["result_mean"], expected["result_mean"], rtol=1e-8, equal_nan=True)
            np.testing.assert_allclose(actual["extinction_mean"], expected["extinction_mean"], rtol=1e-12, equal_nan=True)
