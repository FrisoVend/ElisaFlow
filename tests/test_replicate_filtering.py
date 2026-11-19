import json
import unittest
from dataclasses import FrozenInstanceError, replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
import pymupdf
from fixtures import assay_metadata, borrelia_settings, reference_parameters, samples_from_records
from table_helpers import sample_table, well_table

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.calibration.logistic import FourPLModel
from elisaflow.configuration.assay_settings import AssaySettings, FormulaOverrideSettings, PatientSettings, ReplicateFilterSettings, SampleSettings
from elisaflow.models import AssayData, ResponseScale, SampleType, ValidationMetric, ValidationStatus
from elisaflow.reporting.pdf_report import PdfReportWriter


class ReplicateFilterSettingsTests(unittest.TestCase):
    def test_every_configured_filter_requires_a_valid_explicit_response_scale(self):
        for collection in ("patients", "controls", "standards"):
            for scale in (None, True, "", "automatic", [], {}):
                with self.subTest(collection=collection, scale=scale):
                    rule = {"max_deviation_percent": 15}
                    if scale is not None:
                        rule["response_scale"] = scale
                    config = {
                        "lot_number": "lot",
                        "calculation": {"method": "formula", "expression": "sample", "response_scale": "raw"},
                        collection: {"replicate_filter": rule} if collection == "patients" else {
                            "QC" if collection == "controls" else "1": {"replicate_filter": rule},
                        },
                    }
                    with self.assertRaisesRegex(ValueError, "response_scale"):
                        AssaySettings.from_dict("test.skax", config)

    def test_patient_settings_round_trip_and_reject_unrelated_fields(self) -> None:
        config = {
            "lot_number": "lot",
            "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample"},
            "patients": {"replicate_filter": {"response_scale": "blank_corrected", "max_deviation_percent": 15, "max_exclusions": 2}},
        }
        settings = AssaySettings.from_dict("test.skax", config)
        self.assertEqual(settings.patients, PatientSettings(ReplicateFilterSettings(15, 2, response_scale="blank_corrected")))
        self.assertEqual(settings.to_dict()["patients"], config["patients"])
        self.assertEqual(AssaySettings.from_dict("test.skax", settings.to_dict()), settings)
        self.assertEqual(PatientSettings().to_dict(), {})
        with self.assertRaises(FrozenInstanceError):
            settings.patients.replicate_filter = None
        for field in ("lot_number", "preprocessing", "controls"):
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, f"Onbekende velden.*{field}"):
                PatientSettings.from_dict({field: "unexpected"})
        with self.assertRaisesRegex(ValueError, "ReplicateFilterSettings"):
            PatientSettings(replicate_filter={})
        with self.assertRaisesRegex(ValueError, "patients"):
            replace(settings, patients={})

    def test_filter_round_trips_through_assay_settings_and_is_immutable(self) -> None:
        rule = {"response_scale": "blank_corrected", "max_deviation_percent": 15, "max_exclusions": 1}
        config = {
            "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / cal", "references": {"cal": "Cal"}},
            "lot_number": "lot",
            "controls": {"Cal": {"replicate_filter": rule}},
            "standards": {"10": {"replicate_filter": rule}},
        }
        settings = AssaySettings.from_dict("test.skax", config)
        parsed = settings.controls["Cal"].replicate_filter
        self.assertEqual(parsed, ReplicateFilterSettings(15, 1, response_scale="blank_corrected"))
        exported = json.loads(json.dumps(settings.to_dict(), allow_nan=False))
        self.assertEqual(exported["controls"]["Cal"]["replicate_filter"], rule)
        self.assertEqual(AssaySettings.from_dict("test.skax", exported), settings)
        rule["max_deviation_percent"] = 99
        self.assertEqual(parsed.max_deviation_percent, 15)
        with self.assertRaises(FrozenInstanceError):
            parsed.max_deviation_percent = 20
        self.assertEqual(ReplicateFilterSettings.from_dict({"response_scale": "blank_corrected", "max_deviation_percent": 0}), ReplicateFilterSettings(0, response_scale="blank_corrected"))
        self.assertNotIn("replicate_filter", SampleSettings().to_dict())

    def test_invalid_filter_settings_are_rejected_at_both_boundaries(self) -> None:
        for name, values in (
            ("max_deviation_percent", [-1, True, "15", None, float("nan"), float("inf")]),
            ("max_exclusions", [-1, True, 1.0, "1", None]),
        ):
            for value in values:
                with self.subTest(field=name, value=value):
                    config = {"response_scale": "blank_corrected", "max_deviation_percent": 15} | {name: value}
                    with self.assertRaisesRegex(ValueError, name):
                        ReplicateFilterSettings.from_dict(config)
                    with self.assertRaisesRegex(ValueError, name):
                        ReplicateFilterSettings(**config)
        for config in ({}, None, [], {"response_scale": "blank_corrected", "max_deviation_percent": 15, "metric": "result"}):
            with self.subTest(config=config), self.assertRaises(ValueError):
                SampleSettings.from_dict({"replicate_filter": config})
        for name in ("expected_count", "min_remaining"):
            with self.subTest(field=name), self.assertRaisesRegex(ValueError, f"Onbekende velden.*{name}"):
                ReplicateFilterSettings.from_dict({"response_scale": "blank_corrected", "max_deviation_percent": 15, name: 3})
        with self.assertRaisesRegex(ValueError, "ReplicateFilterSettings"):
            SampleSettings(replicate_filter={})
        with self.assertRaisesRegex(ValueError, 'Replicaatfiltering van blanco’s'):
            AssaySettings.from_dict(
                "test.skax",
                {
                    "lot_number": "lot",
                    "blank": {"replicate_filter": {"response_scale": "blank_corrected", "max_deviation_percent": 15}},
                    "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / cal", "references": {"cal": "Cal"}},
                    "result": {"label": "Ratio"},
                },
            )


class PatientReplicateFilteringTests(unittest.TestCase):
    def setUp(self) -> None:
        self.settings = AssaySettings.from_dict(
            "test.skax",
            {
                "lot_number": "lot",
                "patients": {"replicate_filter": {"response_scale": "blank_corrected", "max_deviation_percent": 15, "max_exclusions": 1}},
                "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample * 2"},
            },
        )
        self.frame = pd.DataFrame(
            {
                "sample_id": ["blank"] + ["a"] * 3 + ["b"] * 3 + ["c"] * 3 + ["qc"] * 3,
                "sample_name": ["Blank"] + ["Shared"] * 6 + ["QC"] * 6,
                "sample_type": ["Blank"] + ["Unknown"] * 9 + ["Control"] * 3,
                "source_extinction": np.array([0, 1, 1.02, 1.3, 10, 10.2, 13, 2, 2.04, 2.6, 1, 1.02, 1.3]) + 0.2,
                "well": [f"A{index + 1}" for index in range(13)],
            },
        )

    def analyze(self, frame=None, settings=None):
        return AssayAnalyzer(
            AssayData(samples_from_records((self.frame if frame is None else frame).to_dict('records')), assay_metadata()), settings or self.settings,
        ).analyze()

    def test_patient_groups_are_independent_and_filtered_before_calculation(self) -> None:
        frame = self.frame.copy()
        frame.index = [0, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]
        original = frame.copy(deep=True)
        result = self.analyze(frame)
        np.testing.assert_array_equal(np.flatnonzero(~well_table(result.samples)["included_in_analysis"]), [3, 6, 9])
        self.assertTrue(well_table(result.samples).iloc[[3, 6, 9]]["result"].isna().all())
        for identity, mean, wells in (("a", 2.02, "A2, A3"), ("b", 20.2, "A5, A6"), ("c", 4.04, "A8, A9")):
            with self.subTest(patient=identity):
                patient = sample_table(result.samples).set_index("sample_id").loc[identity]
                self.assertEqual(patient["count"], 2)
                self.assertAlmostEqual(patient["result_mean"], mean)
                self.assertEqual(patient["used_wells"], wells)
        self.assertEqual(sample_table(result.samples).set_index("sample_id").loc["qc", "count"], 3)
        pd.testing.assert_frame_equal(frame, original)
        np.testing.assert_array_equal(well_table(result.samples)["source_extinction"], original["source_extinction"])

    def test_disabled_patient_filter_retains_every_replicate(self) -> None:
        result = self.analyze(settings=replace(self.settings, patients=PatientSettings()))
        self.assertTrue(well_table(result.samples)["included_in_analysis"].all())
        self.assertTrue(well_table(result.samples)["replicate_mean"].isna().all())

    def test_dilutions_of_one_patient_are_filtered_separately(self) -> None:
        config = self.settings.to_dict()
        config["calculation"] = {"response_scale": "blank_corrected", "method": "calibration", "fit_model": "4pl", "parameters": vars(reference_parameters())}
        settings = AssaySettings.from_dict("test.skax", config)
        frame = self.frame.iloc[:7].copy()
        frame.loc[frame["sample_type"].eq("Unknown"), "sample_id"] = "patient"
        frame["sample_dilution"] = [1, 1, 1, 1, 10, 10, 10]
        frame["source_extinction"] = np.array([0, 0.5, 0.51, 0.65, 1, 1.02, 1.3]) + 0.2
        result = self.analyze(frame, settings)
        np.testing.assert_array_equal(np.flatnonzero(~well_table(result.samples)["included_in_analysis"]), [3, 6])
        patients = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")].set_index("sample_dilution")
        for dilution, extinctions in ((1, [0.5, 0.51]), (10, [1, 1.02])):
            with self.subTest(dilution=dilution):
                expected = np.mean(FourPLModel().calculate_concentrations(np.array(extinctions), reference_parameters())) * dilution
                self.assertEqual(patients.loc[dilution, "count"], 2)
                self.assertAlmostEqual(patients.loc[dilution, "result_mean"], expected)

    def test_patient_filter_enforces_limits_and_rejects_unusable_values(self) -> None:
        for values, message in (
            ([1, 1, 1.8], "maximaal 1"),
            ([0, 0, 0], "ongelijk aan nul"),
            ([1, np.nan, 1], "eindig"),
        ):
            with self.subTest(values=values):
                frame = self.frame.copy()
                frame.loc[[1, 2, 3], "source_extinction"] = np.array(values) + 0.2
                with self.assertRaisesRegex(ValueError, f"Patiëntmonster 'Shared'.*{message}"):
                    self.analyze(frame)

    def test_instrument_exclusions_do_not_count_as_patient_outliers(self) -> None:
        frame = self.frame.copy()
        frame["disabled"] = False
        frame.loc[[2, 3], "disabled"] = True
        result = self.analyze(frame)
        patient = sample_table(result.samples).set_index("sample_id").loc["a"]
        self.assertEqual(patient["count"], 1)
        self.assertAlmostEqual(patient["result_mean"], 2)
        self.assertEqual(patient["used_wells"], "A2")

    def test_patient_exclusions_are_listed_in_pdf_with_original_readings(self) -> None:
        result = self.analyze()
        with TemporaryDirectory() as directory:
            output = PdfReportWriter(result, analyst="operator").write(Path(directory) / "report.pdf")
            with pymupdf.open(output) as document:
                text = "\n".join(page.get_text() for page in document)
        self.assertIn("Uitgesloten replicaten", text)
        self.assertIn("A4", text)
        self.assertIn("1.500", text)
        self.assertIn("15%", text)

    def test_raw_filter_report_keeps_the_applied_scale_and_original_mean(self):
        rule = replace(self.settings.patients.replicate_filter, response_scale="raw", max_deviation_percent=14)
        result = self.analyze(settings=replace(self.settings, patients=PatientSettings(rule)))
        patient = next(sample for sample in result.samples if sample.sample_id == "a")
        excluded = next(well for well in patient.measurements if not well.included)
        self.assertAlmostEqual(excluded.replicate_mean, np.mean([1.2, 1.22, 1.5]))
        self.assertIn("vóór blancoaftrek", excluded.exclusion_reason)
        with TemporaryDirectory() as directory:
            path = PdfReportWriter(result, analyst="operator").write(Path(directory) / "report.pdf")
            with pymupdf.open(path) as document:
                text = " ".join(" ".join(page.get_text() for page in document).split())
        self.assertIn("vóór blancoaftrek", text)


class ReplicateFilterScaleTests(unittest.TestCase):
    def test_scale_controls_filtering_after_reference_correction_and_before_calculation(self):
        for kind in (SampleType.UNKNOWN, SampleType.CONTROL, SampleType.STANDARD):
            for has_blank in (False, True):
                for scale in ResponseScale:
                    with self.subTest(kind=kind, has_blank=has_blank, scale=scale):
                        records = []
                        for index, response in enumerate((0.3, 0.32, 0.4)):
                            reference = 0.05 * (index + 1)
                            record = {
                                "sample_id": "replicates", "sample_name": "QC", "sample_type": kind.value,
                                "sample_concentration": 1 if kind == SampleType.STANDARD else None, "well": f"A{index + 1}",
                            }
                            records.extend([
                                record | {"wavelength_nm": 450, "source_extinction": response + reference},
                                record | {"wavelength_nm": 620, "source_extinction": reference},
                            ])
                        if has_blank:
                            blank = {"sample_id": "blank", "sample_name": "Blanco", "sample_type": "Blank", "well": "B1"}
                            records.extend([
                                blank | {"wavelength_nm": 450, "source_extinction": 0.25},
                                blank | {"wavelength_nm": 620, "source_extinction": 0.05},
                            ])
                        rule = {"max_deviation_percent": 30, "response_scale": scale.value}
                        group = {"replicate_filter": rule}
                        collection = {SampleType.UNKNOWN: "patients", SampleType.CONTROL: "controls", SampleType.STANDARD: "standards"}[kind]
                        config = {
                            "lot_number": "lot",
                            "measurement": {"reference_wavelength_nm": 620},
                            "calculation": {"method": "formula", "expression": "sample", "response_scale": "blank_corrected"},
                            collection: group if kind == SampleType.UNKNOWN else {"QC" if kind == SampleType.CONTROL else "1": group},
                        }
                        settings = AssaySettings.from_dict("test.skax", config)
                        data = AssayData(samples_from_records(records), assay_metadata())
                        result = AssayAnalyzer(data, settings).analyze()
                        sample = next(sample for sample in result.samples if sample.sample_type == kind)
                        should_exclude = has_blank and scale == ResponseScale.BLANK_CORRECTED
                        self.assertEqual([well.included for well in sample.measurements], [True, True, not should_exclude])
                        inputs = np.array([0.3, 0.32, 0.4])
                        if has_blank and scale == ResponseScale.BLANK_CORRECTED:
                            inputs -= 0.2
                        np.testing.assert_allclose([well.replicate_mean for well in sample.measurements], inputs.mean())
                        corrected = np.array([0.3, 0.32, 0.4]) - (0.2 if has_blank else 0)
                        self.assertAlmostEqual(sample.result.mean, np.mean(corrected[:2] if should_exclude else corrected))


class ReplicateFilteringTests(unittest.TestCase):
    def setUp(self) -> None:
        self.metadata = assay_metadata(session_time="2026-10-09T10:00:00+02:00")
        self.settings = borrelia_settings()
        self.frame = pd.DataFrame(
            {
                "sample_id": ["blank", "negative", "cal", "cal", "cal", "positive", "patient"],
                "sample_name": ["Blank", "Negatief", "Kalibratie", "Kalibratie", "Kalibratie", "Positief", "Kalibratie"],
                "sample_type": ["Blank", "Control", "Control", "Control", "Control", "Control", "Unknown"],
                "well": ["A1", "B1", "C1", "D1", "E1", "F1", "G1"],
                "source_extinction": [0.1, 0.3, 1.1, 1.12, 1.4, 1.6, 1.11],
            }
        )

    def analyze(self, frame=None, settings=None):
        return AssayAnalyzer(AssayData(samples_from_records((self.frame if frame is None else frame).to_dict('records')), self.metadata), settings or self.settings).analyze()

    def test_borrelia_exclusion_updates_cutoff_ratios_qc_and_used_wells(self) -> None:
        original = self.frame.copy(deep=True)
        analyzer = AssayAnalyzer(AssayData(samples_from_records(self.frame.to_dict('records')), self.metadata), self.settings)
        result = analyzer.analyze()

        self.assertAlmostEqual(well_table(result.samples).loc[1, "result"], 0.2 / 1.01)
        self.assertAlmostEqual(well_table(result.samples).loc[5, "result"], 1.5 / 1.01)
        self.assertAlmostEqual(well_table(result.samples).loc[6, "result"], 1 / 0.430)
        self.assertEqual(result.validation_status, ValidationStatus.PASSED)
        calibrator = sample_table(result.samples).loc[
            sample_table(result.samples)["sample_type"].eq("Control") & sample_table(result.samples)["sample_name"].eq("Kalibratie")
        ].iloc[0]
        self.assertEqual(calibrator["count"], 2)
        self.assertEqual(calibrator["used_wells"], "C1, D1")
        self.assertAlmostEqual(calibrator["extinction_mean"], 1.01)
        self.assertAlmostEqual(calibrator["extinction_cv_percent"], 100 * np.std([1, 1.02], ddof=1) / 1.01)
        excluded = well_table(result.samples).loc[4]
        self.assertFalse(excluded["included_in_analysis"])
        self.assertAlmostEqual(excluded["source_extinction"], 1.4)
        self.assertAlmostEqual(excluded["extinction"], 1.3)
        self.assertAlmostEqual(excluded["replicate_mean"], np.mean([1, 1.02, 1.3]))
        self.assertGreater(excluded["replicate_deviation_percent"], 15)
        self.assertIn("15%", excluded["exclusion_reason"])
        self.assertTrue(pd.isna(excluded["result"]))
        pd.testing.assert_frame_equal(self.frame, original)
        pd.testing.assert_frame_equal(well_table(analyzer.analyze().samples), well_table(result.samples))

    def test_filter_uses_blank_corrected_extinction_before_control_overrides(self) -> None:
        frame = self.frame.copy()
        frame["source_extinction"] += 10
        calculation = replace(
            self.settings.calculation,
            control_overrides={"Positief": FormulaOverrideSettings("sample * kalibratie", response_scale="blank_corrected")},
        )
        result = self.analyze(frame, replace(self.settings, calculation=calculation))
        self.assertFalse(well_table(result.samples).loc[4, "included_in_analysis"])
        self.assertAlmostEqual(well_table(result.samples).loc[5, "result"], 1.5 * 1.01)

    def test_exclusion_matches_well_rows_with_duplicate_dataframe_indexes(self) -> None:
        frame = self.frame.copy()
        frame.index = [0, 1, 2, 3, 4, 5, 4]
        result = self.analyze(frame)
        self.assertFalse(well_table(result.samples).iloc[4]["included_in_analysis"])
        self.assertTrue(well_table(result.samples).iloc[6]["included_in_analysis"])
        self.assertAlmostEqual(well_table(result.samples).iloc[6]["result"], 1 / 0.430)
        self.assertEqual(len(sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")]), 1)

    def test_values_at_fifteen_percent_are_retained(self) -> None:
        frame = self.frame.copy()
        frame.loc[[2, 3, 4], "source_extinction"] = [0.95, 1.1, 1.25]
        result = self.analyze(frame)
        self.assertTrue(well_table(result.samples)["included_in_analysis"].all())

    def test_unusable_replicate_groups_stop_analysis(self) -> None:
        for values, rule, message in (
            ([1, 1, 1.8], ReplicateFilterSettings(15, 1, response_scale="blank_corrected"), "maximaal 1"),
            ([1, 1.02, 1.3], ReplicateFilterSettings(15, 0, response_scale="blank_corrected"), "maximaal 0"),
            ([1, 2, 5], ReplicateFilterSettings(15, 3, response_scale="blank_corrected"), "alle metingen"),
            ([0, 0, 0], ReplicateFilterSettings(15, response_scale="blank_corrected"), "ongelijk aan nul"),
            ([1, np.nan, 1], ReplicateFilterSettings(15, response_scale="blank_corrected"), 'eindig'),
        ):
            with self.subTest(values=values, rule=rule):
                frame = self.frame.copy()
                frame.loc[[2, 3, 4], "source_extinction"] = np.array(values) + 0.1
                controls = dict(self.settings.controls)
                controls["Kalibratie"] = replace(controls["Kalibratie"], replicate_filter=rule)
                with self.assertRaisesRegex(ValueError, message):
                    self.analyze(frame, replace(self.settings, controls=controls))
        with self.assertRaisesRegex(ValueError, 'geen metingen'):
            self.analyze(self.frame.drop(index=[2, 3, 4]))

    def test_filter_accepts_any_replicate_count(self) -> None:
        for count in (1, 2, 4, 7):
            with self.subTest(count=count):
                replicates = pd.concat([self.frame.loc[[2]]] * count, ignore_index=True)
                replicates["well"] = [f"C{index + 1}" for index in range(count)]
                frame = pd.concat([self.frame.drop(index=[2, 3, 4]), replicates], ignore_index=True)
                result = self.analyze(frame)
                calibrator = sample_table(result.samples).loc[
                    sample_table(result.samples)["sample_type"].eq("Control") & sample_table(result.samples)["sample_name"].eq("Kalibratie")
                ].iloc[0]
                self.assertEqual(calibrator["count"], count)
                self.assertTrue(well_table(result.samples)["included_in_analysis"].all())

    def test_filter_can_retain_a_single_measurement(self) -> None:
        frame = self.frame.copy()
        frame.loc[[2, 3, 4], "source_extinction"] = np.array([0.1, 1, 1.9]) + 0.1
        controls = dict(self.settings.controls)
        controls["Kalibratie"] = replace(controls["Kalibratie"], replicate_filter=ReplicateFilterSettings(15, 2, response_scale="blank_corrected"))
        result = self.analyze(frame, replace(self.settings, controls=controls))
        calibrator = sample_table(result.samples).loc[
            sample_table(result.samples)["sample_type"].eq("Control") & sample_table(result.samples)["sample_name"].eq("Kalibratie")
        ].iloc[0]
        self.assertEqual(calibrator["count"], 1)
        self.assertEqual(calibrator["used_wells"], "D1")
        self.assertAlmostEqual(calibrator["extinction_mean"], 1)

    def test_an_assay_without_a_rule_keeps_all_replicates(self) -> None:
        controls = dict(self.settings.controls)
        controls["Kalibratie"] = replace(controls["Kalibratie"], replicate_filter=None)
        result = self.analyze(settings=replace(self.settings, controls=controls))
        self.assertTrue(well_table(result.samples)["included_in_analysis"].all())

    def test_standard_filtering_removes_outlier_from_curve_fit_and_recovery(self) -> None:
        model = FourPLModel()
        parameters = reference_parameters()
        nominal = np.repeat([0, 1, 5, 25, 100], 3)
        extinction = np.asarray(model.calculate_extinctions(nominal, parameters))
        extinction[5] *= 1.3
        frame = pd.DataFrame(
            {
                "sample_id": [f"standard-{x}" for x in nominal],
                "sample_name": "Patient",
                "sample_type": "Standard",
                "sample_concentration": nominal,
                "source_extinction": extinction,
                "well": [f"A{x + 1}" for x in range(len(nominal))],
            }
        )
        settings = AssaySettings.from_dict(
            "test.skax",
            {
                "lot_number": "lot",
                "standards": {"1": {"replicate_filter": {"response_scale": "blank_corrected", "max_deviation_percent": 15}}},
                "calculation": {"response_scale": "raw", "method": "calibration", "fit_model": "4pl"},
                "result": {"label": "Conc"},
            },
        )
        result = self.analyze(frame, settings)
        standard = sample_table(result.samples).loc[sample_table(result.samples)["sample_concentration"].eq(1)].iloc[0]
        self.assertEqual(standard["count"], 2)
        self.assertEqual(standard["used_wells"], "A4, A5")
        self.assertAlmostEqual(standard["recovery_percent"], 100, places=4)
        self.assertTrue(pd.isna(well_table(result.samples).loc[5, "result"]))
        self.assertTrue(pd.isna(well_table(result.samples).loc[5, "back_calculated_concentration"]))
        self.assertEqual(len(result.calculation_details.curve.diagnostics.back_calculated), len(nominal) - 1)

    def test_pdf_preserves_raw_reading_and_identifies_exclusion(self) -> None:
        result = self.analyze()
        with TemporaryDirectory() as directory:
            output = PdfReportWriter(result, analyst="operator").write(Path(directory) / "report.pdf")
            with pymupdf.open(output) as document:
                text = "\n".join(page.get_text() for page in document)
        self.assertIn("1.400", text)
        self.assertIn("Uitgesloten replicaten", text)
        self.assertIn("E1", text)
        self.assertIn("17.5", text)
        self.assertIn("15%", text)

    def test_qc_cv_uses_only_retained_wells_across_sample_ids(self) -> None:
        frame = self.frame.copy()
        frame.loc[[2, 3, 4], "sample_id"] = ["cal-1", "cal-2", "cal-2"]
        result = self.analyze(frame)
        check = next(
            check for check in result.validation_results
            if check.sample_key == (SampleType.CONTROL, "Kalibratie") and check.metric == ValidationMetric.EXTINCTION
        )
        self.assertAlmostEqual(check.cv_percent, 100 * np.std([1, 1.02], ddof=1) / 1.01)


if __name__ == "__main__":
    unittest.main()
