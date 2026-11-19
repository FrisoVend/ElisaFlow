import json
import unittest
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Event
from unittest.mock import patch
from zipfile import BadZipFile

from skax_fixtures import basic_records, write_skax

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.assays.registry import AssayRegistry
from elisaflow.configuration.application import ApplicationSettings
from elisaflow.configuration.assay_settings import AssaySettings
from elisaflow.runs.experiments import load_experiment, prepare_live_experiment, prepare_reanalysis
from elisaflow.runs.workflow import RunCancelled, create_report, run_workflow


class SavedExperimentTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.folder = self.directory / "saved experiment"
        self.folder.mkdir()
        self.config = ApplicationSettings(
            self.directory / "missing-protocols", self.directory / "configs", self.directory / "results", self.directory / "editor.hash",
        )
        self.config.configs_dir.mkdir()
        (self.config.configs_dir / "current.json").write_text("invalid JSON", encoding="utf-8")
        self.settings = AssaySettings.from_dict("Original.skax", {
            "lot_number": "original-lot", "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / cutoff", "references": {"cutoff": "Cutoff"}},
            "result": {"label": "Ratio"}, "validation": {"result": {"controls": {"Cutoff": {"result": {"min": 0.9, "max": 1.1}}}}},
        })
        self.session = write_skax(self.folder / "archived_session.skax", name="Original")
        self.snapshot = self.folder / "archived_settings.json"
        self.snapshot.write_text(json.dumps({self.settings.protocol_name: self.settings.to_dict()}), encoding="utf-8")
        (self.folder / "archived_report.pdf").write_bytes(b"%PDF-original")
        for target in ("select_protocol", "open_protocol", "wait_for_measurement", "AssayRegistry"):
            mocked = patch(f"elisaflow.runs.workflow.{target}", side_effect=AssertionError(f"Replay called {target}"))
            mocked.start()
            self.addCleanup(mocked.stop)

    def source_contents(self) -> dict:
        return {path.name: (path.read_bytes(), path.stat().st_mtime_ns) for path in self.folder.iterdir()}

    def test_preparation_preserves_inputs_and_returns_paths_before_analysis(self):
        original = self.source_contents()
        output = self.directory / "prepared" / "comparison.pdf"
        prepared = prepare_reanalysis(load_experiment(self.folder), self.config.results_dir, analyst="operator", output=output)
        self.assertEqual(prepared.directory, output.parent)
        self.assertEqual(prepared.report, output)
        self.assertEqual(prepared.settings, self.settings)
        self.assertFalse(prepared.report.exists())
        self.assertEqual(load_experiment(prepared.directory).settings, self.settings)
        self.assertEqual(prepared.session.read_bytes(), self.session.read_bytes())
        self.assertEqual(self.source_contents(), original)

    def test_custom_live_preparation_keeps_session_in_run_folder_and_snapshots_beside_both_paths(self):
        output = self.directory / "custom-live" / "comparison.pdf"
        prepared = prepare_live_experiment(self.session, self.settings, self.config.results_dir, analyst="operator", output=output)
        self.assertEqual(prepared.report, output)
        self.assertEqual(load_experiment(prepared.directory).settings, self.settings)
        self.assertEqual(prepared.session.parent, prepared.directory)
        self.assertEqual(list(output.parent.glob("*.skax")), [])
        snapshot, = output.parent.glob("*_settings.json")
        self.assertEqual(json.loads(snapshot.read_text(encoding="utf-8")), {self.settings.protocol_name: self.settings.to_dict()})

    def test_live_run_produces_a_complete_experiment_that_can_be_reanalyzed(self):
        original = self.source_contents()
        templates, configs = self.directory / "templates", self.directory / "live-configs"
        templates.mkdir()
        configs.mkdir()
        template = write_skax(templates / self.settings.protocol_name, completed=False)
        (configs / "assay.json").write_text(json.dumps({template.name: self.settings.to_dict()}), encoding="utf-8")
        for output in (None, self.directory / "custom" / "report.pdf"):
            with self.subTest(output=output):
                with (
                    patch("elisaflow.runs.workflow.AssayRegistry", AssayRegistry), patch("elisaflow.runs.workflow.open_protocol") as launch,
                    patch("elisaflow.runs.workflow.wait_for_measurement", side_effect=lambda session, **kwargs: write_skax(session, name="Original")),
                ):
                    live = run_workflow(replace(self.config, protocols_dir=templates, configs_dir=configs), protocol_template=template, output=output, open_report=False, analyst="operator")
                launch.assert_called_once()
                live_folder = launch.call_args.args[0].parent
                saved = load_experiment(live_folder)
                self.assertEqual(saved.settings, self.settings)
                self.assertEqual(saved.session.suffix, ".skax")
                self.assertFalse(list(live_folder.glob("*.xml")))
                rerun = run_workflow(self.config, experiment_dir=live_folder, open_report=False, analyst="operator")
                self.assertNotEqual(rerun.parent, live_folder)
                self.assertEqual(load_experiment(rerun.parent).settings, self.settings)
                self.assertTrue(live.read_bytes().startswith(b"%PDF-"))
                self.assertTrue(rerun.read_bytes().startswith(b"%PDF-"))
                self.assertEqual(len(list(rerun.parent.iterdir())), 3)
        self.assertEqual(self.source_contents(), original)

    def test_reanalysis_uses_archived_settings_and_real_inputs_without_modifying_the_source(self) -> None:
        original = self.source_contents()
        with patch("elisaflow.runs.workflow.create_report", wraps=create_report) as report_writer:
            report = run_workflow(self.config, experiment_dir=self.folder, open_report=False, analyst="operator")
        data, settings, output = report_writer.call_args.args
        self.assertEqual(settings, self.settings)
        result = AssayAnalyzer(data, settings).analyze()
        patient = next(sample for sample in result.samples if sample.name == "Patient")
        self.assertAlmostEqual(patient.result.mean, 2)
        self.assertEqual(result.validation_status.value, "passed")
        self.assertEqual(report, output)
        self.assertTrue(report.read_bytes().startswith(b"%PDF-"))
        self.assertTrue(report.parent.is_relative_to(self.config.results_dir))
        self.assertTrue(report.parent.name.startswith("Original_"))
        self.assertEqual(self.source_contents(), original)

        saved = load_experiment(report.parent)
        self.assertEqual(saved.session.read_bytes(), self.session.read_bytes())
        self.assertEqual(saved.settings, self.settings)
        self.assertEqual(
            {path.name for path in report.parent.iterdir()},
            {report.name, saved.session.name, f"{report.parent.name}_settings.json"},
        )

    def test_custom_output_is_self_contained_and_can_be_reanalyzed(self) -> None:
        output = self.directory / "comparison" / "new_report.pdf"
        with patch("elisaflow.runs.workflow.create_report", side_effect=lambda _data, _settings, path, **kwargs: path):
            with patch("elisaflow.runs.workflow.open_document") as open_report:
                self.assertEqual(run_workflow(self.config, experiment_dir=self.folder, output=output, analyst="operator"), output)
            open_report.assert_called_once_with(output)
            saved = load_experiment(output.parent)
            self.assertEqual(saved.settings, self.settings)
            replay = run_workflow(self.config, experiment_dir=output.parent, open_report=False, analyst="operator")
        self.assertEqual(load_experiment(replay.parent).settings, self.settings)
        self.assertTrue(replay.parent.name.startswith("Original_"))

    def test_missing_or_ambiguous_inputs_fail_before_creating_output(self) -> None:
        for source in (self.session, self.snapshot):
            contents = source.read_bytes()
            with self.subTest(missing=source.name):
                source.unlink()
                with self.assertRaisesRegex(ValueError, "precies één"):
                    run_workflow(self.config, experiment_dir=self.folder, open_report=False, analyst="operator")
                source.write_bytes(contents)
            with self.subTest(ambiguous=source.name):
                extra = source.with_name(f"extra_{source.name}")
                extra.write_bytes(contents)
                with self.assertRaisesRegex(ValueError, "precies één"):
                    run_workflow(self.config, experiment_dir=self.folder, open_report=False, analyst="operator")
                extra.unlink()
        self.assertFalse(self.config.results_dir.exists())

    def test_invalid_or_multiple_assay_settings_fail_without_using_the_library(self) -> None:
        invalid_settings = dict(self.settings.to_dict(), calculation={"method": "unknown"})
        for snapshot in ("not JSON", "[]", "{}", '{"first": {}, "second": {}}', json.dumps({"Original.skax": invalid_settings})):
            with self.subTest(snapshot=snapshot):
                self.snapshot.write_text(snapshot, encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "instellingen|precies één assay"):
                    run_workflow(self.config, experiment_dir=self.folder, open_report=False, analyst="operator")
        self.assertFalse(self.config.results_dir.exists())

    def test_settings_snapshot_accepts_utf8_bom(self) -> None:
        self.snapshot.write_bytes(b"\xef\xbb\xbf" + self.snapshot.read_bytes())
        self.assertEqual(load_experiment(self.folder).settings, self.settings)

    def test_missing_experiment_directory_is_reported(self) -> None:
        with self.assertRaisesRegex(FileNotFoundError, "experimentmap"):
            run_workflow(self.config, experiment_dir=self.folder / "missing", open_report=False, analyst="operator")
        self.assertFalse(self.config.results_dir.exists())

    def test_live_protocol_cannot_be_combined_with_an_experiment(self) -> None:
        with self.assertRaisesRegex(ValueError, "Kies een experimentmap"):
            run_workflow(self.config, experiment_dir=self.folder, protocol_template=self.session, analyst="operator")
        self.assertFalse(self.config.results_dir.exists())

    def test_output_cannot_be_written_inside_the_original_experiment(self) -> None:
        original = self.source_contents()
        for output in (self.folder / "archived_report.pdf", self.folder / "new.pdf", self.folder / "nested" / "new.pdf"):
            with self.subTest(output=output), self.assertRaisesRegex(ValueError, "buiten de oorspronkelijke"):
                run_workflow(self.config, experiment_dir=self.folder, output=output, open_report=False, analyst="operator")
        for root in (self.folder, self.folder / "nested"):
            with self.subTest(results=root), self.assertRaisesRegex(ValueError, "buiten de oorspronkelijke"):
                run_workflow(replace(self.config, results_dir=root), experiment_dir=self.folder, open_report=False, analyst="operator")
        self.assertEqual(self.source_contents(), original)

    def test_existing_output_artifacts_are_never_overwritten(self) -> None:
        output = self.directory / "comparison" / "new_report.pdf"
        output.parent.mkdir()
        for suffix in ("report.pdf", "session.skax", "settings.json"):
            existing = output.with_name(f"new_{suffix}")
            existing.write_bytes(b"keep existing")
            with self.subTest(suffix=suffix), self.assertRaisesRegex(FileExistsError, "nieuw uitvoerpad"):
                run_workflow(self.config, experiment_dir=self.folder, output=output, open_report=False, analyst="operator")
            self.assertEqual(list(output.parent.iterdir()), [existing])
            self.assertEqual(existing.read_bytes(), b"keep existing")
            existing.unlink()

    def test_cancelled_reanalysis_creates_no_output(self) -> None:
        cancellation = Event()
        cancellation.set()
        with self.assertRaises(RunCancelled):
            run_workflow(self.config, experiment_dir=self.folder, cancel_event=cancellation, analyst="operator")
        self.assertFalse(self.config.results_dir.exists())

    def test_failed_session_loading_preserves_source_and_copied_inputs(self):
        self.session.write_bytes(b"invalid archive")
        original = self.source_contents()
        with self.assertRaises(BadZipFile):
            run_workflow(self.config, experiment_dir=self.folder, open_report=False, analyst="operator")
        snapshot, = self.config.results_dir.rglob("*_settings.json")
        saved = load_experiment(snapshot.parent)
        self.assertEqual(saved.settings, self.settings)
        self.assertEqual(saved.session.read_bytes(), self.session.read_bytes())
        self.assertFalse(list(saved.directory.glob("*.pdf")))
        self.assertEqual(self.source_contents(), original)

    def test_missing_qc_record_stops_reanalysis_before_report_creation(self):
        layout = basic_records()
        write_skax(self.session, layout[1:], layout=layout)
        original = self.source_contents()
        with patch("elisaflow.runs.workflow.create_report") as report, self.assertRaisesRegex(ValueError, "A1 bij 450 nm"):
            run_workflow(self.config, experiment_dir=self.folder, open_report=False, analyst="operator")
        report.assert_not_called()
        self.assertFalse(list(self.config.results_dir.rglob("*.pdf")))
        self.assertEqual(self.source_contents(), original)

if __name__ == "__main__":
    unittest.main()
