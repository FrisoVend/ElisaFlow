import json
import stat
import unittest
from dataclasses import replace
from datetime import datetime
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Event
from unittest.mock import patch

from skax_fixtures import replace_member, write_skax

from elisaflow.assays.registry import AssayRegistry
from elisaflow.configuration.application import ApplicationSettings
from elisaflow.configuration.assay_settings import AssaySettings, ReportSettings
from elisaflow.readers.assay_loader import AssayLoader
from elisaflow.runs.experiments import load_experiment, prepare_live_experiment
from elisaflow.runs.workflow import RunCancelled, _file_is_stable, run_workflow, wait_for_measurement


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.config = ApplicationSettings(
            self.directory / "protocols", self.directory / "configs", self.directory / "results", self.directory / "editor.hash",
        )
        self.config.protocols_dir.mkdir()
        self.config.configs_dir.mkdir()
        self.template = write_skax(self.config.protocols_dir / "session.skax", completed=False)
        self.settings = AssaySettings.from_dict(self.template.name, {
            "lot_number": "lot", "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / cutoff", "references": {"cutoff": "Cutoff"}},
        })
        self.settings_path = self.config.configs_dir / "session.json"
        self.settings_path.write_text(json.dumps({self.template.name: self.settings.to_dict()}), encoding="utf-8")


    @staticmethod
    def complete_measurement(session: Path, **kwargs):
        return write_skax(session)

    def test_live_run_waits_for_the_copied_skax_and_archives_every_input(self):
        original = self.template.read_bytes()
        with (
            patch("elisaflow.runs.workflow.open_protocol") as launch,
            patch("elisaflow.runs.workflow.wait_for_measurement", side_effect=self.complete_measurement) as wait,
        ):
            report = run_workflow(self.config, protocol_template=self.template, open_report=False, analyst="operator")
        saved = load_experiment(report.parent)
        self.assertEqual(self.template.read_bytes(), original)
        self.assertEqual(saved.settings, self.settings)
        self.assertEqual(saved.session.name, f"{report.parent.name}_session.skax")
        launch.assert_called_once_with(saved.session)
        wait.assert_called_once_with(saved.session, cancel_event=None)
        self.assertTrue(report.read_bytes().startswith(b"%PDF"))
        self.assertEqual({path.suffix for path in report.parent.iterdir()}, {".skax", ".json", ".pdf"})

    def test_selector_enters_the_live_workflow(self):
        with (
            patch("elisaflow.runs.workflow.select_protocol", return_value=self.template) as select,
            patch("elisaflow.runs.workflow.open_protocol"),
            patch("elisaflow.runs.workflow.wait_for_measurement", side_effect=self.complete_measurement),
        ):
            report = run_workflow(self.config, open_report=False, analyst="operator")
        select.assert_called_once_with(self.config.protocols_dir)
        self.assertTrue(report.is_file())

    def test_each_run_uses_the_entered_name_for_its_files_and_report(self):
        for analyst in ("First Analyst", "Second Analyst"):
            with (
                self.subTest(analyst=analyst),
                patch("elisaflow.runs.workflow.open_protocol"),
                patch("elisaflow.runs.workflow.wait_for_measurement", side_effect=self.complete_measurement),
                patch("elisaflow.reporting.pdf_report.PdfReportWriter") as writer,
            ):
                writer.return_value.write.side_effect = lambda path: path
                report = run_workflow(self.config, protocol_template=self.template, open_report=False, analyst=analyst)
                self.assertEqual(writer.call_args.kwargs["analyst"], analyst)
                self.assertTrue(report.parent.name.endswith(f"_{analyst}"))

    def test_custom_pdf_destination_keeps_the_run_folder_reanalyzable(self):
        output = self.directory / "custom" / "comparison.pdf"
        with patch("elisaflow.runs.workflow.open_protocol") as launch, patch("elisaflow.runs.workflow.wait_for_measurement", side_effect=self.complete_measurement):
            self.assertEqual(run_workflow(self.config, protocol_template=self.template, output=output, open_report=False, analyst="operator"), output)
        self.assertEqual(load_experiment(launch.call_args.args[0].parent).settings, self.settings)
        self.assertEqual(AssayRegistry(output.parent).get(self.template), self.settings)

    def test_cancellation_before_start_creates_no_files(self):
        cancellation = Event()
        cancellation.set()
        with self.assertRaises(RunCancelled):
            run_workflow(self.config, protocol_template=self.template, cancel_event=cancellation, analyst="operator")
        self.assertFalse(self.config.results_dir.exists())

    def test_cancelling_selection_creates_no_run(self):
        with patch("elisaflow.runs.workflow.select_protocol", side_effect=RuntimeError("Cancelled")), self.assertRaisesRegex(RuntimeError, "Cancelled"):
            run_workflow(self.config, analyst="operator")
        self.assertFalse(self.config.results_dir.exists())

    def test_live_wait_can_be_cancelled_and_keeps_its_snapshots(self):
        cancellation = Event()
        with patch("elisaflow.runs.workflow.open_protocol", side_effect=lambda _: cancellation.set()), self.assertRaises(RunCancelled):
            run_workflow(self.config, protocol_template=self.template, cancel_event=cancellation, analyst="operator")
        snapshot, = self.config.results_dir.rglob("*_settings.json")
        self.assertEqual(AssayRegistry(snapshot.parent).get(self.template), self.settings)
        self.assertTrue(list(snapshot.parent.glob("*.skax")))
        self.assertFalse(list(snapshot.parent.glob("*.pdf")))

    def test_wait_cancels_without_waiting_for_the_timeout(self):
        cancellation = Event()
        with patch.object(cancellation, "wait", side_effect=lambda _: cancellation.set()), self.assertRaises(RunCancelled):
            wait_for_measurement(self.template, timeout=3600, poll_interval=0, cancel_event=cancellation)

    def test_cancellation_after_loading_prevents_report_generation(self):
        cancellation = Event()
        original_load = AssayLoader.load
        def load_and_cancel(*args, **kwargs):
            data = original_load(*args, **kwargs)
            cancellation.set()
            return data
        with (
            patch("elisaflow.runs.workflow.open_protocol"), patch("elisaflow.runs.workflow.wait_for_measurement", side_effect=self.complete_measurement),
            patch("elisaflow.runs.workflow.AssayLoader.load", side_effect=load_and_cancel), patch("elisaflow.runs.workflow.create_report") as report,
            self.assertRaises(RunCancelled),
        ):
            run_workflow(self.config, protocol_template=self.template, cancel_event=cancellation, analyst="operator")
        report.assert_not_called()

    def test_library_changes_during_measurement_do_not_change_the_run_settings(self):
        def complete_after_edit(session, **kwargs):
            self.settings_path.chmod(self.settings_path.stat().st_mode | stat.S_IWUSR)
            self.settings_path.write_text(json.dumps({self.template.name: self.settings.to_dict() | {"lot_number": "new lot"}}), encoding="utf-8")
            return self.complete_measurement(session)
        with patch("elisaflow.runs.workflow.open_protocol"), patch("elisaflow.runs.workflow.wait_for_measurement", side_effect=complete_after_edit):
            report = run_workflow(self.config, protocol_template=self.template, open_report=False, analyst="operator")
        self.assertEqual(load_experiment(report.parent).settings, self.settings)
        self.assertEqual(AssayRegistry(self.config.configs_dir).get(self.template).lot_number, "new lot")

    def test_only_a_stable_completed_archive_with_valid_hashes_is_ready(self):
        with patch("elisaflow.runs.workflow.time.sleep"):
            self.assertFalse(_file_is_stable(self.template))
            write_skax(self.template)
            self.assertTrue(_file_is_stable(self.template))
            replace_member(self.template, "Measurements.data", b"corrupt", update_hashes=False)
            self.assertFalse(_file_is_stable(self.template))
            self.template.write_bytes(b"not a zip archive")
            self.assertFalse(_file_is_stable(self.template))
        with self.assertRaisesRegex(TimeoutError, "Sla de gemeten sessie"):
            wait_for_measurement(self.template, timeout=0)

    def test_run_directories_use_iso_year_and_week_at_year_boundaries(self) -> None:
        for started, year, week in (
            (datetime(2026, 10, 6, 10, 30), "2026", "week_41"),
            (datetime(2024, 12, 30, 10, 30), "2025", "week_01"),
            (datetime(2021, 1, 1, 10, 30), "2020", "week_53"),
        ):
            with self.subTest(started=started), patch("elisaflow.runs.experiments.datetime") as clock:
                clock.now.return_value = started
                directory = prepare_live_experiment(self.template, self.settings, self.config.results_dir, analyst="operator").directory
                self.assertEqual(directory.parent, self.config.results_dir / year / week)
                self.assertEqual(directory.name, f"session_{started:%Y-%m-%d_%H-%M-%S}_operator")
                self.assertTrue(directory.is_dir())


    def test_same_second_runs_keep_distinct_folders_without_overwriting(self) -> None:
        with patch("elisaflow.runs.experiments.datetime") as clock:
            clock.now.return_value = datetime(2026, 10, 22, 14, 6, 12)
            first = prepare_live_experiment(self.template, self.settings, self.config.results_dir, analyst="operator").directory
            marker = first / "keep.txt"
            marker.write_text("previous run", encoding="utf-8")
            second = prepare_live_experiment(self.template, self.settings, self.config.results_dir, analyst="operator").directory
        self.assertEqual(first.name, "session_2026-10-22_14-06-12_operator")
        self.assertEqual(second.name, f"{first.name}_2")
        self.assertEqual(marker.read_text(encoding="utf-8"), "previous run")


    def test_run_directories_preserve_names_and_sanitize_filename_characters(self):
        for index, (analyst, expected) in enumerate((
            ("Friso van Dijk", "Friso van Dijk"),
            ("Friso / van Dijk", "Friso _ van Dijk"),
            (r"LAB\operator:name?", "LAB_operator_name_"),
            ("  Operator  ", "Operator"),
        )):
            with self.subTest(analyst=analyst):
                directory = prepare_live_experiment(self.template, self.settings, self.config.results_dir / str(index), analyst=analyst).directory
                self.assertTrue(directory.name.endswith(f"_{expected}"))


    def test_existing_snapshots_are_preserved_and_reports_have_distinct_snapshots(self) -> None:
        settings = AssayRegistry(self.config.configs_dir).get(self.template)
        output = self.directory / "custom" / "first.pdf"
        prepare_live_experiment(self.template, settings, self.config.results_dir, analyst="operator", output=output)
        snapshot = output.parent / "first_settings.json"
        original_bytes = snapshot.read_bytes()
        original_time = snapshot.stat().st_mtime_ns
        prepare_live_experiment(self.template, settings, self.config.results_dir, analyst="operator", output=output)
        self.assertEqual(snapshot.stat().st_mtime_ns, original_time)

        changed = AssaySettings.from_dict(self.template.name, settings.to_dict() | {"lot_number": "new lot"})
        with self.assertRaisesRegex(FileExistsError, 'andere instellingen'):
            prepare_live_experiment(self.template, changed, self.config.results_dir, analyst="operator", output=output)
        self.assertEqual(snapshot.read_bytes(), original_bytes)
        prepare_live_experiment(self.template, changed, self.config.results_dir, analyst="operator", output=output.with_name("second.pdf"))
        second = output.parent / "second_settings.json"
        self.assertNotEqual(snapshot, second)
        self.assertEqual(
            AssaySettings.from_dict(self.template.name, json.loads(second.read_text(encoding="utf-8"))[self.template.name]),
            changed,
        )


    def test_snapshot_preserves_effective_report_options_and_rejects_changed_options(self) -> None:
        defaults = AssayRegistry(self.config.configs_dir).get(self.template)
        output = self.directory / "report.pdf"
        prepare_live_experiment(self.template, defaults, self.config.results_dir, analyst="operator", output=output)
        snapshot = output.parent / "report_settings.json"
        contents = json.loads(snapshot.read_text(encoding="utf-8"))[self.template.name]
        self.assertEqual(contents["report"], ReportSettings().to_dict())
        custom = replace(defaults, report=ReportSettings(calibration_decimals=3, cv_warning_percent=12.5, recovery_warning_percent=15))
        with self.assertRaisesRegex(FileExistsError, 'andere instellingen'):
            prepare_live_experiment(self.template, custom, self.config.results_dir, analyst="operator", output=output)
        prepare_live_experiment(self.template, custom, self.config.results_dir, analyst="operator", output=output.with_name("custom.pdf"))
        snapshot = output.parent / "custom_settings.json"
        contents = json.loads(snapshot.read_text(encoding="utf-8"))[self.template.name]
        self.assertEqual(AssaySettings.from_dict(self.template.name, contents), custom)


    def test_settings_snapshots_normalize_newlines(self) -> None:
        settings = AssayRegistry(self.config.configs_dir).get(self.template)
        output = self.directory / "snapshots" / "report.pdf"
        prepare_live_experiment(self.template, settings, self.config.results_dir, analyst="operator", output=output)
        snapshot = output.parent / "report_settings.json"
        snapshot.write_bytes(snapshot.read_bytes().replace(b"\n", b"\r\n"))
        original_time = snapshot.stat().st_mtime_ns
        prepare_live_experiment(self.template, settings, self.config.results_dir, analyst="operator", output=output)
        self.assertEqual(snapshot.stat().st_mtime_ns, original_time)
