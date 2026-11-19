from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Event, get_ident
from unittest.mock import Mock, patch

from skax_fixtures import write_skax

from elisaflow.assays.access_control import AccessControl
from elisaflow.assays.management import AssayManagement
from elisaflow.configuration.application import ApplicationSettings
from elisaflow.runs.experiments import load_experiment
from elisaflow.runs.workflow import RunCancelled
from elisaflow.ui.desktop import ProtocolSelectionCancelled

try:
    from PySide6.QtCore import Qt, QTimer
    from PySide6.QtGui import QCloseEvent
    from PySide6.QtTest import QTest
    from PySide6.QtWidgets import QApplication

    from elisaflow.ui.main_window import MainWindow
except ModuleNotFoundError as error:
    if error.name != "PySide6":
        raise
    QApplication = None


@unittest.skipIf(QApplication is None, "Install the optional PySide6 UI dependency")
class MainWindowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        environment = patch.dict(os.environ, {"QT_QPA_PLATFORM": "offscreen"})
        environment.start()
        cls.addClassCleanup(environment.stop)
        cls.application = QApplication.instance() or QApplication([])

    def setUp(self) -> None:
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        directory = Path(temporary.name).resolve()
        self.config = ApplicationSettings(*(directory / name for name in ("protocols", "configs", "results", "editor.hash")))
        self.assays = Mock(spec=AssayManagement)
        self.assays.list_assays.return_value = ("assay.skax",)
        self.access = Mock(spec=AccessControl)
        self.access.can_edit_assays.return_value = False
        name_prompt = patch("elisaflow.ui.main_window.prompt_analyst_name", return_value="operator")
        self.name_prompt = name_prompt.start()
        self.addCleanup(name_prompt.stop)
        self.window = MainWindow(self.config, self.assays, self.access)
        self.release = Event()
        self.addCleanup(self.clean_up_window)
        self.protocol = self.config.protocols_dir / "assay.skax"
        self.report = self.config.results_dir / "report.pdf"

    def clean_up_window(self) -> None:
        self.release.set()
        worker = self.window._worker
        if worker is not None:
            worker.cancel()
            worker.wait(5000)
        self.application.processEvents()
        self.window.deleteLater()

    def wait_until(self, predicate, timeout: float = 5) -> None:
        deadline = time.monotonic() + timeout
        while not predicate() and time.monotonic() < deadline:
            self.application.processEvents()
            time.sleep(0.01)
        self.assertTrue(predicate(), "Desktop operation did not complete")

    def saved_experiment(self) -> Path:
        folder = self.config.results_dir / "2026" / "week_45" / "Opgeslagen meting"
        folder.mkdir(parents=True)
        write_skax(folder / "original_session.skax", name="Original")
        settings = {
            "Original.skax": {
                "lot_number": "saved-lot",
                "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / cutoff", "references": {"cutoff": "Cutoff"}},
                "result": {"label": "Ratio"},
            }
        }
        (folder / "original_settings.json").write_text(json.dumps(settings), encoding="utf-8")
        (folder / "original_report.pdf").write_bytes(b"%PDF-original")
        return folder

    def test_reanalysis_selects_from_configured_results_and_preserves_the_saved_experiment(self):
        folder = self.saved_experiment()
        original = {path.name: (path.read_bytes(), path.stat().st_mtime_ns) for path in folder.iterdir()}
        with (
            patch("elisaflow.ui.main_window.QFileDialog.getExistingDirectory", return_value=str(folder)) as picker,
            patch("elisaflow.ui.main_window.select_protocol", side_effect=AssertionError("Live protocol selection")),
            patch("elisaflow.runs.workflow.open_protocol", side_effect=AssertionError("SkanIt launched")),
            patch("elisaflow.runs.workflow.wait_for_measurement", side_effect=AssertionError("Waiting for measurement")),
            patch("elisaflow.runs.workflow.AssayRegistry", side_effect=AssertionError("Current assay library used")),
            patch("elisaflow.ui.main_window.open_document") as open_report,
            patch("elisaflow.ui.main_window.QMessageBox.warning") as warning,
        ):
            self.window.reanalyze_button.click()
            self.assertFalse(self.window.start_button.isEnabled())
            self.assertFalse(self.window.reanalyze_button.isEnabled())
            self.wait_until(lambda: self.window._worker is None, timeout=20)
        picker.assert_called_once_with(self.window, "Kies een experimentmap voor heranalyse", str(self.config.results_dir))
        self.name_prompt.assert_called_once_with(parent=self.window)
        warning.assert_not_called()
        report = self.window._report
        self.assertTrue(report.read_bytes().startswith(b"%PDF-"))
        self.assertTrue(report.is_relative_to(self.config.results_dir))
        self.assertFalse(report.is_relative_to(folder))
        self.assertEqual(load_experiment(report.parent).settings.lot_number, "saved-lot")
        self.assertEqual(original, {path.name: (path.read_bytes(), path.stat().st_mtime_ns) for path in folder.iterdir()})
        open_report.assert_called_once_with(report)
        self.assertTrue(self.window.start_button.isEnabled())
        self.assertTrue(self.window.reanalyze_button.isEnabled())

    def test_cancelled_reanalysis_folder_selection_keeps_the_previous_report(self):
        self.config.results_dir.mkdir()
        self.window._report = self.report
        with (
            patch("elisaflow.ui.main_window.QFileDialog.getExistingDirectory", return_value=""),
            patch("elisaflow.ui.main_window.RunWorker") as worker,
            patch("elisaflow.ui.main_window.QMessageBox.warning") as warning,
        ):
            self.window.reanalyze_button.click()
        self.name_prompt.assert_not_called()
        worker.assert_not_called()
        warning.assert_not_called()
        self.assertEqual(self.window._report, self.report)
        self.assertEqual(list(self.config.results_dir.iterdir()), [])

    def test_invalid_reanalysis_folder_is_reported_before_asking_for_a_name(self):
        self.config.results_dir.mkdir()
        self.window._report = self.report
        with (
            patch("elisaflow.ui.main_window.QFileDialog.getExistingDirectory", return_value=str(self.config.results_dir)),
            patch("elisaflow.ui.main_window.RunWorker") as worker,
            patch("elisaflow.ui.main_window.QMessageBox.warning") as warning,
        ):
            self.window.reanalyze_button.click()
        worker.assert_not_called()
        self.name_prompt.assert_not_called()
        warning.assert_called_once()
        self.assertEqual(warning.call_args.args[1], "Heranalyse kon niet gestart worden")
        self.assertIn("precies", warning.call_args.args[2])
        self.assertEqual(self.window._report, self.report)
        self.assertTrue(self.window.reanalyze_button.isEnabled())

    def test_missing_results_directory_is_reported_without_opening_the_picker(self):
        with (
            patch("elisaflow.ui.main_window.QFileDialog.getExistingDirectory") as picker,
            patch("elisaflow.ui.main_window.QMessageBox.warning") as warning,
        ):
            self.window.reanalyze_button.click()
        picker.assert_not_called()
        self.name_prompt.assert_not_called()
        warning.assert_called_once()
        self.assertIn("resultatenmap bestaat niet", warning.call_args.args[2])
        self.assertFalse(self.config.results_dir.exists())

    def test_cancelling_the_reanalysis_name_does_not_start_processing(self):
        folder = self.saved_experiment()
        original = set(self.config.results_dir.rglob("*"))
        self.name_prompt.return_value = None
        self.window._report = self.report
        with (
            patch("elisaflow.ui.main_window.QFileDialog.getExistingDirectory", return_value=str(folder)),
            patch("elisaflow.ui.main_window.RunWorker") as worker,
        ):
            self.window.reanalyze_button.click()
        worker.assert_not_called()
        self.assertEqual(self.window._report, self.report)
        self.assertEqual(set(self.config.results_dir.rglob("*")), original)

    def test_running_reanalysis_can_be_cancelled_and_blocks_a_second_run(self):
        folder = self.saved_experiment()
        started = Event()

        def process(_config, **kwargs):
            self.assertEqual(kwargs["experiment_dir"], folder)
            self.assertIsNone(kwargs["protocol_template"])
            self.assertEqual(kwargs["analyst"], "operator")
            started.set()
            kwargs["cancel_event"].wait(5)
            raise RunCancelled()

        with (
            patch("elisaflow.ui.main_window.QFileDialog.getExistingDirectory", return_value=str(folder)) as picker,
            patch("elisaflow.ui.run_worker.run_workflow", side_effect=process),
            patch("elisaflow.ui.main_window.QMessageBox.warning") as warning,
        ):
            self.window.reanalyze_button.click()
            self.wait_until(started.is_set)
            self.assertFalse(self.window.start_button.isEnabled())
            self.assertFalse(self.window.reanalyze_button.isEnabled())
            self.window._start_reanalysis()
            self.window._start_run()
            picker.assert_called_once()
            self.window.cancel_button.click()
            self.wait_until(lambda: self.window._worker is None)
        warning.assert_not_called()
        self.assertIn("geannuleerd", self.window.status.text())
        self.assertTrue(self.window.start_button.isEnabled())
        self.assertTrue(self.window.reanalyze_button.isEnabled())

    def test_one_window_contains_runs_and_settings_with_existing_permissions(self) -> None:
        self.assertEqual(self.window.tabs.count(), 2)
        self.assertEqual(self.window.tabs.tabText(1), "Assay beheer")
        self.assertEqual(self.window.tabs.currentIndex(), 0)
        self.window.tabs.setCurrentWidget(self.window.assay_settings)
        self.assertIs(self.window.tabs.currentWidget(), self.window.assay_settings)
        self.assertEqual(self.window.assay_settings.assay_list.count(), 1)
        self.assertTrue(self.window.assay_settings.new_button.isEnabled())
        self.window.assay_settings.assay_list.setCurrentRow(0)
        self.assertTrue(self.window.assay_settings.view_button.isEnabled())
        self.assertTrue(self.window.assay_settings.edit_button.isEnabled())

    def test_run_actions_are_reachable_at_minimum_window_size(self):
        self.config.results_dir.mkdir()
        self.window.resize(self.window.minimumSize())
        self.window.show()
        page = self.window.tabs.widget(0)
        with (
            patch("elisaflow.ui.main_window.select_protocol", side_effect=ProtocolSelectionCancelled()) as select,
            patch("elisaflow.ui.main_window.QFileDialog.getExistingDirectory", return_value="") as picker,
            patch("elisaflow.ui.main_window.RunWorker") as worker,
        ):
            for button in (self.window.start_button, self.window.reanalyze_button):
                page.ensureWidgetVisible(button)
                self.application.processEvents()
                self.assertTrue(page.viewport().rect().contains(button.mapTo(page.viewport(), button.rect().center())))
                QTest.mouseClick(button, Qt.MouseButton.LeftButton)
            select.assert_called_once()
            picker.assert_called_once_with(
                self.window, "Kies een experimentmap voor heranalyse", str(self.config.results_dir),
            )
            worker.assert_not_called()
            self.name_prompt.assert_not_called()

    def test_documentation_button_is_available_on_both_tabs_without_starting_a_run(self):
        self.window.show()
        with (
            patch("elisaflow.ui.main_window.open_documentation") as open_manual,
            patch("elisaflow.ui.main_window.RunWorker") as worker,
        ):
            for size in (self.window.minimumSize(), self.window.size()):
                self.window.resize(size)
                for index in range(self.window.tabs.count()):
                    self.window.tabs.setCurrentIndex(index)
                    self.application.processEvents()
                    button = self.window.documentation_button
                    self.assertTrue(button.isVisible())
                    content = self.window.centralWidget()
                    self.assertTrue(content.rect().contains(button.mapTo(content, button.rect().bottomRight())))
                    self.assertGreater(
                        button.mapTo(content, button.rect().topLeft()).y(),
                        self.window.tabs.mapTo(content, self.window.tabs.rect().bottomRight()).y(),
                    )
                    QTest.mouseClick(button, Qt.MouseButton.LeftButton)
            self.assertEqual(open_manual.call_count, 4)
            worker.assert_not_called()

    def test_documentation_error_does_not_interrupt_measurement_state(self):
        with (
            patch("elisaflow.ui.main_window.open_documentation", side_effect=FileNotFoundError("handleiding ontbreekt")),
            patch("elisaflow.ui.main_window.QMessageBox.warning") as warning,
        ):
            self.window.documentation_button.click()
        warning.assert_called_once_with(self.window, "Documentatie kon niet geopend worden", "handleiding ontbreekt")
        self.assertTrue(self.window.start_button.isEnabled())
        self.assertIsNone(self.window._worker)

    def test_cancelling_selection_does_not_start_a_worker_or_show_an_error(self) -> None:
        with (
            patch("elisaflow.ui.main_window.select_protocol", side_effect=ProtocolSelectionCancelled()),
            patch("elisaflow.ui.main_window.RunWorker") as worker,
            patch("elisaflow.ui.main_window.QMessageBox.warning") as warning,
        ):
            self.window.start_button.click()
        worker.assert_not_called()
        warning.assert_not_called()
        self.assertTrue(self.window.start_button.isEnabled())
        self.assertFalse(self.config.results_dir.exists())

    def test_run_is_in_a_worker_and_the_desktop_remains_responsive(self) -> None:
        self.window.show()
        started = Event()
        thread_ids = []
        arguments = []

        def process(config, **kwargs):
            thread_ids.append(get_ident())
            arguments.append((config, kwargs))
            started.set()
            self.release.wait(5)
            return self.report

        heartbeat = []
        report_threads = []
        with (
            patch("elisaflow.ui.main_window.select_protocol", return_value=self.protocol) as select,
            patch("elisaflow.ui.run_worker.run_workflow", side_effect=process) as run,
            patch("elisaflow.ui.main_window.open_document", side_effect=lambda _path: report_threads.append(get_ident())) as open_report,
            patch("elisaflow.ui.main_window.open_documentation") as open_manual,
        ):
            self.window.start_button.click()
            self.wait_until(started.is_set)
            self.assertFalse(self.window.start_button.isEnabled())
            self.assertFalse(self.window.reanalyze_button.isEnabled())
            self.assertFalse(self.window.progress.isHidden())
            self.window._start_run()
            select.assert_called_once()
            QTimer.singleShot(0, lambda: heartbeat.append(get_ident()))
            self.window.tabs.setCurrentWidget(self.window.assay_settings)
            self.wait_until(lambda: bool(heartbeat))
            QTest.mouseClick(self.window.documentation_button, Qt.MouseButton.LeftButton)
            open_manual.assert_called_once()
            self.assertIsNotNone(self.window._worker)
            self.assertNotEqual(thread_ids[0], get_ident())
            self.assertEqual(heartbeat, [get_ident()])
            self.assertIs(arguments[0][0], self.config)
            self.assertEqual(arguments[0][1]["protocol_template"], self.protocol)
            self.assertFalse(arguments[0][1]["open_report"])
            self.assertEqual(arguments[0][1]["analyst"], "operator")
            self.name_prompt.assert_called_once_with(parent=self.window)
            self.assertIsInstance(arguments[0][1]["cancel_event"], Event)
            self.release.set()
            self.wait_until(lambda: self.window._worker is None)
            run.assert_called_once()
            open_report.assert_called_once_with(self.report)
            self.assertEqual(report_threads, [get_ident()])
        self.assertTrue(self.window.start_button.isEnabled())
        self.assertTrue(self.window.reanalyze_button.isEnabled())
        self.assertTrue(self.window.progress.isHidden())
        self.assertEqual(self.window.status_title.text(), "Rapport gereed")
        self.assertFalse(self.window.open_report_button.isHidden())
        self.assertEqual(self.window.report_path.text(), str(self.report))
        self.assertTrue(self.window.report_path.isReadOnly())

    def test_failure_is_reported_and_the_operator_can_retry(self) -> None:
        with (
            patch("elisaflow.ui.main_window.select_protocol", return_value=self.protocol),
            patch("elisaflow.ui.run_worker.run_workflow", side_effect=[ValueError("Invalid measurement"), self.report]) as run,
            patch("elisaflow.ui.main_window.QMessageBox.warning") as warning,
            patch("elisaflow.ui.main_window.open_document") as open_report,
            self.assertLogs("elisaflow.ui.run_worker", level="ERROR"),
        ):
            self.window.start_button.click()
            self.wait_until(lambda: self.window._worker is None)
            warning.assert_called_once()
            self.assertIn("Invalid measurement", self.window.status.text())
            self.assertTrue(self.window.start_button.isEnabled())
            self.window.start_button.click()
            self.wait_until(lambda: self.window._worker is None)
            self.assertEqual(run.call_count, 2)
            open_report.assert_called_once_with(self.report)

    def test_start_prompts_for_the_name_on_the_gui_thread_before_running(self):
        prompt_threads = []
        def enter_name(*, parent):
            self.assertIs(parent, self.window)
            prompt_threads.append(get_ident())
            return "Friso van Dijk"
        with (
            patch("elisaflow.ui.main_window.prompt_analyst_name", side_effect=enter_name) as prompt,
            patch("elisaflow.ui.main_window.select_protocol", return_value=self.protocol),
            patch("elisaflow.ui.run_worker.run_workflow", return_value=self.report) as run,
            patch("elisaflow.ui.main_window.open_document"),
        ):
            self.window.start_button.click()
            self.wait_until(lambda: self.window._worker is None)
            self.assertEqual(prompt_threads, [get_ident()])
            prompt.assert_called_once_with(parent=self.window)
            self.assertEqual(run.call_args.kwargs["analyst"], "Friso van Dijk")

    def test_cancelling_name_entry_keeps_the_previous_report_and_creates_no_run(self):
        self.window._report = self.report
        self.name_prompt.return_value = None
        with (
            patch("elisaflow.ui.main_window.select_protocol", return_value=self.protocol),
            patch("elisaflow.ui.main_window.RunWorker") as worker,
            patch("elisaflow.ui.main_window.QMessageBox.warning") as warning,
        ):
            self.window.start_button.click()
        worker.assert_not_called()
        warning.assert_not_called()
        self.assertEqual(self.window._report, self.report)
        self.assertTrue(self.window.start_button.isEnabled())
        self.assertFalse(self.config.results_dir.exists())

    def test_cancellation_keeps_the_worker_alive_until_it_stops_safely(self) -> None:
        started = Event()

        def process(_config, **kwargs):
            started.set()
            kwargs["cancel_event"].wait(5)
            raise RunCancelled()

        with (
            patch("elisaflow.ui.main_window.select_protocol", return_value=self.protocol),
            patch("elisaflow.ui.run_worker.run_workflow", side_effect=process),
            patch("elisaflow.ui.main_window.QMessageBox.information") as information,
            patch("elisaflow.ui.main_window.QMessageBox.warning") as warning,
            patch("elisaflow.ui.main_window.open_document") as open_report,
        ):
            self.window.start_button.click()
            self.wait_until(started.is_set)
            closing = QCloseEvent()
            self.window.closeEvent(closing)
            self.assertFalse(closing.isAccepted())
            information.assert_called_once()
            self.window.cancel_button.click()
            self.assertFalse(self.window.cancel_button.isEnabled())
            self.wait_until(lambda: self.window._worker is None)
            warning.assert_not_called()
            open_report.assert_not_called()
        self.assertIn("geannuleerd", self.window.status.text())
        self.assertTrue(self.window.start_button.isEnabled())

    def test_report_opening_errors_keep_the_successful_report_available(self) -> None:
        self.window._report = self.report
        with (
            patch("elisaflow.ui.main_window.open_document", side_effect=OSError("No PDF viewer")),
            patch("elisaflow.ui.main_window.QMessageBox.warning") as warning,
        ):
            self.window._run_succeeded(self.report)
        warning.assert_called_once()
        self.assertEqual(self.window._report, self.report)
        self.assertFalse(self.window.open_report_button.isHidden())
        self.assertFalse(self.window.report_path.isHidden())

    def test_desktop_entry_points_start_and_close_in_a_fresh_process(self) -> None:
        script = """
import json
import runpy
import sys
from contextlib import chdir
from pathlib import Path
from tempfile import TemporaryDirectory
from elisaflow.cli import main

selected_tabs = []

# Close each real window after it is shown, without importing Qt before startup.
def close_after_show(frame, event, function):
    if event == 'c_return' and frame.f_code.co_name == 'launch' and getattr(function, '__name__', '') == 'show':
        assert not any(name in sys.modules for name in ('scipy.optimize', 'matplotlib', 'reportlab', 'elisaflow.reporting.pdf_report'))
        from PySide6.QtCore import QTimer
        window = frame.f_locals['window']
        assert window.tabs.count() == 2
        selected_tabs.append(window.tabs.currentIndex())
        QTimer.singleShot(0, window.close)

with TemporaryDirectory() as temporary:
    directory = Path(temporary)
    (directory / 'configs').mkdir()
    (directory / 'configs' / 'assay.json').write_text(json.dumps({
        'assay.skax': {'lot_number': 'lot', 'calculation': {'method': 'formula', 'expression': 'sample', 'response_scale': 'blank_corrected'}}
    }))
    config = directory / 'app.toml'
    config.write_text('[paths]\\nconfigs = "configs"\\n')
    sys.setprofile(close_after_show)
    try:
        with chdir(directory):
            assert main([]) == 0
            assert main(['settings']) == 0
            sys.argv = ['elisaflow']
            try:
                runpy.run_module('elisaflow', run_name='__main__')
            except SystemExit as error:
                assert error.code == 0, error.code
    finally:
        sys.setprofile(None)
    assert selected_tabs == [0, 1, 0], selected_tabs
"""
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_first_report_loads_the_solvers_and_renderer_in_a_worker_after_startup(self) -> None:
        script = """
import sys
from contextlib import chdir
from pathlib import Path
from tempfile import TemporaryDirectory
from elisaflow.cli import main

sys.path.insert(0, sys.argv[1])
workers, errors = [], []

def start_report_after_show(frame, event, function):
    if event != 'c_return' or frame.f_code.co_name != 'launch' or getattr(function, '__name__', '') != 'show':
        return
    assert not any(name in sys.modules for name in ('scipy.optimize', 'matplotlib', 'reportlab'))
    from PySide6.QtCore import QThread, QTimer
    window = frame.f_locals['window']
    class ReportWorker(QThread):
        def run(self):
            try:
                from test_multi_wavelength import dual_data, dual_settings
                from elisaflow.runs.workflow import create_report
                create_report(dual_data(), dual_settings(), report, analyst="operator")
            except Exception as error:
                errors.append(error)
    worker = ReportWorker()
    workers.append(worker)
    worker.finished.connect(window.close)
    QTimer.singleShot(0, worker.start)

with TemporaryDirectory() as temporary:
    directory = Path(temporary)
    (directory / 'configs').mkdir()
    config = directory / 'app.toml'
    config.write_text('[paths]\\n')
    report = directory / 'first_report.pdf'
    sys.setprofile(start_report_after_show)
    try:
        with chdir(directory):
            assert main([]) == 0
    finally:
        sys.setprofile(None)
    assert len(workers) == 1
    workers[0].wait()
    if errors:
        raise errors[0]
    assert report.read_bytes().startswith(b'%PDF-')
"""
        result = subprocess.run(
            [sys.executable, "-c", script, str(Path(__file__).resolve().parent)],
            capture_output=True, text=True, timeout=30,
        )
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
