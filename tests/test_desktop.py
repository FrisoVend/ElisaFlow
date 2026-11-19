from __future__ import annotations

import builtins
import os
import subprocess
import sys
import unittest
from contextlib import chdir
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from elisaflow.ui.desktop import open_documentation, prompt_analyst_name, select_protocol

try:
    from PySide6.QtCore import QTimer
    from PySide6.QtWidgets import QApplication, QInputDialog

    from elisaflow.ui.run_selection import RunSelectionDialog, prompt_run_analyst, select_run_protocol
except ModuleNotFoundError as error:
    if error.name != "PySide6":
        raise
    QApplication = None


class ProtocolDirectoryFixture:
    def setUp(self) -> None:
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.first = self.directory / "Alpha.skax"
        self.second = self.directory / "beta.SKAX"
        self.first.write_bytes(b"protocol one")
        self.second.write_bytes(b"protocol two")


class DocumentationTests(unittest.TestCase):
    def test_documentation_opens_from_source_and_bundle_independent_of_working_directory(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source project" / "src" / "elisaflow" / "ui" / "desktop.py"
            executable = root / "USB installatie" / "ElisaFlow.exe"
            for frozen, directory in ((False, root / "source project" / "build"), (True, executable.parent)):
                with self.subTest(frozen=frozen):
                    document = directory / "docs" / "index.html"
                    document.parent.mkdir(parents=True)
                    document.write_text("<html>Handleiding</html>", encoding="utf-8")
                    with (
                        chdir(root),
                        patch("elisaflow.ui.desktop.__file__", str(source)),
                        patch("elisaflow.ui.desktop.sys.executable", str(executable)),
                        patch("elisaflow.ui.desktop.sys.frozen", frozen, create=True),
                        patch("elisaflow.ui.desktop.webbrowser.open", return_value=True) as browser,
                    ):
                        open_documentation()
                    browser.assert_called_once_with(document.resolve().as_uri())

    def test_missing_bundled_documentation_does_not_open_a_browser(self):
        with TemporaryDirectory() as temporary:
            executable = Path(temporary) / "ElisaFlow.exe"
            with (
                patch("elisaflow.ui.desktop.sys.executable", str(executable)),
                patch("elisaflow.ui.desktop.sys.frozen", True, create=True),
                patch("elisaflow.ui.desktop.webbrowser.open") as browser,
                self.assertRaisesRegex(FileNotFoundError, "handleiding ontbreekt"),
            ):
                open_documentation()
            browser.assert_not_called()

    def test_browser_failure_is_reported(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            document = root / "docs" / "index.html"
            document.parent.mkdir()
            document.write_text("<html>Handleiding</html>", encoding="utf-8")
            with (
                patch("elisaflow.ui.desktop.sys.executable", str(root / "ElisaFlow.exe")),
                patch("elisaflow.ui.desktop.sys.frozen", True, create=True),
                patch("elisaflow.ui.desktop.webbrowser.open", return_value=False),
                self.assertRaisesRegex(OSError, "standaardbrowser"),
            ):
                open_documentation()


class DesktopSelectionTests(ProtocolDirectoryFixture, unittest.TestCase):
    def test_missing_or_empty_protocol_directory_is_reported_before_loading_qt(self) -> None:
        with self.assertRaisesRegex(FileNotFoundError, 'De protocolmap bestaat niet'):
            select_protocol(self.directory / "missing")
        empty = self.directory / "empty"
        empty.mkdir()
        (empty / "ignored.txt").write_text("not a protocol", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, 'Geen .skax-protocollen gevonden'):
            select_protocol(empty)

    def test_missing_qt_explains_optional_dependency_and_command_line_alternative(self) -> None:
        original_import = builtins.__import__

        def without_qt(name: str, *args, **kwargs):
            if name == "elisaflow.ui.run_selection":
                raise ModuleNotFoundError("No module named 'PySide6'", name="PySide6")
            return original_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=without_qt):
            with self.assertRaises(RuntimeError) as caught:
                select_protocol(self.directory)
        self.assertIn('python -m pip install -e ".[ui]"', str(caught.exception))
        self.assertIn("--protocol", str(caught.exception))

    def test_command_line_help_does_not_import_qt(self) -> None:
        script = """
import runpy
import sys

class BlockQt:
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in ('PySide6', 'shiboken6', 'tkinter'):
            raise RuntimeError('Desktop toolkit imported during command-line help')

sys.meta_path.insert(0, BlockQt())
sys.argv = ['elisaflow', '--help']
runpy.run_module('elisaflow', run_name='__main__')
"""
        result = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True,
            text=True,
            timeout=30,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--protocol", result.stdout)


@unittest.skipIf(QApplication is None, "Install the optional PySide6 UI dependency")
class RunSelectionTests(ProtocolDirectoryFixture, unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        environment = patch.dict(os.environ, {"QT_QPA_PLATFORM": "offscreen"})
        environment.start()
        cls.addClassCleanup(environment.stop)
        cls.application = QApplication.instance() or QApplication([])

    def test_dropdown_requires_selection_and_returns_the_selected_protocol(self) -> None:
        dialog = RunSelectionDialog((self.first, self.second))
        self.addCleanup(dialog.deleteLater)
        self.assertEqual(dialog.protocol.currentIndex(), -1)
        self.assertFalse(dialog.protocol.isEditable())
        self.assertFalse(dialog.start_button.isEnabled())
        with self.assertRaisesRegex(RuntimeError, "Geen SkanIt protocol geselecteerd"):
            dialog.selected_protocol()
        dialog.protocol.setCurrentIndex(1)
        self.assertTrue(dialog.start_button.isEnabled())
        self.assertEqual(dialog.selected_protocol(), self.second)
        dialog.start_button.click()
        self.assertEqual(dialog.result(), dialog.DialogCode.Accepted)

    def test_name_prompt_accepts_a_full_name_and_reuses_the_application(self):
        def enter_name():
            dialog = self.application.activeModalWidget()
            self.assertIsInstance(dialog, QInputDialog)
            dialog.setTextValue("  Friso van Dijk  ")
            dialog.accept()
        QTimer.singleShot(0, enter_name)
        self.assertEqual(prompt_run_analyst(), "Friso van Dijk")
        self.assertIs(QApplication.instance(), self.application)

    def test_name_prompt_retries_whitespace_and_can_be_cancelled(self):
        with (
            patch("elisaflow.ui.run_selection.QInputDialog.getText", side_effect=[("  ", True), ("", False)]) as enter,
            patch("elisaflow.ui.run_selection.QMessageBox.warning") as warning,
        ):
            self.assertIsNone(prompt_run_analyst())
        self.assertEqual(enter.call_count, 2)
        warning.assert_called_once()

    def test_name_entry_uses_a_dialog_for_the_desktop_and_for_consoleless_runs(self):
        parent = RunSelectionDialog((self.first,))
        self.addCleanup(parent.deleteLater)
        for owner in (parent, None):
            with (
                self.subTest(parent=owner), patch("elisaflow.ui.desktop.sys.stdin", None),
                patch("elisaflow.ui.run_selection.QInputDialog.getText", return_value=("Friso van Dijk", True)) as enter,
                patch("builtins.input") as console,
            ):
                self.assertEqual(prompt_analyst_name(parent=owner), "Friso van Dijk")
            self.assertIs(enter.call_args.args[0], owner)
            console.assert_not_called()

    def test_selector_uses_a_modal_dialog_and_reuses_existing_application(self) -> None:
        execute = RunSelectionDialog.exec

        def start_run(dialog: RunSelectionDialog) -> int:
            def choose() -> None:
                dialog.protocol.setCurrentIndex(1)
                dialog.start_button.click()

            QTimer.singleShot(0, choose)
            return execute(dialog)

        with patch.object(RunSelectionDialog, "exec", start_run):
            selected = select_run_protocol((self.first, self.second))
        self.assertEqual(selected, self.second)
        self.assertIs(QApplication.instance(), self.application)

    def test_cancel_closes_modal_selector_without_selecting_a_protocol(self) -> None:
        execute = RunSelectionDialog.exec

        def cancel(dialog: RunSelectionDialog) -> int:
            QTimer.singleShot(0, dialog.reject)
            return execute(dialog)

        with patch.object(RunSelectionDialog, "exec", cancel):
            with self.assertRaisesRegex(RuntimeError, "Geen SkanIt protocol geselecteerd"):
                select_run_protocol((self.first, self.second))

    def test_directory_listing_filters_and_sorts_protocol_files(self) -> None:
        (self.directory / "ignored.txt").write_text("not a protocol", encoding="utf-8")
        (self.directory / "folder.skax").mkdir()
        with patch("elisaflow.ui.run_selection.select_run_protocol", return_value=self.second) as select:
            chosen = select_protocol(self.directory)
        select.assert_called_once_with((self.first, self.second))
        self.assertEqual(chosen, self.second.resolve())

    def test_selected_protocol_must_still_exist_in_the_configured_directory(self) -> None:
        outside = self.directory / "buiten"
        outside.mkdir()
        protocol = outside / "other.skax"
        protocol.write_bytes(b"outside protocol")
        with patch("elisaflow.ui.run_selection.select_run_protocol", return_value=protocol):
            with self.assertRaisesRegex(ValueError, 'Selecteer een .skax-bestand uit'):
                select_protocol(self.directory)

        def remove_before_return(_protocols) -> Path:
            self.first.unlink()
            return self.first

        with patch("elisaflow.ui.run_selection.select_run_protocol", side_effect=remove_before_return):
            with self.assertRaises(FileNotFoundError):
                select_protocol(self.directory)

    def test_selector_can_create_and_reuse_application_in_a_standalone_process(self) -> None:
        script = """
from pathlib import Path
from unittest.mock import patch
from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QApplication
from elisaflow.ui.run_selection import RunSelectionDialog, select_run_protocol

assert QApplication.instance() is None
execute = RunSelectionDialog.exec

def start(dialog):
    def choose():
        dialog.protocol.setCurrentIndex(0)
        dialog.start_button.click()
    QTimer.singleShot(0, choose)
    return execute(dialog)

with patch.object(RunSelectionDialog, 'exec', start):
    assert select_run_protocol((Path('test.skax'),)) == Path('test.skax')
    application = QApplication.instance()
    assert application is not None
    assert select_run_protocol((Path('second.skax'),)) == Path('second.skax')
    assert QApplication.instance() is application
"""
        result = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True,
            text=True,
            timeout=30,
        )
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
