import builtins
import contextlib
import io
import subprocess
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

from elisaflow.cli import main
from elisaflow.configuration.application import ApplicationSettings, default_config_path


class EntryPointTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = ApplicationSettings(*(Path(name) for name in ("protocols", "configs", "results", "editor.hash")))
        load = patch("elisaflow.cli.ApplicationSettings.from_file", return_value=self.config)
        self.load = load.start()
        self.addCleanup(load.stop)
        name_prompt = patch("elisaflow.cli.prompt_analyst_name", return_value="operator")
        self.name_prompt = name_prompt.start()
        self.addCleanup(name_prompt.stop)

    def test_default_command_opens_the_shared_desktop(self) -> None:
        with patch("elisaflow.ui.application.launch", return_value=0) as launch:
            self.assertEqual(main([]), 0)
        launch.assert_called_once_with(self.config, show_settings=False)
        self.name_prompt.assert_not_called()

    def test_settings_shortcut_opens_the_same_desktop_on_its_settings_tab(self) -> None:
        with patch("elisaflow.ui.application.launch", return_value=0) as launch:
            self.assertEqual(main(["settings"]), 0)
        self.load.assert_called_once_with(default_config_path())
        launch.assert_called_once_with(self.config, show_settings=True)

    def test_run_command_passes_workflow_arguments(self) -> None:
        options = ["--protocol", "assay.skax", "--output", "report.pdf"]
        with patch("elisaflow.runs.workflow.run_workflow") as run:
            self.assertEqual(main(["run", *options]), 0)
        run.assert_called_once_with(
            self.config, analyst="operator", protocol_template=Path("assay.skax"), experiment_dir=None, output=Path("report.pdf"),
        )
        self.name_prompt.assert_called_once_with()

    def test_experiment_option_uses_saved_inputs_without_a_protocol_argument(self) -> None:
        with patch("elisaflow.runs.workflow.run_workflow") as run:
            self.assertEqual(main(["run", "--experiment", "saved run", "--output", "new/report.pdf"]), 0)
        run.assert_called_once_with(
            self.config, analyst="operator", protocol_template=None,
            experiment_dir=Path("saved run"), output=Path("new/report.pdf"),
        )

    def test_run_without_a_protocol_still_uses_the_existing_selector(self) -> None:
        with patch("elisaflow.runs.workflow.run_workflow") as run:
            self.assertEqual(main(["run"]), 0)
        self.assertIsNone(run.call_args.kwargs["protocol_template"])

    def test_cancelling_name_entry_does_not_start_the_workflow(self):
        self.name_prompt.return_value = None
        with patch("elisaflow.runs.workflow.run_workflow") as run:
            self.assertEqual(main(["run", "--protocol", "assay.skax"]), 0)
        run.assert_not_called()

    def test_invalid_argument_combinations_are_rejected_before_loading_settings(self) -> None:
        for options in (
            ["--config", "custom.toml"],
            ["settings", "--config", "custom.toml"],
            ["run", "--config", "custom.toml"],
            ["settings", "--protocol", "assay.skax"],
            ["--protocol", "assay.skax"],
            ["--measurement", "run.xml"],
            ["--output", "report.pdf"],
            ["--no-open"],
            ["run", "--no-open"],
            ["run", "--measurement", "run.xml"],
            ["run", "--protocol", "assay.skax", "--measurement", "run.xml"],
            ["run", "--protocol", "assay.skax", "--clipboard"],
            ["--clipboard"],
            ["settings", "--clipboard"],
            ["run", "--protocol", "assay.skax", "--clipboard", "--measurement", "run.xml"],
            ["--experiment", "saved-run"],
            ["settings", "--experiment", "saved-run"],
            ["run", "--experiment", "saved-run", "--protocol", "assay.skax"],
            ["run", "--experiment", "saved-run", "--measurement", "run.xml"],
            ["run", "--experiment", "saved-run", "--clipboard"],
        ):
            with self.subTest(options=options), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as caught:
                    main(options)
                self.assertEqual(caught.exception.code, 2)
        self.load.assert_not_called()

    def test_configuration_and_processing_errors_return_failure(self) -> None:
        self.load.side_effect = FileNotFoundError("Missing application settings")
        with self.assertLogs("elisaflow.cli", level="ERROR"):
            self.assertEqual(main([]), 1)
        self.load.side_effect = None
        with patch("elisaflow.runs.workflow.run_workflow", side_effect=ValueError("Invalid measurement")):
            with self.assertLogs("elisaflow.cli", level="ERROR"):
                self.assertEqual(main(["run", "--protocol", "assay.skax"]), 1)

    def test_direct_runs_work_without_importing_a_desktop_toolkit(self) -> None:
        script = """
import sys
from unittest.mock import patch

class BlockQt:
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in ('PySide6', 'shiboken6', 'tkinter'):
            raise RuntimeError('Desktop toolkit imported for a direct run')

sys.meta_path.insert(0, BlockQt())
from elisaflow.cli import main
with patch('elisaflow.cli.ApplicationSettings.from_file'), patch('elisaflow.runs.workflow.run_workflow') as run, patch('builtins.input', return_value='operator'):
    assert main(['run', '--protocol', 'assay.skax']) == 0
    run.assert_called_once()
    run.reset_mock()
    assert main(['run', '--experiment', 'saved-run']) == 0
    run.assert_called_once()
"""
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_missing_qt_explains_how_to_install_the_desktop_dependency(self) -> None:
        from elisaflow.ui.application import launch

        original_import = builtins.__import__

        def without_qt(name, *args, **kwargs):
            if name == "elisaflow.ui.main_window":
                raise ModuleNotFoundError("No module named 'PySide6'", name="PySide6")
            return original_import(name, *args, **kwargs)

        with patch("elisaflow.ui.application.AssayManagement"), patch("builtins.__import__", side_effect=without_qt):
            with self.assertRaisesRegex(RuntimeError, 'Installeer eerst de UI-afhankelijkheid'):
                launch(self.config)


if __name__ == "__main__":
    unittest.main()
