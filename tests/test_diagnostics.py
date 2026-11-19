import io
import logging
import os
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from elisaflow.cli import main
from elisaflow.diagnostics import configure_logging


class WindowedDiagnosticsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        root = logging.getLogger()
        previous_handlers, previous_level = root.handlers[:], root.level
        root.handlers.clear()

        def restore_logging() -> None:
            for handler in root.handlers[:]:
                root.removeHandler(handler)
                handler.close()
            root.handlers[:] = previous_handlers
            root.setLevel(previous_level)

        self.addCleanup(restore_logging)
        for patcher in (
            patch.dict(os.environ, {"LOCALAPPDATA": str(self.directory)}),
            patch.object(sys, "stderr", None),
            patch.object(sys, "stdout", None),
            patch.object(sys, "excepthook"),
        ):
            patcher.start()
            self.addCleanup(patcher.stop)

    def test_startup_failure_is_logged_and_shown_without_a_console(self) -> None:
        missing_config = self.directory / "app.toml"
        with (
            patch("elisaflow.cli.default_config_path", return_value=missing_config),
            patch("elisaflow.cli.show_error") as show_error,
        ):
            self.assertEqual(main([]), 1)
        log_path = self.directory / "ElisaFlow" / "ElisaFlow.log"
        contents = log_path.read_text(encoding="utf-8")
        self.assertIn("app.toml", contents)
        self.assertIn("FileNotFoundError", contents)
        self.assertEqual(show_error.call_args.args[1], log_path)
        self.assertIn("app.toml", show_error.call_args.args[0])

    def test_unexpected_error_keeps_a_traceback_and_shows_the_log_path(self) -> None:
        log_path = configure_logging()
        try:
            raise LookupError("Unexpected startup failure")
        except LookupError:
            with patch("elisaflow.diagnostics.show_error") as show_error:
                sys.excepthook(*sys.exc_info())
        show_error.assert_called_once_with("Unexpected startup failure", log_path)
        contents = log_path.read_text(encoding="utf-8")
        self.assertIn("Traceback", contents)
        self.assertIn("LookupError: Unexpected startup failure", contents)

    def test_packaged_application_logs_to_file_with_inherited_output_handles(self) -> None:
        inherited_output = io.StringIO()
        with patch.object(sys, "frozen", True, create=True), patch.object(sys, "stderr", inherited_output):
            log_path = configure_logging()
            logging.getLogger("elisaflow").info("Processing a measurement")
        self.assertIn("Processing a measurement", log_path.read_text(encoding="utf-8"))
        self.assertEqual(inherited_output.getvalue(), "")


if __name__ == "__main__":
    unittest.main()
