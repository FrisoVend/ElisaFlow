"""Keep diagnostics accessible when Windows launches ElisaFlow without a console."""

import ctypes
import logging
import os
import sys
from logging.handlers import RotatingFileHandler
from pathlib import Path

from elisaflow.branding import APP_NAME


def show_error(message: str, log_path: Path | None = None) -> None:
    if log_path is not None:
        message += f"\n\nLogbestand: {log_path}"
    ctypes.windll.user32.MessageBoxW(None, message, APP_NAME, 0x10)


def configure_logging() -> Path | None:
    if not getattr(sys, "frozen", False) and sys.stderr is not None:
        logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
        return None

    directory = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local")) / APP_NAME
    directory.mkdir(parents=True, exist_ok=True)
    if getattr(sys, "frozen", False):
        # The folder build has stable font paths, so reuse the font cache rather
        # than rebuilding it in PyInstaller's temporary Matplotlib directory.
        os.environ["MPLCONFIGDIR"] = str(directory / "matplotlib")
    log_path = directory / f"{APP_NAME}.log"
    handler = RotatingFileHandler(log_path, maxBytes=1_048_576, backupCount=3, encoding="utf-8")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        handlers=[handler],
        force=True,
    )

    def report_exception(exception_type, exception, traceback) -> None:
        logging.getLogger(__name__).error("Unexpected error", exc_info=(exception_type, exception, traceback))
        show_error(str(exception), log_path)

    sys.excepthook = report_exception
    return log_path
