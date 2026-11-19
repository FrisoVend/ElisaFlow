"""Execute a live measurement or reanalysis without blocking the desktop."""

import logging
from pathlib import Path
from threading import Event

from elisaflow.configuration.application import ApplicationSettings
from elisaflow.runs.workflow import RunCancelled, run_workflow

# Load calculation dependencies before Qt installs its import hooks.
# isort: split
from PySide6.QtCore import QThread, Signal

logger = logging.getLogger(__name__)


class RunWorker(QThread):
    succeeded = Signal(object)
    failed = Signal(str)
    cancelled = Signal()

    def __init__(
        self, config: ApplicationSettings, analyst: str, parent=None,
        *, protocol_template: Path | None = None, experiment_dir: Path | None = None,
    ) -> None:
        super().__init__(parent)
        self._config = config
        self._protocol = protocol_template
        self._experiment_dir = experiment_dir
        self._cancel_event = Event()
        self._analyst = analyst

    def cancel(self) -> None:
        self._cancel_event.set()

    def run(self) -> None:
        try:
            report = run_workflow(
                self._config,
                protocol_template=self._protocol,
                experiment_dir=self._experiment_dir,
                open_report=False,
                cancel_event=self._cancel_event,
                analyst=self._analyst,
            )
        except RunCancelled:
            self.cancelled.emit()
        except Exception as error:
            logger.exception("Analysis failed")
            self.failed.emit(str(error))
        else:
            self.succeeded.emit(report)
