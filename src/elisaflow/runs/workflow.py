import logging
import time
from pathlib import Path
from threading import Event
from xml.etree import ElementTree as ET
from zipfile import BadZipFile

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.assays.protocol_validation import validate_protocol_measurement
from elisaflow.assays.registry import AssayRegistry
from elisaflow.configuration.application import ApplicationSettings
from elisaflow.configuration.assay_settings import AssaySettings
from elisaflow.models import AssayData
from elisaflow.readers.assay_loader import AssayLoader
from elisaflow.readers.skax_document import SkaxDocument
from elisaflow.runs.analyst import validate_analyst_name
from elisaflow.runs.experiments import load_experiment, prepare_live_experiment, prepare_reanalysis
from elisaflow.ui.desktop import open_document, open_protocol, select_protocol

logger = logging.getLogger(__name__)


class RunCancelled(RuntimeError):
    """Processing was stopped at a safe boundary; existing files are retained."""


def _check_cancelled(cancel_event: Event | None) -> None:
    if cancel_event is not None and cancel_event.is_set():
        raise RunCancelled("Meting geannuleerd.")


def create_report(data: AssayData, settings: AssaySettings, output: Path, *, analyst: str) -> Path:
    """Analyze loaded measurements and write their PDF report."""
    from elisaflow.reporting.pdf_report import PdfReportWriter

    result = AssayAnalyzer(data, settings).analyze()
    return PdfReportWriter(result, analyst=analyst).write(output)


def wait_for_measurement(
    session: Path, timeout: float = 3600, poll_interval: float = 1.0, *, cancel_event: Event | None = None,
) -> Path:
    """Wait until SkanIt saves the copied session with completed measurement data."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        _check_cancelled(cancel_event)
        if _file_is_stable(session, poll_interval):
            _check_cancelled(cancel_event)
            return session
        if cancel_event is not None:
            cancel_event.wait(poll_interval)
        else:
            time.sleep(poll_interval)
    _check_cancelled(cancel_event)
    raise TimeoutError(f"Geen opgeslagen, voltooide SKAX-meting binnen {timeout} seconden. Sla de gemeten sessie in SkanIt op: {session}")


def _file_is_stable(path: Path, delay: float = 1.0) -> bool:
    """Require a stable, readable archive with matching hashes and a finished run."""
    try:
        before = path.stat()
        time.sleep(delay)
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            return False
        document = SkaxDocument.from_file(path)
        document.validate_hashes()
        return document.completed
    except (FileNotFoundError, PermissionError, BadZipFile, ET.ParseError, ValueError):
        return False


def run_workflow(
    config: ApplicationSettings, *, analyst: str, protocol_template: Path | None = None, experiment_dir: Path | None = None,
    output: Path | None = None, open_report: bool = True, cancel_event: Event | None = None,
) -> Path:
    """Analyze a completed SKAX from a live run or a saved experiment folder."""
    _check_cancelled(cancel_event)
    analyst = validate_analyst_name(analyst)
    if experiment_dir is not None and protocol_template is not None:
        raise ValueError("Kies een experimentmap of een protocol voor een live meting.")
    experiment = load_experiment(experiment_dir) if experiment_dir is not None else None
    if experiment is not None:
        template, settings = experiment.session, experiment.settings
    else:
        template = protocol_template.resolve() if protocol_template else select_protocol(config.protocols_dir)
        settings = AssayRegistry(config.configs_dir).get(template)
    _check_cancelled(cancel_event)
    if experiment is not None:
        prepared = prepare_reanalysis(experiment, config.results_dir, output=output, analyst=analyst)
    else:
        validate_protocol_measurement(template, settings.measurement)
        _check_cancelled(cancel_event)
        prepared = prepare_live_experiment(template, settings, config.results_dir, output=output, analyst=analyst)
    session, output = prepared.session, prepared.report
    settings = prepared.settings
    if experiment is None:
        _check_cancelled(cancel_event)
        open_protocol(session)
        logger.info("Waiting for the completed SKAX session: %s", session)
        wait_for_measurement(session, cancel_event=cancel_event)
    _check_cancelled(cancel_event)
    data = AssayLoader.load(session, measurement_settings=settings.measurement)
    _check_cancelled(cancel_event)
    report = create_report(data, settings, output, analyst=analyst)
    _check_cancelled(cancel_event)
    logger.info("Report saved to %s", report)
    if open_report:
        open_document(report)
    return report
