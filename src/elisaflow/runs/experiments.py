"""Own experiment discovery, preserved inputs and output-path conventions."""

import json
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from shutil import copyfile

from elisaflow.configuration.assay_settings import AssaySettings
from elisaflow.configuration.serialization import serialize_assay_configs
from elisaflow.runs.analyst import validate_analyst_name


@dataclass(frozen=True)
class SavedExperiment:
    directory: Path
    session: Path
    settings: AssaySettings


@dataclass(frozen=True)
class PreparedExperiment:
    """Preserved session, chosen report path and the settings for this run."""

    directory: Path
    session: Path
    report: Path
    settings: AssaySettings


def load_experiment(folder: Path) -> SavedExperiment:
    """Find one SKAX and parse its snapshot; measurement validity is checked later."""
    directory = Path(folder).expanduser().resolve()
    if not directory.is_dir():
        raise FileNotFoundError(f"De experimentmap bestaat niet: {directory}")

    def find_one(pattern: str) -> Path:
        matches = sorted(path for path in directory.glob(pattern) if path.is_file())
        if len(matches) != 1:
            raise ValueError(f"Experimentmap moet precies één {pattern}-bestand bevatten; {len(matches)} gevonden: {directory}")
        return matches[0]

    session = find_one("*.skax")
    settings_path = find_one("*_settings.json")
    try:
        snapshot = json.loads(settings_path.read_text(encoding="utf-8-sig"))
    except (ValueError, UnicodeError) as error:
        raise ValueError(f"Ongeldige assay-instellingen in {settings_path}: {error}") from error
    if not isinstance(snapshot, dict) or len(snapshot) != 1:
        raise ValueError(f"Instellingensnapshot moet precies één assay bevatten: {settings_path}")
    protocol_name, config = next(iter(snapshot.items()))
    try:
        settings = AssaySettings.from_dict(protocol_name, config)
    except ValueError as error:
        raise ValueError(f"Ongeldige assay-instellingen in {settings_path}: {error}") from error
    return SavedExperiment(directory, session, settings)


def _create_result_directory(template: Path, results_dir: Path, *, analyst: str) -> Path:
    """Group unique run folders by ISO week-year and week of processing."""
    started = datetime.now()
    year, week, _ = started.isocalendar()
    filename_name = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", validate_analyst_name(analyst)).strip(" .") or "_"
    name = f"{template.stem}_{started:%Y-%m-%d_%H-%M-%S}_{filename_name}"
    week_dir = results_dir / str(year) / f"week_{week:02d}"
    sequence = 1
    while True:
        run_dir = week_dir / (name if sequence == 1 else f"{name}_{sequence}")
        try:
            run_dir.mkdir(parents=True)
        except FileExistsError:
            sequence += 1
        else:
            return run_dir


def _result_file(output: Path, suffix: str) -> Path:
    """Give every saved artifact the same prefix as its report."""
    return output.with_name(f"{output.stem.removesuffix('_report')}_{suffix}")


def _write_snapshot(path: Path, contents: str | bytes, conflict: str) -> Path:
    """Create a snapshot once, accepting an existing copy only when it matches."""
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open("xb") as snapshot:
            snapshot.write(contents.encode("utf-8") if isinstance(contents, str) else contents)
    except FileExistsError:
        existing = path.read_text(encoding="utf-8") if isinstance(contents, str) else path.read_bytes()
        if existing != contents:
            raise FileExistsError(f"{conflict}: {path}") from None
    return path


def _export_assay_settings(settings: AssaySettings, output: Path) -> Path:
    """Preserve the run's settings beside its PDF, without replacing an older snapshot."""
    contents = serialize_assay_configs({settings.protocol_name: settings.to_dict()})
    return _write_snapshot(_result_file(output, "settings.json"), contents, "Er bestaat al een kopie van de assay-instellingen met andere instellingen")


def prepare_live_experiment(
    template: Path, settings: AssaySettings, results_dir: Path, *, analyst: str, output: Path | None = None,
) -> PreparedExperiment:
    """Copy a live template and settings, including snapshots for a custom PDF path."""
    analyst = validate_analyst_name(analyst)
    directory = _create_result_directory(template, results_dir, analyst=analyst)
    archived_report = directory / f"{directory.name}_report.pdf"
    output = output or archived_report
    session = directory / f"{directory.name}_session.skax"
    copyfile(template, session)
    if output.parent.resolve() != directory.resolve():
        _export_assay_settings(settings, archived_report)
    _export_assay_settings(settings, output)
    return PreparedExperiment(directory, session, output, settings)


def prepare_reanalysis(
    experiment: SavedExperiment, results_dir: Path, *, analyst: str, output: Path | None = None,
) -> PreparedExperiment:
    """Preserve reanalysis inputs outside their source; retained copies survive later failures."""
    analyst = validate_analyst_name(analyst)
    destination = (output.parent if output is not None else results_dir).expanduser().resolve()
    if destination.is_relative_to(experiment.directory) or (output is not None and output.resolve().is_relative_to(experiment.directory)):
        raise ValueError("Heranalyse moet buiten de oorspronkelijke experimentmap worden opgeslagen.")
    if output is None:
        directory = _create_result_directory(Path(experiment.settings.protocol_name), destination, analyst=analyst)
        output = directory / f"{directory.name}_report.pdf"
    else:
        output = output.expanduser().resolve()
    session = _result_file(output, "session.skax")
    for path in (output, session, _result_file(output, "settings.json")):
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"Kies een nieuw uitvoerpad voor heranalyse; bestand bestaat al: {path}")
    _write_snapshot(session, experiment.session.read_bytes(), "Er bestaat al een kopie van de SKAX-sessie met andere inhoud")
    _export_assay_settings(experiment.settings, output)
    return PreparedExperiment(session.parent, session, output, experiment.settings)
