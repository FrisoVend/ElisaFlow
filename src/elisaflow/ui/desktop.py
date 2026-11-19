"""Desktop interactions kept outside readers, analysis, and PDF rendering."""

import os
import sys
import webbrowser
from pathlib import Path

from elisaflow.runs.analyst import validate_analyst_name
from elisaflow.xml_utils import require_file


class ProtocolSelectionCancelled(RuntimeError):
    """The operator closed the protocol selector without starting a run."""


def prompt_analyst_name(parent=None) -> str | None:
    """Ask for a fresh analyst name on the caller's UI or console thread."""
    if parent is not None or sys.stdin is None:
        try:
            from elisaflow.ui.run_selection import prompt_run_analyst
        except ModuleNotFoundError as error:
            if error.name and error.name.split(".")[0] in {"PySide6", "shiboken6"}:
                raise RuntimeError('Naam invoeren in de desktop vereist PySide6. Installeer dit met python -m pip install -e ".[ui]".') from error
            raise
        return prompt_run_analyst(parent)
    while True:
        try:
            value = input("Vul uw naam in: ")
        except (EOFError, KeyboardInterrupt):
            return None
        try:
            return validate_analyst_name(value)
        except ValueError as error:
            print(error)


def select_protocol(protocols_dir: Path) -> Path:
    """Let the operator choose one protocol from the configured library."""
    protocols_dir = Path(protocols_dir).expanduser().resolve()
    if not protocols_dir.is_dir():
        raise FileNotFoundError(f"De protocolmap bestaat niet: {protocols_dir}")
    protocols = tuple(
        sorted(
            (path for path in protocols_dir.iterdir() if path.is_file() and path.suffix.lower() == ".skax"),
            key=lambda path: path.name.casefold(),
        )
    )
    if not protocols:
        raise ValueError(f"Geen .skax-protocollen gevonden in {protocols_dir}.")
    try:
        from elisaflow.ui.run_selection import select_run_protocol
    except ModuleNotFoundError as error:
        if error.name and error.name.split(".")[0] in {"PySide6", "shiboken6"}:
            raise RuntimeError(
                'Protocolselectie vereist PySide6. Installeer dit met python -m pip install -e ".[ui]" of geef --protocol op.'
            ) from error
        raise
    protocol = require_file(select_run_protocol(protocols), ".skax", "SkanIt protocol")
    if protocol.parent != protocols_dir:
        raise ValueError(f"Selecteer een .skax-bestand uit {protocols_dir}.")
    return protocol


def open_protocol(path: Path) -> None:
    if os.name != "nt":
        raise RuntimeError("SkanIt openen vereist Windows; gebruik --experiment om een opgeslagen experiment opnieuw te verwerken.")
    os.startfile(path)


def open_document(path: Path) -> None:
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    if os.name == "nt":
        os.startfile(path)
    else:
        webbrowser.open(path.as_uri())


def open_documentation() -> None:
    """Open the bundled manual, independent of the current working directory."""
    directory = Path(sys.executable).resolve().parent if getattr(sys, "frozen", False) else Path(__file__).resolve().parents[3] / "build"
    document = directory / "docs" / "index.html"
    if not document.is_file():
        raise FileNotFoundError("De handleiding ontbreekt. Gebruik de volledige applicatiemap met de map 'docs'.")
    if not webbrowser.open(document.as_uri()):
        raise OSError("De handleiding kon niet in de standaardbrowser worden geopend.")
