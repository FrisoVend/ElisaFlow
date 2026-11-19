"""Read the shared archive and protocol structure of a SkanIt session."""

from base64 import b64encode
from dataclasses import dataclass
from datetime import datetime
from hashlib import sha256
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from elisaflow.xml_utils import child_text, require_file


@dataclass(frozen=True)
class SkaxDocument:
    path: Path
    header: ET.Element
    session: ET.Element
    measurements: bytes | None
    disabled: bytes | None
    session_bytes: bytes

    @classmethod
    def from_file(cls, path: str | Path) -> SkaxDocument:
        path = require_file(path, ".skax", "SkanIt session")
        with ZipFile(path) as archive:
            members = set(archive.namelist())
            if len(members) != len(archive.namelist()):
                raise ValueError("Het SKAX-archief bevat dubbele onderdelen.")
            if missing := {"Header.xml", "Session.xml"} - members:
                raise ValueError(f"Vereiste bestanden ontbreken in het SKAX-archief: {', '.join(sorted(missing))}")
            return cls(
                path,
                ET.fromstring(archive.read("Header.xml")),
                ET.fromstring(archive.read("Session.xml")),
                archive.read("Measurements.data") if "Measurements.data" in members else None,
                archive.read("Disableds.data") if "Disableds.data" in members else None,
                archive.read("Session.xml"),
            )

    def validate_hashes(self) -> None:
        for name, contents, field in (
            ("Session.xml", self.session_bytes, "SessionHash"),
            ("Measurements.data", self.measurements, "MeasurementHash"),
            ("Disableds.data", self.disabled, "DisabledHash"),
        ):
            if contents is not None and child_text(self.header, field) != b64encode(sha256(contents).digest()).decode("ascii"):
                raise ValueError(f"{name} komt niet overeen met de SKAX-{field}.")

    @property
    def completed(self) -> bool:
        start, end = (child_text(self.session, name) for name in ("ExecutionStart", "ExecutionEnd"))
        if not self.measurements or not start or not end:
            return False
        try:
            started, finished = datetime.fromisoformat(start), datetime.fromisoformat(end)
            return started.year > 1 and finished >= started
        except (TypeError, ValueError):
            return False

    @property
    def protocol_steps(self) -> list[ET.Element]:
        return self.session.findall("{*}Protocol/{*}ProtocolSteps/*")

    @property
    def measurement_steps(self) -> list[ET.Element]:
        return sorted((step for step in self.protocol_steps if child_text(step, "Type") == "PhotometricFilterMeasurementStep"), key=_step_order)

    @property
    def wavelengths_nm(self) -> tuple[int, ...]:
        wavelengths = []
        filters = self.instrument_filters
        for step in self.measurement_steps:
            photometric_filters = step.findall("{*}Filters/{*}PhotometricFilter")
            if len(photometric_filters) != 1:
                raise ValueError("Elke absorbantiestap moet precies één fotometrisch filter gebruiken.")
            wavelength = child_text(photometric_filters[0], "Wavelength")
            position = child_text(photometric_filters[0], "FilterPosition")
            if not wavelength or not wavelength.isdecimal() or int(wavelength) <= 0:
                raise ValueError("Absorbantiestappen vereisen een positieve golflengte in de filtermetadata.")
            if position not in filters or filters[position] != wavelength:
                raise ValueError(f"De filterpositie {position!r} van de stap bij {wavelength} nm komt niet overeen met de filtermetadata van het instrument.")
            wavelengths.append(int(wavelength))
        if not wavelengths or len(set(wavelengths)) != len(wavelengths):
            raise ValueError("SKAX moet verschillende ruwe absorbantiegolflengten bevatten.")
        return tuple(wavelengths)

    @property
    def instrument_filters(self) -> dict[str | None, str | None]:
        return {child_text(node, "Position"): child_text(node, "Wavelength") for node in self.session.findall("{*}Instrument//{*}Filters/{*}Filter")}


def _step_order(step: ET.Element) -> int:
    return int(child_text(step, "OrderNumber") or 0)
