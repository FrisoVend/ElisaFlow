"""Read protocol layouts, recorded measurements and metadata from SkanIt sessions."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, TypedDict
from xml.etree import ElementTree as ET

import pandas as pd

from elisaflow.readers.measurement_data import read_measurements
from elisaflow.readers.skax_document import SkaxDocument
from elisaflow.xml_utils import child_element, child_text, xml_attribute


class ExecutionMetadata(TypedDict):
    start: str | None
    end: str | None


class InstrumentMetadata(TypedDict):
    name: str | None


class SoftwareMetadata(TypedDict):
    name: str | None
    version: str | None


class ProtocolSessionMetadata(TypedDict):
    id: str | None
    name: str | None
    description: str | None
    execution: ExecutionMetadata
    instrument: InstrumentMetadata
    software: SoftwareMetadata


class ProtocolMetadata(TypedDict):
    source_path: str
    session: ProtocolSessionMetadata


@dataclass(frozen=True)
class SkaxReader:
    """A fully read SkanIt SKAX session or protocol, constructed with from_file."""

    path: Path
    metadata: ProtocolMetadata
    wells: pd.DataFrame
    measurements: pd.DataFrame
    completed: bool
    wavelengths_nm: tuple[int, ...]

    @classmethod
    def from_file(cls, path: str | Path) -> SkaxReader:
        document = SkaxDocument.from_file(path)
        wavelengths_nm = document.wavelengths_nm
        header, session = document.header, document.session
        references = cls._build_reference_index(session)
        plate_indexes = cls._plate_indexes(header)
        wells = cls._read_wells(session, references, plate_indexes)
        wells = wells.sort_values(["plate_index", "well_x", "well_y"]).reset_index(drop=True)
        metadata = cls._build_metadata(document.path, session, references)

        measurements = pd.DataFrame()
        if document.measurements is not None:
            document.validate_hashes()
            measurements = cls._read_measurements(document, wells)
        return cls(document.path, metadata, wells, measurements, document.completed, wavelengths_nm)

    @classmethod
    def _read_measurements(cls, document: SkaxDocument, wells: pd.DataFrame) -> pd.DataFrame:
        readings = read_measurements(document.measurements, document.disabled)
        steps = dict(zip((child_text(step, "Id") for step in document.measurement_steps), document.wavelengths_nm, strict=True))
        wavelengths = {
            int(cls._required_text(pair, "Key")): steps[child_text(pair, "Value")]
            for pair in document.header.findall(".//{*}Steps/*") if child_text(pair, "Value") in steps
        }
        expected = readings["step_index"].map(wavelengths)
        if expected.isna().any() or readings["wavelength_nm"].ne(expected).any():
            raise ValueError("De ruwe meetgolflengten moeten overeenkomen met de SKAX-meetstappen en filtermetadata.")
        if readings.duplicated(["plate_index", "well_x", "well_y", "step_index"]).any():
            raise ValueError("Dubbele wellmetingen bij dezelfde golflengte.")
        columns = ["plate_index", "well_x", "well_y", "wavelength_nm", "source_extinction", "saturated", "disabled", "step_index"]
        combined = readings[columns].merge(wells, on=["plate_index", "well_x", "well_y"], how="left", validate="many_to_one", indicator=True)
        if combined["_merge"].ne("both").any():
            raise ValueError("Een SKAX-meetrecord verwijst naar een well buiten de protocolindeling.")
        return combined.drop(columns="_merge").sort_values(["step_index", "plate_index", "well_x", "well_y"], kind="stable").reset_index(drop=True)

    @classmethod
    def _read_wells(
        cls,
        session: ET.Element,
        references: dict[str, ET.Element],
        plate_indexes: dict[str, int],
    ) -> pd.DataFrame:
        rows: list[dict[str, Any]] = []
        plates = session.find("{*}PlateLayout/{*}Plates")
        if plates is None:
            raise ValueError("Session.xml bevat geen PlateLayout/Plates")

        for plate_node in plates:
            plate = cls._resolve(plate_node, references)
            plate_id = child_text(plate, "Id")
            plate_index = plate_indexes.get(plate_id) if plate_id is not None else None
            plate_name = cls._value(plate, "Name", references)
            wells = child_element(plate, "Wells")
            if wells is None:
                continue

            for well_node in wells:
                rows.append(cls._read_well(well_node, references, plate_index, plate_name))

        frame = pd.DataFrame(rows)
        frame["sample_concentration"] = pd.to_numeric(frame["sample_concentration"], errors="coerce")
        # Do not turn a malformed dilution into an apparently unspecified factor.
        try:
            frame["sample_dilution"] = pd.to_numeric(frame["sample_dilution"], errors="raise")
        except (TypeError, ValueError) as error:
            raise ValueError("Verdunningsfactoren moeten numeriek zijn.") from error
        return frame

    @classmethod
    def _read_well(
        cls,
        well_node: ET.Element,
        references: dict[str, ET.Element],
        plate_index: int | None,
        plate_name: str | None,
    ) -> dict[str, Any]:
        well = cls._resolve(well_node, references)
        x = int(cls._required_text(well, "WellLocationX"))
        y = int(cls._required_text(well, "WellLocationY"))

        replicate = cls._resolve(child_element(well, "Replicate"), references)
        contents = child_element(replicate, "Contents")
        sample_node = next(iter(contents), None) if contents is not None else None
        sample = cls._resolve(sample_node, references)

        return {
            "plate_index": plate_index,
            "plate_name": plate_name,
            "well": f"{cls._row_label(y)}{x + 1}",
            "well_x": x,
            "well_y": y,
            "sample_name": cls._value(sample, "Name", references),
            "sample_id": child_text(sample, "Id"),
            "sample_type": child_text(sample, "SampleType"),
            "sample_concentration": child_text(sample, "Concentration"),
            "sample_dilution": child_text(sample, "Dilution") or "1",
        }

    @classmethod
    def _build_metadata(
        cls,
        path: Path,
        session: ET.Element,
        references: dict[str, ET.Element],
    ) -> ProtocolMetadata:
        instrument_node = child_element(session, "Instrument")
        instrument = cls._resolve(instrument_node, references) if instrument_node is not None else None

        return {
            "source_path": str(path),
            "session": {
                "id": child_text(session, "Id"),
                "name": child_text(session, "Name"),
                "description": child_text(session, "Description"),
                "execution": {
                    "start": child_text(session, "ExecutionStart"),
                    "end": child_text(session, "ExecutionEnd"),
                },
                "instrument": {
                    "name": cls._value(instrument, "Name", references),
                },
                "software": cls._software(session),
            },
        }

    @classmethod
    def _plate_indexes(cls, header: ET.Element) -> dict[str, int]:
        return {cls._required_text(pair, "Value"): int(cls._required_text(pair, "Key")) for pair in header.findall(".//{*}Plates/*")}

    @staticmethod
    def _software(session: ET.Element) -> SoftwareMetadata:
        version_xml = child_text(session, "Version")
        if not version_xml:
            return {"name": None, "version": None}

        version_info = ET.fromstring(version_xml)
        version = child_element(version_info, "SkanItVersion")
        parts = [child_text(version, part) for part in ("_Major", "_Minor", "_Build", "_Revision")]
        version_number = ".".join(part for part in parts if part is not None) or None
        return {
            "name": child_text(version_info, "ApplicationName"),
            "version": version_number,
        }

    @staticmethod
    def _required_text(element: ET.Element, name: str) -> str:
        value = child_text(element, name)
        if value is None:
            raise ValueError(f"{name!r} verwacht in Session.xml")
        return value

    @staticmethod
    def _resolve(
        element: ET.Element | None,
        references: dict[str, ET.Element],
    ) -> ET.Element:
        if element is None:
            raise ValueError("Session.xml bevat een verwijzing naar een ontbrekend element")
        return references.get(xml_attribute(element, "Ref"), element)

    @classmethod
    def _value(
        cls,
        element: ET.Element | None,
        name: str,
        references: dict[str, ET.Element],
    ) -> str | None:
        if element is None:
            return None
        child = child_element(cls._resolve(element, references), name)
        return cls._resolve(child, references).text if child is not None else None

    @staticmethod
    def _build_reference_index(session: ET.Element) -> dict[str, ET.Element]:
        return {reference: element for element in session.iter() if (reference := xml_attribute(element, "Id")) is not None}

    @staticmethod
    def _row_label(index: int) -> str:
        label = ""
        number = index + 1
        while number:
            number, remainder = divmod(number - 1, 26)
            label = chr(65 + remainder) + label
        return label
