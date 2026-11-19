"""Independent SkanIt 8 archives for reader and workflow integration tests."""

import csv
from base64 import b64encode
from hashlib import sha256
from io import StringIO
from math import isfinite
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from fixtures import add_measurement_steps


def basic_records(wavelengths=(450,)) -> list[dict]:
    return [
        dict(well=f"A{index + 1}", sample_id=name, sample_name=name, sample_type=kind, source_extinction=value, wavelength_nm=wavelength)
        for wavelength in wavelengths
        for index, (name, kind, value) in enumerate((("Blank", "Blank", 0.1), ("Cutoff", "Control", 0.6), ("Patient", "Unknown", 1.1)))
    ]


def write_skax(path: Path, records=None, *, completed=True, wavelengths=None, layout=None, disabled_indices=None, name="test.skax") -> Path:
    records = basic_records() if records is None else list(records)
    layout = records if layout is None else list(layout)
    wavelengths = tuple(dict.fromkeys(row["wavelength_nm"] for row in layout)) if wavelengths is None else tuple(wavelengths)
    session = ET.Element("Session")
    ET.SubElement(session, "Name").text = name
    ET.SubElement(session, "ExecutionStart").text = "2026-10-06T10:00:00+02:00" if completed else ""
    ET.SubElement(session, "ExecutionEnd").text = "2026-10-06T10:00:02+02:00" if completed else ""
    add_measurement_steps(session, wavelengths)
    header = ET.Element("Header")
    step_indexes = ET.SubElement(header, "Steps")
    for index, step in enumerate(session.findall("Protocol/ProtocolSteps/*"), 1):
        ET.SubElement(step, "Id").text = f"step-{index}"
        pair = ET.SubElement(step_indexes, "Pair")
        ET.SubElement(pair, "Key").text = str(index)
        ET.SubElement(pair, "Value").text = f"step-{index}"
    plate_indexes = ET.SubElement(header, "Plates")
    plates = ET.SubElement(ET.SubElement(session, "PlateLayout"), "Plates")
    plate_nodes, seen = {}, set()
    for record in layout:
        plate_index = int(record.get("plate_index", 1))
        if plate_index not in plate_nodes:
            pair = ET.SubElement(plate_indexes, "Pair")
            ET.SubElement(pair, "Key").text = str(plate_index)
            ET.SubElement(pair, "Value").text = f"plate-{plate_index}"
            plate = ET.SubElement(plates, "Plate")
            ET.SubElement(plate, "Id").text = f"plate-{plate_index}"
            ET.SubElement(plate, "Name").text = record.get("plate_name") or f"plate-{plate_index}"
            plate_nodes[plate_index] = ET.SubElement(plate, "Wells")
        coordinate = record["well"]
        if (plate_index, coordinate) in seen:
            continue
        seen.add((plate_index, coordinate))
        well = ET.SubElement(plate_nodes[plate_index], "Well")
        ET.SubElement(well, "WellLocationX").text = str(int(coordinate[1:]) - 1)
        ET.SubElement(well, "WellLocationY").text = str(ord(coordinate[0]) - ord("A"))
        sample = ET.SubElement(ET.SubElement(ET.SubElement(well, "Replicate"), "Contents"), "Sample")
        for field, value in (
            ("Id", record["sample_id"]), ("Name", record["sample_name"]), ("SampleType", record["sample_type"]),
            ("Dilution", record.get("sample_dilution", 1)),
        ):
            ET.SubElement(sample, field).text = str(value)
        concentration = record.get("sample_concentration")
        if concentration is not None and isfinite(float(concentration)):
            ET.SubElement(sample, "Concentration").text = str(concentration)
    members = {"Session.xml": ET.tostring(session, encoding="utf-8")}
    if completed:
        output = StringIO(newline="")
        writer = csv.writer(output, lineterminator="\r\n")
        disabled = []
        for index, record in enumerate(records):
            coordinate, wavelength = record["well"], record["wavelength_nm"]
            writer.writerow([
                record.get("plate_index", 1), int(coordinate[1:]) - 1, ord(coordinate[0]) - ord("A"),
                record.get("well_point_x", 0), record.get("well_point_y", 0), record.get("step_index", wavelengths.index(wavelength) + 1),
                record.get("iteration", 0), int(record.get("saturated", False)), record["source_extinction"], 1000 + index,
                wavelength, 0, 0, 0, 0, 0,
            ])
            if record.get("disabled", False):
                disabled.append(index)
        members["Measurements.data"] = output.getvalue().encode("utf-8")
        if disabled_indices is not None:
            disabled = disabled_indices
        if disabled:
            members["Disableds.data"] = "".join(f"{index}\r\n" for index in disabled).encode("utf-8")
    members["Header.xml"] = ET.tostring(header)
    write_members(path, members)
    return path


def write_members(path: Path, members: dict[str, bytes], *, update_hashes=True) -> None:
    if update_hashes:
        header = ET.fromstring(members["Header.xml"])
        for name, field in (("Session.xml", "SessionHash"), ("Measurements.data", "MeasurementHash"), ("Disableds.data", "DisabledHash")):
            if name in members:
                node = header.find(field)
                if node is None:
                    node = ET.SubElement(header, field)
                node.text = b64encode(sha256(members[name]).digest()).decode("ascii")
        members = members | {"Header.xml": ET.tostring(header)}
    with ZipFile(path, "w") as archive:
        for name, contents in members.items():
            archive.writestr(name, contents)


def replace_member(path: Path, name: str, contents: bytes, *, update_hashes=True) -> None:
    with ZipFile(path) as archive:
        members = {entry.filename: archive.read(entry) for entry in archive.infolist()}
    write_members(path, members | {name: contents}, update_hashes=update_hashes)
