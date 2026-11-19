"""Read SkanIt 8 measurement records and their disabled-record indices."""

import csv
from io import StringIO

import pandas as pd

MEASUREMENT_COLUMNS = (
    "plate_index", "well_x", "well_y", "well_point_x", "well_point_y", "step_index", "iteration", "saturated",
    "source_extinction", "relative_measurement_time", "wavelength_nm", "emission_wavelength", "no_filter_used", "decay_time", "polarisation", "g_factor",
)


def read_measurements(contents: bytes, disabled_contents: bytes | None) -> pd.DataFrame:
    records = []
    try:
        for row in csv.reader(StringIO(contents.decode("utf-8-sig"), newline="")):
            if not row:
                continue
            if len(row) != len(MEASUREMENT_COLUMNS):
                raise ValueError("Elk meetrecord in Measurements.data moet precies 16 velden bevatten.")
            records.append([float(value) if index in (8, 15) else int(value) for index, value in enumerate(row)])
        disabled = {int(line) for line in disabled_contents.decode("utf-8-sig").splitlines() if line.strip()} if disabled_contents is not None else set()
    except (UnicodeError, ValueError, csv.Error) as error:
        raise ValueError(f"Ongeldige SKAX-meetgegevens: {error}") from error
    if not records:
        raise ValueError("Measurements.data bevat geen metingen.")
    if any(index < 0 or index >= len(records) for index in disabled):
        raise ValueError("Disableds.data verwijst naar een niet-bestaand meetrecord.")
    frame = pd.DataFrame(records, columns=MEASUREMENT_COLUMNS)
    if not frame[["saturated", "no_filter_used"]].isin((0, 1)).all().all():
        raise ValueError("SKAX-instrumentmarkeringen moeten 0 of 1 zijn.")
    if frame[["well_point_x", "well_point_y", "iteration"]].ne(0).any().any():
        raise ValueError("Analyse ondersteunt alleen enkelvoudige absorbantiemetingen, zonder kinetiek of multipointmetingen.")
    frame["saturated"] = frame["saturated"].astype(bool)
    frame["disabled"] = [index in disabled for index in range(len(frame))]
    return frame
