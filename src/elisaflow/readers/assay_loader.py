"""Load an executed SKAX session as typed assay measurements."""

from collections import defaultdict
from dataclasses import replace
from math import isnan
from pathlib import Path

import pandas as pd

from elisaflow.configuration.assay_settings import MeasurementSettings
from elisaflow.models import AssayData, AssayMetaData, Sample, SampleType, WellMeasurement
from elisaflow.readers.skax_reader import SkaxReader


class AssayLoader:
    @staticmethod
    def load(session_path: Path, *, measurement_settings: MeasurementSettings) -> AssayData:
        reader = SkaxReader.from_file(session_path)
        if not reader.completed:
            raise ValueError("Analyse vereist een opgeslagen, voltooide SKAX-meting met meetgegevens en uitvoeringstijden.")
        measurement_settings.analysis_wavelengths(reader.wavelengths_nm)
        readings = reader.measurements
        if readings["plate_index"].nunique(dropna=False) != 1:
            raise ValueError("Analyse vereist precies een gemeten plaat.")
        wavelengths = tuple(dict.fromkeys(readings["wavelength_nm"]))
        if set(wavelengths) != set(reader.wavelengths_nm):
            raise ValueError("De opgeslagen meting moet alle ruwe protocolgolflengten bevatten, inclusief een eventuele referentie.")
        recorded = set(readings[["plate_index", "well", "wavelength_nm"]].itertuples(index=False, name=None))
        required = reader.wells.loc[reader.wells["sample_type"].isin(("Blank", "Control", "Standard"))]
        for well in required.itertuples(index=False):
            for wavelength in reader.wavelengths_nm:
                if (well.plate_index, well.well, wavelength) not in recorded:
                    raise ValueError(
                        f"Meetrecord ontbreekt: {well.well} bij {wavelength} nm (monster {well.sample_name!r}, plaat {well.plate_index})."
                    )
        session = reader.metadata["session"]
        instrument, software = session["instrument"], session["software"]
        metadata = AssayMetaData(
            assay_name=session["name"],
            session_time=session["execution"]["start"],
            result_step_name=" / ".join(f"{wavelength} nm" for wavelength in wavelengths),
            assay_description=session["description"],
            instrument_name=instrument["name"],
            software_name=software["name"], software_version=software["version"],
            measurement_wavelengths=tuple(str(wavelength) for wavelength in wavelengths),
        )
        return AssayData(samples=_samples_from_measurements(readings), metadata=metadata)


def _samples_from_measurements(readings: pd.DataFrame) -> tuple[Sample, ...]:
    groups = {}
    measurements = defaultdict(list)
    for order, record in enumerate(readings.itertuples(index=False)):
        sample = Sample(
            sample_id=record.sample_id,
            name=record.sample_name,
            sample_type=SampleType(record.sample_type),
            measurements=(),
            dilution=record.sample_dilution,
            concentration=None if isnan(record.sample_concentration) else record.sample_concentration,
        )
        groups.setdefault(sample.key, sample)
        measurements[sample.key].append(WellMeasurement(
            source_extinction=record.source_extinction,
            well=record.well,
            wavelength_nm=record.wavelength_nm,
            plate_name="" if pd.isna(record.plate_name) else record.plate_name,
            plate_index=record.plate_index,
            well_x=record.well_x,
            well_y=record.well_y,
            order=order,
            saturated=record.saturated,
            disabled=record.disabled,
        ))
    return tuple(replace(sample, measurements=tuple(measurements[key])) for key, sample in groups.items())
