"""Validate raw filter metadata before starting a live SkanIt measurement."""

from pathlib import Path

from elisaflow.configuration.assay_settings import MeasurementSettings
from elisaflow.readers.skax_document import SkaxDocument


def validate_protocol_measurement(path: Path, settings: MeasurementSettings) -> None:
    document = SkaxDocument.from_file(path)
    document.validate_hashes()
    if document.measurements is not None:
        raise ValueError("Gebruik voor een live meting een protocol zonder meetgegevens; heranalyseer een bestaande meting via de experimentmap.")
    settings.analysis_wavelengths(document.wavelengths_nm)
