"""Sample identities supplied by a SkanIt protocol, separate from QC settings."""

from collections.abc import Mapping
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
from types import MappingProxyType
from xml.etree import ElementTree as ET
from zipfile import BadZipFile

from elisaflow.configuration.assay_settings import AssaySettings, FormulaSettings, WavelengthCalibrationSettings
from elisaflow.models import SampleType
from elisaflow.readers.skax_reader import SkaxReader


@dataclass(frozen=True)
class ProtocolSamples:
    controls: tuple[str, ...]
    standards: Mapping[float, tuple[str, ...]]
    wavelengths_nm: tuple[int, ...]
    has_blank: bool = False

    @classmethod
    def from_file(cls, path: Path) -> ProtocolSamples:
        try:
            protocol = SkaxReader.from_file(path)
            wells = protocol.wells
        except (BadZipFile, ET.ParseError, KeyError) as error:
            raise ValueError(f"Het SKAX-protocol {path.name} kon niet worden gelezen. Controleer of het bestand een geldig SkanIt-protocol is.") from error

        controls = {}
        standards: dict[float, list[str]] = {}
        samples = wells.loc[
            wells["sample_type"].isin((SampleType.CONTROL.value, SampleType.STANDARD.value)),
            ["sample_type", "sample_name", "sample_concentration"],
        ].drop_duplicates()
        for sample_type, name, concentration in samples.itertuples(index=False, name=None):
            if not isinstance(name, str) or not name.strip():
                raise ValueError(f"Namen van SKAX-monsters ({sample_type.lower()}) moeten niet-lege teksten zijn.")
            if sample_type == SampleType.CONTROL.value:
                controls[name] = None
            else:
                concentration = float(concentration)
                if not isfinite(concentration) or concentration < 0:
                    raise ValueError(f"SKAX-standaard {name!r} vereist een eindige, niet-negatieve concentratie.")
                standards.setdefault(concentration, []).append(name)
        return cls(
            tuple(controls),
            MappingProxyType({value: tuple(names) for value, names in standards.items()}),
            protocol.wavelengths_nm,
            bool(wells["sample_type"].eq(SampleType.BLANK.value).any()),
        )

    def validate_measurements(self, settings: AssaySettings) -> None:
        analysis = settings.measurement.analysis_wavelengths(self.wavelengths_nm)
        if isinstance(settings.calculation, WavelengthCalibrationSettings) and set(settings.calculation.wavelengths) != set(analysis):
            raise ValueError("De kalibratiegolflengten moeten overeenkomen met de SKAX-analysegolflengten.")
        if set(settings.validation.wavelengths) - set(analysis):
            raise ValueError("De validatiegolflengten moeten bij de SKAX-analysegolflengten horen.")

    def validate_samples(self, settings: AssaySettings) -> None:
        """Do not silently discard QC settings for identities absent from SKAX."""
        configured_controls = set(settings.controls)
        configured_standards = set(settings.standards)
        for rules in (settings.validation.result, *settings.validation.wavelengths.values()):
            configured_controls.update(rules.controls)
            configured_standards.update(rules.standards)
        if isinstance(settings.calculation, FormulaSettings):
            configured_controls.update(settings.calculation.control_overrides)
            configured_controls.update(settings.calculation.control_results)
        missing_controls = configured_controls - set(self.controls)
        if missing_controls:
            names = ", ".join(repr(name) for name in sorted(missing_controls))
            raise ValueError(f"Ingestelde controles ontbreken in het SKAX-protocol: {names}.")
        missing_standards = configured_standards - set(self.standards)
        if missing_standards:
            values = ", ".join(str(value) for value in sorted(missing_standards))
            raise ValueError(f"Ingestelde standaardconcentraties ontbreken in het SKAX-protocol: {values}.")

    def validate_settings(self, settings: AssaySettings) -> None:
        """Every control reference must resolve to a control in the protocol."""
        self.validate_samples(settings)
        self.validate_measurements(settings)
        if isinstance(settings.calculation, FormulaSettings):
            formulas = [(settings.calculation.formula, "Resultaatformule")]
        elif isinstance(settings.calculation, WavelengthCalibrationSettings):
            formulas = [
                (curve.preprocessing, f"Kalibratievoorbewerking bij {wavelength} nm")
                for wavelength, curve in settings.calculation.wavelengths.items()
            ]
        else:
            formulas = [(settings.calculation.preprocessing, "Kalibratievoorbewerking")]
        for formula, label in formulas:
            if formula is not None:
                for name in formula.references.values():
                    if name not in self.controls:
                        raise ValueError(f"{label} verwijst naar een controle die ontbreekt in het SKAX-protocol: {name!r}.")
