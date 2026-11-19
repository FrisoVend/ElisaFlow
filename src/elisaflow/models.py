from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from enum import Enum, StrEnum
from math import isfinite, isnan
from types import MappingProxyType
from typing import TYPE_CHECKING, Literal

import numpy as np

if TYPE_CHECKING:
    from elisaflow.calibration.curves import CalibrationCurve, CalibrationParameters
    from elisaflow.configuration.assay_settings import AssaySettings, CalibrationSettings, FormulaSettings


class CalculationMethod(Enum):
    FORMULA = "formula"
    CALIBRATION = "calibration"


class ResponseScale(StrEnum):
    """Shared extinction input, after reference and saturation handling."""

    RAW = "raw"
    BLANK_CORRECTED = "blank_corrected"

    @classmethod
    def parse(cls, value: object) -> ResponseScale:
        if isinstance(value, str):
            try:
                return cls(value)
            except ValueError:
                pass
        raise ValueError("De responsschaal (response_scale) moet 'raw' of 'blank_corrected' zijn.")

    @classmethod
    def default_for_blank(cls, has_blank: bool) -> ResponseScale:
        return cls.BLANK_CORRECTED if has_blank else cls.RAW

    def value_for(self, well: WellMeasurement) -> float:
        return well.analysis_extinction if self == self.RAW else well.extinction

    def statistics(self, sample: Sample) -> Statistics:
        return Statistics.calculate(self.value_for(well) for well in sample.retained)


class SampleType(Enum):
    BLANK = "Blank"
    UNKNOWN = "Unknown"
    STANDARD = "Standard"
    CONTROL = "Control"


class ValidationMetric(Enum):
    EXTINCTION = "extinction"
    RESULT = "result"


class ValidationStatus(Enum):
    PASSED = "passed"
    FAILED = "failed"
    NOT_EVALUATED = "not_evaluated"


class ResultStatus(Enum):
    IN_RANGE = "in_range"
    BELOW_RANGE = "below_range"
    ABOVE_RANGE = "above_range"
    UNAVAILABLE = "unavailable"


@dataclass(frozen=True)
class AssayMetaData:
    assay_name: str
    session_time: str
    result_step_name: str
    assay_description: str | None
    instrument_name: str | None
    software_name: str | None
    software_version: str | None
    measurement_wavelengths: tuple[str, ...]

    def __post_init__(self) -> None:
        for name in ("assay_name", "session_time", "result_step_name"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"Assaymetadataveld {name!r} moet een niet-lege tekst zijn.")
        object.__setattr__(self, "measurement_wavelengths", tuple(self.measurement_wavelengths))


@dataclass(frozen=True)
class Statistics:
    mean: float = float("nan")
    sd: float = float("nan")
    cv_percent: float = float("nan")

    @classmethod
    def calculate(cls, values: Iterable[float | None]) -> Statistics:
        values = np.asarray(tuple(values), dtype=float)
        if not values.size or not np.isfinite(values).all():
            return cls()
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            mean = float(np.mean(values))
            sd = float(np.std(values, ddof=1)) if values.size > 1 else float("nan")
            cv = 100 * sd / abs(mean) if mean else float("nan")
        return cls(mean, sd, cv)


@dataclass(frozen=True)
class ResultQualification:
    """Consistent numeric, bounded or unavailable meaning shared by result callers."""

    status: ResultStatus = ResultStatus.UNAVAILABLE
    limit: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.status, ResultStatus):
            raise ValueError("Een resultaatkwalificatie vereist een ResultStatus.")
        if self.status in (ResultStatus.BELOW_RANGE, ResultStatus.ABOVE_RANGE):
            if self.limit is None or not isfinite(self.limit):
                raise ValueError("Een resultaatgrens moet eindig zijn.")
        elif self.limit is not None:
            raise ValueError("Alleen een resultaat buiten bereik heeft een resultaatgrens.")

    @classmethod
    def from_measurements(cls, measurements: Iterable[WellMeasurement]) -> ResultQualification:
        qualifications = []
        for well in measurements:
            if well.result_status == ResultStatus.IN_RANGE and isfinite(well.result) and isnan(well.result_limit):
                qualifications.append(cls(ResultStatus.IN_RANGE))
            elif (
                well.result_status in (ResultStatus.BELOW_RANGE, ResultStatus.ABOVE_RANGE)
                and isnan(well.result) and isfinite(well.result_limit)
            ):
                qualifications.append(cls(well.result_status, well.result_limit))
            else:
                return cls()
        return cls.combine(qualifications)

    @classmethod
    def combine(cls, qualifications: Iterable[ResultQualification]) -> ResultQualification:
        """A group agrees only when every result has the same qualification and bound."""
        qualifications = set(qualifications)
        return next(iter(qualifications)) if len(qualifications) == 1 else cls()


@dataclass(frozen=True, kw_only=True)
class WellMeasurement:
    """One well at one wavelength; original readings survive every processing step."""

    source_extinction: float
    well: str
    wavelength_nm: int
    plate_name: str = ""
    plate_index: int = 0
    well_x: int | None = None
    well_y: int | None = None
    order: int = 0
    saturated: bool = False
    disabled: bool = False
    reference_invalid: bool = False
    raw_extinction: float | None = None
    reference_extinction: float | None = None
    reference_saturated: bool = False
    reference_disabled: bool = False
    converted_extinction: float | None = None
    conversion_source_wavelength_nm: int | None = None
    included: bool = True
    exclusion_reason: str = ""
    extinction: float = float("nan")
    replicate_mean: float = float("nan")
    replicate_deviation_percent: float = float("nan")
    result: float = float("nan")
    result_status: ResultStatus = ResultStatus.UNAVAILABLE
    result_limit: float = float("nan")
    # Calibration input after optional preprocessing for fixed curves.
    curve_extinction: float | None = None
    back_calculated_concentration: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.well, str) or not self.well.strip():
            raise ValueError("Elke meting vereist een wellcoördinaat.")
        for name in ("saturated", "disabled", "reference_invalid", "reference_saturated", "reference_disabled", "included"):
            if not isinstance(getattr(self, name), (bool, np.bool_)):
                raise ValueError(f"Instrumentmarkering {name!r} mag alleen waar of onwaar bevatten.")
        if (
            isinstance(self.wavelength_nm, bool) or not isinstance(self.wavelength_nm, (int, np.integer)) or self.wavelength_nm <= 0
        ):
            raise ValueError("Metingen vereisen positieve gehele getallen voor wavelength_nm.")
        if not isinstance(self.result_status, ResultStatus):
            raise ValueError("Een meetresultaat vereist een ResultStatus.")

    @property
    def key(self) -> tuple[int, str]:
        return self.plate_index, self.well

    @property
    def analysis_extinction(self) -> float:
        return self.converted_extinction if self.converted_extinction is not None else self.source_extinction

    @property
    def saturation_unresolved(self) -> bool:
        return self.saturated and self.converted_extinction is None


@dataclass(frozen=True)
class Sample:
    """A sample identity and its replicate readings, across one or more wavelengths."""

    sample_id: str
    name: str
    sample_type: SampleType
    measurements: tuple[WellMeasurement, ...]
    dilution: float = 1.0
    concentration: float | None = None
    selected_wavelength_nm: int | None = None
    selection_reason: str = ""
    calculation_input_scale: ExtinctionScale | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.sample_type, SampleType):
            raise ValueError("Een monster vereist een SampleType.")
        readings = tuple(self.measurements)
        if not all(isinstance(reading, WellMeasurement) for reading in readings):
            raise ValueError("Monstermetingen moeten WellMeasurement-objecten zijn.")
        if self.calculation_input_scale is not None and not isinstance(self.calculation_input_scale, ExtinctionScale):
            raise ValueError("De gebruikte rekeninvoer vereist een ExtinctionScale.")
        object.__setattr__(self, "measurements", readings)

    @property
    def key(self) -> tuple[str, str, SampleType, float]:
        return self.sample_id, self.name, self.sample_type, self.dilution

    @property
    def retained(self) -> tuple[WellMeasurement, ...]:
        return tuple(well for well in self.measurements if well.included)

    @property
    def result_measurements(self) -> tuple[WellMeasurement, ...]:
        """Keep curve saturation bounds and conflicts with surviving replicates.

        A stored curve input (including NaN) identifies a completed calibration;
        formula assays continue to use only their retained readings.
        """
        return tuple(
            well for well in self.measurements if well.included or (
                self.sample_type == SampleType.UNKNOWN and well.saturation_unresolved and not well.disabled
                and well.curve_extinction is not None
            )
        )

    @property
    def reportable(self) -> bool:
        return (
            bool(self.retained)
            or self.sample_type == SampleType.UNKNOWN
            and any(not well.disabled and (well.saturated or well.reference_invalid) for well in self.measurements)
        )

    @property
    def result(self) -> Statistics:
        if self.result_qualification.status != ResultStatus.IN_RANGE:
            return Statistics()
        return Statistics.calculate(well.result for well in self.result_measurements)

    @property
    def used_wells(self) -> tuple[str, ...]:
        return tuple(well.well for well in self.retained if well.well)

    @property
    def result_status(self) -> ResultStatus:
        return self.result_qualification.status

    @property
    def result_limit(self) -> float | None:
        return self.result_qualification.limit

    @property
    def result_qualification(self) -> ResultQualification:
        return ResultQualification.from_measurements(self.result_measurements)

    @property
    def curve_extinction(self) -> Statistics:
        return Statistics.calculate(well.curve_extinction for well in self.retained)

    @property
    def back_calculated_concentration(self) -> Statistics:
        return Statistics.calculate(well.back_calculated_concentration for well in self.retained)

    @property
    def recovery_percent(self) -> float:
        if self.concentration is None or self.concentration <= 0:
            return float("nan")
        return 100 * self.back_calculated_concentration.mean / self.concentration


@dataclass(frozen=True)
class ExtinctionScale:
    """Mathematical meaning of the input actually used by a calculation."""

    stage: Literal["formula_input", "curve_input"]
    wavelength_nm: int
    reference_wavelength_nm: int | None
    blank_subtracted: float | None
    preprocessing_expression: str | None
    response_scale: ResponseScale = ResponseScale.BLANK_CORRECTED

    def __post_init__(self) -> None:
        if self.stage not in ("formula_input", "curve_input"):
            raise ValueError("Onbekende schaal voor de rekeninvoer.")
        object.__setattr__(self, "response_scale", ResponseScale.parse(self.response_scale))
        if self.response_scale == ResponseScale.RAW and self.blank_subtracted is not None:
            raise ValueError("Ruwe rekeninvoer gebruikt geen blancoaftrek.")
        if self.blank_subtracted is not None and not isfinite(self.blank_subtracted):
            raise ValueError("Het afgetrokken blancogemiddelde moet eindig zijn.")
        if self.stage == "formula_input" and self.preprocessing_expression is not None:
            raise ValueError("Kalibratievoorbewerking hoort alleen bij curve-invoer.")

    @classmethod
    def for_calculation(
        cls, calculation: FormulaSettings | CalibrationSettings, *, wavelength_nm: int,
        reference_wavelength_nm: int | None, blank_value: float | None, control_name: str | None = None,
    ) -> ExtinctionScale:
        """Resolve input semantics once, including an independent control formula scale."""
        calibration = calculation.method == CalculationMethod.CALIBRATION
        response_scale = (
            calculation.response_scale_for_control(control_name)
            if not calibration and control_name is not None else calculation.response_scale
        )
        return cls(
            stage="curve_input" if calibration else "formula_input",
            wavelength_nm=wavelength_nm,
            reference_wavelength_nm=reference_wavelength_nm,
            blank_subtracted=blank_value if response_scale == ResponseScale.BLANK_CORRECTED else None,
            preprocessing_expression=calculation.preprocessing.expression if calibration and calculation.preprocessing is not None else None,
            response_scale=response_scale,
        )

    def statistics(self, sample: Sample) -> Statistics:
        """Use retained wells on this exact scale, without substituting other values."""
        return sample.curve_extinction if self.stage == "curve_input" else self.response_scale.statistics(sample)

    def summarize(self, sample: Sample) -> ScaledStatistics | None:
        """Missing executed curve input stays missing; partial input yields unavailable statistics."""
        if self.stage == "curve_input" and not any(well.curve_extinction is not None for well in sample.retained):
            return None
        if any(well.wavelength_nm != self.wavelength_nm for well in sample.retained):
            return None
        return ScaledStatistics(self.statistics(sample), self)


@dataclass(frozen=True)
class ScaledStatistics:
    statistics: Statistics
    scale: ExtinctionScale


@dataclass(frozen=True)
class AssayData:
    samples: tuple[Sample, ...]
    metadata: AssayMetaData

    def __post_init__(self) -> None:
        samples = tuple(self.samples)
        if not all(isinstance(sample, Sample) for sample in samples):
            raise ValueError("Assaygegevens moeten Sample-objecten bevatten.")
        object.__setattr__(self, "samples", samples)


@dataclass(frozen=True)
class CalculationDetails:
    """Executed formula or calibration, with its original curve and settings."""

    formula_settings: FormulaSettings | None = None
    curve: CalibrationCurve[CalibrationParameters] | None = None
    calibration_settings: CalibrationSettings | None = None
    curve_r_squared: float | None = None
    calibration_error: str | None = None


@dataclass(frozen=True)
class ValidationResult:
    """The outcome of one configured sample/metric acceptance rule."""

    sample_type: SampleType
    sample_name: str
    metric: ValidationMetric
    value: float | None
    minimum: float | None
    maximum: float | None
    status: ValidationStatus
    message: str
    result_status: ResultStatus | None = None
    result_limit: float | None = None
    wavelength_nm: int | None = None
    nominal_concentration: float | None = None
    cv_percent: float | None = None
    response_scale: ResponseScale | None = None

    @property
    def sample_key(self) -> tuple[SampleType, str | float]:
        identity = self.nominal_concentration if self.sample_type == SampleType.STANDARD else self.sample_name
        return self.sample_type, identity


@dataclass(frozen=True, kw_only=True)
class ResultData:
    """Immutable processed samples and QC checks for one processing scope."""

    samples: tuple[Sample, ...]
    calculation_details: CalculationDetails
    blank_value: float | None = None
    validation_results: tuple[ValidationResult, ...] = ()

    def __post_init__(self) -> None:
        samples = tuple(self.samples)
        if not all(isinstance(sample, Sample) for sample in samples):
            raise ValueError("Analyseresultaten moeten Sample-objecten bevatten.")
        object.__setattr__(self, "samples", samples)
        object.__setattr__(self, "validation_results", tuple(self.validation_results))


@dataclass(frozen=True)
class MeasurementResult(ResultData):
    """One independently processed measurement; never a nested assay run."""

    wavelength_nm: int
    label: str
    calculation_input_scale: ExtinctionScale | None = None

    def calculation_extinction(self, sample: Sample) -> ScaledStatistics | None:
        """Summarize this sample's actual calculation input on its recorded scale."""
        scale = sample.calculation_input_scale
        if scale is None or self.calculation_input_scale is None:
            return None
        if scale.wavelength_nm != self.wavelength_nm or scale.stage != self.calculation_input_scale.stage:
            return None
        return scale.summarize(sample)

    @property
    def usable_for_results(self) -> bool:
        """Calibration availability is independent of QC acceptance."""
        return self.calculation_details.calibration_error is None


@dataclass(frozen=True)
class AnalysisResult:
    """Compose selected results and independent measurements for one assay."""

    metadata: AssayMetaData
    settings: AssaySettings
    selected: ResultData
    measurement_results: Mapping[int, MeasurementResult]

    def __post_init__(self) -> None:
        if not self.measurement_results:
            raise ValueError("Een analyseresultaat vereist ten minste één meting.")
        object.__setattr__(self, "measurement_results", MappingProxyType(dict(self.measurement_results)))

    @property
    def samples(self) -> tuple[Sample, ...]:
        return self.selected.samples

    def patient_extinction(self, sample: Sample) -> ScaledStatistics | None:
        """Resolve the selected scale while using the supplied sample's current wells."""
        measurement = (
            next(iter(self.measurement_results.values()))
            if len(self.measurement_results) == 1
            else self.measurement_results.get(sample.selected_wavelength_nm)
        )
        return measurement.calculation_extinction(sample) if measurement is not None else None

    @property
    def calculation_details(self) -> CalculationDetails:
        return self.selected.calculation_details

    @property
    def blank_value(self) -> float | None:
        return self.selected.blank_value

    @property
    def validation_results(self) -> tuple[ValidationResult, ...]:
        return (
            *self.selected.validation_results,
            *(check for measurement in self.measurement_results.values() for check in measurement.validation_results),
        )

    @property
    def validation_status(self) -> ValidationStatus:
        """Apply all configured QC independently of curve selection."""
        statuses = {result.status for result in self.validation_results}
        if ValidationStatus.FAILED in statuses:
            return ValidationStatus.FAILED
        if not statuses or ValidationStatus.NOT_EVALUATED in statuses:
            return ValidationStatus.NOT_EVALUATED
        if self.settings.measurement.switching and not any(measurement.usable_for_results for measurement in self.measurement_results.values()):
            return ValidationStatus.NOT_EVALUATED
        return ValidationStatus.PASSED
