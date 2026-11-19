"""Coordinate preparation, calculation, selection and QC for typed samples."""

from dataclasses import replace
from math import isfinite, isinf

import numpy as np

from elisaflow.analysis.calculation import calculate_samples
from elisaflow.analysis.sample_updates import clear_result, map_wells
from elisaflow.analysis.wavelengths import measurement_inputs, select_measurement_results
from elisaflow.calibration.errors import CalibrationError
from elisaflow.configuration.assay_settings import AssaySettings, ReplicateFilterSettings, ValidationScope, WavelengthCalibrationSettings
from elisaflow.models import (
    AnalysisResult,
    AssayData,
    CalculationDetails,
    CalculationMethod,
    MeasurementResult,
    ResponseScale,
    ResultData,
    Sample,
    SampleType,
)
from elisaflow.quality_control.validation import AssayValidator

NAN = float("nan")
_EMPTY_WELL_THRESHOLD = 0.009


class AssayAnalyzer:
    def __init__(self, data: AssayData, settings: AssaySettings):
        self.data = data
        self.settings = settings

    def analyze(self) -> AnalysisResult:
        inputs = measurement_inputs(self.data, self.settings.measurement)
        calculation = self.settings.calculation
        if isinstance(calculation, WavelengthCalibrationSettings) and set(calculation.wavelengths) != set(inputs):
            raise ValueError("De kalibratiegolflengten moeten overeenkomen met de gemeten golflengten in calculation.wavelengths.")
        validation = self.settings.validation
        unknown = set(validation.wavelengths) - set(inputs)
        if unknown:
            raise ValueError(f"De golflengten voor validatie ontbreken in de metingen: {', '.join(map(str, sorted(unknown)))} nm.")
        prepared = self._prepare_measurements(inputs)
        measurements = {}
        for key, measurement in prepared.items():
            try:
                measurements[key] = self._analyze_measurement(measurement, key)
            except CalibrationError as error:
                raise type(error)(f"Kalibratie/analyse bij {key} nm: {error}") from error
        samples = select_measurement_results(measurements, self.settings.measurement)
        preferred = self.settings.measurement.preferred_wavelength_nm if self.settings.measurement.switching else next(iter(measurements))
        primary = measurements[preferred]
        checks = AssayValidator(validation.result).validate(samples, blank_value=primary.blank_value)
        selected = ResultData(
            samples=samples,
            calculation_details=CalculationDetails() if self.settings.measurement.switching else primary.calculation_details,
            blank_value=primary.blank_value,
            validation_results=checks,
        )
        return AnalysisResult(self.data.metadata, self.settings, selected, measurements)

    def _prepare_measurements(self, inputs: dict[int, AssayData]) -> dict[int, AssayData]:
        """Prepare flags before empty-well selection, always before blank correction."""
        prepared = {key: replace(data, samples=self._prepare_samples(data.samples)) for key, data in inputs.items()}
        settings = self.settings.measurement
        preferred = settings.preferred_wavelength_nm if settings.switching else next(iter(prepared))
        alternatives = {
            sample.key: sample
            for wavelength, data in prepared.items() if wavelength != preferred
            for sample in data.samples
        }
        retained = set()
        for sample in prepared[preferred].samples:

            def empty(well):
                return (
                    sample.sample_type == SampleType.UNKNOWN
                    and well.included
                    and isfinite(well.analysis_extinction)
                    and well.analysis_extinction < _EMPTY_WELL_THRESHOLD
                )

            if settings.switching:
                if all(empty(well) for well in (*sample.measurements, *alternatives[sample.key].measurements)):
                    continue
                retained.update(well.key for well in sample.measurements)
            else:
                retained.update(well.key for well in sample.measurements if not empty(well))
        for key, data in prepared.items():
            samples = tuple(
                replace(sample, measurements=tuple(well for well in sample.measurements if well.key in retained)) for sample in data.samples
            )
            prepared[key] = replace(data, samples=tuple(sample for sample in samples if sample.measurements))
        return prepared

    def _prepare_samples(self, samples: tuple[Sample, ...]) -> tuple[Sample, ...]:
        prepared = []
        for sample in samples:
            if isinstance(sample.dilution, (bool, np.bool_)):
                raise ValueError("Verdunningsfactoren moeten getallen zijn, geen waarden waar of onwaar.")
            try:
                dilution = float(sample.dilution)
            except (TypeError, ValueError) as error:
                raise ValueError("Verdunningsfactoren moeten numeriek zijn.") from error
            patient = sample.sample_type == SampleType.UNKNOWN
            if not isfinite(dilution) or patient and dilution < 1:
                raise ValueError("Verdunningsfactoren voor patiëntmonsters moeten eindig en ten minste 1 zijn.")
            if not patient and dilution not in (0, 1):
                raise ValueError("Verdunning wordt alleen ondersteund voor patiëntmonsters (Unknown).")
            if patient and self.settings.calculation.method == CalculationMethod.FORMULA and dilution != 1:
                raise ValueError("Formuleanalyse vereist een patiëntverdunning van 1; er is geen assayspecifieke verdunningscorrectie ingesteld.")
            wells = []
            for well in sample.measurements:
                reason = ""
                for flag, message in (
                    (well.saturation_unresolved, "Meting verzadigd in SkanIt, geen conversie uit andere golflengte mogelijk"),
                    (well.disabled, "Well uitgeschakeld in SkanIt."),
                    (well.reference_invalid, "Achtergrond-golflengte mist, is verzadigd of is uitgeschakeld."),
                ):
                    if flag:
                        reason = message
                wells.append(replace(
                    clear_result(well),
                    included=not bool(reason), exclusion_reason=reason,
                    extinction=NAN, replicate_mean=NAN, replicate_deviation_percent=NAN,
                    curve_extinction=None, back_calculated_concentration=None,
                ))
            prepared.append(replace(sample, dilution=dilution if patient else 1.0, measurements=tuple(wells), calculation_input_scale=None))
        return tuple(prepared)

    def _analyze_measurement(self, data: AssayData, wavelength_nm: int) -> MeasurementResult:
        samples, blank_value = self._correct_blank(data.samples)
        samples = self._filter_replicates(samples)
        calculation = calculate_samples(
            samples, self.settings.calculation_for_wavelength(wavelength_nm),
            wavelength_nm=wavelength_nm,
            reference_wavelength_nm=self.settings.measurement.reference_wavelength_nm,
            blank_value=blank_value,
            allow_fit_failure=self.settings.measurement.switching,
        )
        scope = self.settings.validation.wavelengths.get(wavelength_nm, ValidationScope())
        checks = AssayValidator(scope).validate(calculation.samples, blank_value=blank_value, wavelength_nm=wavelength_nm)
        return MeasurementResult(
            samples=calculation.samples,
            wavelength_nm=wavelength_nm,
            label=data.metadata.result_step_name,
            calculation_details=calculation.details,
            blank_value=blank_value,
            validation_results=checks,
            calculation_input_scale=calculation.input_scale,
        )

    @staticmethod
    def _correct_blank(samples: tuple[Sample, ...]) -> tuple[tuple[Sample, ...], float | None]:
        blanks = tuple(well for sample in samples if sample.sample_type == SampleType.BLANK for well in sample.retained)
        if any(isinf(well.analysis_extinction) for sample in samples for well in sample.retained):
            raise ValueError("Bronmetingen mogen niet oneindig zijn.")
        blank_value = None
        if any(sample.sample_type == SampleType.BLANK for sample in samples):
            if not blanks:
                raise ValueError("Er zijn geen bruikbare blancometingen over.")
            if not all(isfinite(well.analysis_extinction) for well in blanks):
                raise ValueError("Blancometingen moeten eindig zijn.")
            blank_value = float(np.mean([well.analysis_extinction for well in blanks]))
        corrected = map_wells(samples, lambda _, well: replace(well, extinction=well.analysis_extinction - (blank_value or 0) if well.included else NAN))
        return corrected, blank_value

    def _filter_replicates(self, samples: tuple[Sample, ...]) -> tuple[Sample, ...]:
        patient_rule = self.settings.patients.replicate_filter
        if patient_rule is not None:
            samples = tuple(
                _filter_replicate_group(
                    (sample,), f"Patiëntmonster {sample.name!r} (ID {sample.sample_id!r}, verdunning {sample.dilution:g})", patient_rule,
                )[0] if sample.sample_type == SampleType.UNKNOWN else sample
                for sample in samples
            )
        for kind, groups in ((SampleType.CONTROL, self.settings.controls), (SampleType.STANDARD, self.settings.standards)):
            for identity, settings in groups.items():
                if settings.replicate_filter is None:
                    continue
                selected = tuple(
                    sample for sample in samples
                    if sample.sample_type == kind and (sample.name == identity if kind == SampleType.CONTROL else sample.concentration == identity)
                )
                label = f"Controle {identity!r}" if kind == SampleType.CONTROL else f"Standaard {identity:g}"
                updates = {sample.key: sample for sample in _filter_replicate_group(selected, label, settings.replicate_filter)}
                samples = tuple(updates.get(sample.key, sample) for sample in samples)
        return samples


def _filter_replicate_group(samples: tuple[Sample, ...], label: str, settings: ReplicateFilterSettings) -> tuple[Sample, ...]:
    """Compare retained readings once with their original group mean."""
    wells = tuple(well for sample in samples for well in sample.retained)
    if samples and not wells:
        return samples
    if not wells:
        raise ValueError(f"{label}: geen metingen gevonden voor replicaatfiltering.")
    values = np.array([settings.response_scale.value_for(well) for well in wells])
    if not np.isfinite(values).all():
        raise ValueError(f"{label}: extincties van replicaten moeten eindig zijn.")
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        mean = float(np.mean(values))
        deviations = 100 * np.abs(values - mean) / abs(mean)
    if not isfinite(mean) or mean == 0:
        raise ValueError(f"{label}: replicaatfiltering vereist een eindige gemiddelde extinctie ongelijk aan nul.")
    if not np.isfinite(deviations).all():
        raise ValueError(f"{label}: afwijkingen van replicaten moeten eindig zijn.")
    excluded = (deviations > settings.max_deviation_percent) & ~np.isclose(deviations, settings.max_deviation_percent, rtol=0, atol=1e-12)
    if excluded.sum() > settings.max_exclusions:
        raise ValueError(
            f"{label}: {excluded.sum()} replicaten overschrijden de afwijkingsgrens; er mogen maximaal {settings.max_exclusions} worden uitgesloten."
        )
    if excluded.all():
        raise ValueError(f"{label}: replicaatfiltering mag niet alle metingen uitsluiten.")
    scale_label = "vóór blancoaftrek" if settings.response_scale == ResponseScale.RAW else "na blancoaftrek"
    reason = f"Extinctie wijkt meer dan {settings.max_deviation_percent:g}% af van het oorspronkelijke replicaatgemiddelde ({mean:g}; {scale_label})."
    updates = {
        well.key: replace(
            well, replicate_mean=mean, replicate_deviation_percent=float(deviation), included=not bool(exclude), exclusion_reason=reason if exclude else "",
        )
        for well, deviation, exclude in zip(wells, deviations, excluded, strict=True)
    }
    return map_wells(samples, lambda _, well: updates[well.key] if well.included else well)
