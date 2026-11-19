"""Align wavelength readings and choose a complete result for each sample."""

from dataclasses import replace
from math import isfinite

from elisaflow.analysis.sample_updates import clear_result, map_wells
from elisaflow.configuration.assay_settings import MeasurementSettings, SaturationSettings
from elisaflow.models import AssayData, MeasurementResult, ResultQualification, ResultStatus, Sample, SampleType, WellMeasurement

_WellMap = dict[tuple[int, str], tuple[Sample, WellMeasurement]]


def _well_map(samples: tuple[Sample, ...]) -> _WellMap:
    readings = [(well.key, (sample, well)) for sample in samples for well in sample.measurements]
    mapping = dict(readings)
    if len(mapping) != len(readings):
        raise ValueError("Dubbele wellmetingen bij dezelfde golflengte.")
    return mapping


def _aligned(primary: _WellMap, other: _WellMap, *, reference: bool = False) -> None:
    if primary.keys() != other.keys():
        if reference:
            raise ValueError("De analyse- en referentiegolflengten moeten dezelfde wells bevatten.")
        raise ValueError("Beide golflengten moeten dezelfde protocolwells bevatten.")
    for key, (sample, well) in primary.items():
        other_sample, other_well = other[key]
        if (sample.key, sample.concentration, well.plate_name) != (other_sample.key, other_sample.concentration, other_well.plate_name):
            target = "overeenkomen met de referentiegolflengte" if reference else "tussen golflengten overeenkomen"
            raise ValueError(f"Monsteridentiteiten, verdunningen en standaardconcentraties moeten {target}.")


def _convert_saturation(channel: AssayData, source: _WellMap, rule: SaturationSettings) -> AssayData:
    """Convert unsaturated sources after reference correction and before blanking."""
    def convert(_sample: Sample, well: WellMeasurement) -> WellMeasurement:
        if well.saturated and not well.disabled and not well.reference_invalid:
            donor = source[well.key][1]
            if not donor.saturated and not donor.disabled and not donor.reference_invalid and isfinite(donor.source_extinction):
                converted = rule.intercept + rule.slope * donor.source_extinction
                if isfinite(converted):
                    return replace(well, converted_extinction=converted, conversion_source_wavelength_nm=rule.source_wavelength_nm)
        return well

    return replace(channel, samples=map_wells(channel.samples, convert))


def _subtract_reference(channel: AssayData, references: _WellMap) -> AssayData:
    def subtract(_sample: Sample, well: WellMeasurement) -> WellMeasurement:
        ref = references[well.key][1]
        invalid = not isfinite(ref.source_extinction) or ref.saturated or ref.disabled
        return replace(
            well,
            raw_extinction=well.source_extinction,
            reference_extinction=ref.source_extinction,
            reference_saturated=ref.saturated,
            reference_disabled=ref.disabled,
            reference_invalid=bool(invalid),
            source_extinction=well.source_extinction - ref.source_extinction if not invalid else float("nan"),
        )

    return replace(channel, samples=map_wells(channel.samples, subtract))


def measurement_inputs(data: AssayData, settings: MeasurementSettings) -> dict[int, AssayData]:
    wavelengths = tuple(dict.fromkeys(well.wavelength_nm for sample in data.samples for well in sample.measurements))
    settings.analysis_wavelengths(wavelengths)
    channels = {wavelength: [] for wavelength in wavelengths}
    for sample in data.samples:
        grouped = {}
        for well in sample.measurements:
            grouped.setdefault(well.wavelength_nm, []).append(well)
        for wavelength, wells in grouped.items():
            channels[wavelength].append(replace(sample, measurements=tuple(wells)))
    inputs = {
        wavelength: AssayData(tuple(samples), replace(
            data.metadata, result_step_name=settings.wavelength_label(wavelength), measurement_wavelengths=(str(wavelength),),
        ))
        for wavelength, samples in channels.items()
    }
    reference = settings.reference_wavelength_nm
    # Align all raw channels before reference correction or saturation conversion.
    originals = {wavelength: _well_map(channel.samples) for wavelength, channel in inputs.items()}
    primary = originals[next(wavelength for wavelength in inputs if wavelength != reference)]
    for wavelength, wells in originals.items():
        if wells is not primary:
            _aligned(primary, wells, reference=wavelength == reference)
    if reference is not None:
        inputs.pop(reference)
        inputs = {wavelength: _subtract_reference(channel, originals[reference]) for wavelength, channel in inputs.items()}
    # Snapshot unconverted, reference-corrected sources to prevent conversion chains.
    sources = {wavelength: _well_map(channel.samples) for wavelength, channel in inputs.items()} if reference is not None else originals
    for wavelength, rule in settings.saturation.items():
        if rule.mode == "calculate":
            inputs[wavelength] = _convert_saturation(inputs[wavelength], sources[rule.source_wavelength_nm], rule)
    return inputs


def select_measurement_results(results: dict[int, MeasurementResult], settings: MeasurementSettings) -> tuple[Sample, ...]:
    """Choose complete sample results after independently analyzing each wavelength."""
    if not settings.switching:
        return next(iter(results.values())).samples
    preferred = settings.preferred_wavelength_nm
    alternate = next(wavelength for wavelength in results if wavelength != preferred)
    primary, other = results[preferred], results[alternate]
    alternatives = {sample.key: sample for sample in other.samples}
    primary_usable, alternate_usable = primary.usable_for_results, other.usable_for_results
    selected = []
    for sample in primary.samples:
        chosen, reason = _select_sample(
            sample, alternatives[sample.key], settings.cutoff_concentration,
            primary_usable=primary_usable, alternate_usable=alternate_usable,
        )
        if chosen is None:
            sample = replace(
                sample,
                measurements=tuple(map(clear_result, sample.measurements)),
                selected_wavelength_nm=None,
                selection_reason="no_usable_result",
            )
        else:
            wavelength = preferred if chosen is sample else alternate
            order = {well.key: well.order for well in sample.measurements}
            sample = replace(
                chosen,
                selected_wavelength_nm=wavelength,
                selection_reason=reason,
                measurements=tuple(replace(well, order=order[well.key]) for well in chosen.measurements),
            )
        selected.append(sample)
    return tuple(selected)


def _select_sample(
    preferred: Sample, alternate: Sample, cutoff: float, *, primary_usable: bool, alternate_usable: bool,
) -> tuple[Sample | None, str]:
    """Apply the concentration cutoff first, then handle unusable results and bounds."""
    outcome = _channel_result(preferred) if primary_usable else ResultQualification()
    other_outcome = _channel_result(alternate) if alternate_usable else ResultQualification()
    status, limit = outcome.status, outcome.limit
    other_status, other_limit = other_outcome.status, other_outcome.limit
    dilution = preferred.dilution if preferred.sample_type == SampleType.UNKNOWN else 1.0
    if any(isfinite(well.result) and well.result / dilution >= cutoff for well in preferred.retained):
        return (alternate if other_status != ResultStatus.UNAVAILABLE else None), "concentration_cutoff"
    if status == ResultStatus.IN_RANGE:
        return preferred, "preferred"
    if other_status == ResultStatus.IN_RANGE:
        if not primary_usable:
            reason = "preferred_curve_unusable"
        elif status in (ResultStatus.BELOW_RANGE, ResultStatus.ABOVE_RANGE):
            reason = "preferred_out_of_range"
        else:
            reason = "preferred_measurement_unavailable"
        return alternate, reason
    if status == ResultStatus.UNAVAILABLE:
        return (alternate if other_status != ResultStatus.UNAVAILABLE else None), other_status.value
    if other_status == ResultStatus.UNAVAILABLE:
        return preferred, status.value
    if status == other_status:
        use_alternate = other_limit > limit if status == ResultStatus.ABOVE_RANGE else other_limit < limit
        return (alternate if use_alternate else preferred), status.value
    return None, "no_usable_result"


def _channel_result(sample: Sample) -> ResultQualification:
    if any(
        (well.reference_invalid or well.saturation_unresolved and well.result_status == ResultStatus.UNAVAILABLE) and not well.disabled
        for well in sample.measurements
    ):
        return ResultQualification()
    return sample.result_qualification
