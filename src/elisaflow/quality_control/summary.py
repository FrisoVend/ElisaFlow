"""Unrounded QC statistics and plate-ordered groups for reporting."""

from collections.abc import Iterable
from dataclasses import dataclass, replace
from math import isfinite

from elisaflow.configuration.assay_settings import ValidationCriterion
from elisaflow.models import (
    AnalysisResult,
    MeasurementResult,
    ResultQualification,
    Sample,
    SampleType,
    Statistics,
    ValidationMetric,
    ValidationResult,
)
from elisaflow.sample_selection import select_samples


@dataclass(frozen=True)
class QualityControlSummary:
    """Unrounded acceptance mean, pooled replicate CV and shared result qualification."""

    mean: float | None
    cv_percent: float | None = None
    message: str = ""
    result_qualification: ResultQualification | None = None


def summarize_quality_control(
    samples: Iterable[Sample],
    criterion: ValidationCriterion,
    *,
    blank_value: float | None = None,
) -> QualityControlSummary:
    """Select and summarize the criterion's target on its validated response scale.

    Acceptance averages sample means equally; CV pools retained replicate values.
    Blank extinction uses the stored usable blank mean, including zero, and can
    still be evaluated when the caller supplies no blank sample records.
    """
    samples = criterion.select(samples)
    extinction = criterion.metric == ValidationMetric.EXTINCTION
    scale = criterion.limits.response_scale
    qualification = ResultQualification.combine(sample.result_qualification for sample in samples) if not extinction and samples else None

    def summary(mean: float | None = None, message: str = "") -> QualityControlSummary:
        return QualityControlSummary(mean, message=message, result_qualification=qualification)

    statistics = tuple(scale.statistics(sample) if extinction else sample.result for sample in samples)
    if criterion.uses_blank_value:
        result = summary(
            blank_value if blank_value is not None and isfinite(blank_value) else None,
            "Geen metingen gevonden voor dit monster.",
        )
    elif not samples:
        return summary(message="Geen metingen gevonden voor dit monster.")
    elif not all(isfinite(value.mean) for value in statistics):
        return summary(message="Elk monster moet een eindige gemiddelde meetwaarde voor deze grootheid hebben.")
    else:
        mean = Statistics.calculate(value.mean for value in statistics).mean
        result = summary(mean) if isfinite(mean) else summary(message="De gemiddelde meetwaarde is niet eindig.")
    if not samples or result.mean is None:
        return result
    wavelengths = {sample.selected_wavelength_nm for sample in samples if sample.selected_wavelength_nm is not None}
    if len(wavelengths) > 1 and extinction:
        return summary(message="Extincties van verschillende golflengten kunnen niet worden samengevoegd.")
    if len(samples) == 1 and not criterion.uses_blank_value:
        cv = statistics[0].cv_percent
    else:
        values = (
            scale.value_for(well) if extinction else well.result
            for sample in samples for well in sample.retained
        )
        cv = Statistics.calculate(values).cv_percent
    return replace(result, cv_percent=float(cv) if isfinite(cv) else None)


@dataclass(frozen=True)
class QualityControlGroup:
    name: str
    checks: tuple[ValidationResult, ...]
    coordinates: tuple[str, ...] = ()


def quality_control_groups(
    context: AnalysisResult | MeasurementResult,
    *,
    blank_only: bool = False,
) -> tuple[QualityControlGroup, ...]:
    groups, labels = {}, {}
    for check in context.validation_results:
        groups.setdefault(check.sample_key, []).append(check)
        labels[check.sample_key] = check.sample_name
    for sample in select_samples(context.samples, SampleType.CONTROL):
        key = SampleType.CONTROL, sample.name
        groups.setdefault(key, [])
        labels[key] = sample.name
    groups = {key: checks for key, checks in groups.items() if (key[0] == SampleType.BLANK) == blank_only}
    readings = [(sample, well) for sample in context.samples for well in sample.measurements]
    readings.sort(
        key=lambda item: (
            (item[1].plate_index, item[1].well_x, item[1].well_y)
            if item[1].well_x is not None and item[1].well_y is not None
            else (item[1].plate_index, item[1].order, 0)
        )
    )
    positions, coordinates = {}, {}
    for position, (sample, well) in enumerate(readings):
        identity = (
            sample.concentration if sample.sample_type == SampleType.STANDARD else "Blanco" if sample.sample_type == SampleType.BLANK else sample.name
        )
        key = sample.sample_type, identity
        positions.setdefault(key, position)
        if well.well.strip():
            coordinates.setdefault(key, {})[well.well] = None
    return tuple(
        QualityControlGroup(
            labels[key],
            tuple(groups[key]),
            tuple(coordinates.get(key, ())),
        )
        for key in sorted(groups, key=lambda key: positions.get(key, len(readings)))
    )
