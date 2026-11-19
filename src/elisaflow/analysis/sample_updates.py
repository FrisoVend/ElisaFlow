"""Immutable updates shared by measurement preparation and analysis."""

from collections.abc import Callable, Iterable
from dataclasses import replace

from elisaflow.models import ResultStatus, Sample, WellMeasurement


def map_wells(samples: Iterable[Sample], update: Callable[[Sample, WellMeasurement], WellMeasurement]) -> tuple[Sample, ...]:
    """Update wells with their sample as context, preserving identities and order."""
    return tuple(replace(sample, measurements=tuple(update(sample, well) for well in sample.measurements)) for sample in samples)


def clear_result(well: WellMeasurement) -> WellMeasurement:
    return replace(well, result=float("nan"), result_status=ResultStatus.UNAVAILABLE, result_limit=float("nan"))
