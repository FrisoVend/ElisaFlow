"""Shared sample selection."""

from collections.abc import Iterable

from elisaflow.models import Sample, SampleType


def select_samples(
    samples: Iterable[Sample],
    sample_type: SampleType,
    *,
    name: str | None = None,
    concentration: float | None = None,
) -> tuple[Sample, ...]:
    return tuple(
        sample
        for sample in samples
        if sample.reportable
        and sample.sample_type == sample_type
        and (name is None or sample.name == name)
        and (concentration is None or sample.concentration == concentration)
    )
