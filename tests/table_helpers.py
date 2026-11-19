"""Independent tabular views for test fixtures and assertions."""

from collections.abc import Iterable
from dataclasses import asdict
from math import isnan

import pandas as pd

from elisaflow.models import ResponseScale, Sample


def well_table(samples: Iterable[Sample]) -> pd.DataFrame:
    rows = []
    for sample in samples:
        for well in sample.measurements:
            row = asdict(well)
            row["result_status"] = well.result_status.value
            row["included_in_analysis"] = row.pop("included")
            row.update(
                sample_id=sample.sample_id,
                sample_name=sample.name,
                sample_type=sample.sample_type.value,
                sample_dilution=sample.dilution,
                sample_concentration=sample.concentration,
                selected_wavelength_nm=sample.selected_wavelength_nm,
                selection_reason=sample.selection_reason,
            )
            rows.append(row)
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame = frame.sort_values("order", kind="stable").set_index("order")
    frame.index.name = None
    for name in ("curve_extinction", "back_calculated_concentration", "raw_extinction", "reference_extinction", "sample_concentration",
                 "converted_extinction", "conversion_source_wavelength_nm"):
        if frame[name].isna().all():
            frame = frame.drop(columns=name)
    return frame


def sample_table(samples: Iterable[Sample]) -> pd.DataFrame:
    rows = []
    for sample in samples:
        if not sample.reportable:
            continue
        row = {
            "sample_id": sample.sample_id,
            "sample_name": sample.name,
            "sample_type": sample.sample_type.value,
            "sample_dilution": sample.dilution,
            "sample_concentration": sample.concentration if sample.concentration is not None else float("nan"),
            "count": sum(not isnan(well.extinction) for well in sample.retained),
            "used_wells": ", ".join(sample.used_wells),
            "result_status": sample.result_status.value,
            "result_limit": sample.result_limit if sample.result_limit is not None else float("nan"),
            "selected_wavelength_nm": sample.selected_wavelength_nm if sample.selected_wavelength_nm is not None else float("nan"),
            "selection_reason": sample.selection_reason,
        }
        for name in ("extinction", "result", "curve_extinction", "back_calculated_concentration"):
            statistics = ResponseScale.BLANK_CORRECTED.statistics(sample) if name == "extinction" else getattr(sample, name)
            row.update({f"{name}_{metric}": value for metric, value in asdict(statistics).items()})
        row["recovery_percent"] = sample.recovery_percent
        rows.append(row)
    columns = [
        "sample_id",
        "sample_name",
        "sample_type",
        "sample_dilution",
        "sample_concentration",
        "count",
        "used_wells",
        "result_status",
        "result_limit",
        "selected_wavelength_nm",
        "selection_reason",
        *(
            f"{metric}_{statistic}"
            for metric in ("extinction", "result", "curve_extinction", "back_calculated_concentration")
            for statistic in ("mean", "sd", "cv_percent")
        ),
        "recovery_percent",
    ]
    return pd.DataFrame(rows, columns=columns)
