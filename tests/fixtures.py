"""Fresh, explicit reference data independent of the editable assay library."""

from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import fields, replace
from math import isfinite, isnan, sqrt
from xml.etree import ElementTree as ET

from elisaflow.calibration.logistic import FourPLParameters
from elisaflow.configuration.assay_settings import AssaySettings
from elisaflow.models import AssayMetaData, ResultStatus, Sample, SampleType, WellMeasurement


def _measurement_fields(values: dict, order: int) -> dict:
    """Give compact test records a physical well and an explicit 450 nm default."""
    defaults = {"well": f"{chr(65 + order % 8)}{order // 8 % 12 + 1}", "wavelength_nm": 450}
    if "well" not in values:
        defaults["plate_index"] = order // 96
    return defaults | values


def samples_from_records(records: Iterable[Mapping]) -> tuple[Sample, ...]:
    """Build analysis inputs, including optional correction fields, for tests."""
    groups = {}
    measurements = defaultdict(list)
    well_fields = {
        "source_extinction",
        "well",
        "wavelength_nm",
        "plate_name",
        "plate_index",
        "well_x",
        "well_y",
        "saturated",
        "disabled",
        "raw_extinction",
        "reference_extinction",
        "reference_invalid",
        "reference_saturated",
        "reference_disabled",
        "converted_extinction",
        "conversion_source_wavelength_nm",
    }
    for order, record in enumerate(records):
        missing = {"sample_id", "sample_name", "sample_type", "source_extinction"} - record.keys()
        if missing:
            raise ValueError(f"Ontbrekende monsterkolommen: {', '.join(sorted(missing))}.")
        sample = Sample(
            sample_id=record["sample_id"],
            name=record["sample_name"],
            sample_type=SampleType(record["sample_type"]),
            measurements=(),
            dilution=record.get("sample_dilution", 1.0),
            concentration=_optional_number(record.get("sample_concentration")),
        )
        values = {name: value for name, value in record.items() if name in well_fields}
        values["source_extinction"] = float(values["source_extinction"])
        values["order"] = order
        for name in ("well_x", "well_y", "converted_extinction", "conversion_source_wavelength_nm"):
            if name in values and _optional_number(values[name]) is None:
                values.pop(name)
        for name in ("plate_index", "plate_name"):
            if name in values and (values[name] is None or isinstance(values[name], float) and isnan(values[name])):
                values.pop(name)
        groups.setdefault(sample.key, sample)
        measurements[sample.key].append(WellMeasurement(**_measurement_fields(values, order)))
    return tuple(replace(sample, measurements=tuple(measurements[key])) for key, sample in groups.items())


def _optional_number(value) -> float | None:
    return None if value is None or isnan(float(value)) else float(value)


def processed_samples(records) -> tuple[Sample, ...]:
    """Build processed samples with explicit retained and excluded well values."""
    groups, readings = {}, defaultdict(list)
    allowed = {field.name for field in fields(WellMeasurement)}
    for index, row in enumerate(records):
        sample = Sample(
            str(row.get("sample_id", index)),
            row.get("sample_name", "Sample"),
            SampleType(row["sample_type"]),
            (),
            concentration=row.get("sample_concentration"),
            dilution=row.get("sample_dilution", 1),
        )
        values = {key: value for key, value in row.items() if key in allowed}
        values["source_extinction"] = row.get("source_extinction", row.get("extinction", float("nan")))
        values["included"] = row.get("included_in_analysis", True)
        values["order"] = index
        values["result_status"] = ResultStatus(row.get("result_status", "in_range" if isfinite(row.get("result", float("nan"))) else "unavailable"))
        readings[sample.key].append(WellMeasurement(**_measurement_fields(values, index)))
        groups[sample.key] = sample
    return tuple(replace(sample, measurements=tuple(readings[key])) for key, sample in groups.items())


def samples_with_statistics(records, *, source_samples=()) -> tuple[Sample, ...]:
    """Make concrete replicate readings for report and QC fixtures."""
    scales = {(sample.sample_id, sample.sample_type, sample.selected_wavelength_nm): sample.calculation_input_scale for sample in source_samples}
    samples = []
    order = 0
    for index, row in enumerate(records):

        def number(name, default=float("nan")):
            try:
                return float(row.get(name, default))
            except (TypeError, ValueError):
                return float("nan")

        metrics = ("extinction", "result", "curve_extinction", "back_calculated_concentration")
        count = int(number("count", 2)) or 1
        if not any(f"{metric}_sd" in row or f"{metric}_cv_percent" in row for metric in metrics):
            count = 1
        for metric in metrics:
            if isfinite(number(f"{metric}_mean")):
                spread = f"{metric}_cv_percent" if f"{metric}_cv_percent" in row else f"{metric}_sd"
                if spread in row and not isfinite(number(spread)):
                    count = 1
        values = {}
        for metric in metrics:
            mean = number(f"{metric}_mean")
            sd = abs(mean) * number(f"{metric}_cv_percent") / 100 if f"{metric}_cv_percent" in row else number(f"{metric}_sd", 0)
            offset = sd * sqrt((count - 1) / 2) if count > 1 and isfinite(sd) else 0
            values[metric] = [mean - offset, mean + offset, *([mean] * max(0, count - 2))] if count > 1 else [mean]
        wells = []
        coordinates = str(row.get("used_wells", "")).split(", ")
        wavelength = number("selected_wavelength_nm", 450)
        for replicate in range(count):
            status = row.get("result_status", "in_range" if isfinite(values["result"][replicate]) else "unavailable")
            identity = _measurement_fields({}, order)
            if replicate < len(coordinates) and coordinates[replicate]:
                identity["well"] = coordinates[replicate]
            identity["wavelength_nm"] = int(wavelength) if isfinite(wavelength) else 450
            wells.append(
                WellMeasurement(
                    source_extinction=number("source_extinction", values["extinction"][replicate]),
                    **identity,
                    order=order,
                    extinction=values["extinction"][replicate],
                    result=values["result"][replicate],
                    result_status=ResultStatus(status),
                    result_limit=number("result_limit"),
                    curve_extinction=values["curve_extinction"][replicate] if "curve_extinction_mean" in row else None,
                    back_calculated_concentration=values["back_calculated_concentration"][replicate]
                    if "back_calculated_concentration_mean" in row
                    else None,
                )
            )
            order += 1
        concentration = number("sample_concentration")
        wavelength = number("selected_wavelength_nm")
        sample_id = str(row.get("sample_id", index))
        sample_type = SampleType(row["sample_type"])
        selected_wavelength = int(wavelength) if isfinite(wavelength) else None
        samples.append(
            Sample(
                sample_id,
                str(row.get("sample_name", f"Sample {index}")),
                sample_type,
                tuple(wells),
                dilution=number("sample_dilution", 1),
                concentration=concentration if isfinite(concentration) else None,
                selected_wavelength_nm=selected_wavelength,
                selection_reason=row.get("selection_reason", ""),
                calculation_input_scale=scales.get((sample_id, sample_type, selected_wavelength)),
            )
        )
    return tuple(samples)


def add_measurement_steps(session: ET.Element, wavelengths=(405,)) -> None:
    protocol = ET.SubElement(session, "Protocol")
    steps = ET.SubElement(protocol, "ProtocolSteps")
    ET.SubElement(protocol, "ResultSteps")
    hardware = ET.SubElement(ET.SubElement(session, "Instrument"), "Filters")
    for index, wavelength in enumerate(wavelengths, 1):
        step = ET.SubElement(steps, "PhotometricFilterMeasurementStep")
        ET.SubElement(step, "Type").text = "PhotometricFilterMeasurementStep"
        ET.SubElement(step, "Name").text = f"Misleading name {index}"
        ET.SubElement(step, "OrderNumber").text = str(index)
        filter_node = ET.SubElement(ET.SubElement(step, "Filters"), "PhotometricFilter")
        ET.SubElement(filter_node, "Wavelength").text = str(wavelength)
        ET.SubElement(filter_node, "FilterPosition").text = str(index)
        device_filter = ET.SubElement(hardware, "Filter")
        ET.SubElement(device_filter, "Wavelength").text = str(wavelength)
        ET.SubElement(device_filter, "Position").text = str(index)


def assay_metadata(*, name: str = "test.skax", session_time: str = "2026-10-06T10:00:00+02:00") -> AssayMetaData:
    return AssayMetaData(
        assay_name=name,
        session_time=session_time,
        result_step_name="405nm-620nm",
        assay_description=None,
        instrument_name="reader",
        software_name="SkanIt",
        software_version="1",
        measurement_wavelengths=("405", "620"),
    )


def reference_parameters() -> FourPLParameters:
    return FourPLParameters(2.2, -1.3, 2, 0.2)


def borrelia_settings() -> AssaySettings:
    return AssaySettings.from_dict(
        "Borellia_C10_as.skax",
        {
            "lot_number": "25070204",
            "controls": {"Negatief": {}, "Kalibratie": {"replicate_filter": {"response_scale": "blank_corrected", "max_deviation_percent": 15, "max_exclusions": 1}}, "Positief": {}},
            "validation": {
                "result": {
                    "controls": {
                        "Negatief": {"extinction": {"response_scale": "blank_corrected", "max": 0.25}, "result": {"max": 0.9}},
                        "Kalibratie": {"extinction": {"response_scale": "blank_corrected", "min": 0.3}},
                        "Positief": {"extinction": {"response_scale": "blank_corrected", "min": 0.5}, "result": {"min": 1.25}},
                    }
                }
            },
            "calculation": {
                "response_scale": "blank_corrected",
                "method": "formula",
                "expression": "sample / (kalibratie * 0.43)",
                "references": {"kalibratie": "Kalibratie"},
                "control_overrides": {"Negatief": {"expression": "sample / kalibratie", "response_scale": "blank_corrected"}, "Positief": {"expression": "sample / kalibratie", "response_scale": "blank_corrected"}},
            },
            "result": {"label": "Ratio"},
        },
    )


def with_validation(settings, **targets):
    """Replace selected acceptance criteria while retaining other assay settings."""
    from dataclasses import replace

    from elisaflow.configuration.assay_settings import ValidationScope

    scope = ValidationScope.from_dict(targets)
    return replace(
        settings,
        validation=replace(
            settings.validation,
            result=replace(settings.validation.result, **{key: getattr(scope, key) for key in targets}),
        ),
    )
