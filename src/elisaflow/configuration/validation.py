"""Validate object fields and finite numbers, independently of assay QC."""

from math import isfinite
from numbers import Real
from typing import Any


def finite_number(value: object, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"'{field_name}' moet een eindig getal zijn.")
    try:
        number = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"'{field_name}' moet een eindig getal zijn.") from error
    if not isfinite(number):
        raise ValueError(f"'{field_name}' moet een eindig getal zijn.")
    return number

def validate_config(value: object, allowed_fields: set[str], name: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{name} moet een object zijn.")
    unknown = set(value) - allowed_fields
    if unknown:
        fields = ", ".join(sorted(map(str, unknown)))
        raise ValueError(f"Onbekende velden in {name}: {fields}.")
    return value
