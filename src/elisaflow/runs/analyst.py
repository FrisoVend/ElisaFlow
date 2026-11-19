"""Validate the analyst name entered for a measurement or reanalysis."""


def validate_analyst_name(value: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("Vul uw naam in.")
    if any(ord(character) < 32 or ord(character) == 127 for character in value):
        raise ValueError("De naam mag geen regeleinden of besturingstekens bevatten.")
    return value.strip()
