"""Small XML and input-file helpers shared across ELISA modules."""

from pathlib import Path
from xml.etree import ElementTree as ET


def require_file(path: str | Path, suffix: str, description: str) -> Path:
    path = Path(path).expanduser().resolve()
    if path.suffix.lower() != suffix:
        raise ValueError(f"{description} moet een {suffix}-bestand zijn: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"{description} bestaat niet: {path}")
    return path


def element_name(element: ET.Element) -> str:
    """Return an ElementTree tag without its optional namespace."""
    return element.tag.rsplit("}", 1)[-1]


def child_element(element: ET.Element | None, name: str) -> ET.Element | None:
    """Return the first matching direct ElementTree child."""
    if element is None:
        return None
    return next((child for child in element if element_name(child) == name), None)


def child_text(element: ET.Element | None, name: str) -> str | None:
    """Return the text of the first matching direct ElementTree child."""
    child = child_element(element, name)
    return child.text if child is not None else None


def xml_attribute(element: ET.Element, name: str) -> str | None:
    """Return an ElementTree attribute by namespace-neutral name."""
    return next(
        (value for key, value in element.attrib.items() if key.rsplit("}", 1)[-1] == name),
        None,
    )
