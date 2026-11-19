"""Sphinx configuration for the ElisaFlow documentation."""

from __future__ import annotations

import sys
import tomllib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from elisaflow.branding import ACCENT, BRAND, ICON_PATH, LOGO_PATH  # noqa: E402

project = "Reinier Diagnostiek ElisaFlow"
copyright = "2026, Reinier Diagnostiek"
author = "Friso Venderbosch"
release = tomllib.loads((Path(__file__).resolve().parents[1] / "pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
language = "nl"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx-mathjax-offline",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "furo"
html_static_path = ["_static"]
html_css_files = ["elisaflow.css"]
html_logo = str(LOGO_PATH)
html_favicon = str(ICON_PATH)
html_title = f"ElisaFlow {release}"
html_theme_options = {
    "light_css_variables": {
        "color-brand-primary": BRAND,
        "color-brand-content": BRAND,
    },
    "dark_css_variables": {
        "color-brand-primary": ACCENT,
        "color-brand-content": ACCENT,
    },
}

autodoc_typehints = "description"
autodoc_member_order = "bysource"
