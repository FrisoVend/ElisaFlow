"""Application paths and access control, separate from assay settings."""

from __future__ import annotations

import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ApplicationSettings:
    protocols_dir: Path
    configs_dir: Path
    results_dir: Path
    password_file: Path

    @classmethod
    def from_file(cls, path: str | Path) -> ApplicationSettings:
        path = Path(path).expanduser().resolve()
        with path.open("rb") as source:
            config = tomllib.load(source)
        paths = config.get("paths", {})
        access_control = config.get("access_control", {})
        if not all(isinstance(section, dict) for section in (paths, access_control)):
            raise ValueError("De toepassingsinstellingen 'paths' en 'access_control' moeten TOML-tabellen zijn.")
        resolved = {}
        for name, default in {"protocols": "protocols", "configs": "configs", "results": "resultaten"}.items():
            value = paths.get(name, default)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"Het toepassingspad '{name}' moet een niet-lege tekst zijn.")
            resolved[f"{name}_dir"] = (path.parent / Path(value).expanduser()).resolve()
        password_file = access_control.get("password_file", ".secrets/assay-editor.hash")
        if not isinstance(password_file, str) or not password_file.strip():
            raise ValueError("'access_control.password_file' moet een niet-lege tekst zijn.")
        return cls(
            **resolved,
            password_file=(path.parent / Path(password_file).expanduser()).resolve(),
        )


def default_config_path() -> Path:
    """Use app.toml beside the executable or in the source working directory."""
    directory = Path(sys.executable).resolve().parent if getattr(sys, "frozen", False) else Path.cwd()
    return directory / "app.toml"
