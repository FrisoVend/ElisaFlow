"""Load assay configurations and look them up by protocol filename."""

import json
from pathlib import Path

from elisaflow.configuration.assay_settings import AssaySettings
from elisaflow.configuration.protection import make_read_only


class AssayRegistry:
    """Protect and read the JSON library, creating one settings object per protocol."""

    def __init__(self, configs_dir: str | Path) -> None:
        configs_dir = Path(configs_dir).expanduser().resolve()
        if not configs_dir.is_dir():
            raise FileNotFoundError(f"De configuratiemap bestaat niet: {configs_dir}")

        self._settings: dict[str, AssaySettings] = {}
        self._sources: dict[str, tuple[Path, bytes]] = {}
        for path in sorted(configs_dir.glob("*.json")):
            make_read_only(path)
            contents = path.read_bytes()
            try:
                configs = json.loads(contents.decode("utf-8-sig"))
            except json.JSONDecodeError as error:
                raise ValueError(f"Ongeldige JSON in {path} op regel {error.lineno}, kolom {error.colno}.") from error
            if not isinstance(configs, dict):
                raise ValueError(f"{path} moet een object met protocolbestandsnamen als sleutels bevatten.")

            for protocol_name, config in configs.items():
                if protocol_name in self._settings:
                    raise ValueError(f"Dubbel protocol {protocol_name!r} in {path}.")
                try:
                    self._settings[protocol_name] = AssaySettings.from_dict(protocol_name, config)
                    self._sources[protocol_name] = (path, contents)
                except ValueError as error:
                    raise ValueError(f"Ongeldige instellingen voor {protocol_name!r} in {path}: {error}") from error

    def protocol_names(self) -> tuple[str, ...]:
        """Return the configured protocol filenames for an assay selection list."""
        return tuple(sorted(self._settings))

    def get(self, protocol_path: str | Path) -> AssaySettings:
        """Return settings for a protocol filename or path; fail if it is unknown."""
        protocol_name = Path(protocol_path).name
        try:
            return self._settings[protocol_name]
        except KeyError:
            raise KeyError(f"Geen assay-instellingen geregistreerd voor {protocol_name!r}.") from None

    def source_snapshot(self, protocol_path: str | Path) -> tuple[Path, bytes]:
        """Return the source file and its exact contents when settings were loaded."""
        protocol_name = self.get(protocol_path).protocol_name
        return self._sources[protocol_name]
