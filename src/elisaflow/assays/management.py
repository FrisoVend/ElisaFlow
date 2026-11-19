"""Assay loading, validation and authorized persistence, without widget imports."""

import json
import os
import stat
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path
from tempfile import NamedTemporaryFile

if os.name == "nt":
    import msvcrt
else:
    import fcntl

from elisaflow.assays.access_control import AccessControl
from elisaflow.assays.protocol_samples import ProtocolSamples
from elisaflow.assays.registry import AssayRegistry
from elisaflow.configuration.assay_settings import AssaySettings
from elisaflow.configuration.protection import make_read_only
from elisaflow.configuration.serialization import serialize_assay_configs
from elisaflow.xml_utils import require_file


@contextmanager
def _assay_lock(source: Path) -> Iterator[None]:
    """Coordinate saves using an OS lock, released even if an editor crashes."""
    lock_path = source.with_suffix(source.suffix + ".lock")
    flags = os.O_RDWR | os.O_CREAT
    if os.name == "nt":
        # Windows deletes the sidecar after its last handle closes, including
        # after a crash. Do not unlink after unlocking: that could race a save.
        flags |= os.O_BINARY | os.O_TEMPORARY
    # POSIX sidecars must remain: unlinking could create two lockable inodes.
    try:
        lock = os.fdopen(os.open(lock_path, flags, 0o666), "r+b")
    except OSError as error:
        raise OSError(f"Het bestand {source.name} kon niet worden vergrendeld voor opslaan. Probeer het opnieuw.") from error
    with lock:
        try:
            if os.name == "nt":
                msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            raise OSError(f"Het bestand {source.name} kon niet worden vergrendeld voor opslaan. Probeer het opnieuw.") from error
        try:
            # Initialize only after acquisition; another writer may already
            # have locked the first byte of a newly created sidecar.
            if lock.seek(0, os.SEEK_END) == 0:
                lock.write(b"\0")
                lock.flush()
            yield
        finally:
            lock.seek(0)
            if os.name == "nt":
                msvcrt.locking(lock.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def _write_file(directory: Path, contents: bytes, *, prefix: str, suffix: str) -> Path:
    """Write and flush a unique file, cleaning it up if writing fails."""
    path = None
    try:
        with NamedTemporaryFile(dir=directory, prefix=prefix, suffix=suffix, delete=False) as file:
            path = Path(file.name)
            file.write(contents)
            file.flush()
            os.fsync(file.fileno())
        return path
    except BaseException:
        if path is not None:
            path.unlink(missing_ok=True)
        raise


def _remove_staged_file(path: Path) -> None:
    if os.name == "nt" and path.exists():
        path.chmod(stat.S_IREAD | stat.S_IWRITE)
    path.unlink(missing_ok=True)


def _publish_new_file(path: Path, contents: bytes, *, read_only: bool = False) -> None:
    """Publish a complete file atomically, refusing to replace an existing path."""
    temporary = _write_file(path.parent, contents, prefix=f".{path.stem}.", suffix=".tmp")
    try:
        if read_only:
            make_read_only(temporary)
        if os.name == "nt":
            # Windows rename fails when the destination already exists.
            os.rename(temporary, path)
        else:
            os.link(temporary, path)
    finally:
        _remove_staged_file(temporary)


class AssayManagement:
    """Manage the assay library against its SKAX protocols and editing permissions."""

    def __init__(
        self,
        configs_dir: Path,
        access: AccessControl,
        *,
        protocols_dir: Path,
    ) -> None:
        self._configs_dir = Path(configs_dir).expanduser().resolve()
        self._protocols_dir = Path(protocols_dir).expanduser().resolve()
        self._registry = AssayRegistry(self._configs_dir)
        self._access = access

    def list_assays(self) -> tuple[str, ...]:
        return self._registry.protocol_names()

    def load_assay(self, protocol_name: str) -> AssaySettings:
        # Reopening an editor after a conflict must show the latest disk settings.
        self._registry = AssayRegistry(self._configs_dir)
        return self._registry.get(protocol_name)

    def validate_assay(
        self,
        protocol_name: str,
        config: dict,
        *,
        protocol_path: Path | None = None,
    ) -> AssaySettings:
        """Validate settings against the selected or registered SKAX protocol."""
        settings = AssaySettings.from_dict(protocol_name, config)
        if protocol_path is not None:
            path = require_file(protocol_path, ".skax", "SkanIt protocol")
            if path.name != protocol_name:
                raise ValueError("De assay-instellingen moeten de bestandsnaam van het geselecteerde protocol gebruiken.")
            ProtocolSamples.from_file(path).validate_settings(settings)
        else:
            self.load_protocol_samples(protocol_name).validate_settings(settings)
        return settings

    def _protocol_file(self, protocol_name: str) -> Path:
        path = (self._protocols_dir / protocol_name).resolve()
        if Path(protocol_name).name != protocol_name or path.parent != self._protocols_dir:
            raise ValueError("Assayprotocollen moeten bestandsnamen binnen de protocolmap zijn.")
        return path

    def load_protocol_samples(self, protocol_name: str) -> ProtocolSamples:
        """Read the assay's controls and standard concentrations from SKAX."""
        settings = self._registry.get(protocol_name)
        path = self._protocol_file(settings.protocol_name)
        return ProtocolSamples.from_file(path)

    def _new_assay_paths(self, source: Path) -> tuple[Path, Path]:
        """Refresh the library and reject collisions before importing a protocol."""
        self._registry = AssayRegistry(self._configs_dir)
        if source.name.casefold() in {name.casefold() for name in self.list_assays()}:
            raise ValueError(f"Er is al een assay geregistreerd voor {source.name!r}. Open deze om de instellingen te bewerken.")
        config = self._configs_dir / source.with_suffix(".json").name
        if config.exists():
            raise FileExistsError(f"Er bestaat al een assayconfiguratie: {config.name}")
        protocol = self._protocol_file(source.name)
        if protocol.exists() and protocol.read_bytes() != source.read_bytes():
            raise FileExistsError(f"Er bestaat al een ander protocol: {protocol.name}")
        return config, protocol

    def prepare_new_assay(self, protocol_path: Path) -> ProtocolSamples:
        """Inspect an authorized user's selected protocol without writing files."""
        self._access.require_assay_editor()
        source = require_file(protocol_path, ".skax", "SkanIt protocol")
        self._new_assay_paths(source)
        samples = ProtocolSamples.from_file(source)
        if not samples.controls and not samples.standards:
            raise ValueError("Een assayprotocol moet controles of standaarden bevatten.")
        return samples

    def create_assay(self, settings: AssaySettings, *, protocol_path: Path) -> None:
        """Import a new protocol and register its settings only when Save is used."""
        self._access.require_assay_editor()
        source = require_file(protocol_path, ".skax", "SkanIt protocol")
        if source.name != settings.protocol_name:
            raise ValueError("De assay-instellingen moeten de bestandsnaam van het geselecteerde protocol gebruiken.")
        protocol_contents = source.read_bytes()
        samples = self.prepare_new_assay(source)
        settings = AssaySettings.from_dict(settings.protocol_name, settings.to_dict())
        samples.validate_settings(settings)
        config_contents = serialize_assay_configs({settings.protocol_name: settings.to_dict()}).encode("utf-8")
        config_path = self._configs_dir / source.with_suffix(".json").name
        with _assay_lock(config_path):
            config_path, protocol_path = self._new_assay_paths(source)
            if source.read_bytes() != protocol_contents:
                raise ValueError("Het geselecteerde protocol is gewijzigd. Selecteer het opnieuw voordat je opslaat.")
            protocol_path.parent.mkdir(parents=True, exist_ok=True)
            imported = False
            try:
                if not protocol_path.exists():
                    _publish_new_file(protocol_path, protocol_contents)
                    imported = True
                _publish_new_file(config_path, config_contents, read_only=True)
            except BaseException:
                if imported and protocol_path.exists() and protocol_path.read_bytes() == protocol_contents:
                    protocol_path.unlink()
                raise
            self._registry = AssayRegistry(self._configs_dir)

    def save_assay(
        self,
        settings: AssaySettings,
        *,
        expected_settings: AssaySettings,
    ) -> None:
        """Back up and atomically replace an existing assay's unchanged source file.

        Editors pass their original settings to detect stale drafts even when
        another editor in the same application has refreshed the registry.
        """
        self._access.require_assay_editor()
        if settings.protocol_name not in self._registry.protocol_names():
            raise ValueError(f"Een niet-geregistreerde assay kan niet worden opgeslagen: {settings.protocol_name!r}.")
        settings = self.validate_assay(settings.protocol_name, settings.to_dict())
        loaded = self._registry.get(settings.protocol_name)
        conflict = "De assay-instellingen zijn gewijzigd. Sluit de editor en open deze opnieuw voordat je opslaat."
        if expected_settings.protocol_name != loaded.protocol_name or expected_settings.to_dict() != loaded.to_dict():
            raise ValueError(conflict)
        source, expected_contents = self._registry.source_snapshot(settings.protocol_name)
        with _assay_lock(source):
            original = source.read_bytes()
            if original != expected_contents:
                raise ValueError(conflict)
            if settings.to_dict() == loaded.to_dict():
                return
            configs = json.loads(original.decode("utf-8-sig"))
            configs[settings.protocol_name] = settings.to_dict()
            contents = serialize_assay_configs(configs).encode("utf-8")
            temporary = _write_file(source.parent, contents, prefix=f".{source.stem}.", suffix=".tmp")
            try:
                mode = stat.S_IMODE(source.stat().st_mode)
                os.chmod(temporary, mode & ~0o222)
                history = source.parent / ".history" / source.name
                history.mkdir(parents=True, exist_ok=True)
                backup = _write_file(
                    history,
                    original,
                    prefix=f"{datetime.now(UTC):%Y%m%dT%H%M%S.%fZ}-",
                    suffix=".json",
                )
                make_read_only(backup)
                if source.read_bytes() != original:
                    raise ValueError(conflict)
                # Windows cannot replace a read-only destination. Unlock only
                # for the atomic replacement; the staged file is already protected.
                source.chmod(mode | stat.S_IWUSR)
                try:
                    os.replace(temporary, source)
                finally:
                    make_read_only(source)
            finally:
                _remove_staged_file(temporary)
            self._registry = AssayRegistry(self._configs_dir)
