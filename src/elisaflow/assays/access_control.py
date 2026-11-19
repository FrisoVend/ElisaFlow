"""Password-based assay editing permissions, independent of the GUI toolkit."""

from __future__ import annotations

import hashlib
import hmac
import os
import re
import secrets
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from tempfile import NamedTemporaryFile

if os.name == "nt":
    import msvcrt
else:
    import fcntl

MIN_PASSWORD_LENGTH = 6
MAX_PASSWORD_LENGTH = 1024
_HASH_PREFIX = "scrypt$131072$8$1$"
_HASH_PATTERN = re.compile(re.escape(_HASH_PREFIX) + r"[0-9a-f]{32}\$[0-9a-f]{64}")


@contextmanager
def _password_lock(path: Path) -> Iterator[None]:
    """Serialize password changes across application instances."""
    flags = os.O_RDWR | os.O_CREAT
    if os.name == "nt":
        flags |= os.O_BINARY | os.O_TEMPORARY
    with os.fdopen(os.open(path.with_suffix(path.suffix + ".lock"), flags, 0o600), "r+b") as lock:
        try:
            if os.name == "nt":
                msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            raise OSError("Het wachtwoord wordt elders gewijzigd. Probeer het opnieuw.") from error
        try:
            yield
        finally:
            lock.seek(0)
            if os.name == "nt":
                msvcrt.locking(lock.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def _derive_key(password: str, salt: bytes) -> bytes:
    # OWASP's scrypt baseline: N=2**17, r=8, p=1 (about 128 MiB).
    return hashlib.scrypt(password.encode("utf-8"), salt=salt, n=2**17, r=8, p=1, maxmem=256 * 1024 * 1024, dklen=32)


def _matches_password(password: str, stored_hash: str) -> bool:
    if not isinstance(password, str) or not password or len(password) > MAX_PASSWORD_LENGTH:
        return False
    salt, expected_key = stored_hash.rsplit("$", 2)[1:]
    return hmac.compare_digest(_derive_key(password, bytes.fromhex(salt)), bytes.fromhex(expected_key))


class AccessControl:
    """One shared password, authorizing one editing operation at a time.

    First setup requires write access to the password directory. Subsequent
    changes always require the current password, including during editing.
    """

    def __init__(self, password_file: str | Path) -> None:
        self.password_file = Path(password_file).expanduser().resolve()
        self._authenticated_hash: str | None = None

    def has_password(self) -> bool:
        """Whether a password file exists; an invalid file still counts as set."""
        try:
            self.password_file.stat()
        except FileNotFoundError:
            return False
        return True

    def _read_hash(self) -> str | None:
        try:
            with self.password_file.open("rb") as source:
                contents = source.read(512)
        except FileNotFoundError:
            return None
        try:
            stored_hash = contents.decode("ascii").strip()
        except UnicodeDecodeError:
            stored_hash = ""
        if len(contents) == 512 or _HASH_PATTERN.fullmatch(stored_hash) is None:
            raise ValueError("Het wachtwoordbestand is ongeldig. Laat de beheerder dit herstellen.")
        return stored_hash

    def verify_password(self, entered_password: str) -> bool:
        """Check the password without granting edit permission."""
        stored_hash = self._read_hash()
        return stored_hash is not None and _matches_password(entered_password, stored_hash)

    @contextmanager
    def edit_assays(self, entered_password: str) -> Iterator[None]:
        """Authorize one editor; always revoke permission when it closes.

        Every call verifies the password again. Cancellation and exceptions
        also end permission, so subsequent editors require a new password.
        """
        self.lock()
        stored_hash = self._read_hash()
        if stored_hash is None or not _matches_password(entered_password, stored_hash):
            raise PermissionError("Het wachtwoord is onjuist.")
        self._authenticated_hash = stored_hash
        try:
            yield
        finally:
            self.lock()

    def lock(self) -> None:
        """Revoke permission for the active editor."""
        self._authenticated_hash = None

    def can_edit_assays(self) -> bool:
        """Check the active editing operation and detect changed credentials."""
        if self._authenticated_hash is None:
            return False
        try:
            stored_hash = self._read_hash()
        except (OSError, ValueError):
            stored_hash = None
        if stored_hash != self._authenticated_hash:
            self.lock()
            return False
        return True

    def require_assay_editor(self) -> None:
        """Call this from write services, even when the UI disables editing."""
        if not self.can_edit_assays():
            raise PermissionError("Je hebt geen rechten om assays te bewerken. Voer het wachtwoord in om een assay te bewerken.")

    def set_new_password(self, plain_password: str, *, current_password: str | None = None) -> None:
        """Set up or change the password atomically, then lock this instance.

        An existing password can only be replaced with the correct current
        password. Missing, unreadable or corrupt credentials never grant access.
        """
        self.password_file.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        with _password_lock(self.password_file):
            stored_hash = self._read_hash()
            if stored_hash is not None and (current_password is None or not _matches_password(current_password, stored_hash)):
                raise PermissionError("Het huidige wachtwoord is onjuist.")
            if not isinstance(plain_password, str) or not MIN_PASSWORD_LENGTH <= len(plain_password) <= MAX_PASSWORD_LENGTH:
                raise ValueError(f"Gebruik een wachtwoord van {MIN_PASSWORD_LENGTH} tot {MAX_PASSWORD_LENGTH} tekens.")
            if not plain_password.strip():
                raise ValueError("Het wachtwoord mag niet alleen uit spaties bestaan.")
            salt = secrets.token_bytes(16)
            new_hash = _HASH_PREFIX + salt.hex() + "$" + _derive_key(plain_password, salt).hex()
            temporary = None
            try:
                with NamedTemporaryFile(dir=self.password_file.parent, prefix=f".{self.password_file.name}.", suffix=".tmp", delete=False) as source:
                    temporary = Path(source.name)
                    source.write((new_hash + "\n").encode("ascii"))
                    source.flush()
                    os.fsync(source.fileno())
                os.replace(temporary, self.password_file)
            finally:
                if temporary is not None:
                    temporary.unlink(missing_ok=True)
            self.lock()
