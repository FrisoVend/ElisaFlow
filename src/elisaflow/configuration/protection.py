"""Read-only file attributes to prevent accidental configuration edits."""

import stat
from pathlib import Path


def make_read_only(path: Path) -> None:
    """Remove write permission without changing other file permissions."""
    mode = stat.S_IMODE(path.stat().st_mode)
    if mode & 0o222:
        path.chmod(mode & ~0o222)
