"""Consistent JSON encoding for assay libraries and per-run snapshots."""

import json
from collections.abc import Mapping


def serialize_assay_configs(configs: Mapping[str, object]) -> str:
    """Preserve Unicode, reject non-finite values, and finish with one newline."""
    return json.dumps(dict(configs), indent=4, ensure_ascii=False, allow_nan=False) + "\n"
