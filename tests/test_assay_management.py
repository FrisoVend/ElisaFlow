import json
import os
import stat
import subprocess
import sys
import unittest
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from skax_fixtures import basic_records, write_skax

from elisaflow.assays.access_control import AccessControl
from elisaflow.assays.management import AssayManagement, _assay_lock
from elisaflow.assays.registry import AssayRegistry


class TestAccess(AccessControl):
    def __init__(self, editable: bool = True) -> None:
        self.editable = editable

    def can_edit_assays(self) -> bool:
        return self.editable


class AssayManagementTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name) / "configs"
        self.directory.mkdir()
        self.protocols = Path(temporary.name) / "protocols"
        self.protocols.mkdir()
        self.source = self.directory / "library.json"
        self.configs = {
            "first.skax": {
                "lot_number": "old lot",
                "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / grens", "references": {"grens": "Grens"}},
                "result": {"label": "Ratio"},
            },
            "second.skax": {
                "lot_number": "other lot",
                "controls": {"Positiëf": {}},
                "validation": {"result": {"controls": {"Positiëf": {"result": {"min": 1.2}}}}},
                "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / grens", "references": {"grens": "Grens"}},
                "result": {"label": "Ratio"},
            },
        }
        self.source.write_text(json.dumps(self.configs), encoding="utf-8-sig")
        self.original = self.source.read_bytes()
        records = basic_records()
        records[1]["sample_name"] = "Grens"
        records.append(records[1] | {"well": "B1", "sample_id": "positive", "sample_name": "Positiëf"})
        for name in self.configs:
            write_skax(self.protocols / name, records, completed=False)
        self.access = TestAccess()
        self.assays = AssayManagement(self.directory, self.access, protocols_dir=self.protocols)
        self.settings = self.assays.load_assay("first.skax")
        self.updated = replace(self.settings, lot_number="new lot")

    def test_loaded_library_is_read_only(self) -> None:
        self.assertFalse(self.source.stat().st_mode & stat.S_IWUSR)
        if os.name == "nt":
            with self.assertRaises(PermissionError):
                self.source.write_text("accidental edit", encoding="utf-8")
        self.assertEqual(self.source.read_bytes(), self.original)

    def test_save_preserves_other_assays_and_backs_up_exact_source(self) -> None:
        self.assays.save_assay(self.updated, expected_settings=self.settings)
        saved = json.loads(self.source.read_text(encoding="utf-8"))
        self.assertEqual(saved["first.skax"], self.updated.to_dict())
        self.assertEqual(saved["second.skax"], self.configs["second.skax"])
        self.assertEqual(self.assays._registry.get("first.skax"), self.updated)
        self.assertEqual(AssayRegistry(self.directory).get("first.skax"), self.updated)
        backups = list((self.directory / ".history" / self.source.name).glob("*.json"))
        self.assertEqual(len(backups), 1)
        self.assertEqual(backups[0].read_bytes(), self.original)
        self.assertFalse(self.source.stat().st_mode & stat.S_IWUSR)
        self.assertFalse(backups[0].stat().st_mode & stat.S_IWUSR)

    def test_denied_save_does_not_write_anything(self) -> None:
        self.access.editable = False
        with self.assertRaises(PermissionError):
            self.assays.save_assay(self.updated, expected_settings=self.settings)
        self.assertEqual(self.source.read_bytes(), self.original)
        self.assertEqual(list(self.directory.iterdir()), [self.source])
        self.assertFalse(self.source.stat().st_mode & stat.S_IWUSR)

    def test_external_edit_is_preserved_and_reopening_refreshes_settings(self) -> None:
        self.configs["second.skax"]["lot_number"] = "external lot"
        self.source.chmod(self.source.stat().st_mode | stat.S_IWUSR)
        self.source.write_text(json.dumps(self.configs), encoding="utf-8")
        external = self.source.read_bytes()
        with self.assertRaisesRegex(ValueError, 'gewijzigd'):
            self.assays.save_assay(self.updated, expected_settings=self.settings)
        self.assertEqual(self.source.read_bytes(), external)
        self.assertFalse((self.directory / ".history").exists())
        refreshed = self.assays.load_assay("first.skax")
        self.assays.save_assay(self.updated, expected_settings=refreshed)
        self.assertEqual(
            json.loads(self.source.read_text())["second.skax"]["lot_number"],
            "external lot",
        )

    def test_stale_editor_in_same_application_cannot_overwrite_saved_settings(self) -> None:
        second_editor = self.assays.load_assay("first.skax")
        self.assays.save_assay(self.updated, expected_settings=self.settings)
        stale_update = replace(second_editor, lot_number="stale lot")
        with self.assertRaisesRegex(ValueError, 'gewijzigd'):
            self.assays.save_assay(stale_update, expected_settings=second_editor)
        self.assertEqual(self.assays.load_assay("first.skax"), self.updated)

    def test_failed_atomic_replace_preserves_source_and_allows_retry(self) -> None:
        def fail_replace(staged: Path, source: Path) -> None:
            self.assertTrue(source.stat().st_mode & stat.S_IWUSR)
            self.assertFalse(staged.stat().st_mode & stat.S_IWUSR)
            raise OSError("disk error")

        with patch("elisaflow.assays.management.os.replace", side_effect=fail_replace):
            with self.assertRaisesRegex(OSError, "disk error"):
                self.assays.save_assay(self.updated, expected_settings=self.settings)
        self.assertEqual(self.source.read_bytes(), self.original)
        self.assertEqual(list(self.directory.glob("*.tmp")), [])
        self.assertFalse(self.source.stat().st_mode & stat.S_IWUSR)
        self.assays.save_assay(self.updated, expected_settings=self.settings)
        self.assertEqual(self.assays.load_assay("first.skax"), self.updated)

    def test_unchanged_save_preserves_original_format_without_backup(self) -> None:
        self.assays.save_assay(self.settings, expected_settings=self.settings)
        self.assertEqual(self.source.read_bytes(), self.original)
        self.assertFalse((self.directory / ".history").exists())
        self.assertFalse(self.source.stat().st_mode & stat.S_IWUSR)

    def test_failed_backup_keeps_source_read_only_and_removes_staged_file(self) -> None:
        with patch("elisaflow.assays.management.Path.mkdir", side_effect=OSError("backup error")):
            with self.assertRaisesRegex(OSError, "backup error"):
                self.assays.save_assay(self.updated, expected_settings=self.settings)
        self.assertEqual(self.source.read_bytes(), self.original)
        self.assertFalse(self.source.stat().st_mode & stat.S_IWUSR)
        self.assertEqual(list(self.directory.glob("*.tmp")), [])

    def test_unknown_assay_is_not_created(self) -> None:
        with self.assertRaisesRegex(ValueError, 'niet-geregistreerde'):
            self.assays.save_assay(replace(self.updated, protocol_name="unknown.skax"), expected_settings=self.settings)
        self.assertEqual(list(self.directory.iterdir()), [self.source])

    def test_simultaneous_save_is_rejected_until_lock_is_released(self) -> None:
        with _assay_lock(self.source):
            with self.assertRaises(OSError):
                self.assays.save_assay(self.updated, expected_settings=self.settings)
            self.assertEqual(self.source.read_bytes(), self.original)
        self.assays.save_assay(self.updated, expected_settings=self.settings)
        self.assertEqual(self.assays.load_assay("first.skax"), self.updated)

    @unittest.skipUnless(os.name == "nt", "Windows deletes sidecars on close")
    def test_completed_and_unchanged_saves_remove_the_lock_sidecar(self) -> None:
        lock_path = self.source.with_suffix(".json.lock")
        # Reuse and clean up a sidecar left by an older application version.
        lock_path.write_bytes(b"\0")
        self.assays.save_assay(self.updated, expected_settings=self.settings)
        self.assertFalse(lock_path.exists())
        self.assays.save_assay(self.updated, expected_settings=self.updated)
        self.assertFalse(lock_path.exists())

    @unittest.skipUnless(os.name == "nt", "Windows deletes sidecars on close")
    def test_failed_save_removes_the_lock_sidecar_and_allows_retry(self) -> None:
        lock_path = self.source.with_suffix(".json.lock")
        with patch("elisaflow.assays.management.os.replace", side_effect=OSError("disk error")):
            with self.assertRaisesRegex(OSError, "disk error"):
                self.assays.save_assay(self.updated, expected_settings=self.settings)
        self.assertFalse(lock_path.exists())
        self.assertEqual(self.source.read_bytes(), self.original)
        self.assays.save_assay(self.updated, expected_settings=self.settings)
        self.assertFalse(lock_path.exists())

    @unittest.skipUnless(os.name == "nt", "Windows deletes sidecars on close")
    def test_contending_process_does_not_remove_an_active_lock(self) -> None:
        script = """
import sys
from pathlib import Path
from elisaflow.assays.management import _assay_lock

try:
    with _assay_lock(Path(sys.argv[1])):
        raise AssertionError("Acquired another process's lock")
except OSError:
    pass
"""
        lock_path = self.source.with_suffix(".json.lock")
        with _assay_lock(self.source):
            for _ in range(2):
                result = subprocess.run([sys.executable, "-c", script, str(self.source)], capture_output=True, text=True, timeout=20)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertTrue(lock_path.exists())
            with self.assertRaises(OSError):
                self.assays.save_assay(self.updated, expected_settings=self.settings)
            self.assertEqual(self.source.read_bytes(), self.original)
        self.assertFalse(lock_path.exists())
        self.assays.save_assay(self.updated, expected_settings=self.settings)
        self.assertFalse(lock_path.exists())

    @unittest.skipUnless(os.name == "nt", "Windows deletes sidecars on close")
    def test_process_exit_without_finally_still_removes_the_lock(self) -> None:
        script = """
import os
import sys
from pathlib import Path
from elisaflow.assays.management import _assay_lock

source = Path(sys.argv[1])
with _assay_lock(source):
    assert source.with_suffix('.json.lock').exists()
    os._exit(0)
"""
        result = subprocess.run([sys.executable, "-c", script, str(self.source)], capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(self.source.with_suffix(".json.lock").exists())
        self.assays.save_assay(self.updated, expected_settings=self.settings)

    def test_edit_during_staging_is_detected_before_replacement(self) -> None:
        self.configs["first.skax"]["lot_number"] = "external lot"
        external = json.dumps(self.configs).encode("utf-8")
        fsync = os.fsync
        edited = False

        def edit_during_flush(descriptor: int) -> None:
            nonlocal edited
            if not edited:
                self.source.chmod(self.source.stat().st_mode | stat.S_IWUSR)
                self.source.write_bytes(external)
                edited = True
            fsync(descriptor)

        with patch("elisaflow.assays.management.os.fsync", side_effect=edit_during_flush):
            with self.assertRaisesRegex(ValueError, 'gewijzigd'):
                self.assays.save_assay(self.updated, expected_settings=self.settings)
        self.assertEqual(self.source.read_bytes(), external)
        self.assertEqual(list(self.directory.glob("*.tmp")), [])


if __name__ == "__main__":
    unittest.main()
