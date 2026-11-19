import json
import os
import stat
import unittest
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from skax_fixtures import write_skax

from elisaflow.assays.access_control import MAX_PASSWORD_LENGTH, AccessControl, _password_lock
from elisaflow.assays.management import AssayManagement

PASSWORD = "a long editor passphrase"
NEW_PASSWORD = "a different editor passphrase"


class AccessControlTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.password_file = Path(temporary.name) / ".secrets" / "assay-editor.hash"
        self.access = AccessControl(self.password_file)

    def test_missing_credentials_never_grant_access_or_create_a_file(self) -> None:
        self.assertFalse(self.access.has_password())
        self.assertFalse(self.access.verify_password(PASSWORD))
        with self.assertRaises(PermissionError), self.access.edit_assays(PASSWORD):
            self.fail("Missing credentials must not authorize editing")
        self.assertFalse(self.access.can_edit_assays())
        with self.assertRaises(PermissionError):
            self.access.require_assay_editor()
        self.assertFalse(self.password_file.exists())

    def test_setup_persists_only_a_salted_hash_and_starts_locked(self) -> None:
        self.access.set_new_password(PASSWORD)
        stored = self.password_file.read_text(encoding="ascii")
        self.assertTrue(stored.startswith("scrypt$131072$8$1$"))
        self.assertNotIn(PASSWORD, stored)
        self.assertTrue(self.access.has_password())
        self.assertFalse(self.access.can_edit_assays())
        other = AccessControl(self.password_file)
        self.assertTrue(other.verify_password(PASSWORD))
        self.assertFalse(other.can_edit_assays())

    def test_permission_lasts_only_for_one_editing_operation(self) -> None:
        self.access.set_new_password(PASSWORD)
        self.assertTrue(self.access.verify_password(PASSWORD))
        self.assertFalse(self.access.can_edit_assays())
        with self.access.edit_assays(PASSWORD):
            self.assertIs(self.access.can_edit_assays(), True)
            self.access.require_assay_editor()
        self.assertFalse(self.access.can_edit_assays())

    def test_each_edit_requires_a_valid_password_even_after_a_successful_edit(self) -> None:
        self.access.set_new_password(PASSWORD)
        with self.access.edit_assays(PASSWORD):
            self.assertTrue(self.access.can_edit_assays())
        with self.assertRaises(PermissionError), self.access.edit_assays("incorrect"):
            self.fail("Previous authorization must not carry over")
        self.assertFalse(self.access.can_edit_assays())

    def test_editor_exceptions_revoke_permission(self) -> None:
        self.access.set_new_password(PASSWORD)
        with self.assertRaisesRegex(RuntimeError, "Editor failed"), self.access.edit_assays(PASSWORD):
            self.assertTrue(self.access.can_edit_assays())
            raise RuntimeError("Editor failed")
        self.assertFalse(self.access.can_edit_assays())

    def test_changes_always_require_the_current_password(self) -> None:
        self.access.set_new_password(PASSWORD)
        original = self.password_file.read_bytes()
        with self.access.edit_assays(PASSWORD):
            for current in (None, "incorrect"):
                with self.subTest(current=current), self.assertRaises(PermissionError):
                    self.access.set_new_password(NEW_PASSWORD, current_password=current)
                self.assertEqual(self.password_file.read_bytes(), original)
        self.assertTrue(self.access.verify_password(PASSWORD))

    def test_change_rejects_the_old_password_and_revokes_active_editors(self) -> None:
        self.access.set_new_password(PASSWORD)
        other = AccessControl(self.password_file)
        with self.access.edit_assays(PASSWORD), other.edit_assays(PASSWORD):
            other.set_new_password(NEW_PASSWORD, current_password=PASSWORD)
            self.assertFalse(other.can_edit_assays())
            self.assertFalse(self.access.can_edit_assays())
        self.assertFalse(self.access.verify_password(PASSWORD))
        with self.access.edit_assays(NEW_PASSWORD):
            self.assertTrue(self.access.can_edit_assays())

    def test_deleting_credentials_revokes_access(self) -> None:
        self.access.set_new_password(PASSWORD)
        with self.access.edit_assays(PASSWORD):
            self.password_file.unlink()
            self.assertFalse(self.access.can_edit_assays())
            with self.assertRaises(PermissionError):
                self.access.require_assay_editor()

    def test_corrupt_credentials_cannot_be_used_or_reset(self) -> None:
        self.password_file.parent.mkdir()
        for contents in (b"", b"not a hash", b"\xff", b"scrypt$999999999$8$1$" + b"0" * 32 + b"$" + b"0" * 64):
            with self.subTest(contents=contents):
                self.password_file.write_bytes(contents)
                self.assertTrue(self.access.has_password())
                with self.assertRaises(ValueError), self.access.edit_assays(PASSWORD):
                    self.fail("Corrupt credentials must not authorize editing")
                self.assertFalse(self.access.can_edit_assays())
                with self.assertRaises(ValueError):
                    self.access.set_new_password(NEW_PASSWORD)
                self.assertEqual(self.password_file.read_bytes(), contents)

    def test_corruption_or_read_failure_revokes_the_active_editor(self) -> None:
        self.access.set_new_password(PASSWORD)
        with self.access.edit_assays(PASSWORD):
            with patch.object(Path, "open", side_effect=PermissionError("Access denied")):
                self.assertFalse(self.access.can_edit_assays())
            self.assertFalse(self.access.can_edit_assays())
        with self.access.edit_assays(PASSWORD):
            self.password_file.write_bytes(b"invalid")
            self.assertFalse(self.access.can_edit_assays())

    def test_oversized_credentials_cannot_hide_trailing_data(self) -> None:
        self.access.set_new_password(PASSWORD)
        self.password_file.write_bytes(self.password_file.read_bytes() + b" " * 512 + b"extra")
        with self.assertRaises(ValueError), self.access.edit_assays(PASSWORD):
            self.fail("Oversized credentials must not authorize editing")
        self.assertFalse(self.access.can_edit_assays())

    def test_assay_service_checks_real_permissions_before_saving(self) -> None:
        configs_dir = self.password_file.parent.parent
        config = configs_dir / "assay.json"
        config.write_text(
            json.dumps({"test.skax": {"lot_number": "original lot", "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample"}}}),
            encoding="utf-8",
        )
        protocols_dir = configs_dir / "protocols"
        protocols_dir.mkdir()
        write_skax(protocols_dir / "test.skax", completed=False)
        assays = AssayManagement(configs_dir, self.access, protocols_dir=protocols_dir)
        settings = assays.load_assay("test.skax")
        updated = replace(settings, lot_number="new lot")
        with self.assertRaises(PermissionError):
            assays.save_assay(updated, expected_settings=settings)
        self.access.set_new_password(PASSWORD)
        with self.access.edit_assays(PASSWORD):
            assays.save_assay(updated, expected_settings=settings)
            self.assertEqual(assays.load_assay("test.skax").lot_number, "new lot")
            other = AccessControl(self.password_file)
            other.set_new_password(NEW_PASSWORD, current_password=PASSWORD)
            original = config.read_bytes()
            with self.assertRaises(PermissionError):
                assays.save_assay(settings, expected_settings=updated)
        with self.assertRaises(PermissionError):
            assays.save_assay(settings, expected_settings=updated)
        self.assertEqual(config.read_bytes(), original)

    def test_invalid_new_passwords_preserve_the_existing_hash(self) -> None:
        self.access.set_new_password(PASSWORD)
        original = self.password_file.read_bytes()
        for password in ("", "short", " " * 12, "x" * (MAX_PASSWORD_LENGTH + 1)):
            with self.subTest(length=len(password)), self.assertRaises(ValueError):
                self.access.set_new_password(password, current_password=PASSWORD)
            self.assertEqual(self.password_file.read_bytes(), original)

    def test_unicode_and_long_passwords_are_not_truncated(self) -> None:
        password = "\N{LATIN SMALL LETTER E WITH ACUTE}" * 80 + "original"
        self.access.set_new_password(password)
        self.assertTrue(self.access.verify_password(password))
        self.assertFalse(self.access.verify_password("\N{LATIN SMALL LETTER E WITH ACUTE}" * 80 + "changed"))

    def test_each_setup_uses_a_new_random_salt(self) -> None:
        self.access.set_new_password(PASSWORD)
        first = self.password_file.read_bytes()
        self.access.set_new_password(PASSWORD, current_password=PASSWORD)
        self.assertNotEqual(self.password_file.read_bytes(), first)
        self.assertTrue(self.access.verify_password(PASSWORD))

    def test_failed_write_preserves_credentials_and_cleans_up_temporary_files(self) -> None:
        self.access.set_new_password(PASSWORD)
        original = self.password_file.read_bytes()
        with patch("elisaflow.assays.access_control.os.replace", side_effect=OSError("Disk error")):
            with self.assertRaises(OSError):
                self.access.set_new_password(NEW_PASSWORD, current_password=PASSWORD)
        self.assertEqual(self.password_file.read_bytes(), original)
        self.assertEqual(list(self.password_file.parent.glob("*.tmp")), [])
        self.assertTrue(self.access.verify_password(PASSWORD))

    def test_simultaneous_password_changes_are_rejected_without_overwriting(self) -> None:
        self.access.set_new_password(PASSWORD)
        original = self.password_file.read_bytes()
        with _password_lock(self.password_file):
            with self.assertRaises(OSError):
                self.access.set_new_password(NEW_PASSWORD, current_password=PASSWORD)
        self.assertEqual(self.password_file.read_bytes(), original)
        self.access.set_new_password(NEW_PASSWORD, current_password=PASSWORD)
        self.assertTrue(self.access.verify_password(NEW_PASSWORD))

    @unittest.skipIf(os.name == "nt", "Windows permissions are inherited from the directory ACL")
    def test_posix_credentials_are_owner_only(self) -> None:
        self.access.set_new_password(PASSWORD)
        self.assertEqual(stat.S_IMODE(self.password_file.stat().st_mode), 0o600)
        self.assertEqual(stat.S_IMODE(self.password_file.parent.stat().st_mode), 0o700)


if __name__ == "__main__":
    unittest.main()
