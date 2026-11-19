import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch

from elisaflow.assays.access_control import AccessControl
from elisaflow.assays.management import AssayManagement

try:
    from PySide6.QtWidgets import QApplication, QDialog, QLineEdit

    from elisaflow.ui.assay_settings import AssaySettingsWidget
except ModuleNotFoundError as error:
    if error.name != "PySide6":
        raise
    QApplication = None

PASSWORD = "a long editor passphrase"
NEW_PASSWORD = "a different editor passphrase"


@unittest.skipIf(QApplication is None, "Install the optional PySide6 UI dependency")
class AccessControlWidgetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        environment = patch.dict(os.environ, {"QT_QPA_PLATFORM": "offscreen"})
        environment.start()
        cls.addClassCleanup(environment.stop)
        cls.application = QApplication.instance() or QApplication([])

    def setUp(self) -> None:
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.access = AccessControl(Path(temporary.name) / "editor.hash")
        self.assays = Mock(spec=AssayManagement)
        self.assays.list_assays.return_value = ("test.skax",)
        self.widget = AssaySettingsWidget(self.assays, self.access)
        self.addCleanup(self.widget.deleteLater)
        self.widget.assay_list.setCurrentRow(0)
        warning = patch("elisaflow.ui.assay_settings.QMessageBox.warning")
        self.warning = warning.start()
        self.addCleanup(warning.stop)
        information = patch("elisaflow.ui.assay_settings.QMessageBox.information")
        self.information = information.start()
        self.addCleanup(information.stop)

    def test_first_setup_requires_confirmation_and_keeps_editing_locked(self) -> None:
        password = "edit42"
        with patch("elisaflow.ui.assay_settings.QInputDialog.getText", side_effect=[(password, True), (password, True)]) as prompt:
            self.widget.password_button.click()
        self.assertTrue(self.access.verify_password(password))
        self.assertFalse(self.access.can_edit_assays())
        self.assertTrue(self.widget.view_button.isEnabled())
        self.assertTrue(all(call.args[-1] == QLineEdit.EchoMode.Password for call in prompt.call_args_list))
        self.warning.assert_not_called()

    def test_search_hides_nonmatching_assays_and_clears_an_invisible_selection(self) -> None:
        self.widget._populate_assay_list(("Alpha.skax", "Beta.skax"))
        self.widget._filter_assays()
        self.widget.assay_list.setCurrentRow(0)
        self.widget.search.setText(" BETA ")
        self.assertTrue(self.widget.assay_list.item(0).isHidden())
        self.assertFalse(self.widget.assay_list.item(1).isHidden())
        self.assertIsNone(self.widget.assay_list.currentItem())
        self.assertFalse(self.widget.view_button.isEnabled())
        self.assertFalse(self.widget.edit_button.isEnabled())
        self.widget.assay_list.setCurrentRow(1)
        with patch("elisaflow.ui.assay_settings.AssayEditor") as editor:
            self.widget.view_button.click()
        self.assays.load_assay.assert_called_once_with("Beta.skax")
        self.assertTrue(editor.call_args.kwargs["read_only"])

    def test_search_empty_state_recovers_when_the_query_is_cleared(self) -> None:
        self.widget.search.setText("missing assay")
        self.assertTrue(self.widget.assay_list.isHidden())
        self.assertFalse(self.widget.empty_state.isHidden())
        self.assertTrue(self.widget.new_button.isEnabled())
        self.widget.search.clear()
        self.assertFalse(self.widget.assay_list.isHidden())
        self.assertTrue(self.widget.empty_state.isHidden())
        self.assertFalse(self.widget.assay_list.item(0).isHidden())

    def test_each_edit_prompts_again_and_permission_ends_when_the_editor_closes(self) -> None:
        self.access.set_new_password(PASSWORD)

        def edit():
            self.assertTrue(self.access.can_edit_assays())
            self.access.require_assay_editor()
            return QDialog.DialogCode.Accepted

        with (
            patch("elisaflow.ui.assay_settings.QInputDialog.getText", return_value=(PASSWORD, True)) as prompt,
            patch("elisaflow.ui.assay_settings.AssayEditor") as editor,
        ):
            editor.return_value.exec.side_effect = edit
            for count in (1, 2):
                self.widget.edit_button.click()
                self.assertEqual(prompt.call_count, count)
                self.assertEqual(editor.call_count, count)
                self.assertFalse(self.access.can_edit_assays())
        self.assertEqual(prompt.call_args.args[-1], QLineEdit.EchoMode.Password)
        self.warning.assert_not_called()

    def test_incorrect_or_cancelled_password_prevents_edit_and_create(self) -> None:
        self.access.set_new_password(PASSWORD)
        with (
            patch("elisaflow.ui.assay_settings.AssayEditor") as editor,
            patch("elisaflow.ui.assay_settings.QFileDialog.getOpenFileName") as picker,
        ):
            for button in (self.widget.edit_button, self.widget.new_button):
                for reply in (("incorrect", True), (PASSWORD, False)):
                    with self.subTest(action=button.text(), reply=reply), patch("elisaflow.ui.assay_settings.QInputDialog.getText", return_value=reply):
                        button.click()
                        self.assertFalse(self.access.can_edit_assays())
            editor.assert_not_called()
            picker.assert_not_called()
        self.assays.load_assay.assert_not_called()
        self.assertEqual(self.warning.call_count, 2)

    def test_new_assay_always_prompts_and_picker_cancellation_revokes_permission(self) -> None:
        self.access.set_new_password(PASSWORD)

        def cancel_picker(*args):
            self.assertTrue(self.access.can_edit_assays())
            return "", ""

        with (
            patch("elisaflow.ui.assay_settings.QInputDialog.getText", return_value=(PASSWORD, True)) as prompt,
            patch("elisaflow.ui.assay_settings.QFileDialog.getOpenFileName", side_effect=cancel_picker) as picker,
            patch("elisaflow.ui.assay_settings.AssayEditor") as editor,
        ):
            for count in (1, 2):
                self.widget.new_button.click()
                self.assertEqual(prompt.call_count, count)
                self.assertEqual(picker.call_count, count)
                self.assertFalse(self.access.can_edit_assays())
            editor.assert_not_called()

    def test_cancelling_a_new_editor_revokes_permission(self) -> None:
        self.access.set_new_password(PASSWORD)

        def cancel_editor():
            self.assertTrue(self.access.can_edit_assays())
            return QDialog.DialogCode.Rejected

        with (
            patch("elisaflow.ui.assay_settings.QInputDialog.getText", return_value=(PASSWORD, True)),
            patch("elisaflow.ui.assay_settings.QFileDialog.getOpenFileName", return_value=("new.skax", "")),
            patch("elisaflow.ui.assay_settings.AssayEditor") as editor,
        ):
            editor.return_value.exec.side_effect = cancel_editor
            self.widget.new_button.click()
        editor.assert_called_once()
        self.assertFalse(self.access.can_edit_assays())

    def test_editor_errors_revoke_permission(self) -> None:
        self.access.set_new_password(PASSWORD)
        with (
            patch("elisaflow.ui.assay_settings.QInputDialog.getText", return_value=(PASSWORD, True)),
            patch("elisaflow.ui.assay_settings.AssayEditor") as editor,
        ):
            editor.return_value.exec.side_effect = ValueError("Editor error")
            self.widget.edit_button.click()
        self.warning.assert_called_once()
        self.assertFalse(self.access.can_edit_assays())

    def test_view_and_double_click_never_prompt_for_a_password(self) -> None:
        with (
            patch("elisaflow.ui.assay_settings.QInputDialog.getText") as prompt,
            patch("elisaflow.ui.assay_settings.AssayEditor") as editor,
        ):
            self.widget.view_button.click()
            self.widget.assay_list.itemDoubleClicked.emit(self.widget.assay_list.currentItem())
        prompt.assert_not_called()
        self.assertEqual(editor.call_count, 2)
        self.assertTrue(all(call.kwargs["read_only"] for call in editor.call_args_list))
        self.assertFalse(self.access.can_edit_assays())

    def test_missing_credentials_guide_the_user_to_password_setup(self) -> None:
        with patch("elisaflow.ui.assay_settings.QInputDialog.getText") as prompt:
            self.widget.edit_button.click()
            self.widget.new_button.click()
        prompt.assert_not_called()
        self.assertEqual(self.information.call_count, 2)
        self.assertFalse(self.access.has_password())

    def test_mismatched_confirmation_or_cancellation_does_not_save(self) -> None:
        for replies in (
            [(PASSWORD, False)],
            [(PASSWORD, True), (PASSWORD, False)],
            [(PASSWORD, True), (NEW_PASSWORD, True)],
        ):
            with self.subTest(replies=replies), patch("elisaflow.ui.assay_settings.QInputDialog.getText", side_effect=replies):
                self.widget.password_button.click()
                self.assertFalse(self.access.has_password())
        self.warning.assert_called_once()

    def test_change_requires_current_password_and_locks_editing(self) -> None:
        self.access.set_new_password(PASSWORD)
        with patch("elisaflow.ui.assay_settings.QInputDialog.getText", side_effect=[(PASSWORD, True), (NEW_PASSWORD, True), (NEW_PASSWORD, True)]):
            self.widget.password_button.click()
        self.assertTrue(self.access.verify_password(NEW_PASSWORD))
        self.assertFalse(self.access.verify_password(PASSWORD))
        self.assertFalse(self.access.can_edit_assays())
        self.warning.assert_not_called()

    def test_wrong_current_password_stops_the_change(self) -> None:
        self.access.set_new_password(PASSWORD)
        original = self.access.password_file.read_bytes()
        with patch("elisaflow.ui.assay_settings.QInputDialog.getText", return_value=("incorrect", True)) as prompt:
            self.widget.password_button.click()
        prompt.assert_called_once()
        self.warning.assert_called_once()
        self.assertEqual(self.access.password_file.read_bytes(), original)

    def test_corrupt_credentials_show_an_error_and_keep_editing_locked(self) -> None:
        self.access.password_file.write_bytes(b"invalid")
        with patch("elisaflow.ui.assay_settings.QInputDialog.getText", return_value=(PASSWORD, True)):
            self.widget.edit_button.click()
        self.warning.assert_called_once()
        self.assertFalse(self.access.can_edit_assays())


if __name__ == "__main__":
    unittest.main()
