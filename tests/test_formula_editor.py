"""Controls are inserted and bound without user-managed aliases."""

import os
import unittest
from unittest.mock import patch

from elisaflow.formatting import format_formula
from elisaflow.formulas.expressions import Formula

try:
    from PySide6.QtWidgets import QApplication

    from elisaflow.ui.editors.formula import FormulaDialog, FormulaEditor, FormulaInput
except ModuleNotFoundError as error:
    if error.name != "PySide6":
        raise
    QApplication = None


@unittest.skipIf(QApplication is None, "Install the optional PySide6 UI dependency")
class FormulaEditorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        environment = patch.dict(os.environ, {"QT_QPA_PLATFORM": "offscreen"})
        environment.start()
        cls.addClassCleanup(environment.stop)
        cls.application = QApplication.instance() or QApplication([])

    def editor(self, controls, config=None, *, sample_name="Patient"):
        editor = FormulaInput(
            (config or {}).get("expression", "sample"), tuple(controls), dict((config or {}).get("references", {})), sample_name=sample_name,
        )
        self.addCleanup(editor.deleteLater)
        return editor

    def insert(self, editor, control):
        editor.operand.setCurrentIndex(editor.operand.findData(control))
        editor.insert_button.click()

    def test_selected_control_replaces_selected_text_and_evaluates_after_reopening(self):
        editor = self.editor(("Kalibratie", "Negatief "))
        self.assertFalse(editor.insert_button.isEnabled())
        editor.expression.setText("Patient / vervang * 2")
        editor.expression.setSelection(len("Patient / "), len("vervang"))
        self.insert(editor, "Kalibratie")
        self.assertEqual(editor.expression.text(), "Patient / Kalibratie * 2")
        config = editor.to_config()
        self.assertEqual(config["references"], {"kalibratie": "Kalibratie"})
        reopened = self.editor(("Kalibratie", "Negatief "), config)
        self.assertEqual(reopened.to_config(), config)
        self.assertEqual(Formula.from_dict(config).evaluate(6, {"Kalibratie": 3}), 4)

    def test_edit_save_and_reopen_keep_the_same_names_as_the_summary_and_report(self):
        for sample_name in ("Patient", "Monster", "Negatief"):
            with self.subTest(sample_name=sample_name):
                config = {"expression": "sample / (kalibratie * 0.40)", "references": {"kalibratie": "Kalibratie"}}
                editor = FormulaEditor(config, ("Kalibratie",), sample_name=sample_name)
                self.addCleanup(editor.deleteLater)
                expected = f"{sample_name} / (Kalibratie * 0.40)"
                self.assertEqual(editor.summary.text(), expected)

                def save(dialog):
                    self.assertEqual(dialog.input.expression.text(), expected)
                    dialog.input.expression.insert(" + 1")
                    dialog._save()
                    return dialog.result()

                with patch.object(FormulaDialog, "exec", save):
                    editor._edit()
                saved = editor.to_config()
                formula = Formula.from_dict(saved)
                self.assertAlmostEqual(formula.evaluate(6, {"Kalibratie": 3}), 6)
                self.assertEqual(editor.summary.text(), expected + " + 1")
                self.assertEqual(format_formula(formula, sample_name), editor.summary.text())
                reopened = self.editor(("Kalibratie",), saved, sample_name=sample_name)
                self.assertEqual(reopened.expression.text(), editor.summary.text())

    def test_named_sample_inputs_preserve_arithmetic_and_control_references(self):
        config = {"expression": "sample / kalibratie + sample * 2", "references": {"kalibratie": "Kalibratie"}}
        for sample_name, symbol in (
            ("Patient", "Patient"), ("Monster", "Monster"), ("Kalibratie", "Kalibratie"),
            ("Negatief", "Negatief"), ("Positief", "Positief"), ("Negatief ", "Negatief"),
            ("QC 1", "'QC 1'"), ("QC-1", "'QC-1'"), ("class", "'class'"),
            ("sample", "'sample'"), ("x", "'x'"), ("kalibratie", "kalibratie"),
        ):
            with self.subTest(sample_name=sample_name):
                editor = self.editor(("Kalibratie",), config, sample_name=sample_name)
                denominator = "'Kalibratie (controlegemiddelde)'" if sample_name == "Kalibratie" else "Kalibratie"
                self.assertEqual(editor.expression.text(), f"{symbol} / {denominator} + {symbol} * 2")
                self.assertEqual(editor.expression.text(), format_formula(Formula.from_dict(config), sample_name))
                self.assertEqual(editor.expression.placeholderText(), "")
                self.assertEqual(editor.operand.itemText(0), f"{sample_name.strip()} (meetwaarde)")
                self.assertEqual(editor.operand.itemText(1), "Kalibratie (controlegemiddelde)")
                editor.expression.clear()
                self.insert(editor, None)
                editor.expression.insert(" / ")
                self.insert(editor, "Kalibratie")
                editor.expression.insert(" + ")
                self.insert(editor, None)
                editor.expression.insert(" * 2")
                saved = editor.to_config()
                self.assertEqual(saved, config)
                self.assertEqual(Formula.from_dict(saved).evaluate(6, {"Kalibratie": 3}), 14)
                reopened = self.editor(("Kalibratie",), saved, sample_name=sample_name)
                self.assertEqual(reopened.expression.text(), editor.expression.text())

    def test_repeated_insertion_reuses_binding_and_deleting_formula_removes_it(self):
        editor = self.editor(("Grens",))
        editor.expression.setText("Patient / ")
        self.insert(editor, "Grens")
        editor.expression.insert(" + ")
        self.insert(editor, "Grens")
        self.assertEqual(editor.to_config(), {"expression": "sample / grens + grens", "references": {"grens": "Grens"}})
        editor.expression.setText("Patient * 2")
        self.assertEqual(editor.to_config(), {"expression": "sample * 2"})

    def test_insertion_uses_the_cursor_position(self):
        editor = self.editor(("Positief",))
        editor.expression.setText("Patient / ( * 0.4)")
        editor.expression.setCursorPosition(len("Patient / ("))
        self.insert(editor, "Positief")
        self.assertEqual(editor.expression.text(), "Patient / (Positief * 0.4)")
        self.assertEqual(Formula.from_dict(editor.to_config()).evaluate(4, {"Positief": 2}), 5)

    def test_similar_and_reserved_control_names_keep_exact_distinct_bindings(self):
        controls = ("QC-1", "QC 1", "QC_1", "x", "sample", "Patient", "class", "123", "Negatief ", "µ controle")
        editor = self.editor(controls)
        editor.expression.clear()
        for index, control in enumerate(controls):
            if index:
                editor.expression.insert(" + ")
            self.insert(editor, control)
        config = editor.to_config()
        self.assertEqual(set(config["references"].values()), set(controls))
        self.assertEqual(len(config["references"]), len(controls))
        self.assertNotIn("x", config["references"])
        self.assertNotIn("sample", config["references"])
        self.assertNotIn("Patient", config["references"])
        self.assertEqual(Formula.from_dict(config).evaluate(0, dict.fromkeys(controls, 1)), len(controls))
        editor.expression.insert(" + ")
        self.insert(editor, None)
        self.assertEqual(Formula.from_dict(editor.to_config()).evaluate(2, dict.fromkeys(controls, 1)), len(controls) + 2)

    def test_unknown_names_and_missing_controls_fail_on_save(self):
        editor = self.editor(("Grens",))
        editor.expression.setText("Patient / onbekend")
        with self.assertRaisesRegex(ValueError, "Onbekende formulereferenties"):
            editor.to_config()
        editor = self.editor(("Grens",), {"expression": "sample / ref", "references": {"ref": "Ontbreekt"}})
        with self.assertRaisesRegex(ValueError, "Ontbreekt.*SKAX"):
            editor.to_config()
        editor.expression.setText("Patient / ")
        self.insert(editor, "Grens")
        self.assertEqual(editor.to_config()["references"], {"grens": "Grens"})

    def test_no_controls_still_allows_arithmetic(self):
        editor = self.editor(())
        self.assertTrue(editor.operand.isEnabled())
        self.assertEqual(editor.operand.count(), 1)
        self.assertFalse(editor.insert_button.isEnabled())
        editor.expression.clear()
        self.insert(editor, None)
        editor.expression.insert(" * 2")
        self.assertEqual(editor.to_config(), {"expression": "sample * 2"})

    def test_view_only_editor_cannot_insert_controls(self):
        editor = self.editor(("Grens",))
        editor.setEnabled(False)
        editor.operand.setCurrentIndex(0)
        editor.insert_button.click()
        self.assertEqual(editor.expression.text(), "Patient")
        self.assertEqual(editor.to_config(), {"expression": "sample"})

    def test_cancel_discards_equation_and_control_bindings(self):
        editor = FormulaEditor({"expression": "sample * 2"}, ("Grens",), sample_name="Patient")
        self.addCleanup(editor.deleteLater)

        def cancel(dialog):
            dialog.input.expression.setText("Patient / ")
            self.insert(dialog.input, "Grens")
            dialog.reject()
            return dialog.result()

        with patch.object(FormulaDialog, "exec", cancel):
            editor._edit()
        self.assertEqual(editor.to_config(), {"expression": "sample * 2"})
        self.assertEqual(editor.references, {})
        self.assertEqual(editor.summary.text(), "Patient * 2")

    def test_dialog_save_rejects_invalid_syntax_and_unknown_controls(self):
        dialog = FormulaDialog("Controleformule", None, ("Grens",), {}, sample_name="Patient")
        self.addCleanup(dialog.deleteLater)
        for expression in ("Patient / (", "Patient / onbekend", "x / grens"):
            dialog.input.expression.setText(expression)
            with patch("elisaflow.ui.editors.formula.QMessageBox.warning") as warning:
                dialog._save()
            warning.assert_called_once()
            self.assertNotEqual(dialog.result(), dialog.DialogCode.Accepted)
            with self.assertRaises(ValueError):
                dialog.to_config()
        dialog.input.expression.setText("Patient / ")
        self.insert(dialog.input, "Grens")
        dialog._save()
        self.assertEqual(dialog.result(), dialog.DialogCode.Accepted)
        self.assertEqual(dialog.to_config(), {"expression": "sample / grens", "references": {"grens": "Grens"}})

    def test_reset_action_is_independent_of_empty_summary_text(self):
        for empty_text in ("Patientformule", "Volg de patiëntformule"):
            with self.subTest(empty_text=empty_text):
                editor = FormulaEditor(
                    {"expression": "sample * 2"}, (), sample_name="Controle", optional=True,
                    empty_text=empty_text, clear_text="Herstellen",
                )
                self.addCleanup(editor.deleteLater)
                self.assertEqual(editor.clear_button.text(), "Herstellen")
                editor.clear_button.click()
                self.assertIsNone(editor.to_config())
                self.assertEqual(editor.summary.text(), empty_text)
                self.assertFalse(editor.clear_button.isEnabled())

    def test_saved_optional_equation_is_committed_and_can_be_removed(self):
        editor = FormulaEditor(None, ("Grens",), sample_name="Patient", optional=True)
        self.addCleanup(editor.deleteLater)
        self.assertIsNone(editor.to_config())
        self.assertEqual(editor.clear_button.text(), "Verwijderen")

        def save(dialog):
            self.assertEqual(dialog.input.expression.text(), "")
            self.assertEqual(dialog.input.expression.placeholderText(), "")
            self.insert(dialog.input, None)
            dialog.input.expression.insert(" / ")
            self.insert(dialog.input, "Grens")
            dialog._save()
            return dialog.result()

        with patch.object(FormulaDialog, "exec", save):
            editor._edit()
        self.assertEqual(editor.to_config(), {"expression": "sample / grens", "references": {"grens": "Grens"}})
        self.assertEqual(editor.references, {"grens": "Grens"})
        editor.clear_button.click()
        self.assertIsNone(editor.to_config())
        self.assertEqual(editor.summary.text(), "Geen voorbewerking")
