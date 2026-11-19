"""Visible formula names must survive editing without changing calculations."""

import unittest

from elisaflow.formatting import format_formula
from elisaflow.formulas.expressions import Formula
from elisaflow.formulas.names import FormulaNames


class FormulaNamesTests(unittest.TestCase):
    def test_actual_names_round_trip_with_the_same_notation_in_reports(self):
        for control, visible in (
            ("Kalibratie", "Kalibratie"), ("Negatief ", "'Negatief '"),
            ("QC 1", "'QC 1'"), ("QC-1", "'QC-1'"), ("class", "'class'"),
            ("123", "'123'"), ("sample", "'sample'"), ("x", "'x'"),
            ("Patient", "'Patient (controlegemiddelde)'"),
        ):
            with self.subTest(control=control):
                formula = Formula("sample / (ref * 0.40)", {"ref": control})
                names = FormulaNames("Patient", formula.references)
                display = f"Patient / ({visible} * 0.40)"
                self.assertEqual(names.format(formula.expression), display)
                self.assertEqual(format_formula(formula, "Patient"), display)
                parsed = names.parse(display)
                self.assertEqual(parsed, formula)
                self.assertAlmostEqual(parsed.evaluate(6, {control: 3}), 5)

    def test_inserted_controls_use_actual_names_and_reuse_existing_bindings(self):
        names = FormulaNames("Patient", {"cutoff_alias": "Kalibratie"})
        self.assertEqual(names.control_symbol("Kalibratie"), "Kalibratie")
        self.assertEqual(names.control_symbol("Positief"), "Positief")
        formula = names.parse("Patient / Kalibratie + Positief")
        self.assertEqual(formula.expression, "sample / cutoff_alias + positief")
        self.assertEqual(formula.evaluate(6, {"Kalibratie": 3, "Positief": 4}), 6)
        self.assertEqual(format_formula(formula, "Patient"), "Patient / Kalibratie + Positief")

    def test_exact_control_names_with_identical_trimmed_labels_stay_distinct(self):
        formula = Formula("sample - first + second", {"first": "Negatief", "second": "Negatief "})
        names = FormulaNames("Patient", formula.references)
        display = "Patient - Negatief + 'Negatief '"
        self.assertEqual(names.format(formula.expression), display)
        self.assertEqual(names.parse(display).evaluate(2, {"Negatief": 3, "Negatief ": 8}), 7)

    def test_complex_sample_name_and_its_own_control_mean_are_distinct(self):
        formula = Formula("sample / mean", {"mean": "QC-1"})
        names = FormulaNames("QC-1", formula.references)
        display = "'QC-1' / 'QC-1 (controlegemiddelde)'"
        self.assertEqual(names.format(formula.expression), display)
        self.assertEqual(names.parse(display).evaluate(6, {"QC-1": 3}), 2)

    def test_new_control_names_cannot_change_the_meaning_of_existing_symbols(self):
        for first, second in (("Patient", "Patient (controlegemiddelde)"), ("Patient (controlegemiddelde)", "Patient")):
            with self.subTest(first=first):
                names = FormulaNames("Patient", {})
                first_symbol = names.control_symbol(first)
                second_symbol = names.control_symbol(second)
                self.assertNotEqual(first_symbol, second_symbol)
                self.assertEqual(names.control_symbol(first), first_symbol)
                display = f"Patient / {first_symbol} + {second_symbol}"
                formula = names.parse(display)
                self.assertEqual(formula.evaluate(6, {first: 3, second: 7}), 9)
                self.assertEqual(format_formula(formula, "Patient"), display)

    def test_quoted_names_keep_operators_and_quotes_out_of_the_arithmetic(self):
        control = 'QC "A" / 2'
        formula = Formula("sample / ref", {"ref": control})
        names = FormulaNames("Patient", formula.references)
        display = names.format(formula.expression)
        self.assertEqual(names.parse(display), formula)
        self.assertEqual(names.parse(display).evaluate(6, {control: 3}), 2)
        self.assertEqual(names.parse('Patient / "QC \\"A\\" / 2"'), formula)

    def test_internal_aliases_and_legacy_inputs_are_not_accepted_as_display_names(self):
        names = FormulaNames("Patient", {"kalibratie": "Kalibratie"})
        for expression in ("Patient / kalibratie", "sample / Kalibratie", "x / Kalibratie", "Patient / Missing"):
            with self.subTest(expression=expression), self.assertRaisesRegex(ValueError, "Onbekende formulereferenties"):
                names.parse(expression)
        with self.assertRaisesRegex(ValueError, "Ongeldige rekenkundige formule"):
            names.parse("Patient / (")
