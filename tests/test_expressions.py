import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
from fixtures import assay_metadata, samples_from_records
from table_helpers import sample_table, well_table

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.assays.registry import AssayRegistry
from elisaflow.configuration.assay_settings import FormulaOverrideSettings, FormulaSettings, ValidationPlan
from elisaflow.formulas.expressions import Formula
from elisaflow.models import AssayData, ValidationStatus


class FormulaTests(unittest.TestCase):
    def test_arithmetic_preserves_scalar_and_array_inputs(self):
        formula = Formula("-(sample - reference) * +2 / 4 + 3", {"reference": "Control"})
        self.assertEqual(formula.evaluate(5, {"Control": 1}), 1)
        np.testing.assert_allclose(formula.evaluate([1, 5], {"Control": 1}), [3, 1])
        sample = np.array([1, 5, np.nan])
        original = sample.copy()
        np.testing.assert_allclose(formula.evaluate(sample, {"Control": 1}), [3, 1, np.nan])
        np.testing.assert_allclose(sample, original)
        self.assertTrue(np.isnan(Formula("100").evaluate(np.nan)))

    def test_precedence_and_parentheses(self):
        self.assertEqual(Formula("sample + 2 * 3").evaluate(1), 7)
        self.assertEqual(Formula("(sample + 2) * 3").evaluate(1), 9)
        self.assertEqual(Formula("sample / 2 / 4").evaluate(16), 2)

    def test_non_arithmetic_syntax_and_unknown_names_are_rejected(self):
        for expression in (
            "",
            " ",
            "sample +",
            "unknown",
            "sample.__class__",
            "sample[0]",
            "abs(sample)",
            "__import__('os')",
            "sample ** 2",
            "sample // 2",
            "sample % 2",
            "sample > 0",
            "[sample]",
            "sample if sample else 0",
            "True",
            "None",
            "'sample'",
            "1j",
            "1e1000",
            "lambda: sample",
            " + ".join(["sample"] * 100),
            "(" * 1000 + "sample" + ")" * 1000,
        ):
            with self.subTest(expression=expression), self.assertRaises(ValueError):
                Formula(expression)

    def test_references_and_numerical_failures(self):
        for alias in ("sample", "x", "_private", "a b", "1a", "class", 1):
            with self.subTest(alias=alias), self.assertRaises(ValueError):
                Formula("sample", {alias: "Control"})
        for value in (None, "", 1):
            with self.subTest(control=value), self.assertRaises(ValueError):
                Formula("sample", {"ref": value})
        formula = Formula("sample / (positive - negative)", {"positive": "P", "negative": "N"})
        for controls in ({}, {"P": 1, "N": 1}, {"P": np.nan, "N": 1}, {"P": np.inf, "N": 1}, {"P": True, "N": 1}):
            with self.subTest(controls=controls), self.assertRaises(ValueError):
                formula.evaluate(1, controls)
        for sample in (np.inf, [[1, 2]], "invalid"):
            with self.subTest(sample=sample), self.assertRaises(ValueError):
                Formula("sample").evaluate(sample)
        with self.assertRaisesRegex(ValueError, 'door nul delen'):
            Formula("1 / sample").evaluate([1, 0])
        with self.assertRaisesRegex(ValueError, 'oneindig'):
            Formula("sample * 1e308").evaluate(10)
        self.assertEqual(Formula("sample - ref", {"ref": "Zero"}).evaluate(2, {"Zero": 0}), 2)


class FormulaAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.settings = AssayRegistry(Path(__file__).resolve().parents[1] / "configs").get("AP50.skax")
        self.frame = pd.DataFrame(
            {
                "sample_id": ["blank", "p", "p", "n", "n", "x", "x", "low", "high"],
                "sample_name": ["Blank", "Positief", "Positief", "Negatief", "Negatief", "Patient", "Patient", "Low", "High"],
                "sample_type": ["Blank", "Control", "Control", "Control", "Control", "Unknown", "Unknown", "Unknown", "Unknown"],
                "source_extinction": [0.1, 1.5, 1.7, 0.5, 0.7, 1.0, 1.2, 0.3, 1.9],
            }
        )

    def analyze(self, settings=None):
        return AssayAnalyzer(AssayData(samples_from_records(self.frame.to_dict('records')), assay_metadata(name='AP50.skax')), settings or self.settings).analyze()

    def test_ap50_uses_control_means_and_does_not_clip_scores(self):
        original = self.frame.copy(deep=True)
        result = self.analyze()
        np.testing.assert_allclose(well_table(result.samples)["result"].iloc[1:], [90, 110, -10, 10, 40, 60, -30, 130])
        patients = sample_table(result.samples).set_index("sample_id")
        self.assertAlmostEqual(patients.loc["x", "result_mean"], 50)
        self.assertAlmostEqual(patients.loc["p", "result_mean"], 100)
        self.assertAlmostEqual(patients.loc["n", "result_mean"], 0)
        self.assertEqual(self.settings.result.to_dict(), {"label": "Score", "unit": "%"})
        self.assertFalse(self.settings.report.show_result_cv)
        pd.testing.assert_frame_equal(self.frame, original)

    def test_result_validation_and_control_overrides_share_reference_semantics(self):
        settings = replace(
            self.settings,
            calculation=replace(self.settings.calculation, control_overrides={
                "Positief": FormulaOverrideSettings("100 * (sample / negative - negative) / (positive - negative)", response_scale="blank_corrected"),
                "Negatief": FormulaOverrideSettings("100 * (sample * 10 - negative) / (positive - negative)", response_scale="blank_corrected"),
            }),
            validation=ValidationPlan.from_dict({"result": {"controls": {"Positief": {"result": {"min": 200, "max": 300}}}}}),
        )
        result = self.analyze(settings)
        # Control overrides never change references used elsewhere.
        np.testing.assert_allclose(well_table(result.samples)["result"].iloc[1:5], [230, 270, 350, 550])
        self.assertAlmostEqual(sample_table(result.samples).set_index("sample_id").loc["x", "result_mean"], 50)
        self.assertEqual(result.validation_status, ValidationStatus.PASSED)

    def test_equal_control_means_and_missing_controls_stop_analysis(self):
        self.frame.loc[self.frame.sample_id.eq("p"), "source_extinction"] = [0.5, 0.7]
        with self.assertRaisesRegex(ValueError, 'door nul delen'):
            self.analyze()
        self.frame = self.frame.loc[~self.frame.sample_id.eq("n")]
        with self.assertRaisesRegex(ValueError, 'Geen metingen'):
            self.analyze()

    def test_control_overrides_replace_the_main_formula_for_those_wells(self):
        self.frame = pd.DataFrame(
            {"sample_id": ["p", "x"], "sample_name": ["Positief", "Patient"], "sample_type": ["Control", "Unknown"], "source_extinction": [1, 2]}
        )
        settings = replace(
            self.settings,
            calculation=FormulaSettings(
                Formula("sample / (sample - 1)"),
                {"Positief": FormulaOverrideSettings("sample", response_scale="blank_corrected")},
                response_scale="blank_corrected",
            ),
        )
        np.testing.assert_allclose(well_table(self.analyze(settings).samples)["result"], [1, 2])

    def test_required_reference_means_are_checked_when_all_readings_saturated(self):
        self.frame = self.frame.loc[~self.frame.sample_type.eq("Blank")].copy()
        self.frame["saturated"] = True
        with self.assertRaisesRegex(ValueError, 'Geen metingen'):
            self.analyze()
