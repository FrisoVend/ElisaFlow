import unittest

import numpy as np
import pandas as pd

from elisaflow.configuration.assay_settings import ReportSettings
from elisaflow.formatting import format_formula, format_number, format_result, format_text, metric_decimals, metric_label
from elisaflow.formulas.expressions import Formula
from elisaflow.models import ResultStatus, ValidationMetric


class ReportFormattingTests(unittest.TestCase):
    def test_formula_display_uses_sample_name_and_keeps_control_names_intact(self) -> None:
        formula = Formula("sample / sample_control", {"sample_control": "Kalibratie"})
        for sample_name, expected in (
            ("Patient", "Patient / Kalibratie"), ("Negatief", "Negatief / Kalibratie"),
            ("Positief", "Positief / Kalibratie"), ("QC sample_control", "'QC sample_control' / Kalibratie"),
        ):
            with self.subTest(sample_name=sample_name):
                self.assertEqual(format_formula(formula, sample_name), expected)

    def test_missing_numeric_values_use_the_selected_marker(self) -> None:
        for value in (None, pd.NA, np.nan, np.inf, -np.inf):
            with self.subTest(value=value):
                self.assertEqual(format_number(value, 2), "-")
                self.assertEqual(format_number(value, 2, missing=""), "")
                self.assertEqual(format_text(value), "-")
        self.assertEqual(format_number(1.23456, 3), "1.235")
        self.assertEqual(format_number(10), "10")
        self.assertEqual(format_text("Patient"), "Patient")

    def test_concentration_bounds_share_numeric_and_missing_value_rules(self) -> None:
        self.assertEqual(format_result(None, ResultStatus.BELOW_RANGE.value, 1.234, 2), "<1.23")
        self.assertEqual(format_result(None, ResultStatus.ABOVE_RANGE.value, 99, 1), ">99.0")
        for limit in (None, pd.NA, np.nan, np.inf):
            with self.subTest(limit=limit):
                self.assertEqual(format_result(2, ResultStatus.BELOW_RANGE.value, limit, 2), "2.00")

    def test_metric_presentation_uses_report_precision_and_concentration_units(self) -> None:
        settings = ReportSettings(extinction_decimals=2, result_decimals=3, calibration_decimals=4)
        for metric, label, decimals in (
            (ValidationMetric.EXTINCTION, "Ext", 2),
            (ValidationMetric.RESULT, "Ratio", 3),
            (ValidationMetric.RESULT, "Score [%]", 3),
        ):
            with self.subTest(metric=metric):
                self.assertEqual(metric_label(metric, "Score", "%") if label == "Score [%]" else metric_label(metric, "Ratio"), label)
                self.assertEqual(metric_decimals(metric, settings), decimals)


if __name__ == "__main__":
    unittest.main()
