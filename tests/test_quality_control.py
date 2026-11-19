import unittest

import numpy as np
import pandas as pd
from fixtures import processed_samples, samples_with_statistics

from elisaflow.configuration.assay_settings import AssaySettings, ValidationCriterion, ValidationLimits
from elisaflow.models import ResultStatus, SampleType, ValidationMetric
from elisaflow.quality_control.summary import summarize_quality_control
from elisaflow.quality_control.validation import AssayValidator


class QualityControlSummaryTests(unittest.TestCase):
    def test_acceptance_averages_sample_means_and_cv_pools_retained_wells(self) -> None:
        samples = pd.DataFrame(
            {
                "sample_id": ["a", "b", "patient"],
                "sample_name": ["QC"] * 3,
                "sample_type": ["Control", "Control", "Unknown"],
                "extinction_mean": [0.5, 1.5, 100],
                "count": [3, 1, 1],
            }
        )
        wells = pd.DataFrame(
            {
                "sample_id": ["a", "a", "a", "b", "b", "patient"],
                "sample_name": ["QC"] * 6,
                "sample_type": ["Control"] * 5 + ["Unknown"],
                "extinction": [0.4, 0.5, 0.6, 1.5, 9, 100],
                "included_in_analysis": [True, True, True, True, False, True],
            }
        )
        criterion = ValidationCriterion(SampleType.CONTROL, "QC", ValidationMetric.EXTINCTION, ValidationLimits(response_scale="blank_corrected"))
        summary = summarize_quality_control(processed_samples(wells.to_dict("records")), criterion)
        self.assertEqual(summary.mean, 1)
        self.assertAlmostEqual(summary.cv_percent, 100 * np.std([0.4, 0.5, 0.6, 1.5], ddof=1) / 0.75)
        settings = AssaySettings.from_dict(
            "test.skax",
            {
                "lot_number": "lot",
                "controls": {"QC": {}},
                "validation": {"result": {"controls": {"QC": {"extinction": {"response_scale": "blank_corrected", "min": 1, "max": 1}}}}},
                "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / qc", "references": {"qc": "QC"}},
                "result": {"label": "Ratio"},
            },
        )
        self.assertEqual(AssayValidator(settings.validation.result).validate(samples_with_statistics(samples.to_dict("records")))[0].value, summary.mean)

    def test_unavailable_mean_does_not_pool_away_an_incomplete_sample(self) -> None:
        samples = pd.DataFrame(
            {
                "sample_id": ["a", "b"],
                "sample_name": ["QC"] * 2,
                "sample_type": ["Control"] * 2,
                "extinction_mean": [1, pd.NA],
            }
        )
        criterion = ValidationCriterion(SampleType.CONTROL, "QC", ValidationMetric.EXTINCTION, ValidationLimits(response_scale="blank_corrected"))
        summary = summarize_quality_control(samples_with_statistics(samples.to_dict("records")), criterion)
        self.assertIsNone(summary.mean)
        self.assertIsNone(summary.cv_percent)
        self.assertIn('Elk monster', summary.message)

    def test_blank_mean_uses_its_criterion_and_keeps_zero_without_sample_records(self):
        criterion = ValidationCriterion(SampleType.BLANK, "Blanco", ValidationMetric.EXTINCTION, ValidationLimits(response_scale="raw"))
        summary = summarize_quality_control((), criterion, blank_value=0)
        self.assertEqual(summary.mean, 0)
        self.assertIsNone(summary.cv_percent)
        self.assertIsNone(summary.result_qualification)

    def test_result_criterion_cannot_use_stored_blank_as_its_mean(self):
        criterion = ValidationCriterion(SampleType.CONTROL, "QC", ValidationMetric.RESULT, ValidationLimits())
        samples = samples_with_statistics([{"sample_name": "QC", "sample_type": "Control", "result_mean": 12}])
        summary = summarize_quality_control(samples, criterion, blank_value=99)
        self.assertEqual(summary.mean, 12)
        self.assertEqual(summary.result_qualification.status, ResultStatus.IN_RANGE)

    def test_standard_criterion_selects_concentration_and_ignores_other_sample_types(self):
        criterion = ValidationCriterion(SampleType.STANDARD, "Standaard 2", ValidationMetric.RESULT, ValidationLimits(), 2)
        samples = samples_with_statistics([
            {"sample_type": "Standard", "sample_concentration": 2, "result_mean": 3},
            {"sample_type": "Standard", "sample_concentration": 5, "result_mean": 99},
            {"sample_type": "Control", "sample_concentration": 2, "result_mean": 99},
        ])
        self.assertEqual(summarize_quality_control(samples, criterion).mean, 3)


if __name__ == "__main__":
    unittest.main()
