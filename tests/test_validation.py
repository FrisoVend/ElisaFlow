import unittest

import numpy as np
import pandas as pd
from fixtures import processed_samples, samples_with_statistics

from elisaflow.configuration.assay_settings import ValidationScope
from elisaflow.models import ResultStatus, SampleType, ValidationMetric, ValidationStatus
from elisaflow.quality_control.validation import AssayValidator


class ValidationTests(unittest.TestCase):
    def test_qc_cvs_use_each_metric_and_exclude_rejected_replicates(self) -> None:
        samples = pd.DataFrame({
            "sample_id": ["a", "b"], "sample_name": ["QC", "QC"], "sample_type": ["Control", "Control"],
            "extinction_mean": [1.5, 4.0], "result_mean": [2.5, 17.0],
        })
        wells = pd.DataFrame({
            "sample_id": ["a", "a", "a", "b", "b"], "sample_name": ["QC"] * 5, "sample_type": ["Control"] * 5,
            "extinction": [1, 2, 9, 3, 5], "result": [1, 4, 81, 9, 25], "included_in_analysis": [True, True, False, True, True],
        })
        scope = self.settings(controls={"QC": {"extinction": {"response_scale": "blank_corrected", "min": 0}, "result": {"min": 0}}})
        checks = AssayValidator(scope).validate(processed_samples(wells.to_dict("records")))
        for check, values in zip(checks, ([1, 2, 3, 5], [1, 4, 9, 25]), strict=True):
            expected = pd.Series(values, dtype=float)
            self.assertAlmostEqual(check.cv_percent, 100 * expected.std() / expected.mean())
            self.assertEqual(check.value, samples[f"{check.metric.value}_mean"].mean())

    @staticmethod
    def settings(**samples: object) -> ValidationScope:
        return ValidationScope.from_dict(samples)

    def test_inclusive_limits_and_one_sided_limits(self) -> None:
        for limits, value, expected in (
            ({"min": 1}, 1, ValidationStatus.PASSED),
            ({"min": 1}, 0.9, ValidationStatus.FAILED),
            ({"max": 2}, 2, ValidationStatus.PASSED),
            ({"max": 2}, 2.1, ValidationStatus.FAILED),
            ({"min": 1, "max": 2}, 1, ValidationStatus.PASSED),
            ({"min": 1, "max": 2}, 2, ValidationStatus.PASSED),
            ({"min": 1, "max": 2}, 1.5, ValidationStatus.PASSED),
            ({"min": 1, "max": 2}, 0.9, ValidationStatus.FAILED),
            ({"min": 1, "max": 2}, 2.1, ValidationStatus.FAILED),
        ):
            with self.subTest(limits=limits, value=value):
                settings = self.settings(controls={"Positive": {"result": limits}})
                samples = pd.DataFrame({"sample_name": ["Positive"], "sample_type": ["Control"], "result_mean": [value]})
                (result,) = AssayValidator(settings).validate(samples_with_statistics(samples.to_dict("records")))
                self.assertEqual(result.status, expected)
                self.assertEqual(result.value, value)
                self.assertEqual(result.minimum, limits.get("min"))
                self.assertEqual(result.maximum, limits.get("max"))
                self.assertTrue(result.message)

    def test_control_concentration_bounds_are_preserved_without_passing_validation(self) -> None:
        settings = self.settings(controls={"Positive": {"result": {"min": 30, "max": 88}}})
        for statuses, bounds, expected_status, expected_limit in (
            (["below_range", "below_range"], [40, 40], ResultStatus.BELOW_RANGE, 40),
            (["above_range", "above_range"], [50, 50], ResultStatus.ABOVE_RANGE, 50),
            (["below_range", "above_range"], [40, 50], ResultStatus.UNAVAILABLE, None),
            (["above_range", "unavailable"], [50, np.nan], ResultStatus.UNAVAILABLE, None),
            (["above_range", "above_range"], [50, 100], ResultStatus.UNAVAILABLE, None),
        ):
            with self.subTest(statuses=statuses, bounds=bounds):
                samples = pd.DataFrame(
                    {
                        "sample_name": ["Positive"] * 2,
                        "sample_type": "Control",
                        "result_mean": np.nan,
                        "result_status": statuses,
                        "result_limit": bounds,
                    }
                )
                (result,) = AssayValidator(settings).validate(samples_with_statistics(samples.to_dict("records")))
                self.assertEqual(result.status, ValidationStatus.NOT_EVALUATED)
                self.assertIsNone(result.value)
                self.assertEqual(result.result_status, expected_status)
                self.assertEqual(result.result_limit, expected_limit)

    def test_concentration_bounds_fail_only_when_they_prove_a_limit_is_violated(self):
        for status, bound, limits, expected in (
            ("above_range", 9, {"min": 2.42, "max": 4.03}, ValidationStatus.FAILED),
            ("below_range", 0.11, {"min": 0.44, "max": 0.92}, ValidationStatus.FAILED),
            ("above_range", 9, {"max": 9}, ValidationStatus.FAILED),
            ("below_range", 0.11, {"min": 0.11}, ValidationStatus.FAILED),
            ("above_range", 9, {"max": 9.01}, ValidationStatus.NOT_EVALUATED),
            ("below_range", 0.11, {"min": 0.10}, ValidationStatus.NOT_EVALUATED),
            ("above_range", 9.004, {"max": 9.003}, ValidationStatus.FAILED),
            ("above_range", 9.003, {"max": 9.004}, ValidationStatus.NOT_EVALUATED),
            ("above_range", 9, {"min": 10}, ValidationStatus.NOT_EVALUATED),
            ("below_range", 0.11, {"max": 0.10}, ValidationStatus.NOT_EVALUATED),
            ("above_range", 9, {"min": 1}, ValidationStatus.NOT_EVALUATED),
            ("below_range", 0.11, {"max": 1}, ValidationStatus.NOT_EVALUATED),
            ("above_range", 9, {}, ValidationStatus.NOT_EVALUATED),
            ("unavailable", np.nan, {"min": 2.42, "max": 4.03}, ValidationStatus.NOT_EVALUATED),
        ):
            with self.subTest(status=status, bound=bound, limits=limits):
                scope = self.settings(controls={"QC": {"result": limits}})
                samples = samples_with_statistics([
                    {"sample_name": "QC", "sample_type": "Control", "result_status": status, "result_limit": bound},
                ])
                (check,) = AssayValidator(scope).validate(samples)
                self.assertEqual(check.status, expected)
                self.assertIsNone(check.value)
                self.assertIsNone(check.cv_percent)
                self.assertTrue(check.message)
                self.assertEqual(check.result_status, ResultStatus(status))
                self.assertEqual(check.result_limit, bound if np.isfinite(bound) else None)

    def test_standard_result_bounds_can_prove_a_failed_acceptance_criterion(self):
        scope = self.settings(standards={"3": {"result": {"min": 2.42, "max": 4.03}}})
        samples = samples_with_statistics([
            {"sample_name": "Standard", "sample_type": "Standard", "sample_concentration": 3, "result_status": "above_range", "result_limit": 9},
        ])
        (check,) = AssayValidator(scope).validate(samples)
        self.assertEqual(check.status, ValidationStatus.FAILED)
        self.assertIsNone(check.value)
        self.assertEqual(check.result_limit, 9)

    def test_control_checks_use_sample_means_and_exclude_other_types(self) -> None:
        settings = self.settings(controls={"Positive": {"extinction": {"response_scale": "blank_corrected", "min": 0.5, "max": 1.5}, "result": {"min": 0.9, "max": 1.1}}})
        samples = pd.DataFrame(
            {
                "sample_id": ["control-a", "control-b", "patient"],
                "sample_name": ["Positive"] * 3,
                "sample_type": ["Control", "Control", "Unknown"],
                "extinction_mean": [0, 2, 100],
                "result_mean": [0, 2, 100],
                "extinction": [100] * 3,
                "result": [100] * 3,
            }
        )
        original = samples.copy(deep=True)
        validator = AssayValidator(settings)
        results = validator.validate(samples_with_statistics(samples.to_dict("records")))
        self.assertEqual(len(results), 2)
        self.assertEqual({result.metric for result in results}, set((ValidationMetric.EXTINCTION, ValidationMetric.RESULT)))
        for result in results:
            self.assertEqual(result.sample_type, SampleType.CONTROL)
            self.assertEqual(result.value, 1)
            self.assertEqual(result.status, ValidationStatus.PASSED)
        pd.testing.assert_frame_equal(samples, original)
        self.assertEqual(validator.validate(samples_with_statistics(samples.to_dict("records"))), results)

    def test_blank_extinction_uses_separately_stored_value(self) -> None:
        settings = self.settings(blank={"extinction": {"response_scale": "raw", "max": 0.25}})
        samples = pd.DataFrame(
            {
                "sample_name": ["Blank A", "Blank B", "Blank"],
                "sample_type": ["Blank", "Blank", "Unknown"],
                "source_extinction": [0.2, 0.4, 0],
                "extinction_mean": [-0.1, 0.1, 0],
            }
        )
        (result,) = AssayValidator(settings).validate(samples_with_statistics(samples.to_dict("records")), blank_value=0.5)
        self.assertEqual(result.sample_type, SampleType.BLANK)
        self.assertAlmostEqual(result.value, 0.5)
        self.assertEqual(result.status, ValidationStatus.FAILED)
        (missing,) = AssayValidator(settings).validate(samples_with_statistics(samples.to_dict("records")))
        self.assertEqual(missing.status, ValidationStatus.NOT_EVALUATED)
        self.assertIsNone(missing.value)
        (stored,) = AssayValidator(settings).validate((), blank_value=0.25)
        self.assertEqual(stored.status, ValidationStatus.PASSED)
        self.assertEqual(stored.value, 0.25)

    def test_standards_match_nominal_concentration_and_type_not_name(self) -> None:
        settings = self.settings(standards={"250": {"extinction": {"response_scale": "blank_corrected", "min": 2}}})
        samples = pd.DataFrame(
            {
                "sample_name": ["Calibrator A", "Calibrator B", "Standaard 250", "C", "D"],
                "sample_type": ["Standard", "Standard", "Unknown", "Control", "Standard"],
                "sample_concentration": [250, "250.0", 250, 250, 50],
                "extinction_mean": [2, 2.2, 100, 100, 0.1],
            }
        )
        (result,) = AssayValidator(settings).validate(samples_with_statistics(samples.to_dict("records")))
        self.assertEqual(result.sample_name, "Standaard 250")
        self.assertEqual(result.sample_type, SampleType.STANDARD)
        self.assertAlmostEqual(result.value, 2.1)
        self.assertEqual(result.status, ValidationStatus.PASSED)

    def test_missing_targets_are_reported(self) -> None:
        settings = self.settings(
            controls={"Positive": {"result": {"min": 1}}}, blank={"extinction": {"response_scale": "raw", "max": 0.25}}, standards={"250": {"extinction": {"response_scale": "blank_corrected", "min": 2}}}
        )
        samples = pd.DataFrame({"sample_name": ["Positive"], "sample_type": ["Unknown"], "extinction_mean": [1], "result_mean": [1]})
        results = AssayValidator(settings).validate(samples_with_statistics(samples.to_dict("records")))
        self.assertEqual(len(results), 3)
        for result in results:
            self.assertEqual(result.status, ValidationStatus.NOT_EVALUATED)
            self.assertIsNone(result.value)
            self.assertIn('Geen metingen', result.message)

    def test_unavailable_metric_does_not_hide_other_results(self) -> None:
        settings = self.settings(controls={"Positive": {"extinction": {"response_scale": "blank_corrected", "min": 0.5}, "result": {"min": 36}}})
        samples = pd.DataFrame({"sample_name": ["Positive"], "sample_type": ["Control"], "extinction_mean": [1]})
        extinction, concentration = AssayValidator(settings).validate(samples_with_statistics(samples.to_dict("records")))
        self.assertEqual(extinction.status, ValidationStatus.PASSED)
        self.assertEqual(concentration.status, ValidationStatus.NOT_EVALUATED)
        self.assertIn("Elk monster", concentration.message)

    def test_incomplete_or_nonfinite_sample_means_are_not_evaluated(self) -> None:
        settings = self.settings(controls={"Positive": {"result": {"min": 1}}})
        for missing in (np.nan, np.inf, -np.inf, None, pd.NA, "invalid"):
            with self.subTest(missing=missing):
                samples = pd.DataFrame({"sample_name": ["Positive"] * 2, "sample_type": ["Control"] * 2, "result_mean": [2, missing]})
                (result,) = AssayValidator(settings).validate(samples_with_statistics(samples.to_dict("records")))
                self.assertEqual(result.status, ValidationStatus.NOT_EVALUATED)
                self.assertIsNone(result.value)
                self.assertIn('Elk monster', result.message)

    def test_standard_without_nominal_concentration_is_not_evaluated(self) -> None:
        settings = self.settings(standards={"250": {"extinction": {"response_scale": "blank_corrected", "min": 2}}})
        samples = pd.DataFrame({"sample_name": ["Standaard 250"], "sample_type": ["Standard"], "extinction_mean": [2.5]})
        (result,) = AssayValidator(settings).validate(samples_with_statistics(samples.to_dict("records")))
        self.assertEqual(result.status, ValidationStatus.NOT_EVALUATED)

    def test_empty_limits_are_not_a_pass(self) -> None:
        settings = self.settings(controls={"Positive": {"result": {}}})
        samples = pd.DataFrame({"sample_name": ["Positive"], "sample_type": ["Control"], "result_mean": [2]})
        (result,) = AssayValidator(settings).validate(samples_with_statistics(samples.to_dict("records")))
        self.assertEqual(result.status, ValidationStatus.NOT_EVALUATED)
        self.assertIn('Geen validatiegrenzen', result.message)

    def test_no_rules_returns_no_checks(self) -> None:
        samples = pd.DataFrame({"sample_name": ["Patient"], "sample_type": ["Unknown"]})
        self.assertEqual(AssayValidator(self.settings()).validate(samples_with_statistics(samples.to_dict("records"))), ())


if __name__ == "__main__":
    unittest.main()
