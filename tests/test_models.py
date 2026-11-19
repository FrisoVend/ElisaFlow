import unittest
from dataclasses import FrozenInstanceError, replace

import numpy as np
from fixtures import assay_metadata
from table_helpers import sample_table, well_table

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.calibration.fitting import FitDiagnostics
from elisaflow.configuration.assay_settings import AssaySettings, FormulaSettings, ValidationScope
from elisaflow.formulas.expressions import Formula
from elisaflow.models import (
    AssayData,
    ResponseScale,
    ResultQualification,
    ResultStatus,
    Sample,
    SampleType,
    Statistics,
    ValidationStatus,
    WellMeasurement,
)
from elisaflow.quality_control.validation import AssayValidator


class ResultOwnershipTests(unittest.TestCase):
    def test_input_and_results_are_immutable_and_tables_are_independent_views(self) -> None:
        readings = [WellMeasurement(wavelength_nm=450, source_extinction=0.5, well="A1")]
        control = Sample("control", "QC", SampleType.CONTROL, readings)
        patients = Sample("patient", "Patient", SampleType.UNKNOWN, (WellMeasurement(wavelength_nm=450, source_extinction=1, well="B1"),))
        samples = [control, patients]
        data = AssayData(samples, assay_metadata())
        readings.clear()
        samples.clear()
        settings = AssaySettings(
            "test.skax", "lot", calculation=FormulaSettings(Formula("sample / qc", {"qc": "QC"}), {}, response_scale="blank_corrected")
        )
        analyzer = AssayAnalyzer(data, settings)
        result = analyzer.analyze()
        self.assertEqual(result.samples[0].result.mean, 1)
        self.assertEqual(result.samples[1].result.mean, 2)
        for obj, attribute, value in ((data, "samples", ()), (control, "name", "Changed"), (control.measurements[0], "source_extinction", 999)):
            with self.subTest(attribute=attribute), self.assertRaises(FrozenInstanceError):
                setattr(obj, attribute, value)
        for table, column in ((sample_table(result.samples), "extinction_mean"), (well_table(result.samples), "extinction")):
            table.iloc[0, table.columns.get_loc(column)] = 999
            self.assertEqual(ResponseScale.BLANK_CORRECTED.statistics(result.samples[0]).mean, 0.5)
        self.assertTrue(np.isnan(data.samples[0].measurements[0].extinction))
        self.assertEqual(ResponseScale.BLANK_CORRECTED.statistics(result.samples[0]).mean, 0.5)
        with self.assertRaises(TypeError):
            result.measurement_results[450] = None

    def test_measurements_require_a_well_and_a_positive_integer_wavelength(self) -> None:
        valid = {"source_extinction": 1, "well": "A1", "wavelength_nm": 450}
        for missing in ("well", "wavelength_nm"):
            with self.subTest(missing=missing), self.assertRaises(TypeError):
                WellMeasurement(**{key: value for key, value in valid.items() if key != missing})
        for coordinate in (None, "", " ", False):
            with self.subTest(coordinate=coordinate), self.assertRaisesRegex(ValueError, "wellcoördinaat"):
                WellMeasurement(**(valid | {"well": coordinate}))
        for wavelength in (None, 0, -1, True, np.bool_(True), 450.0, "450"):
            with self.subTest(wavelength=wavelength), self.assertRaisesRegex(ValueError, "positieve gehele getallen"):
                WellMeasurement(**(valid | {"wavelength_nm": wavelength}))

    def test_raw_frame_arguments_are_not_supported(self) -> None:
        with self.assertRaisesRegex(ValueError, "Sample-objecten"):
            AssayData([{"sample_id": "patient"}], assay_metadata())
        with self.assertRaises(TypeError):
            AssayData(well_data=[], metadata=assay_metadata())

    def test_statistics_are_derived_from_retained_replicates_without_cached_state(self) -> None:
        statistics = ResponseScale.BLANK_CORRECTED.statistics
        sample = Sample(
            "p",
            "Patient",
            SampleType.UNKNOWN,
            (
                WellMeasurement(well="C1", wavelength_nm=450, source_extinction=1, extinction=1),
                WellMeasurement(well="D1", wavelength_nm=450, source_extinction=3, extinction=3),
                WellMeasurement(well="E1", wavelength_nm=450, source_extinction=999, extinction=999, included=False),
            ),
        )
        self.assertEqual(len(sample.retained), 2)
        self.assertEqual(statistics(sample).mean, 2)
        self.assertAlmostEqual(statistics(sample).sd, np.sqrt(2))
        changed = replace(sample, measurements=(replace(sample.measurements[0], extinction=5),))
        self.assertEqual(statistics(changed).mean, 5)
        self.assertTrue(np.isnan(statistics(changed).sd))
        self.assertEqual(statistics(sample).mean, 2)
        self.assertTrue(np.isnan(Statistics.calculate([1, np.nan]).mean))
        self.assertTrue(np.isnan(Statistics.calculate([0, 0]).cv_percent))

    def test_diagnostics_copy_input_arrays_and_reject_in_place_writes(self) -> None:
        values = np.array([1.0, 2.0])
        diagnostics = FitDiagnostics(values, values, values, values, values)
        values[0] = 999
        for name in ("residuals", "predictions", "back_calculated", "relative_errors", "weights"):
            array = getattr(diagnostics, name)
            np.testing.assert_array_equal(array, [1, 2])
            with self.assertRaises(ValueError):
                array[0] = 999

    def test_negative_mean_has_positive_cv_without_changing_mean_or_sd(self) -> None:
        statistics = Statistics.calculate([-0.9, -0.1])
        self.assertAlmostEqual(statistics.mean, -0.5)
        self.assertAlmostEqual(statistics.sd, np.sqrt(0.32))
        self.assertAlmostEqual(statistics.cv_percent, 100 * np.sqrt(0.32) / 0.5)

    def test_contradictory_result_fields_are_unavailable_in_statistics_tables_and_qc(self):
        scope = ValidationScope.from_dict({"controls": {"QC": {"result": {"min": 0, "max": 20}}}})
        for values in (
            {"result": 12, "result_status": ResultStatus.ABOVE_RANGE, "result_limit": 10},
            {"result": 12, "result_status": ResultStatus.UNAVAILABLE},
            {"result": 12, "result_status": ResultStatus.IN_RANGE, "result_limit": 10},
            {"result_status": ResultStatus.IN_RANGE},
            {"result_status": ResultStatus.ABOVE_RANGE},
        ):
            with self.subTest(values=values):
                sample = Sample("c", "QC", SampleType.CONTROL, (WellMeasurement(well="F1", wavelength_nm=450, source_extinction=1, **values),))
                self.assertTrue(np.isnan(sample.result.mean))
                self.assertEqual(sample.result_qualification, ResultQualification())
                self.assertIsNone(sample.result_limit)
                self.assertEqual(sample_table((sample,)).iloc[0]["result_status"], "unavailable")
                (check,) = AssayValidator(scope).validate((sample,))
                self.assertIsNone(check.value)
                self.assertEqual(check.result_status, ResultStatus.UNAVAILABLE)
                self.assertEqual(check.status, ValidationStatus.NOT_EVALUATED)

    def test_excluded_well_does_not_overrule_a_retained_numeric_result(self):
        sample = Sample("p", "Patient", SampleType.UNKNOWN, (
            WellMeasurement(well="G1", wavelength_nm=450, source_extinction=1, result=12, result_status=ResultStatus.IN_RANGE),
            WellMeasurement(well="H1", wavelength_nm=450, source_extinction=1, included=False, result=99, result_status=ResultStatus.ABOVE_RANGE, result_limit=10),
        ))
        self.assertEqual(sample.result_status, ResultStatus.IN_RANGE)
        self.assertEqual(sample.result.mean, 12)

    def test_result_qualification_rejects_missing_nonfinite_or_misplaced_bounds(self):
        for status, limit in (
            (ResultStatus.ABOVE_RANGE, None), (ResultStatus.BELOW_RANGE, float("inf")),
            (ResultStatus.IN_RANGE, 10), (ResultStatus.UNAVAILABLE, 10),
        ):
            with self.subTest(status=status, limit=limit), self.assertRaises(ValueError):
                ResultQualification(status, limit)


if __name__ == "__main__":
    unittest.main()
