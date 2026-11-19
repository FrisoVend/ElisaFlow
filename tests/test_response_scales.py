"""Explicit input scales keep formula, calibration, QC and report values aligned."""

import unittest
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from fixtures import assay_metadata, samples_from_records
from test_multi_wavelength import dual_data, dual_settings

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.assays.protocol_samples import ProtocolSamples
from elisaflow.calibration.fitting import CalibrationDataError
from elisaflow.calibration.logistic import FourPLParameters
from elisaflow.configuration.assay_settings import (
    AssaySettings,
    CalibrationSettings,
    FormulaOverrideSettings,
    FormulaSettings,
    MeasurementSettings,
    ReplicateFilterSettings,
    ValidationPlan,
    ValidationScope,
)
from elisaflow.formulas.expressions import Formula
from elisaflow.models import AssayData, ResponseScale, Sample, SampleType, ValidationMetric, ValidationStatus, WellMeasurement
from elisaflow.quality_control.validation import AssayValidator
from elisaflow.reporting.pdf_report import PdfReportWriter
from elisaflow.reporting.quality_control import validation_table

try:
    from PySide6.QtWidgets import QApplication

    from elisaflow.ui.editors.calculation import FormulaCalculationEditor
    from elisaflow.ui.editors.curve import CurveEditor
    from elisaflow.ui.editors.samples import PatientEditor, SampleEditor
    from elisaflow.ui.editors.validation import ValidationEditor, ValidationScopeEditor
except ModuleNotFoundError as error:
    if error.name != "PySide6":
        raise
    QApplication = None


class ResponseScaleTests(unittest.TestCase):
    def test_direct_settings_require_and_preserve_an_explicit_response_scale(self):
        cases = (
            (FormulaSettings, (Formula("sample"),)),
            (FormulaOverrideSettings, ()),
            (CalibrationSettings, ("4pl",)),
            (CalibrationSettings, ("4pl", FourPLParameters(2.5, -1.3, 2, 0.1))),
            (ReplicateFilterSettings, (15,)),
        )
        for constructor, arguments in cases:
            with self.subTest(settings=constructor.__name__, arguments=arguments):
                with self.assertRaisesRegex(TypeError, "response_scale"):
                    constructor(*arguments)
                for invalid in (None, True, "", "automatic"):
                    with self.subTest(scale=invalid), self.assertRaisesRegex(ValueError, "response_scale"):
                        constructor(*arguments, response_scale=invalid)
                for scale in ResponseScale:
                    settings = constructor(*arguments, response_scale=scale)
                    self.assertIs(settings.response_scale, scale)
                    self.assertEqual(constructor.from_dict(settings.to_dict()), settings)

    def test_formula_settings_require_and_round_trip_explicit_response_scales(self):
        for scale in ResponseScale:
            settings = FormulaSettings.from_dict({"method": "formula", "expression": "sample", "response_scale": scale.value})
            self.assertIs(settings.response_scale, scale)
            self.assertEqual(FormulaSettings.from_dict(settings.to_dict()), settings)
        with self.assertRaisesRegex(ValueError, "response_scale"):
            FormulaSettings.from_dict({"method": "formula", "expression": "sample"})
        for scale in (None, True, [], {}, "automatic", ""):
            with self.subTest(scale=scale), self.assertRaisesRegex(ValueError, "response_scale"):
                FormulaSettings.from_dict({"method": "formula", "expression": "sample", "response_scale": scale})
        for override in ("sample", {"expression": "sample"}, {"expression": "sample", "response_scale": None}):
            with self.subTest(override=override), self.assertRaises(ValueError):
                FormulaSettings.from_dict({
                    "method": "formula", "expression": "sample", "response_scale": "raw", "control_overrides": {"QC": override},
                })

    def test_formula_scale_applies_to_own_values_control_references_overrides_and_cvs(self):
        for blank in (0.1, 0.0, None):
            samples = []
            for sample in self.data().samples:
                if sample.sample_type == SampleType.BLANK:
                    if blank is not None:
                        reference = next(well.source_extinction for well in sample.measurements if well.wavelength_nm == 620)
                        samples.append(replace(sample, measurements=tuple(
                            replace(well, source_extinction=reference + blank) if well.wavelength_nm == 405 else well
                            for well in sample.measurements
                        )))
                    continue
                if sample.sample_type in (SampleType.CONTROL, SampleType.UNKNOWN):
                    extra = tuple(replace(
                        well, well=f"B{1 if sample.sample_type == SampleType.CONTROL else 2}",
                        source_extinction=well.source_extinction + (0.2 if well.wavelength_nm == 405 else 0),
                    ) for well in sample.measurements)
                    rejected = tuple(replace(well, well="C1", source_extinction=99, disabled=True) for well in sample.measurements)
                    sample = replace(
                        sample, dilution=1,
                        measurements=(*sample.measurements, *extra, *(rejected if sample.sample_type == SampleType.CONTROL else ())),
                    )
                samples.append(sample)
            data = AssayData(tuple(samples), assay_metadata())
            for scale in ResponseScale:
                with self.subTest(blank=blank, scale=scale):
                    control_scale = ResponseScale.BLANK_CORRECTED if scale == ResponseScale.RAW else ResponseScale.RAW
                    settings = AssaySettings(
                        "test.skax", "lot",
                        calculation=FormulaSettings(
                            Formula("sample / qc", {"qc": "QC"}),
                            control_overrides={"QC": FormulaOverrideSettings("sample / qc + sample", response_scale=control_scale)}, response_scale=scale,
                        ),
                        measurement=MeasurementSettings(reference_wavelength_nm=620),
                        validation=ValidationPlan(result=ValidationScope.from_dict({"controls": {"QC": {
                            "extinction": {"max": 0.65, "response_scale": "blank_corrected"},
                            "result": {"min": 1, "max": 2},
                        }}})),
                    )
                    settings = AssaySettings.from_dict(settings.protocol_name, settings.to_dict())
                    result = AssayAnalyzer(data, settings).analyze()
                    patient = next(sample for sample in result.samples if sample.sample_type == SampleType.UNKNOWN)
                    control = next(sample for sample in result.samples if sample.sample_type == SampleType.CONTROL)
                    subtracted = (blank or 0) if scale == ResponseScale.BLANK_CORRECTED else 0
                    inputs = np.array([0.95, 1.15]) - subtracted
                    control_mean = 0.69 - subtracted
                    self.assertAlmostEqual(patient.result.mean, inputs.mean() / control_mean)
                    control_subtracted = (blank or 0) if control_scale == ResponseScale.BLANK_CORRECTED else 0
                    self.assertAlmostEqual(control.result.mean, 1 + 0.69 - control_subtracted)
                    actual = result.patient_extinction(patient)
                    self.assertEqual(actual.scale.response_scale, scale)
                    if scale == ResponseScale.BLANK_CORRECTED and blank is not None:
                        self.assertAlmostEqual(actual.scale.blank_subtracted, blank)
                    else:
                        self.assertIsNone(actual.scale.blank_subtracted)
                    self.assertAlmostEqual(actual.statistics.mean, inputs.mean())
                    self.assertAlmostEqual(actual.statistics.cv_percent, 100 * inputs.std(ddof=1) / inputs.mean())
                    measurement = next(iter(result.measurement_results.values()))
                    control_input = measurement.calculation_extinction(control)
                    self.assertEqual(control_input.scale.response_scale, control_scale)
                    self.assertAlmostEqual(control_input.statistics.mean, 0.69 - control_subtracted)
                    check = next(check for check in result.validation_results if check.metric == ValidationMetric.EXTINCTION)
                    self.assertAlmostEqual(check.value, 0.69 - (blank or 0))
                    changed_settings = replace(settings, calculation=replace(settings.calculation, response_scale=ResponseScale.RAW))
                    self.assertEqual(replace(result, settings=changed_settings).patient_extinction(patient), actual)

    def test_control_can_inherit_the_patient_expression_with_an_independent_scale(self):
        data = self.data()
        data = replace(data, samples=tuple(replace(sample, dilution=1) for sample in data.samples))
        calculation = FormulaSettings(
            Formula("sample"), response_scale=ResponseScale.BLANK_CORRECTED,
            control_overrides={"QC": FormulaOverrideSettings(response_scale=ResponseScale.RAW)},
        )
        for expression, factor in (("sample", 1), ("sample * 2", 2)):
            with self.subTest(expression=expression):
                settings = AssaySettings(
                    "test.skax", "lot", calculation=replace(calculation, formula=Formula(expression)),
                    measurement=MeasurementSettings(reference_wavelength_nm=620),
                )
                settings = AssaySettings.from_dict(settings.protocol_name, settings.to_dict())
                result = AssayAnalyzer(data, settings).analyze()
                control = next(sample for sample in result.samples if sample.sample_type == SampleType.CONTROL)
                patient = next(sample for sample in result.samples if sample.sample_type == SampleType.UNKNOWN)
                self.assertAlmostEqual(control.result.mean, factor * 0.59)
                self.assertAlmostEqual(patient.result.mean, factor * 0.85)

    def test_saved_curves_and_extinction_criteria_require_explicit_scales(self):
        for parameters in (None, {"A": 2, "B": -1, "C": 1, "D": 0}):
            with self.subTest(parameters=parameters), self.assertRaisesRegex(ValueError, "response_scale"):
                CalibrationSettings.from_dict({"method": "calibration", "fit_model": "4pl", "parameters": parameters})
        for target in ({"blank": {"extinction": {"max": 1}}}, {"controls": {"QC": {"extinction": {"min": 0}}}}):
            with self.subTest(target=target), self.assertRaisesRegex(ValueError, "response_scale"):
                ValidationScope.from_dict(target)
        for scale in (None, True, [], {}, "automatic", ""):
            with self.subTest(scale=scale), self.assertRaisesRegex(ValueError, "response_scale"):
                ValidationScope.from_dict({"controls": {"QC": {"extinction": {"min": 0, "response_scale": scale}}}})
        with self.assertRaisesRegex(ValueError, "alleen van toepassing"):
            ValidationScope.from_dict({"controls": {"QC": {"result": {"min": 0, "response_scale": "raw"}}}})
        with self.assertRaisesRegex(ValueError, "Blanco-extinctiecriteria"):
            ValidationScope.from_dict({"blank": {"extinction": {"max": 1, "response_scale": "blank_corrected"}}})

    def test_qc_scale_controls_both_mean_and_cv_of_retained_analysis_values(self):
        for kind, target in ((SampleType.CONTROL, "controls"), (SampleType.STANDARD, "standards")):
            samples = (
                Sample("a", "QC", kind, (
                    WellMeasurement(well="A1", wavelength_nm=450, source_extinction=0.1, raw_extinction=99, converted_extinction=0.3, extinction=0.2),
                    WellMeasurement(well="B1", wavelength_nm=450, source_extinction=0.5, extinction=0.4),
                    WellMeasurement(well="C1", wavelength_nm=450, source_extinction=0.7, extinction=0.6),
                    WellMeasurement(well="D1", wavelength_nm=450, source_extinction=99, extinction=98.9, included=False),
                ), concentration=1),
                Sample("b", "QC", kind, (WellMeasurement(well="E1", wavelength_nm=450, source_extinction=0.9, extinction=0.8),), concentration=1),
            )
            for scale, values, mean, status in (
                ("raw", [0.3, 0.5, 0.7, 0.9], 0.7, ValidationStatus.FAILED),
                ("blank_corrected", [0.2, 0.4, 0.6, 0.8], 0.6, ValidationStatus.PASSED),
            ):
                with self.subTest(kind=kind, scale=scale):
                    scope = ValidationScope.from_dict({target: {"QC" if target == "controls" else "1": {
                        "extinction": {"max": 0.65, "response_scale": scale},
                    }}})
                    (check,) = AssayValidator(scope).validate(samples)
                    self.assertAlmostEqual(check.value, mean)
                    self.assertAlmostEqual(check.cv_percent, 100 * np.std(values, ddof=1) / np.mean(values))
                    self.assertEqual(check.status, status)
                    self.assertEqual(check.response_scale, scale)

    @staticmethod
    def data():
        records = []
        entries = [(f"S{x}", SampleType.STANDARD, x, 1) for x in (1, 2, 4, 8) for _ in range(2)]
        entries += [("QC", SampleType.CONTROL, 3, 1), ("Patient", SampleType.UNKNOWN, 5, 2), ("Blanco", SampleType.BLANK, 0, 1)]
        for index, (name, kind, concentration, dilution) in enumerate(entries):
            response = 0.1 if kind == SampleType.BLANK else 0.2 + 0.1 * concentration + 0.01 * concentration ** 2
            reference = 0.01 * (index + 1)
            record = {
                "sample_id": name, "sample_name": name, "sample_type": kind.value,
                "sample_concentration": concentration if kind == SampleType.STANDARD else None,
                "sample_dilution": dilution, "well": f"A{index + 1}",
            }
            records.extend([
                record | {"wavelength_nm": 405, "source_extinction": response + reference},
                record | {"wavelength_nm": 620, "source_extinction": reference},
            ])
        return AssayData(samples_from_records(records), assay_metadata())

    @staticmethod
    def settings(curve_scale, qc_scale):
        return AssaySettings(
            "test.skax", "lot",
            calculation=CalibrationSettings(
                "quadratic", response_scale=curve_scale, weighting_exponent=2,
            ),
            measurement=MeasurementSettings(reference_wavelength_nm=620),
            validation=ValidationPlan(result=ValidationScope.from_dict({
                "controls": {"QC": {
                    "extinction": {"max": 0.55, "response_scale": qc_scale},
                    "result": {"min": 2.99, "max": 3.01},
                }},
                "blank": {"extinction": {"max": 0.15, "response_scale": "raw"}},
                "standards": {"1": {"extinction": {"max": 0.5, "response_scale": "blank_corrected"}}},
            })),
        )

    def test_curve_scale_drives_weights_inversion_and_patient_statistics(self):
        for curve_scale in ("raw", "blank_corrected"):
            for qc_scale in ("raw", "blank_corrected"):
                with self.subTest(curve_scale=curve_scale, qc_scale=qc_scale):
                    settings = self.settings(curve_scale, qc_scale)
                    settings = AssaySettings.from_dict(settings.protocol_name, settings.to_dict())
                    result = AssayAnalyzer(self.data(), settings).analyze()
                    standards = [well for sample in result.samples if sample.sample_type == SampleType.STANDARD for well in sample.retained]
                    blank = 0.1 if curve_scale == "blank_corrected" else 0
                    expected = np.array([well.analysis_extinction - blank for well in standards])
                    diagnostics = result.calculation_details.curve.diagnostics
                    np.testing.assert_allclose(diagnostics.weights, expected ** -2)
                    np.testing.assert_allclose(diagnostics.predictions, expected, atol=1e-10)
                    patient = next(sample for sample in result.samples if sample.sample_type == SampleType.UNKNOWN)
                    self.assertAlmostEqual(patient.result.mean, 10)
                    self.assertAlmostEqual(result.patient_extinction(patient).statistics.mean, 0.95 - blank)
                    applied_blank = result.patient_extinction(patient).scale.blank_subtracted
                    if curve_scale == "blank_corrected":
                        self.assertAlmostEqual(applied_blank, blank)
                    else:
                        self.assertIsNone(applied_blank)
                    check = next(check for check in result.validation_results if check.sample_name == "QC" and check.metric == ValidationMetric.EXTINCTION)
                    self.assertAlmostEqual(check.value, 0.59 if qc_scale == "raw" else 0.49)
                    self.assertEqual(check.status, ValidationStatus.FAILED if qc_scale == "raw" else ValidationStatus.PASSED)
                    self.assertTrue(all(
                        check.status == ValidationStatus.PASSED for check in result.validation_results
                        if check.sample_name != "QC" or check.metric != ValidationMetric.EXTINCTION
                    ))

    def test_no_blank_makes_the_two_scales_equivalent(self):
        data = self.data()
        data = replace(data, samples=tuple(sample for sample in data.samples if sample.sample_type != SampleType.BLANK))
        results = []
        for scale in ("raw", "blank_corrected"):
            result = AssayAnalyzer(data, self.settings(scale, scale)).analyze()
            self.assertIsNone(result.blank_value)
            results.append(result)
        for raw, corrected in zip(results[0].samples, results[1].samples, strict=True):
            self.assertEqual(raw.key, corrected.key)
            np.testing.assert_allclose(
                [[well.extinction, well.curve_extinction, well.result] for well in raw.measurements],
                [[well.extinction, well.curve_extinction, well.result] for well in corrected.measurements], equal_nan=True,
            )
        np.testing.assert_array_equal(
            results[0].calculation_details.curve.diagnostics.weights, results[1].calculation_details.curve.diagnostics.weights,
        )

    def test_wavelengths_keep_independent_curve_and_qc_scales(self):
        config = dual_settings().to_dict()
        config["calculation"]["wavelengths"]["450"]["response_scale"] = "blank_corrected"
        config["validation"]["wavelengths"]["450"]["controls"]["QC"]["extinction"]["response_scale"] = "raw"
        result = AssayAnalyzer(dual_data(), AssaySettings.from_dict("test.skax", config)).analyze()
        for wavelength, measurement in result.measurement_results.items():
            with self.subTest(wavelength=wavelength):
                standards = [well for sample in measurement.samples if sample.sample_type == SampleType.STANDARD for well in sample.retained]
                expected = [well.extinction if wavelength == 450 else well.analysis_extinction for well in standards]
                np.testing.assert_allclose(measurement.calculation_details.curve.diagnostics.weights, np.asarray(expected) ** -2)
        (check,) = [check for check in result.validation_results if check.sample_name == "QC" and check.metric == ValidationMetric.EXTINCTION]
        control = next(sample for sample in result.measurement_results[450].samples if sample.sample_type == SampleType.CONTROL)
        self.assertAlmostEqual(check.value, np.mean([well.analysis_extinction for well in control.retained]))
        self.assertEqual(check.response_scale, "raw")

    def test_nonpositive_corrected_fit_responses_fail_weighted_fitting_without_fallback(self):
        data = self.data()
        data = replace(data, samples=tuple(
            replace(sample, measurements=tuple(
                replace(well, source_extinction=well.source_extinction + 0.3) if well.wavelength_nm == 405 else well
                for well in sample.measurements
            )) if sample.sample_type == SampleType.BLANK else sample
            for sample in data.samples
        ))
        settings = self.settings("blank_corrected", "blank_corrected")
        settings = replace(settings, calculation=replace(settings.calculation, preprocessing=None))
        with self.assertRaisesRegex(CalibrationDataError, "positieve responswaarden"):
            AssayAnalyzer(data, settings).analyze()
        unweighted = replace(settings, calculation=replace(settings.calculation, weighting_exponent=0))
        result = AssayAnalyzer(data, unweighted).analyze()
        self.assertTrue(result.calculation_details.curve.was_fitted)
        np.testing.assert_array_equal(result.calculation_details.curve.diagnostics.weights, 1)

    def test_report_uses_each_qc_criteria_scale_independently_of_curve_input(self):
        for scale, label in (("raw", "Ext"), ("blank_corrected", "ΔExt")):
            with self.subTest(scale=scale), TemporaryDirectory() as directory:
                result = AssayAnalyzer(self.data(), self.settings("raw", scale)).analyze()
                writer = PdfReportWriter(result, analyst="operator")
                writer._initialize_report(Path(directory) / "report.pdf")
                rows = validation_table(writer.context, writer.content_width * 0.60)._cellvalues
                criterion = next(row[2].getPlainText() for row in rows if row[0].getPlainText() == "QC")
                self.assertIn(f"{label}", criterion)
                self.assertEqual("ΔExt" in criterion, scale == "blank_corrected")
                self.assertTrue(writer.write(Path(directory) / "report.pdf").read_bytes().startswith(b"%PDF-"))


@unittest.skipIf(QApplication is None, "Install the optional PySide6 UI dependency")
class ResponseScaleEditorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def test_replicate_filter_defaults_save_and_reopen_independent_choices(self):
        for has_blank in (False, True):
            for kind in ("patients", "controls", "standards"):
                with self.subTest(has_blank=has_blank, kind=kind):
                    def open_editor(config, editable):
                        if kind == "patients":
                            editor = PatientEditor(config, editable, has_blank=has_blank)
                        else:
                            options = {"control": True} if kind == "controls" else {"concentration": 1}
                            editor = SampleEditor(config, editable, has_blank=has_blank, **options)
                        self.addCleanup(editor.deleteLater)
                        return editor

                    editor = open_editor({}, True)
                    rule = editor.replicate_filter
                    expected = "blank_corrected" if has_blank else "raw"
                    self.assertEqual(rule.response_scale.currentData(), expected)
                    self.assertFalse(rule.response_scale.isEnabled())
                    rule.enabled.setChecked(True)
                    self.assertTrue(rule.response_scale.isEnabled())
                    rule.fields["max_deviation_percent"].setText("15")
                    selected = "raw" if has_blank else "blank_corrected"
                    rule.response_scale.setCurrentIndex(rule.response_scale.findData(selected))
                    config = editor.to_config()
                    self.assertEqual(config["replicate_filter"]["response_scale"], selected)
                    reopened = open_editor(config, False)
                    self.assertEqual(reopened.to_config(), config)
                    self.assertFalse(reopened.replicate_filter.response_scale.isEnabled())
                    rule.enabled.setChecked(False)
                    self.assertNotIn("replicate_filter", editor.to_config())

    def test_qc_defaults_follow_protocol_blanks_and_save_concrete_choices(self):
        for has_blank in (True, False):
            with self.subTest(has_blank=has_blank):
                protocol = ProtocolSamples(("QC",), {1.0: ("S1",)}, (405,), has_blank=has_blank)
                editor = ValidationScopeEditor({}, protocol, True)
                self.addCleanup(editor.deleteLater)
                expected = "blank_corrected" if has_blank else "raw"
                self.assertEqual(editor.response_scales["controls", "QC"].currentData(), expected)
                self.assertEqual(editor.response_scales["standards", "1.0"].currentData(), expected)
                self.assertEqual(editor.response_scales["blank", None].currentData(), "raw")
                self.assertFalse(editor.response_scales["blank", None].isEnabled())
                editor.fields["controls", "QC", "extinction"][1].setText("1")
                config = editor.to_config()
                self.assertEqual(config["controls"]["QC"]["extinction"]["response_scale"], expected)
                reopened = ValidationScopeEditor(config, protocol, False)
                self.addCleanup(reopened.deleteLater)
                self.assertEqual(reopened.to_config(), config)
                self.assertFalse(reopened.response_scales["controls", "QC"].isEnabled())

    def test_qc_scale_drafts_survive_switching_wavelength_modes(self):
        protocol = ProtocolSamples(("QC",), {}, (405, 450), has_blank=True)
        editor = ValidationEditor({}, protocol, True)
        self.addCleanup(editor.deleteLater)
        single = MeasurementSettings(reference_wavelength_nm=450)
        editor.set_measurement(single, (405,))
        editor.editors[None].fields["controls", "QC", "extinction"][1].setText("1")
        selector = editor.editors[None].response_scales["controls", "QC"]
        selector.setCurrentIndex(selector.findData("raw"))
        editor.set_measurement(MeasurementSettings("switch_on_concentration", 405, cutoff_concentration=10), (405, 450))
        config = editor.to_config()
        self.assertEqual(config["wavelengths"]["405"]["controls"]["QC"]["extinction"]["response_scale"], "raw")
        other = editor.editors[450].response_scales["controls", "QC"]
        self.assertEqual(other.currentData(), "blank_corrected")
        editor.set_measurement(single, (405,))
        self.assertEqual(editor.to_config()["result"]["controls"]["QC"]["extinction"]["response_scale"], "raw")

    def test_curve_scale_is_visible_for_fitting_and_survives_fixed_mode_changes(self):
        editor = CurveEditor({}, (), editable=True)
        self.addCleanup(editor.deleteLater)
        self.assertFalse(editor.response_scale.isHidden())
        self.assertEqual(editor.to_config()["response_scale"], "raw")
        editor.response_scale.setCurrentIndex(editor.response_scale.findData("blank_corrected"))
        editor.fixed_parameters.setChecked(True)
        for field in editor.parameters.values():
            field.setText("1")
        self.assertEqual(editor.to_config()["response_scale"], "blank_corrected")
        editor.fixed_parameters.setChecked(False)
        self.assertEqual(editor.to_config()["response_scale"], "blank_corrected")
        reopened = CurveEditor(editor.to_config(), (), editable=False)
        self.addCleanup(reopened.deleteLater)
        self.assertEqual(reopened.to_config(), editor.to_config())
        self.assertFalse(reopened.response_scale.isEnabled())

    def test_formula_scale_defaults_follow_blanks_and_reopen_with_control_overrides(self):
        for has_blank in (True, False):
            with self.subTest(has_blank=has_blank):
                initial = {"calculation": {"method": "formula", "expression": "sample", "control_overrides": {"QC": {"expression": "sample * 2", "response_scale": "blank_corrected"}}}}
                editor = FormulaCalculationEditor(initial, ("QC",), True, has_blank=has_blank)
                self.addCleanup(editor.deleteLater)
                self.assertEqual(editor.response_scale.currentData(), "blank_corrected" if has_blank else "raw")
                selected = "raw" if has_blank else "blank_corrected"
                editor.response_scale.setCurrentIndex(editor.response_scale.findData(selected))
                config = editor.to_config()
                settings = FormulaSettings.from_dict(config["calculation"])
                self.assertEqual(settings.response_scale, selected)
                self.assertEqual(settings.control_overrides["QC"].expression, "sample * 2")
                reopened = FormulaCalculationEditor(config, ("QC",), False, has_blank=has_blank)
                self.addCleanup(reopened.deleteLater)
                self.assertEqual(reopened.to_config(), config)
                self.assertFalse(reopened.response_scale.isEnabled())

    def test_control_formula_scale_is_independent_and_linking_tracks_the_patient(self):
        editor = FormulaCalculationEditor({}, ("QC",), True, has_blank=True)
        self.addCleanup(editor.deleteLater)
        card = editor.overrides["QC"]
        self.assertFalse(card.response_scale.isEnabled())
        card.use_sample_result.setChecked(False)
        self.assertTrue(card.response_scale.isEnabled())
        card.response_scale.setCurrentIndex(card.response_scale.findData("raw"))
        self.assertEqual(editor.response_scale.currentData(), "blank_corrected")
        config = editor.to_config()
        self.assertEqual(config["calculation"]["control_overrides"]["QC"], {"response_scale": "raw"})
        reopened = FormulaCalculationEditor(config, ("QC",), True, has_blank=True)
        self.addCleanup(reopened.deleteLater)
        self.assertEqual(reopened.overrides["QC"].response_scale.currentData(), "raw")
        self.assertFalse(reopened.overrides["QC"].use_sample_result.isChecked())
        card.use_sample_result.setChecked(True)
        self.assertEqual(card.response_scale.currentData(), "blank_corrected")
        self.assertFalse(card.response_scale.isEnabled())
        editor.response_scale.setCurrentIndex(editor.response_scale.findData("raw"))
        self.assertEqual(card.response_scale.currentData(), "raw")
        self.assertNotIn("control_overrides", editor.to_config()["calculation"])
        card.use_sample_result.setChecked(False)
        editor.response_scale.setCurrentIndex(editor.response_scale.findData("blank_corrected"))
        self.assertEqual(card.response_scale.currentData(), "raw")
