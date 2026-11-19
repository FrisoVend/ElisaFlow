import json
import unittest
from copy import deepcopy
from dataclasses import FrozenInstanceError, replace
from pathlib import Path

from elisaflow.calibration.curves import create_calibration_curve
from elisaflow.calibration.logistic import FourPLParameters
from elisaflow.configuration.assay_settings import (
    AssaySettings,
    CalibrationSettings,
    FormulaSettings,
    MeasurementSettings,
    ReportSettings,
    ResultSettings,
    SampleSettings,
    ValidationCriterion,
    ValidationLimits,
    ValidationScope,
    WavelengthCalibrationSettings,
)
from elisaflow.formulas.expressions import Formula
from elisaflow.models import SampleType, ValidationMetric


class SettingsTests(unittest.TestCase):
    def test_control_results_round_trip_and_validate_names_and_metadata(self):
        calculation = FormulaSettings(
            Formula("sample"), control_results={"Positief ": ResultSettings(" Controle-index ", " index ")}, response_scale="blank_corrected"
        )
        self.assertEqual(calculation.control_results, {"Positief ": ResultSettings("Controle-index", "index")})
        self.assertEqual(calculation.to_dict()["control_results"], {"Positief ": {"label": "Controle-index", "unit": "index"}})
        self.assertEqual(FormulaSettings.from_dict(calculation.to_dict()), calculation)
        for results in (
            None, [], {"": {"label": "Score"}}, {"Positief": "Controle-index"}, {"Positief": None},
            {"Positief": {"label": " "}}, {"Positief": {"label": 1}}, {"Positief": {"label": "Score", "unit": " "}},
            {"Positief": {"label": "Score", "unit": 1}}, {"Positief": {"label": "Score", "units": "%"}},
        ):
            with self.subTest(results=results), self.assertRaises(ValueError):
                FormulaSettings.from_dict({"response_scale": "blank_corrected", "method": "formula", "expression": "sample", "control_results": results})
        with self.assertRaisesRegex(ValueError, "ResultSettings"):
            FormulaSettings(Formula("sample"), control_results={"Positief": {"label": "Score"}}, response_scale="blank_corrected")

    def test_old_label_only_control_format_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Onbekende velden.*control_labels"):
            FormulaSettings.from_dict({"response_scale": "blank_corrected", "method": "formula", "expression": "sample", "control_labels": {"Positief": "Score"}})

    def test_report_defaults_and_partial_overrides_round_trip(self) -> None:
        config = {
            "lot_number": "lot",
            "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / cal", "references": {"cal": "Cal"}},
            "result": {"label": "Ratio"},
        }
        settings = AssaySettings.from_dict("test.skax", config)
        self.assertEqual(
            settings.report.to_dict(),
            {
                "result_decimals": 1,
                "extinction_decimals": 3,
                "calibration_decimals": 1,
                "percentage_decimals": 1,
                "cv_warning_percent": 20,
                "recovery_warning_percent": 20,
                "show_result_cv": True,
            },
        )
        overrides = {"calibration_decimals": 3, "cv_warning_percent": 12.5, "recovery_warning_percent": 150}
        updated = AssaySettings.from_dict("test.skax", config | {"report": overrides})
        self.assertEqual(updated.report, ReportSettings(calibration_decimals=3, cv_warning_percent=12.5, recovery_warning_percent=150))
        exported = json.loads(json.dumps(updated.to_dict(), allow_nan=False))
        self.assertEqual(exported["report"], updated.report.to_dict())
        self.assertEqual(AssaySettings.from_dict("test.skax", exported), updated)
        with self.assertRaises(FrozenInstanceError):
            updated.report.result_decimals = 3
        with self.assertRaisesRegex(ValueError, 'instellingenobject voor report'):
            replace(settings, report={})

    def test_report_rejects_invalid_precision_and_thresholds_at_both_boundaries(self) -> None:
        for name in ("result_decimals", "extinction_decimals", "calibration_decimals", "percentage_decimals"):
            for value in (-1, 7, 1.5, 1.0, True, "3", None):
                with self.subTest(field=name, value=value):
                    with self.assertRaisesRegex(ValueError, name):
                        ReportSettings.from_dict({name: value})
                    with self.assertRaisesRegex(ValueError, name):
                        replace(ReportSettings(), **{name: value})
            for value in (0, 6):
                self.assertEqual(getattr(ReportSettings.from_dict({name: value}), name), value)
        for name in ("cv_warning_percent", "recovery_warning_percent"):
            for value in (-0.1, float("nan"), float("inf"), True, "20", None):
                with self.subTest(field=name, value=value):
                    with self.assertRaisesRegex(ValueError, name):
                        ReportSettings.from_dict({name: value})
                    with self.assertRaisesRegex(ValueError, name):
                        replace(ReportSettings(), **{name: value})
            self.assertEqual(getattr(ReportSettings.from_dict({name: 0}), name), 0)

    def test_assay_rejects_malformed_report_objects_and_unknown_fields(self) -> None:
        config = {
            "lot_number": "lot",
            "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / cal", "references": {"cal": "Cal"}},
            "result": {"label": "Ratio"},
        }
        for report in (None, [], 20, {"ratio_decimal": 1}):
            with self.subTest(report=report), self.assertRaisesRegex(ValueError, 'Rapportinstellingen'):
                AssaySettings.from_dict("test.skax", config | {"report": report})

    def test_every_bundled_assay_can_be_reloaded_from_serialized_settings(self) -> None:
        directory = Path(__file__).resolve().parents[1] / "configs"
        for path in directory.glob("*.json"):
            for protocol_name, config in json.loads(path.read_text(encoding="utf-8-sig")).items():
                with self.subTest(protocol=protocol_name):
                    settings = AssaySettings.from_dict(protocol_name, config)
                    exported = json.loads(json.dumps(settings.to_dict(), allow_nan=False))
                    self.assertEqual(AssaySettings.from_dict(protocol_name, exported), settings)

    def test_serialization_preserves_optional_fields_and_standard_precision(self):
        preprocessing = {"expression": "sample * 0.86 / cal", "references": {"cal": "Cal"}}
        common = {
            "lot_number": "lot",
            "controls": {"Positive": {"lot_number": "control lot"}},
            "standards": {"0.12345678901234567": {"lot_number": "standard lot"}},
            "blank": {"lot_number": "blank lot"},
            "validation": {"result": {"controls": {"Positive": {"result": {"min": 0, "max": 0}, "extinction": {"response_scale": "blank_corrected"}}}}},
        }
        for calculation in (
            {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / cal", "references": {"cal": "Cal"}, "control_overrides": {"Positive": {"expression": "sample", "response_scale": "blank_corrected"}}},
            {
                "response_scale": "raw",
                "method": "calibration",
                "fit_model": "4pl",
                "weighting_exponent": 1.0,
                "concentration_range": [0, 500],
            },
            {"response_scale": "blank_corrected", "method": "calibration", "fit_model": "4pl", "parameters": {"A": 2, "B": -1, "C": 1, "D": 0}, "preprocessing": preprocessing},
        ):
            with self.subTest(calculation=calculation):
                settings = AssaySettings.from_dict("test.skax", common | {"calculation": calculation, "result": {"label": "Score", "unit": "%"}})
                exported = json.loads(json.dumps(settings.to_dict(), allow_nan=False))
                self.assertEqual(AssaySettings.from_dict("test.skax", exported), settings)

    def test_result_label_and_unit_are_immutable_and_validated(self):
        settings = ResultSettings.from_dict({"label": " Score ", "unit": " % "})
        self.assertEqual(settings.label, "Score")
        self.assertEqual(settings.unit, "%")
        self.assertIsNone(ResultSettings().unit)
        with self.assertRaises(FrozenInstanceError):
            settings.unit = "mg/ml"
        for name in ("label", "unit"):
            for value in ("", " ", 1, True, [], {}):
                with self.subTest(name=name, value=value), self.assertRaises(ValueError):
                    ResultSettings.from_dict({name: value})

    def test_formula_calculation_is_typed_and_detached_from_config(self):
        config = {
            "response_scale": "blank_corrected",
            "method": "formula",
            "expression": "sample / cal",
            "references": {"cal": "Cal"},
            "control_overrides": {"Positive": {"expression": "sample / (cal * 0.5)", "response_scale": "blank_corrected"}},
        }
        settings = FormulaSettings.from_dict(config)
        config["references"]["cal"] = "Other"
        config["control_overrides"].clear()
        self.assertEqual(settings.formula.references["cal"], "Cal")
        self.assertEqual(settings.overrides["Positive"].evaluate(2, {"Cal": 4}), 1)
        with self.assertRaises(TypeError):
            settings.control_overrides["Other"] = "sample"
        with self.assertRaises(FrozenInstanceError):
            settings.formula = Formula("sample")

    def test_old_schema_and_malformed_calculations_are_rejected(self):
        base = {"lot_number": "lot", "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample"}}
        for key in ("type", "cutoff_control", "curve_fitting", "concentration_unit"):
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, 'Onbekende velden'):
                AssaySettings.from_dict("test.skax", base | {key: {}})
        for calculation in (
            None,
            [],
            {},
            {"method": "ratio"},
            {"response_scale": "blank_corrected", "method": "formula"},
            {"response_scale": "blank_corrected", "method": "formula", "expression": "sample", "control_overrides": []},
            {"response_scale": "blank_corrected", "method": "formula", "expression": "sample", "control_overrides": {"P": {"expression": "unknown", "response_scale": "blank_corrected"}}},
        ):
            with self.subTest(calculation=calculation), self.assertRaises(ValueError):
                AssaySettings.from_dict("test.skax", base | {"calculation": calculation})
        with self.assertRaisesRegex(ValueError, 'Onbekende velden'):
            SampleSettings.from_dict({"ratio_reference": "Cal"})

    def test_calibration_preprocessing_is_a_typed_formula(self):
        config = {"expression": "sample * 0.86 / standard", "references": {"standard": "Standaard"}}
        curve = CalibrationSettings.from_dict({
            "response_scale": "raw", "method": "calibration", "fit_model": "4pl",
            "parameters": {"A": 2, "B": -1, "C": 1, "D": 0}, "preprocessing": config,
        })
        config["references"]["standard"] = "Other"
        expected = Formula("sample * 0.86 / standard", {"standard": "Standaard"})
        self.assertEqual(curve.preprocessing, expected)
        self.assertIsNone(CalibrationSettings("4pl", response_scale="raw").preprocessing)
        with self.assertRaises(FrozenInstanceError):
            curve.preprocessing = None

    def test_fitted_curves_reject_preprocessing_in_direct_and_saved_settings(self):
        formula = Formula("sample * 2")
        with self.assertRaisesRegex(ValueError, "Voorbewerking.*vaste parameters"):
            CalibrationSettings("4pl", response_scale="raw", preprocessing=formula)
        with self.assertRaisesRegex(ValueError, "Voorbewerking.*vaste parameters"):
            CalibrationSettings.from_dict({
                "method": "calibration", "fit_model": "4pl", "response_scale": "raw", "preprocessing": formula.to_dict(),
            })

    def test_preprocessing_rejects_old_steps_and_untyped_formulas(self):
        for config in (None, [], {"multiply": 0.43}, [{"operation": "multiply", "value": 2}]):
            with self.subTest(config=config):
                with self.assertRaises(ValueError):
                    CalibrationSettings.from_dict({"response_scale": "raw", "method": "calibration", "fit_model": "4pl", "preprocessing": config})
        for value in ([], {}, "sample", ()):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    CalibrationSettings("4pl", preprocessing=value, response_scale="raw")

    def test_parameters_are_parsed_once_and_detached_from_json(self) -> None:
        parameters = {"A": 2.0, "B": -1.4, "C": 2.5, "D": 0.1}
        settings = CalibrationSettings.from_dict({"response_scale": "blank_corrected", "method": "calibration", "fit_model": "4pl", "parameters": parameters})
        parameters["C"] = 100.0
        curve = create_calibration_curve(settings)
        self.assertIsInstance(settings.parameters, FourPLParameters)
        self.assertIs(curve.parameters, settings.parameters)
        self.assertEqual(settings.parameters.C, 2.5)
        with self.assertRaises(FrozenInstanceError):
            settings.parameters.C = 3.0

    def test_sample_preprocessing_is_rejected_in_every_sample_section(self) -> None:
        steps = {"expression": "sample * 10"}
        base = {
            "lot_number": "lot",
            "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / c", "references": {"c": "c"}},
            "result": {"label": "Ratio"},
        }
        for json_settings in (
            {"blank": {"preprocessing": steps}},
            {"standards": {"1": {"preprocessing": steps}}},
            {"controls": {"c": {"preprocessing": steps}}},
        ):
            with (
                self.subTest(config=json_settings),
                self.assertRaisesRegex(ValueError, "Onbekende velden.*preprocessing"),
            ):
                AssaySettings.from_dict("test.skax", base | json_settings)
        with self.assertRaises(TypeError):
            SampleSettings(preprocessing=Formula("sample * 10"))

    def test_removed_relative_error_setting_is_rejected_without_compatibility(self) -> None:
        single = {"response_scale": "raw", "method": "calibration", "fit_model": "4pl", "maximum_relative_error": 0.1}
        dual = {"method": "calibration", "wavelengths": {"405": {"response_scale": "raw", "fit_model": "4pl"}, "450": {"response_scale": "raw", "fit_model": "4pl"}}, "maximum_relative_error": 0.1}
        for kind, config in ((CalibrationSettings, single), (WavelengthCalibrationSettings, dual)):
            with self.subTest(kind=kind.__name__), self.assertRaisesRegex(ValueError, "Onbekende velden.*maximum_relative_error"):
                kind.from_dict(config)
        with self.assertRaises(TypeError):
            CalibrationSettings("4pl", maximum_relative_error=0.1, response_scale="raw")
        with self.assertRaises(TypeError):
            WavelengthCalibrationSettings(
                {405: CalibrationSettings("4pl", response_scale="raw"), 450: CalibrationSettings("4pl", response_scale="raw")},
                maximum_relative_error=0.1,
            )

    def test_parameters_are_typed_and_obsolete_keys_are_rejected(self) -> None:
        parameters = {"A": 2.0, "B": -1.4, "C": 2.5, "D": 0.1}
        settings = CalibrationSettings.from_dict({"response_scale": "blank_corrected", "method": "calibration", "fit_model": "4pl", "parameters": parameters})
        self.assertEqual(settings.parameters.B, -1.4)
        with self.assertRaisesRegex(ValueError, 'Onbekende velden'):
            CalibrationSettings.from_dict(
                {
                    "response_scale": "blank_corrected",
                    "method": "calibration",
                    "fit_model": "4pl",
                    "parameters": parameters,
                    "fit_parameters": parameters,
                }
            )

    def test_malformed_curve_settings_raise_value_errors(self) -> None:
        cases = [
            [],
            None,
            {},
            {"response_scale": "raw", "fit_model": []},
            {"response_scale": "raw", "fit_model": {}},
            {"response_scale": "raw", "fit_model": "6pl"},
            {"response_scale": "raw", "fit_model": "4pl", "weighting_exponent": {}},
            {"response_scale": "raw", "fit_model": "4pl", "weighting_exponent": "automatic"},
            {"response_scale": "raw", "fit_model": "4pl", "weigthing": "1/y"},
            {"response_scale": "blank_corrected", "fit_model": "4pl", "parameters": []},
            {"response_scale": "blank_corrected", "fit_model": "4pl", "parameters": {}},
        ]
        for config in cases:
            with self.subTest(config=config), self.assertRaises(ValueError):
                CalibrationSettings.from_dict(config)
        for value in ([0], [1, 1], [-1, 2], [0, float("inf")], [False, 2], "0,2"):
            with self.subTest(range=value), self.assertRaises(ValueError):
                CalibrationSettings.from_dict({"response_scale": "raw", "method": "calibration", "fit_model": "4pl", "concentration_range": value})

    def test_fitting_options_cannot_be_used_with_fixed_parameters(self) -> None:
        parameters = FourPLParameters(2.0, -1.4, 2.5, 0.1)
        for options in ({"weighting_exponent": 1.0},):
            with self.subTest(options=options), self.assertRaisesRegex(ValueError, 'alleen van toepassing'):
                CalibrationSettings("4pl", parameters, response_scale="blank_corrected", **options)

    def test_direct_parameters_enforce_numeric_and_model_invariants(self) -> None:
        valid = {"A": 2.0, "B": -1.4, "C": 2.5, "D": 0.1}
        for name in valid:
            for value in (True, "1", float("nan"), float("inf"), 10**1000):
                with self.subTest(parameter=name, value=value), self.assertRaises(ValueError):
                    FourPLParameters(**(valid | {name: value}))
        for invalid in ({"B": 0}, {"C": 0}, {"C": -1}, {"D": 2.0}):
            with self.subTest(parameters=invalid), self.assertRaises(ValueError):
                FourPLParameters(**(valid | invalid))

    def test_public_constructors_reject_untyped_nested_objects(self) -> None:
        for options in (
            {"calculation": {"method": "calibration"}, "result": {"label": "Conc"}},
            {"blank": {}},
            {"controls": {"c": {}}},
            {"standards": {1: {}}},
            {"standards": []},
        ):
            with self.subTest(options=options), self.assertRaises(ValueError):
                AssaySettings("test.skax", "lot", **({"calculation": FormulaSettings(Formula("sample"), response_scale="blank_corrected")} | options))
        with self.assertRaises(ValueError):
            SampleSettings({"extinction": {"response_scale": "blank_corrected", "max": 1}})
        with self.assertRaises(ValueError):
            CalibrationSettings("4pl", parameters={"A": 2, "B": 1, "C": 1, "D": 0}, response_scale="blank_corrected")

    def test_assay_schema_rejects_unknown_fields_and_invalid_nested_settings(self) -> None:
        valid = {
            "lot_number": "lot",
            "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / c", "references": {"c": "c"}},
            "result": {"label": "Ratio"},
        }
        for invalid in (
            {"control": {}},
            {"controls": []},
            {"blank": None},
            {"standards": []},
            {"lot_number": " "},
            {"calculation": {}},
            {"result": []},
        ):
            with self.subTest(settings=invalid), self.assertRaises(ValueError):
                AssaySettings.from_dict("test.skax", valid | invalid)
        for invalid in (
            {"standards": {"1": {}, "1.0": {}}},
            {"validation": {"result": {"standards": {"1": {}, "1.0": {}}}}},
        ):
            with self.subTest(settings=invalid), self.assertRaisesRegex(ValueError, 'Dubbele standaardconcentratie'):
                AssaySettings.from_dict("test.skax", valid | invalid)

    def test_validation_rules_and_standard_keys_are_normalized_and_immutable(self) -> None:
        settings = AssaySettings.from_dict(
            "test.skax",
            {
                "lot_number": "lot",
                "standards": {"1.0": {}},
                "validation": {"result": {"standards": {"1.0": {"extinction": {"response_scale": "blank_corrected", "min": 0, "max": 2}}}}},
                "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / c", "references": {"c": "c"}},
                "result": {"label": "Ratio"},
            },
        )
        limits = settings.validation.result.standards[1][ValidationMetric.EXTINCTION]
        self.assertEqual(limits, ValidationLimits(0.0, 2.0, "blank_corrected"))
        with self.assertRaises(TypeError):
            settings.standards[2] = SampleSettings()
        with self.assertRaises(TypeError):
            settings.validation.result.standards[1][ValidationMetric.RESULT] = ValidationLimits()

    def test_criteria_preserve_identity_order_and_config_without_serializing_derived_state(self):
        config = {
            "controls": {"QC": {"extinction": {"min": 0, "response_scale": "blank_corrected"}, "result": {"max": 20}}},
            "blank": {"extinction": {"max": 0.2, "response_scale": "raw"}},
            "standards": {"1.25": {"result": {"min": 1}}},
        }
        scope = ValidationScope.from_dict(config)
        self.assertEqual([(rule.sample_type, rule.metric) for rule in scope.criteria], [
            (SampleType.CONTROL, ValidationMetric.EXTINCTION), (SampleType.CONTROL, ValidationMetric.RESULT),
            (SampleType.BLANK, ValidationMetric.EXTINCTION), (SampleType.STANDARD, ValidationMetric.RESULT),
        ])
        self.assertEqual(scope.criteria[-1].nominal_concentration, 1.25)
        self.assertEqual(scope.to_dict(), config)
        self.assertEqual(ValidationScope.from_dict(scope.to_dict()), scope)
        with self.assertRaises(FrozenInstanceError):
            scope.criteria[0].metric = ValidationMetric.RESULT

    def test_direct_criteria_reject_invalid_sample_metric_and_scale_combinations(self):
        for kind, metric, scale, concentration in (
            (SampleType.UNKNOWN, ValidationMetric.RESULT, None, None),
            (SampleType.CONTROL, ValidationMetric.EXTINCTION, None, None),
            (SampleType.BLANK, ValidationMetric.EXTINCTION, "blank_corrected", None),
            (SampleType.CONTROL, ValidationMetric.RESULT, "raw", None),
            (SampleType.STANDARD, ValidationMetric.RESULT, None, None),
            (SampleType.CONTROL, ValidationMetric.RESULT, None, 2),
        ):
            with self.subTest(kind=kind, metric=metric, scale=scale, concentration=concentration), self.assertRaises(ValueError):
                ValidationCriterion(kind, "QC", metric, ValidationLimits(response_scale=scale), concentration)


class MeasurementSettingsTests(unittest.TestCase):
    def test_only_processing_selections_are_serialized(self):
        config = {"mode": "switch_on_concentration", "cutoff_concentration": 1000.0, "preferred_wavelength_nm": 450, "reference_wavelength_nm": 620}
        settings = MeasurementSettings.from_dict(config)
        self.assertEqual(settings.to_dict(), config)
        self.assertEqual(settings.analysis_wavelengths((405, 450, 620)), (405, 450))
        self.assertEqual(MeasurementSettings.from_dict(settings.to_dict()), settings)

    def test_json_wavelength_catalogues_are_rejected(self):
        for values in (None, [], [405, 450], "405, 450"):
            with self.subTest(values=values), self.assertRaisesRegex(ValueError, 'Onbekende.*wavelengths_nm'):
                MeasurementSettings.from_dict({"wavelengths_nm": values})

    def test_selections_must_match_the_skax_wavelengths(self):
        settings = MeasurementSettings("switch_on_concentration", 450, 620, cutoff_concentration=1000)
        for values in ((450, 620), (405, 450), (405, 620, 630), (405, 450, 620, 630)):
            with self.subTest(values=values), self.assertRaises(ValueError):
                settings.analysis_wavelengths(values)


class WavelengthCalibrationSettingsTests(unittest.TestCase):
    def setUp(self):
        self.config = {
            "lot_number": "lot",
            "measurement": {"mode": "switch_on_concentration", "cutoff_concentration": 1000.0, "preferred_wavelength_nm": 450},
            "calculation": {
                "method": "calibration",
                "wavelengths": {
                    "405": {"response_scale": "raw", "fit_model": "quadratic", "weighting_exponent": 0.0,
                            "parameters": {"A": 0.03, "B": 0.003, "C": 0.000005}, "preprocessing": {"expression": "sample * 2"}},
                    "450": {"response_scale": "raw", "fit_model": "4pl", "weighting_exponent": 0.0,
                            "parameters": {"A": 2, "B": -1, "C": 1, "D": 0}, "preprocessing": {"expression": "sample / 3"}},
                },
                "concentration_range": [1, 300],
            },
        }

    def test_round_trip_preserves_independent_preprocessing_and_shared_limits(self):
        original = deepcopy(self.config)
        settings = AssaySettings.from_dict("test.skax", self.config)
        self.assertEqual(self.config, original)
        self.assertEqual(settings.to_dict()["calculation"], self.config["calculation"])
        self.assertEqual(AssaySettings.from_dict("test.skax", settings.to_dict()), settings)
        for wavelength, model in ((405, "quadratic"), (450, "4pl")):
            effective = settings.calculation_for_wavelength(wavelength)
            self.assertEqual(effective.fit_model, model)
            self.assertEqual(effective.weighting_exponent, 0)
            self.assertEqual(effective.concentration_range, (1, 300))
            self.assertEqual(effective.preprocessing, Formula("sample * 2" if wavelength == 405 else "sample / 3"))
            self.assertEqual(effective.response_scale, "raw")
        self.config["calculation"]["wavelengths"]["405"]["fit_model"] = "cubic"
        self.assertEqual(settings.calculation.wavelengths[405].fit_model, "quadratic")
        with self.assertRaises(TypeError):
            settings.calculation.wavelengths[405] = CalibrationSettings("cubic", response_scale="raw")
        with self.assertRaises(FrozenInstanceError):
            settings.calculation.wavelengths[405].weighting_exponent = 1.0
        with self.assertRaisesRegex(ValueError, 'Geen kalibratie-instellingen'):
            settings.calculation_for_wavelength(620)

    def test_calibration_and_validation_must_use_the_skax_wavelengths(self):
        from elisaflow.assays.protocol_samples import ProtocolSamples

        config = deepcopy(self.config)
        with self.assertRaisesRegex(ValueError, 'kalibratiegolflengten moeten overeenkomen'):
            ProtocolSamples((), {}, (450, 620)).validate_measurements(AssaySettings.from_dict("test.skax", config))
        config["validation"] = {"wavelengths": {"620": {}}}
        with self.assertRaisesRegex(ValueError, 'validatiegolflengten'):
            ProtocolSamples((), {}, (405, 450)).validate_measurements(AssaySettings.from_dict("test.skax", config))

    def test_rejects_missing_ambiguous_or_malformed_wavelength_calibrations(self):
        fit = {"response_scale": "raw", "fit_model": "4pl", "weighting_exponent": 0.0}
        malformed = [
            None, [], {}, {"450": fit}, {"405": fit, "450": fit, "620": fit},
            {"450": fit, "0450": fit}, {"450": fit, "0": fit}, {"450": fit, "bad": fit},
            {"450": fit, True: fit}, {"450": fit, 405.5: fit},
            {"450": fit, "405": []}, {"450": fit, "405": {}},
            {"450": fit, "405": {"response_scale": "raw", "fit_model": "6pl"}},
            {"450": fit, "405": {"response_scale": "raw", "fit_model": "4pl", "weighting_exponent": "automatic"}},
            {"450": fit, "405": fit | {"weighting_exponent": 2, "parameters": {"A": 2, "B": -1, "C": 1, "D": 0}}},
            {"450": fit, "405": fit | {"concentration_range": [1, 300]}},
        ]
        for curves in malformed:
            config = deepcopy(self.config)
            config["calculation"]["wavelengths"] = curves
            with self.subTest(curves=curves), self.assertRaises(ValueError):
                AssaySettings.from_dict("test.skax", config)
        for name, value in (
            ("fit_model", "4pl"), ("weighting_exponent", "1"), ("parameters", {}), ("response_scale", "raw"),
            ("preprocessing", {"expression": "sample * 2"}),
        ):
            config = deepcopy(self.config)
            config["calculation"][name] = value
            with self.subTest(field=name), self.assertRaisesRegex(ValueError, 'Onbekende velden'):
                AssaySettings.from_dict("test.skax", config)

    def test_switching_requires_peer_calibrations_and_preferred_wavelength(self):
        config = deepcopy(self.config)
        config["calculation"] = {"response_scale": "raw", "method": "calibration", "fit_model": "4pl"}
        with self.assertRaisesRegex(ValueError, 'vereist calculation.wavelengths'):
            AssaySettings.from_dict("test.skax", config)
        config = deepcopy(self.config)
        config["measurement"]["preferred_wavelength_nm"] = 620
        with self.assertRaisesRegex(ValueError, '[Vv]oorkeursgolflengte'):
            AssaySettings.from_dict("test.skax", config)
        config["measurement"] = {"mode": "single"}
        with self.assertRaisesRegex(ValueError, 'vereist meetmethode'):
            AssaySettings.from_dict("test.skax", config)

    def test_direct_construction_validates_typed_curves_and_shared_options(self):
        curves = {405: CalibrationSettings("quadratic", response_scale="raw"), 450: CalibrationSettings("4pl", response_scale="raw")}
        with self.assertRaises(ValueError):
            WavelengthCalibrationSettings(curves, concentration_range=(2, 1))
        with self.assertRaises(TypeError):
            WavelengthCalibrationSettings(curves, preprocessing=Formula("sample * 2"))
        for curve in ({}, CalibrationSettings("4pl", concentration_range=(1, 300), response_scale="raw")):
            with self.subTest(curve=curve), self.assertRaises(ValueError):
                WavelengthCalibrationSettings(curves | {405: curve})

    def test_fixed_wavelength_curves_round_trip_with_independent_response_scales(self):
        for scale in ("raw", "blank_corrected"):
            with self.subTest(scale=scale):
                config = deepcopy(self.config)
                config["calculation"]["wavelengths"]["405"] = {
                    "fit_model": "quadratic", "parameters": {"A": 0.03, "B": 0.003, "C": 0.000005}, "response_scale": scale,
                    "preprocessing": {"expression": "sample * 2"},
                }
                config["calculation"]["wavelengths"]["450"] = {"fit_model": "4pl", "response_scale": "raw", "weighting_exponent": 2.0}
                settings = AssaySettings.from_dict("test.skax", config)
                curve = settings.calculation_for_wavelength(405)
                self.assertEqual(curve.parameters.C, 0.000005)
                self.assertEqual(curve.concentration_range, (1, 300))
                self.assertEqual(curve.preprocessing, Formula("sample * 2"))
                self.assertEqual(curve.response_scale, scale)
                self.assertEqual(curve.weighting_exponent, 0)
                self.assertEqual(settings.calculation_for_wavelength(450).response_scale, "raw")
                self.assertIsNone(settings.calculation_for_wavelength(450).parameters)
                self.assertIsNone(settings.calculation_for_wavelength(450).preprocessing)
                serialized = settings.to_dict()
                self.assertNotIn("concentration_range", serialized["calculation"]["wavelengths"]["405"])
                self.assertEqual(AssaySettings.from_dict("test.skax", serialized), settings)

    def test_fixed_polynomial_wavelength_requires_the_shared_concentration_range(self):
        config = deepcopy(self.config)
        config["calculation"].pop("concentration_range")
        config["calculation"]["wavelengths"]["405"] = {"response_scale": "blank_corrected", "fit_model": "quadratic", "parameters": {"A": 0.03, "B": 0.003, "C": 0.000005}}
        with self.assertRaisesRegex(ValueError, "concentratiebereik"):
            AssaySettings.from_dict("test.skax", config)

    def test_wavelength_preprocessing_is_optional_and_never_inherited(self):
        config = deepcopy(self.config)
        config["calculation"]["wavelengths"]["450"].pop("preprocessing")
        settings = AssaySettings.from_dict("test.skax", config)
        self.assertEqual(settings.calculation_for_wavelength(405).preprocessing, Formula("sample * 2"))
        self.assertIsNone(settings.calculation_for_wavelength(450).preprocessing)
        self.assertNotIn("preprocessing", settings.to_dict()["calculation"])
        self.assertEqual(AssaySettings.from_dict("test.skax", settings.to_dict()), settings)

    def test_invalid_preprocessing_is_rejected_in_each_wavelength(self):
        for wavelength in ("405", "450"):
            for formula in (None, "sample", [], {"expression": "sample ** 2"}, {"expression": "sample / missing"}):
                with self.subTest(wavelength=wavelength, formula=formula), self.assertRaises(ValueError):
                    config = deepcopy(self.config)
                    config["calculation"]["wavelengths"][wavelength]["preprocessing"] = formula
                    AssaySettings.from_dict("test.skax", config)

    def test_each_fitted_wavelength_rejects_preprocessing(self):
        for wavelength in ("405", "450"):
            config = deepcopy(self.config)
            config["calculation"]["wavelengths"][wavelength].pop("parameters")
            with self.subTest(wavelength=wavelength), self.assertRaisesRegex(ValueError, "Voorbewerking.*vaste parameters"):
                AssaySettings.from_dict("test.skax", config)

    def test_direct_fixed_curves_resolve_and_round_trip_the_shared_range(self):
        settings = WavelengthCalibrationSettings({405: CalibrationSettings("4pl", FourPLParameters(2, -1, 1, 0), response_scale="blank_corrected"),
                                                 450: CalibrationSettings("4pl", response_scale="raw")}, concentration_range=(1, 300))
        self.assertEqual(WavelengthCalibrationSettings.from_dict(settings.to_dict()), settings)


if __name__ == "__main__":
    unittest.main()
