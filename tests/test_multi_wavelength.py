import unittest
from dataclasses import replace
from unittest.mock import patch

import numpy as np
import pandas as pd
from fixtures import assay_metadata, samples_from_records
from table_helpers import sample_table, well_table

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.assays.protocol_samples import ProtocolSamples
from elisaflow.calibration.logistic import FourPLModel, FourPLParameters
from elisaflow.calibration.polynomial import QuadraticModel, QuadraticParameters
from elisaflow.configuration.assay_settings import AssaySettings
from elisaflow.models import AssayData, ResultStatus, ValidationStatus


def dual_settings(**changes):
    config = {
        "lot_number": "test",
        "measurement": {"mode": "switch_on_concentration", "cutoff_concentration": 1000.0, "preferred_wavelength_nm": 450},
        "validation": {
            "result": {"controls": {"QC": {"result": {"min": 38, "max": 42}}}},
            "wavelengths": {
                "450": {
                    "blank": {"extinction": {"response_scale": "raw", "max": 0.02}},
                    "controls": {"QC": {"extinction": {"response_scale": "blank_corrected", "min": 2.0}}},
                    "standards": {"300": {"extinction": {"response_scale": "blank_corrected", "min": 1.0}}},
                },
                "405": {"blank": {"extinction": {"response_scale": "raw", "max": 0.01}}},
            },
        },
        "calculation": {
            "method": "calibration",
            "wavelengths": {
                "405": {"response_scale": "raw", "fit_model": "4pl", "weighting_exponent": 2.0},
                "450": {"response_scale": "raw", "fit_model": "4pl", "weighting_exponent": 2.0},
            },
            "concentration_range": [1, 300],
        },
        "result": {"label": "Conc", "unit": "U/mL"},
    }
    config.update(changes)
    return AssaySettings.from_dict("test.skax", config)


def dual_data():
    rows = []
    for nominal in (1, 3, 10, 30, 100, 300):
        for _ in range(2):
            rows.append((f"std-{nominal}", f"Standard {nominal}", "Standard", nominal, 0))
    for identity, name, kind, nominal, dilution in (
        ("p1", "Patient", "Unknown", 20, 1),
        ("p2", "Patient", "Unknown", 2, 1),
        ("qc", "QC", "Control", 40, 0),
        ("blank", "Blank", "Blank", 0, 0),
    ):
        rows.extend([(identity, name, kind, nominal, dilution)] * 2)
    template = pd.DataFrame(rows, columns=["sample_id", "sample_name", "sample_type", "sample_concentration", "sample_dilution"])
    template["well"] = [f"{chr(65 + index % 8)}{index // 8 + 1}" for index in range(len(template))]
    template["plate_name"] = "Plate 1"
    template["saturated"] = False
    template["disabled"] = False
    frames = []
    for wavelength, top, baseline, blank in ((450, 4.5, 0.08, 0.01), (405, 1.4, 0.025, 0.005)):
        frame = template.copy()
        frame["wavelength_nm"] = wavelength
        frame["source_extinction"] = FourPLModel().calculate_extinctions(
            frame["sample_concentration"],
            FourPLParameters(top, -1.2, 30, baseline),
        )
        frame.loc[frame["sample_type"].eq("Blank"), "source_extinction"] = blank
        frames.append(frame)
    return AssayData(samples_from_records(pd.concat(frames, ignore_index=True).to_dict('records')), replace(assay_metadata(), result_step_name='450 nm / 405 nm'))


def saturate(data, wavelength, identity, *, all_replicates=False):
    changed = False
    samples = []
    for sample in data.samples:
        wells = []
        for well in sample.measurements:
            if sample.sample_id == identity and well.wavelength_nm == wavelength and (all_replicates or not changed):
                well = replace(well, saturated=True)
                changed = True
            wells.append(well)
        samples.append(replace(sample, measurements=tuple(wells)))
    return replace(data, samples=tuple(samples))


def mixed_calibration():
    data = dual_data()
    mask = well_table(data.samples)["wavelength_nm"].eq(405) & well_table(data.samples)["sample_type"].ne("Blank")
    _readings = well_table(data.samples)
    _readings.loc[mask, "source_extinction"] = QuadraticModel().calculate_extinctions(
        _readings.loc[mask, "sample_concentration"], QuadraticParameters(0.03, 0.003, 0.000005)
    )
    data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
    config = dual_settings().to_dict()
    config["calculation"]["wavelengths"]["405"] = {"response_scale": "raw", "fit_model": "quadratic", "weighting_exponent": 1.5}
    return data, AssaySettings.from_dict("test.skax", config)


class MultiWavelengthAnalysisTests(unittest.TestCase):
    def analyze(self, data=None, settings=None):
        return AssayAnalyzer(data or dual_data(), settings or dual_settings()).analyze()

    def test_fixed_curve_preprocessing_keeps_independent_fitting_and_wavelength_selection(self):
        for fixed_fallback in (False, True):
            for preprocess_primary in (False, True):
                with self.subTest(fixed_fallback=fixed_fallback, preprocess_primary=preprocess_primary):
                    config = dual_settings().to_dict()
                    config["measurement"]["cutoff_concentration"] = 10
                    curves = config["calculation"]["wavelengths"]
                    if fixed_fallback:
                        curves["405"] = {
                            "fit_model": "4pl", "response_scale": "raw",
                            "parameters": {"A": 2.8, "B": -1.2, "C": 30, "D": 0.05},
                            "preprocessing": {"expression": "sample * 2"},
                        }
                    if preprocess_primary:
                        control_mean = float(FourPLModel().calculate_extinctions(40, FourPLParameters(4.5, -1.2, 30, 0.08)))
                        curves["450"] = {
                            "fit_model": "4pl", "response_scale": "raw",
                            "parameters": {"A": 4.5 / control_mean, "B": -1.2, "C": 30, "D": 0.08 / control_mean},
                            "preprocessing": {"expression": "sample / qc", "references": {"qc": "QC"}},
                        }
                    result = self.analyze(settings=AssaySettings.from_dict("test.skax", config))
                    patients = sample_table(result.samples).set_index("sample_id")
                    self.assertEqual(patients.loc["p1", "selected_wavelength_nm"], 405)
                    self.assertEqual(patients.loc["p2", "selected_wavelength_nm"], 450)
                    np.testing.assert_allclose(patients.loc[["p1", "p2"], "result_mean"], [20, 2], rtol=1e-6)
                    for wavelength, top, baseline, blank in ((405, 1.4, 0.025, 0.005), (450, 4.5, 0.08, 0.01)):
                        measurement = result.measurement_results[wavelength]
                        wells = well_table(measurement.samples)
                        readings = wells.loc[wells["sample_type"].ne("Blank")]
                        control_mean = readings.loc[readings["sample_type"].eq("Control"), "source_extinction"].mean()
                        factor = (2 if fixed_fallback else 1) if wavelength == 405 else 1 / control_mean if preprocess_primary else 1
                        np.testing.assert_allclose(readings["curve_extinction"], readings["source_extinction"] * factor)
                        np.testing.assert_allclose(readings.loc[readings["sample_type"].eq("Control"), "result"], 40, rtol=1e-6)
                        np.testing.assert_allclose(readings["extinction"], readings["source_extinction"] - blank)
                        self.assertAlmostEqual(measurement.calculation_details.curve.parameters.A, top * factor, places=5)
                        self.assertAlmostEqual(measurement.calculation_details.curve.parameters.D, baseline * factor, places=5)
                        self.assertEqual(measurement.calculation_details.curve.was_fitted, not (fixed_fallback if wavelength == 405 else preprocess_primary))

    def test_fixed_and_fitted_curves_use_independent_parameters_and_response_scales(self):
        for scale in ("raw", "blank_corrected"):
            with self.subTest(scale=scale):
                config = dual_settings().to_dict()
                correction = 0.005 if scale == "blank_corrected" else 0
                config["calculation"]["wavelengths"]["405"] = {
                    "fit_model": "4pl", "parameters": {"A": 1.4 - correction, "B": -1.2, "C": 30, "D": 0.025 - correction},
                    "response_scale": scale,
                }
                config["measurement"]["cutoff_concentration"] = 10
                result = self.analyze(settings=AssaySettings.from_dict("test.skax", config))
                patients = sample_table(result.samples).set_index("sample_id")
                self.assertEqual(patients.loc["p1", "selected_wavelength_nm"], 405)
                self.assertEqual(patients.loc["p2", "selected_wavelength_nm"], 450)
                self.assertAlmostEqual(patients.loc["p1", "result_mean"], 20, places=5)
                self.assertAlmostEqual(patients.loc["p2", "result_mean"], 2, places=5)
                fixed = result.measurement_results[405].calculation_details
                self.assertFalse(fixed.curve.was_fitted)
                self.assertIsNone(fixed.curve.diagnostics)
                self.assertAlmostEqual(fixed.curve.parameters.A, 1.4 - correction)
                self.assertTrue(result.measurement_results[450].calculation_details.curve.was_fitted)

    def test_two_fixed_curves_can_switch_without_standards_or_running_a_fitter(self):
        config = dual_settings(validation={}).to_dict()
        config["measurement"]["cutoff_concentration"] = 10
        config["calculation"]["wavelengths"] = {
            str(wavelength): {"fit_model": "4pl", "response_scale": "raw", "parameters": {"A": top, "B": -1.2, "C": 30, "D": baseline}}
            for wavelength, top, baseline in ((450, 4.5, 0.08), (405, 1.4, 0.025))
        }
        data = dual_data()
        data = replace(data, samples=tuple(sample for sample in data.samples if sample.sample_type.value != "Standard"))
        with patch("elisaflow.calibration.curves.fit_curve") as fit:
            result = self.analyze(data, AssaySettings.from_dict("test.skax", config))
        fit.assert_not_called()
        patients = sample_table(result.samples).set_index("sample_id")
        for identity, wavelength, expected in (("p1", 405, 20), ("p2", 450, 2)):
            self.assertEqual(patients.loc[identity, "selected_wavelength_nm"], wavelength)
            self.assertAlmostEqual(patients.loc[identity, "result_mean"], expected, places=5)
        self.assertTrue(all(not measurement.calculation_details.curve.was_fitted for measurement in result.measurement_results.values()))

    def test_fixed_polynomial_uses_the_shared_range_when_selected_instead_of_a_fitted_curve(self):
        data, settings = mixed_calibration()
        config = settings.to_dict()
        config["calculation"]["wavelengths"]["405"] = {
            "fit_model": "quadratic", "parameters": {"A": 0.03, "B": 0.003, "C": 0.000005}, "response_scale": "raw",
        }
        result = self.analyze(saturate(data, 450, "p1"), AssaySettings.from_dict("test.skax", config))
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 405)
        self.assertAlmostEqual(patient["result_mean"], 20, places=5)
        details = result.measurement_results[405].calculation_details
        self.assertIsInstance(details.curve.parameters, QuadraticParameters)
        self.assertFalse(details.curve.was_fitted)
        self.assertEqual(details.curve.concentration_range, (1, 300))

    def test_patient_filter_runs_per_wavelength_and_preserves_saturation_switching(self):
        data = dual_data()
        extra = well_table(data.samples).loc[well_table(data.samples)["sample_id"].eq("p1")].groupby("wavelength_nm").head(1).copy()
        extra["well"] = "H12"
        blank = extra["wavelength_nm"].map({450: 0.01, 405: 0.005})
        extra["source_extinction"] = (extra["source_extinction"] - blank) * 1.3 + blank
        data = replace(data, samples=samples_from_records(pd.concat([well_table(data.samples), extra], ignore_index=True).to_dict("records")))
        _readings = well_table(data.samples)
        _readings.loc[_readings["sample_id"].eq("p1"), "sample_dilution"] = 10
        data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
        data = saturate(data, 450, "p1")
        settings = dual_settings(patients={"replicate_filter": {"response_scale": "blank_corrected", "max_deviation_percent": 15, "max_exclusions": 1}})
        result = self.analyze(data, settings)
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 405)
        self.assertEqual(patient["count"], 2)
        self.assertAlmostEqual(patient["result_mean"], 200, places=4)
        fallback = well_table(result.measurement_results[405].samples)
        excluded = fallback.loc[fallback["well"].eq("H12")].iloc[0]
        self.assertFalse(excluded["included_in_analysis"])
        self.assertGreater(excluded["replicate_deviation_percent"], 15)
        self.assertTrue(pd.isna(excluded["result"]))

    def test_models_weights_and_diagnostics_are_resolved_per_wavelength(self):
        data, settings = mixed_calibration()
        data = saturate(data, 450, "p1")
        data = saturate(data, 450, "std-300", all_replicates=True)
        result = self.analyze(data, settings)
        for wavelength, model, parameter_type, power in (
            (405, "quadratic", QuadraticParameters, 1.5), (450, "4pl", FourPLParameters, 2)
        ):
            with self.subTest(wavelength=wavelength):
                measurement = result.measurement_results[wavelength]
                details = measurement.calculation_details
                self.assertEqual(details.calibration_settings.fit_model, model)
                self.assertIsInstance(details.curve.parameters, parameter_type)
                self.assertAlmostEqual(details.curve_r_squared, 1.0)
                wells = well_table(measurement.samples)
                standards = wells[wells["sample_type"].eq("Standard") & wells["included_in_analysis"]]
                np.testing.assert_allclose(details.curve.diagnostics.weights, 1 / standards["source_extinction"] ** power)
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 405)
        self.assertAlmostEqual(patient["result_mean"], 20, places=5)
        self.assertEqual(result.measurement_results[450].calculation_details.curve.concentration_range, (1, 100))
        self.assertEqual(result.measurement_results[405].calculation_details.curve.concentration_range, (1, 300))
        self.assertIsNone(result.calculation_details.curve)

    def test_configured_calibration_wavelengths_must_match_measurements(self):
        config = dual_settings().to_dict()
        curves = config["calculation"]["wavelengths"]
        curves["620"] = curves.pop("405")
        with self.assertRaisesRegex(ValueError, 'kalibratiegolflengten moeten overeenkomen'):
            self.analyze(settings=AssaySettings.from_dict("test.skax", config))

    def test_one_saturated_duplicate_switches_both_without_mutating_input(self):
        data = dual_data()
        data = saturate(data, 450, "p1")
        original = well_table(data.samples).copy(deep=True)
        result = self.analyze(data)
        samples = sample_table(result.samples).set_index("sample_id")
        self.assertEqual(samples.loc["p1", "selected_wavelength_nm"], 405)
        self.assertEqual(samples.loc["p1", "count"], 2)
        self.assertAlmostEqual(samples.loc["p1", "result_mean"], 20, places=5)
        self.assertEqual(samples.loc["p2", "selected_wavelength_nm"], 450)
        self.assertAlmostEqual(samples.loc["p2", "result_mean"], 2, places=5)
        self.assertEqual(result.validation_status, ValidationStatus.PASSED)
        selected = well_table(result.samples)[well_table(result.samples)["sample_id"].eq("p1")]
        self.assertEqual(selected["selected_wavelength_nm"].tolist(), [405, 405])
        self.assertTrue(selected["included_in_analysis"].all())
        pd.testing.assert_frame_equal(well_table(data.samples), original)

    def test_saturation_is_not_a_global_well_exclusion(self):
        data = dual_data()
        data = saturate(data, 405, "p1", all_replicates=True)
        result = self.analyze(data)
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 450)
        self.assertAlmostEqual(patient["result_mean"], 20, places=5)

    def test_saturation_in_both_channels_does_not_report_a_surviving_duplicate(self):
        data = dual_data()
        data = saturate(data, 450, "p1")
        indices = well_table(data.samples).index[well_table(data.samples)["wavelength_nm"].eq(405) & well_table(data.samples)["sample_id"].eq("p1")]
        _readings = well_table(data.samples)
        _readings.loc[indices[-1], "saturated"] = True
        data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
        result = self.analyze(data)
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertTrue(np.isnan(patient["result_mean"]))
        self.assertEqual(patient["result_status"], ResultStatus.UNAVAILABLE.value)

    def test_missing_fallback_duplicate_keeps_result_unavailable(self):
        data = dual_data()
        data = saturate(data, 450, "p1")
        index = well_table(data.samples).index[well_table(data.samples)["wavelength_nm"].eq(405) & well_table(data.samples)["sample_id"].eq("p1")][0]
        _readings = well_table(data.samples)
        _readings.loc[index, "source_extinction"] = np.nan
        data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
        patient = sample_table(self.analyze(data).samples).set_index("sample_id").loc["p1"]
        self.assertTrue(np.isnan(patient["result_mean"]))

    def test_disabled_saturated_measurement_does_not_trigger_switching(self):
        data = dual_data()
        data = saturate(data, 450, "p1")
        _readings = well_table(data.samples)
        _readings.loc[_readings["saturated"], "disabled"] = True
        data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
        patient = sample_table(self.analyze(data).samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 450)
        self.assertEqual(patient["count"], 1)

    def test_selection_is_independent_of_dilution_correction(self):
        data = dual_data()
        data = saturate(data, 450, "p1")
        _readings = well_table(data.samples)
        _readings.loc[_readings["sample_id"].eq("p1"), "sample_dilution"] = 10
        data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
        patient = sample_table(self.analyze(data).samples).set_index("sample_id").loc["p1"]
        self.assertEqual(patient["selected_wavelength_nm"], 405)
        self.assertAlmostEqual(patient["result_mean"], 200, places=4)

    def test_same_sample_at_different_dilutions_selects_independently(self):
        data = dual_data()
        data = saturate(data, 450, "p1")
        for wavelength in (450, 405):
            index = well_table(data.samples).index[well_table(data.samples)["wavelength_nm"].eq(wavelength) & well_table(data.samples)["sample_id"].eq("p1")][-1]
            _readings = well_table(data.samples)
            _readings.loc[index, "sample_dilution"] = 10
            data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
        samples = sample_table(self.analyze(data).samples)
        patient = samples[samples["sample_id"].eq("p1")].set_index("sample_dilution")
        self.assertEqual(patient.loc[1, "selected_wavelength_nm"], 405)
        self.assertEqual(patient.loc[10, "selected_wavelength_nm"], 450)
        self.assertAlmostEqual(patient.loc[10, "result_mean"], 200, places=4)

    def test_preferred_wavelength_can_be_405(self):
        config = dual_settings().to_dict()
        config["measurement"]["preferred_wavelength_nm"] = 405
        data = dual_data()
        data = saturate(data, 405, "p1")
        result = self.analyze(data, AssaySettings.from_dict("test.skax", config))
        samples = sample_table(result.samples).set_index("sample_id")
        self.assertEqual(samples.loc["p1", "selected_wavelength_nm"], 450)
        self.assertEqual(samples.loc["p2", "selected_wavelength_nm"], 405)

    def test_every_usable_standard_is_fitted_independently_and_ranges_follow_retained_levels(self):
        data = dual_data()
        data = saturate(data, 450, "std-300", all_replicates=True)
        result = self.analyze(data)
        low, high = (result.measurement_results[value].calculation_details for value in (450, 405))
        self.assertEqual(len(low.curve.diagnostics.weights), 10)
        self.assertEqual(len(high.curve.diagnostics.weights), 12)
        self.assertEqual(low.curve.concentration_range, (1, 100))
        self.assertEqual(high.curve.concentration_range, (1, 300))
        rule = next(check for check in result.validation_results if check.sample_name == "Standaard 300")
        self.assertEqual(rule.status, ValidationStatus.NOT_EVALUATED)
        self.assertEqual(rule.wavelength_nm, 450)

    def test_insufficient_preferred_standard_levels_leave_the_other_curve_available(self):
        data = dual_data()
        for identity in ("std-30", "std-100", "std-300"):
            data = saturate(data, 450, identity, all_replicates=True)
        result = self.analyze(data)
        self.assertIn('ten minste 4 verschillende', result.measurement_results[450].calculation_details.calibration_error)
        patients = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")]
        self.assertTrue(patients["selected_wavelength_nm"].eq(405).all())
        self.assertTrue(np.isfinite(patients["result_mean"]).all())

    def test_unflagged_out_of_range_preferred_signal_switches(self):
        data = dual_data()
        data = saturate(data, 450, "std-300", all_replicates=True)
        for wavelength, top, baseline in ((450, 4.5, 0.08), (405, 1.4, 0.025)):
            mask = well_table(data.samples)["wavelength_nm"].eq(wavelength) & well_table(data.samples)["sample_id"].eq("p1")
            _readings = well_table(data.samples)
            _readings.loc[mask, "source_extinction"] = FourPLModel().calculate_extinctions(
                150,
                FourPLParameters(top, -1.2, 30, baseline),
            )
            data = replace(data, samples=samples_from_records(_readings.to_dict("records")))
        config = dual_settings().to_dict()
        config["validation"]["wavelengths"]["450"].pop("standards")
        result = self.analyze(data, AssaySettings.from_dict("test.skax", config))
        patient = sample_table(result.samples).set_index("sample_id").loc["p1"]
        self.assertTrue(result.measurement_results[450].usable_for_results)
        self.assertEqual(patient["selected_wavelength_nm"], 405)
        self.assertEqual(patient["selection_reason"], "preferred_out_of_range")
        self.assertAlmostEqual(patient["result_mean"], 150, places=4)

    def test_scoped_qc_checks_each_wavelength_and_shared_control_checks_selected_result(self):
        data = dual_data()
        data = saturate(data, 450, "qc", all_replicates=True)
        result = self.analyze(data)
        controls = [check for check in result.validation_results if check.sample_name == "QC"]
        shared = next(check for check in controls if check.wavelength_nm is None)
        scoped = next(check for check in controls if check.wavelength_nm == 450)
        self.assertEqual(shared.status, ValidationStatus.PASSED)
        self.assertEqual(scoped.status, ValidationStatus.NOT_EVALUATED)
        self.assertEqual(result.validation_status, ValidationStatus.NOT_EVALUATED)
        patients = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")]
        self.assertTrue(patients["selected_wavelength_nm"].eq(450).all())
        blank_rules = [check for check in result.validation_results if check.sample_name == "Blanco"]
        self.assertEqual({check.wavelength_nm: check.value for check in blank_rules}, {450: 0.01, 405: 0.005})

    def test_unknown_validation_wavelength_and_single_channel_input_are_rejected(self):
        config = dual_settings().to_dict()
        config["validation"]["wavelengths"]["620"] = config["validation"]["wavelengths"].pop("405")
        with self.assertRaisesRegex(ValueError, "620"):
            self.analyze(settings=AssaySettings.from_dict("test.skax", config))
        data = dual_data()
        with self.assertRaisesRegex(ValueError, 'precies 2 analysegolflengte'):
            self.analyze(AssayData(samples_from_records(well_table(data.samples)[well_table(data.samples)['wavelength_nm'].eq(450)].to_dict('records')), data.metadata))


class MultiWavelengthConfigTests(unittest.TestCase):
    def test_round_trip_keeps_all_wavelength_rules(self):
        settings = dual_settings()
        self.assertEqual(AssaySettings.from_dict(settings.protocol_name, settings.to_dict()), settings)
        with self.assertRaises(TypeError):
            settings.validation.wavelengths[450] = None

    def test_rejects_ambiguous_or_unsupported_settings(self):
        for value in (True, 0, -450, 450.5, "450"):
            config = dual_settings().to_dict()
            config["measurement"]["preferred_wavelength_nm"] = value
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'positief.*geheel'):
                AssaySettings.from_dict("test.skax", config)
        config = dual_settings().to_dict()
        config["validation"]["result"]["controls"]["QC"]["extinction"] = {"min": 1, "response_scale": "blank_corrected"}
        with self.assertRaisesRegex(ValueError, 'extinctiecriteria in onder validation.wavelengths'):
            AssaySettings.from_dict("test.skax", config)
        config = dual_settings().to_dict()
        config["validation"]["wavelengths"]["405"]["blank"]["preprocessing"] = [{"operation": "multiply", "value": 2}]
        with self.assertRaisesRegex(ValueError, 'Onbekende velden in Validatiecriteria'):
            AssaySettings.from_dict("test.skax", config)

    def test_protocol_checks_include_wavelength_only_targets(self):
        protocol = ProtocolSamples(("QC",), {1.0: ("Standaard 1",)}, (405, 450))
        with self.assertRaisesRegex(ValueError, "300"):
            protocol.validate_samples(dual_settings())




if __name__ == "__main__":
    unittest.main()
