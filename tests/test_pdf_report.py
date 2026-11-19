import unittest
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pandas as pd
import pymupdf
from fixtures import assay_metadata, reference_parameters, samples_from_records, samples_with_statistics, with_validation
from reportlab.lib import colors
from reportlab.platypus import BaseDocTemplate, Table
from table_helpers import sample_table, well_table
from test_multi_wavelength import dual_data, dual_settings, mixed_calibration, saturate

from elisaflow.analysis.assay import AssayAnalyzer
from elisaflow.assays.registry import AssayRegistry
from elisaflow.calibration.logistic import FivePLModel, FivePLParameters, FourPLModel
from elisaflow.calibration.plot import plot_calibration
from elisaflow.calibration.polynomial import CubicModel, CubicParameters, QuadraticModel, QuadraticParameters
from elisaflow.configuration.assay_settings import (
    AssaySettings,
    CalibrationSettings,
    FormulaSettings,
    MeasurementSettings,
    ReplicateFilterSettings,
    ReportSettings,
    SampleSettings,
)
from elisaflow.formulas.expressions import Formula
from elisaflow.models import AssayData, CalculationMethod, Sample, SampleType, ValidationStatus, WellMeasurement
from elisaflow.readers.skax_reader import SkaxReader
from elisaflow.reporting.calibration import calibration_section, calibration_table
from elisaflow.reporting.patients import patient_pages
from elisaflow.reporting.pdf_report import PdfReportWriter, _format_session_time
from elisaflow.reporting.quality_control import validation_table
from elisaflow.sample_selection import select_samples


def quality_control_table(writer: PdfReportWriter) -> Table | None:
    """Exercise the production QC renderer with an explicit table width."""
    return validation_table(writer.context, writer.content_width * 0.60)


class SessionTimeFormattingTests(unittest.TestCase):
    def test_measurement_times_use_dutch_winter_and_summer_time(self):
        for source, expected in (
            ("2026-11-03T20:39:20.0429605Z", "03-11-2026 21:39:20"),
            ("2026-07-03T20:39:20Z", "03-07-2026 22:39:20"),
            ("2026-11-03T15:39:20-05:00", "03-11-2026 21:39:20"),
            ("2026-11-03T21:39:20+01:00", "03-11-2026 21:39:20"),
            ("2026-07-03T22:39:20+02:00", "03-07-2026 22:39:20"),
            ("2026-11-03T23:30:00Z", "04-11-2026 00:30:00"),
        ):
            with self.subTest(source=source):
                self.assertEqual(_format_session_time(source), expected)

    def test_daylight_saving_transitions_follow_amsterdam_rules(self):
        for source, expected in (
            ("2026-03-29T00:59:59Z", "29-03-2026 01:59:59"),
            ("2026-03-29T01:00:00Z", "29-03-2026 03:00:00"),
            ("2026-10-25T00:59:59Z", "25-10-2026 02:59:59"),
            ("2026-10-25T01:00:00Z", "25-10-2026 02:00:00"),
        ):
            with self.subTest(source=source):
                self.assertEqual(_format_session_time(source), expected)

    def test_timestamps_without_an_offset_are_rejected_instead_of_using_the_computers_timezone(self):
        with self.assertRaisesRegex(ValueError, "UTC-offset"):
            _format_session_time("2026-11-03T20:39:20")


class PdfReportTests(unittest.TestCase):
    @staticmethod
    def patient_table(writer):
        pages = patient_pages(writer.context, writer.report_area, on_page=writer._draw_footer)
        return next(flowable for flowable in pages.flowables if isinstance(flowable, Table))

    def test_negative_mean_cv_is_reported_as_positive_with_warning(self):
        frame = pd.DataFrame({
            "sample_id": ["blank", "patient", "patient"],
            "sample_name": ["Blank", "Patient", "Patient"],
            "sample_type": ["Blank", "Unknown", "Unknown"],
            "source_extinction": [1, 0.1, 0.9],
            "well": ["A1", "A2", "A3"],
        })
        settings = AssaySettings("test.skax", "lot", calculation=FormulaSettings(Formula("sample"), response_scale="blank_corrected"))
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict("records")), assay_metadata()), settings).analyze()
        with TemporaryDirectory() as directory:
            path = PdfReportWriter(result, analyst="operator").write(Path(directory) / "report.pdf")
            with pymupdf.open(path) as document:
                spans = [
                    span for page in document for block in page.get_text("dict")["blocks"]
                    for line in block.get("lines", []) for span in line["spans"]
                ]
        cv_spans = [span for span in spans if span["text"].strip() == "113.1"]
        self.assertEqual(len(cv_spans), 2)
        self.assertTrue(all(span["color"] == 0xB22222 for span in cv_spans))

    def test_report_displays_dutch_measurement_time_without_a_timezone_label(self):
        result = self.make_result(CalculationMethod.FORMULA)
        for source, expected in (
            ("2026-11-03T20:39:20Z", "03-11-2026 21:39:20"),
            ("2026-07-03T20:39:20Z", "03-07-2026 22:39:20"),
        ):
            with self.subTest(source=source), TemporaryDirectory() as directory:
                dated = replace(result, metadata=replace(result.metadata, session_time=source))
                path = PdfReportWriter(dated, analyst="operator").write(Path(directory) / "report.pdf")
                with pymupdf.open(path) as document:
                    text = " ".join(document[0].get_text().split())
                self.assertIn(f"Datum/tijd: {expected} Analist:", text)

    def test_fixed_and_fitted_wavelength_curves_are_reported_with_their_own_settings(self) -> None:
        config = dual_settings().to_dict()
        config["calculation"]["wavelengths"]["405"] = {
            "fit_model": "4pl", "response_scale": "raw", "parameters": {"A": 1.4, "B": -1.2, "C": 30, "D": 0.025},
        }
        result = AssayAnalyzer(dual_data(), AssaySettings.from_dict("test.skax", config)).analyze()
        with TemporaryDirectory() as directory:
            with patch("elisaflow.reporting.calibration.plot_calibration", wraps=plot_calibration) as plot:
                path = PdfReportWriter(result, analyst="operator").write(Path(directory) / "mixed-curves.pdf")
            self.assertEqual(plot.call_count, 2)
            for call, measurement in zip(plot.call_args_list, result.measurement_results.values(), strict=True):
                self.assertEqual(call.args[2], measurement.calculation_details.curve.parameters)
                self.assertEqual(call.kwargs["curve_was_fitted"], measurement.calculation_details.curve.was_fitted)
            with pymupdf.open(path) as document:
                text = " ".join(" ".join(page.get_text().split()) for page in document)
            self.assertIn("Vaste kalibratiecurve (405 nm)", text)
            self.assertNotIn("A=1.4", text)
            self.assertNotIn("D=0.025", text)
            self.assertNotIn("Responsschaal:", text)
            self.assertIn("Kalibratiecurve (450 nm)", text)

    def test_report_is_published_only_after_a_successful_build(self) -> None:
        result = self.make_result(CalculationMethod.FORMULA)
        original_build = BaseDocTemplate.build
        previous_report = b"previous completed report"
        for existing in (False, True):
            with self.subTest(existing=existing), TemporaryDirectory() as directory:
                output = Path(directory) / "nested" / "report.pdf"
                if existing:
                    output.parent.mkdir()
                    output.write_bytes(previous_report)

                def build_document(document, story, *, filename) -> None:
                    temporary_path = Path(filename.name)
                    self.assertNotEqual(temporary_path, output)
                    self.assertEqual(temporary_path.parent, output.parent)
                    self.assertTrue(temporary_path.is_file())
                    original_build(document, story, filename=filename)
                    if existing:
                        self.assertEqual(output.read_bytes(), previous_report)
                    else:
                        self.assertFalse(output.exists())

                with patch("elisaflow.reporting.pdf_report.BaseDocTemplate.build", autospec=True, side_effect=build_document) as build:
                    path = PdfReportWriter(result, analyst="operator").write(output)
                build.assert_called_once()
                self.assertEqual(path, output)
                with pymupdf.open(path) as document:
                    self.assertGreater(document.page_count, 0)
                    self.assertIn("ElisaFlow analyserapport", document[0].get_text())
                self.assertEqual(list(output.parent.iterdir()), [output])

    def test_failed_or_interrupted_build_preserves_the_destination(self) -> None:
        result = self.make_result(CalculationMethod.FORMULA)
        previous_report = b"previous completed report"
        for error in (OSError("Disk full"), KeyboardInterrupt(), SystemExit()):
            for existing in (False, True):
                with self.subTest(error=type(error).__name__, existing=existing), TemporaryDirectory() as directory:
                    output = Path(directory) / "report.pdf"
                    if existing:
                        output.write_bytes(previous_report)

                    def fail_build(document, story, *, filename) -> None:
                        filename.write(b"%PDF-1.4\npartial report")
                        raise error

                    with (
                        patch("elisaflow.reporting.pdf_report.BaseDocTemplate.build", autospec=True, side_effect=fail_build),
                        self.assertRaises(type(error)) as raised,
                    ):
                        PdfReportWriter(result, analyst="operator").write(output)
                    self.assertIs(raised.exception, error)
                    if existing:
                        self.assertEqual(output.read_bytes(), previous_report)
                    else:
                        self.assertFalse(output.exists())
                    self.assertEqual(list(output.parent.iterdir()), [output] if existing else [])

    def test_pdf_serialization_failure_preserves_the_error_and_cleans_up(self) -> None:
        result = self.make_result(CalculationMethod.FORMULA)
        with TemporaryDirectory() as directory:
            output = Path(directory) / "report.pdf"
            output.write_bytes(b"previous completed report")
            error = ValueError("PDF serialization failed")
            with (
                patch("reportlab.pdfbase.pdfdoc.PDFDocument.GetPDFData", side_effect=error),
                self.assertRaises(ValueError) as raised,
            ):
                PdfReportWriter(result, analyst="operator").write(output)
            self.assertIs(raised.exception, error)
            self.assertEqual(output.read_bytes(), b"previous completed report")
            self.assertEqual(list(output.parent.iterdir()), [output])

    def test_report_preparation_failure_removes_the_temporary_file(self) -> None:
        result = self.make_result(CalculationMethod.FORMULA)
        for method in ("_initialize_report", "_write_title"):
            with self.subTest(method=method), TemporaryDirectory() as directory:
                output = Path(directory) / "report.pdf"
                output.write_bytes(b"previous completed report")
                writer = PdfReportWriter(result, analyst="operator")
                error = ValueError("Invalid report content")
                with patch.object(writer, method, side_effect=error), self.assertRaises(ValueError) as raised:
                    writer.write(output)
                self.assertIs(raised.exception, error)
                self.assertEqual(output.read_bytes(), b"previous completed report")
                self.assertEqual(list(output.parent.iterdir()), [output])

    def test_publication_failure_preserves_the_destination_and_removes_the_temporary_file(self) -> None:
        result = self.make_result(CalculationMethod.FORMULA)
        with TemporaryDirectory() as directory:
            output = Path(directory) / "report.pdf"
            output.write_bytes(b"previous completed report")
            error = PermissionError("Report is open in another application")
            with patch.object(Path, "replace", side_effect=error) as publish, self.assertRaises(PermissionError) as raised:
                PdfReportWriter(result, analyst="operator").write(output)
            publish.assert_called_once_with(output)
            self.assertIs(raised.exception, error)
            self.assertEqual(output.read_bytes(), b"previous completed report")
            self.assertEqual(list(output.parent.iterdir()), [output])

    def test_report_preserves_disabled_readings_when_no_patients_are_reportable(self):
        sample = Sample("p", "Patient", SampleType.UNKNOWN, (WellMeasurement(wavelength_nm=450, well="A1", source_extinction=1, disabled=True),))
        settings = AssaySettings.from_dict("test.skax", {"lot_number": "lot", "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample"}})
        result = AssayAnalyzer(AssayData((sample,), assay_metadata()), settings).analyze()
        self.assertFalse(result.samples[0].reportable)
        self.assertTrue(sample_table(result.samples).empty)
        with TemporaryDirectory() as directory:
            path = PdfReportWriter(result, analyst="operator").write(Path(directory) / "disabled.pdf")
            with pymupdf.open(path) as document:
                text = "\n".join(page.get_text() for page in document)
        self.assertIn("UIT", text)
        self.assertIn("Uitgesloten replicaten", text)
        self.assertNotIn("Patiëntresultaten", text)

    def test_failed_qc_remains_reported_without_marking_the_curve_unavailable(self):
        config = dual_settings().to_dict()
        config["validation"]["wavelengths"]["450"]["controls"]["QC"]["extinction"] = {"min": 50, "response_scale": "blank_corrected"}
        result = AssayAnalyzer(dual_data(), AssaySettings.from_dict("test.skax", config)).analyze()
        samples = sample_table(result.samples)
        self.assertTrue(samples.loc[samples["sample_type"].ne("Blank"), "selected_wavelength_nm"].eq(450).all())
        self.assertEqual(result.validation_status, ValidationStatus.FAILED)
        with TemporaryDirectory() as directory:
            path = PdfReportWriter(result, analyst="operator").write(Path(directory) / "failed_qc.pdf")
            with pymupdf.open(path) as document:
                text = " ".join(" ".join(page.get_text().split()) for page in document)
            self.assertIn("Ongeldig", text)
            self.assertIn("50", text)
            self.assertNotIn("Niet bruikbaar:", text)
            self.assertNotIn("QC van gebruikte golflengten", text)

    def test_failed_curve_fits_remain_visible_while_available_curves_are_reported(self):
        for failed in ((450,), (405,), (450, 405)):
            with self.subTest(failed=failed), TemporaryDirectory() as directory:
                data = dual_data()
                for wavelength in failed:
                    for identity in ("std-30", "std-100", "std-300"):
                        data = saturate(data, wavelength, identity, all_replicates=True)
                result = AssayAnalyzer(data, dual_settings()).analyze()
                with patch("elisaflow.reporting.calibration.plot_calibration", wraps=plot_calibration) as plot:
                    report = PdfReportWriter(result, analyst="operator").write(Path(directory) / "failed_curve.pdf")
                self.assertEqual(plot.call_count, 2 - len(failed))
                with pymupdf.open(report) as document:
                    text = " ".join(" ".join(page.get_text().split()) for page in document)
                    self.assertIn("Niet bruikbaar: kalibratie niet beschikbaar", text)
                    self.assertIn('Curvepassing vereist ten minste 4 verschillende concentraties.', text)
                    self.assertNotIn("QC van gebruikte golflengten", text)
                    for wavelength in failed:
                        self.assertIn(f"Curve {wavelength} nm", text)
                    expected_status = ValidationStatus.NOT_EVALUATED if 450 in failed else ValidationStatus.PASSED
                    self.assertEqual(result.validation_status, expected_status)

    def test_qc_value_cv_uses_stored_results_and_respects_result_cv_setting(self):
        base = self.make_result(CalculationMethod.FORMULA)
        settings = with_validation(base.settings, controls={"Positive": {"result": {"min": 1, "max": 2}}})
        result = AssayAnalyzer(AssayData(samples_from_records(well_table(base.samples).to_dict('records')), base.metadata), settings).analyze()
        for visible in (True, False):
            settings = replace(result.settings, report=replace(result.settings.report, show_result_cv=visible))
            writer = PdfReportWriter(replace(result, settings=settings), analyst="operator")
            with (
                TemporaryDirectory() as directory,
                patch("elisaflow.quality_control.summary.summarize_quality_control", side_effect=AssertionError("QC recalculated during report generation")),
            ):
                writer._initialize_report(Path(directory) / "report.pdf")
                table = quality_control_table(writer)
                headers = [cell.getPlainText() for cell in table._cellvalues[0]]
                index = headers.index("Waarde")
                row = self.quality_control_row(writer, "Positive")
                if visible:
                    self.assertEqual(table.repeatRows, 2)
                    subheaders = [cell.getPlainText() for cell in table._cellvalues[1]]
                    self.assertEqual(subheaders[index], "Gem.")
                    self.assertEqual(subheaders[index + 1], "CV%")
                    self.assertEqual(row[index + 1].getPlainText(), f"{result.validation_results[0].cv_percent:.1f}")
                else:
                    self.assertEqual(table.repeatRows, 1)
                    self.assertEqual(headers[index + 1], "Eenheid")

    def test_each_calibration_plot_uses_its_measurements_model(self):
        data, settings = mixed_calibration()
        data = saturate(data, 450, "p1")
        result = AssayAnalyzer(data, settings).analyze()
        with TemporaryDirectory() as directory, patch("elisaflow.reporting.calibration.plot_calibration", wraps=plot_calibration) as plot:
            path = PdfReportWriter(result, analyst="operator").write(Path(directory) / "mixed.pdf")
            self.assertEqual(plot.call_count, 2)
            for call, measurement in zip(plot.call_args_list, result.measurement_results.values(), strict=True):
                model_type = QuadraticModel if measurement.wavelength_nm == 405 else FourPLModel
                self.assertIsInstance(call.args[1], model_type)
                self.assertIs(call.args[1], measurement.calculation_details.curve.model)
                self.assertEqual(call.kwargs["model_name"], measurement.calculation_details.calibration_settings.fit_model)
                self.assertEqual(call.args[2], measurement.calculation_details.curve.parameters)
            with pymupdf.open(path) as document:
                text = " ".join(page.get_text() for page in document)
            self.assertIn("405 nm", text)
            self.assertIn("450 nm", text)

    def test_report_uses_executed_curves_and_settings_after_settings_are_replaced(self):
        config = dual_settings().to_dict()
        config["calculation"]["wavelengths"]["405"] = {
            "fit_model": "4pl", "response_scale": "raw",
            "parameters": {"A": 2.8, "B": -1.2, "C": 30, "D": 0.05}, "preprocessing": {"expression": "sample * 2"},
        }
        result = AssayAnalyzer(dual_data(), AssaySettings.from_dict("test.skax", config)).analyze()
        for wavelength in config["calculation"]["wavelengths"]:
            config["calculation"]["wavelengths"][wavelength] = {
                "response_scale": "blank_corrected",
                "fit_model": "5pl", "parameters": {"A": 8, "B": -1.2, "C": 30, "D": 0.05, "E": 0.8},
                "preprocessing": {"expression": "sample * 9"},
            }
        changed = replace(result, settings=AssaySettings.from_dict("test.skax", config))
        with (
            TemporaryDirectory() as directory,
            patch(
                "elisaflow.calibration.curves.create_calibration_curve", side_effect=AssertionError("Curve reconstructed during reporting"),
            ),
            patch("elisaflow.reporting.calibration.plot_calibration", wraps=plot_calibration) as plot,
        ):
            PdfReportWriter(changed, analyst="operator").write(Path(directory) / "executed.pdf")
        self.assertEqual(plot.call_count, 2)
        for call, measurement in zip(plot.call_args_list, result.measurement_results.values(), strict=True):
            details = measurement.calculation_details
            self.assertIs(call.args[1], details.curve.model)
            self.assertIs(call.args[2], details.curve.parameters)
            self.assertEqual(call.kwargs["model_name"], "4pl")
            self.assertEqual(call.kwargs["weighting_exponent"], details.calibration_settings.weighting_exponent)
            self.assertIs(call.kwargs["preprocessing"], details.calibration_settings.preprocessing)
            self.assertEqual(call.kwargs["r_squared"], details.curve_r_squared)
            self.assertEqual(call.kwargs["curve_was_fitted"], details.curve.was_fitted)
            self.assertEqual(call.kwargs["concentration_range"], None if details.curve.was_fitted else details.curve.concentration_range)

    def test_report_uses_executed_formula_after_formula_settings_are_replaced(self):
        result = self.make_result(CalculationMethod.FORMULA)
        changed = replace(
            result, settings=replace(result.settings, calculation=FormulaSettings(Formula("sample * 99"), response_scale="blank_corrected"))
        )
        with TemporaryDirectory() as directory:
            path = PdfReportWriter(changed, analyst="operator").write(Path(directory) / "report.pdf")
            with pymupdf.open(path) as document:
                text = "\n".join(page.get_text() for page in document)
        self.assertIn("Berekening: Ratio (Monster / Cutoff)", text)
        self.assertNotIn("Monster * 99", text)

    def test_ap50_report_uses_percentage_heading_and_hides_only_result_cv(self):
        project = Path(__file__).resolve().parents[1]
        settings = AssayRegistry(project / "configs").get("AP50.skax")
        frame = pd.DataFrame(
            {
                "sample_id": ["p", "p", "n", "n", "x", "x"],
                "sample_name": ["Positief", "Positief", "Negatief", "Negatief", "Patient", "Patient"],
                "sample_type": ["Control"] * 4 + ["Unknown"] * 2,
                "source_extinction": [1.4, 1.6, 0.1, 0.1, 0.7, 0.9],
                "well": ["A1", "B1", "C1", "D1", "E1", "F1"],
            }
        )
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), assay_metadata(name='AP50.skax')), settings).analyze()
        self.assertTrue(np.isfinite(sample_table(result.samples).iloc[-1]["result_cv_percent"]))
        with TemporaryDirectory() as directory:
            path = PdfReportWriter(result, analyst="operator").write(Path(directory) / "AP50.pdf")
            with pymupdf.open(path) as document:
                text = " ".join(" ".join(page.get_text().split()) for page in document)
            patients = text.split("Pati\u00ebntresultaten")[-1]
            self.assertIn("Score [%]", patients)
            self.assertEqual(patients.count("CV%"), 1)
            self.assertIn("Berekening: Score (100 * (Monster - Negatief) / (Positief - Negatief))", text)
            self.assertNotIn("positive:", text)
            self.assertNotIn("negative:", text)
            fields = [
                "Rungegevens",
                "Datum/tijd:",
                "Analist:",
                "Instrument:",
                "Assaygegevens",
                "Lotnummer kit:",
                "Golflengte:",
                "Berekening:",
            ]
            positions = [text.index(field) for field in fields]
            self.assertEqual(positions, sorted(positions))
            self.assertNotIn("Blanco extinctie", text)

    def test_standard_identity_is_independent_of_its_display_label(self) -> None:
        config = dual_settings().to_dict()
        config["validation"]["result"]["standards"] = {"300": {"result": {"min": 250}}}
        result = AssayAnalyzer(dual_data(), AssaySettings.from_dict("test.skax", config)).analyze()
        original_writer = PdfReportWriter(result, analyst="operator")
        with TemporaryDirectory() as directory:
            original_writer._initialize_report(Path(directory) / "original.pdf")
            original = quality_control_table(original_writer)._cellvalues
        renamed = replace(
            result,
            measurement_results={
                key: replace(
                    measurement,
                    validation_results=tuple(
                        replace(check, sample_name="Calibrator high") if check.sample_type == SampleType.STANDARD else check
                        for check in measurement.validation_results
                    ),
                )
                for key, measurement in result.measurement_results.items()
            },
            selected=replace(
                result.selected,
                validation_results=tuple(
                    replace(check, sample_name="Calibrator high") if check.sample_type == SampleType.STANDARD else check
                    for check in result.selected.validation_results
                ),
            ),
        )
        changed_writer = PdfReportWriter(renamed, analyst="operator")
        with TemporaryDirectory() as directory:
            changed_writer._initialize_report(Path(directory) / "renamed.pdf")
            changed = quality_control_table(changed_writer)._cellvalues
        original_rows = [[cell.getPlainText() for cell in row] for row in original]
        changed_rows = [[cell.getPlainText() for cell in row] for row in changed]
        for row in original_rows:
            if row[0] == "Standaard 300":
                row[0] = "Calibrator high"
        self.assertEqual(changed_rows, original_rows)

    def test_dual_quality_control_combines_wavelengths_and_keeps_rule_provenance(self) -> None:
        data = dual_data()
        data = saturate(data, 450, "qc")
        result = AssayAnalyzer(data, dual_settings()).analyze()
        # Reporting must use stored checks rather than altered measurement summaries.
        measurement_samples = sample_table(result.measurement_results[405].samples)
        measurement_samples["extinction_mean"] = 999
        with (
            TemporaryDirectory() as directory,
            patch("elisaflow.quality_control.summary.summarize_quality_control", side_effect=AssertionError("QC was recalculated")),
        ):
            writer = PdfReportWriter(result, analyst="operator")
            writer._initialize_report(Path(directory) / "report.pdf")
            table = quality_control_table(writer)
            cells = [[cell.getPlainText() for cell in row] for row in table._cellvalues]
            self.assertEqual(cells[0], ["Monster", "Coörd.", "Criterium", "", "Waarde", "", "Eenheid", "Status"])
            self.assertEqual(sum(row[0] == "QC" for row in cells), 1)
            qc = next(index for index, row in enumerate(cells) if row[0] == "QC")
            self.assertEqual(cells[qc][1], "A3, B3")
            self.assertEqual(cells[qc][4], "40.0")
            scoped = next(check for check in result.validation_results if check.sample_name == "QC" and check.wavelength_nm == 450)
            self.assertEqual(cells[qc + 1][4], f"{scoped.value:.3f}")
            self.assertIn("405", cells[qc][2])
            self.assertIn("450", cells[qc + 1][2])
            standard = next(row for row in cells if row[0] == "Standaard 300")
            self.assertEqual(standard[1], "C2, D2")
            standard_check = next(check for check in result.validation_results if check.sample_type == SampleType.STANDARD)
            self.assertEqual(standard[4], f"{standard_check.value:.3f}")
            self.assertEqual(standard[2], "1.000 ≤ ΔExt450")
            self.assertEqual(cells[2][0:2], ["Blanco", "C3, D3"])
            self.assertEqual(cells[2][4:6], ["0.010", "0.0"])
            self.assertEqual(cells[3][4:6], ["0.005", "0.0"])
            self.assertEqual(cells[2][2], "Ext450 ≤ 0.020")
            self.assertEqual(cells[3][2], "Ext405 ≤ 0.010")

    def test_dual_quality_control_stays_beside_metadata_in_one_block(self) -> None:
        data = dual_data()
        data = saturate(data, 450, "qc")
        result = AssayAnalyzer(data, dual_settings()).analyze()
        # One channel can have unavailable CVs while the other retains its duplicates.
        measurement = result.measurement_results[405]
        checks = tuple(
            check if check.sample_type == SampleType.BLANK else replace(check, cv_percent=None)
            for check in measurement.validation_results
        )
        measurements = dict(result.measurement_results)
        measurements[405] = replace(measurement, validation_results=checks)
        result = replace(result, measurement_results=measurements)
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(result, analyst="operator")
            with patch("elisaflow.reporting.calibration.plot_calibration", wraps=plot_calibration) as plot:
                path = writer.write(Path(directory) / "report.pdf")
            self.assertEqual(plot.call_count, 2)
            for call, measurement in zip(plot.call_args_list, result.measurement_results.values(), strict=True):
                self.assertEqual(call.kwargs["model_name"], "4pl")
                self.assertEqual(call.args[2], measurement.calculation_details.curve.parameters)
            headers = [cell.getPlainText() for cell in quality_control_table(writer)._cellvalues[0]]
            self.assertEqual(headers, ["Monster", "Coörd.", "Criterium", "", "Waarde", "", "Eenheid", "Status"])
            with pymupdf.open(path) as document:
                first = document[0]
                metadata = first.search_for("Rungegevens")[0]
                validation = first.search_for("Kwaliteitscontrole")
                self.assertEqual(len(validation), 1)
                self.assertAlmostEqual(metadata.y0, validation[0].y0, places=2)
                self.assertGreater(validation[0].x0, metadata.x1)
                words = first.get_text("words")
                self.assertEqual(sum(word[4] == "Waarde" for word in words), 1)
                self.assertEqual(sum(word[4] == "Eenheid" for word in words), 1)
                self.assertEqual(sum(word[4] == "CV%" for word in words), 1)
                text = " ".join(" ".join(page.get_text().split()) for page in document).replace("∆", "Δ")
                self.assertNotIn("Parameters", text)
                self.assertNotIn("Curvemodel", text)
                self.assertIn("Assaygegevens", text)
                self.assertNotIn("Blanco extinctie", text)
                self.assertNotIn("Kwaliteitscontrole geselecteerde resultaten", text)
                self.assertNotIn("Kwaliteitscontrole 450 nm", text)
                self.assertNotIn("Kwaliteitscontrole 405 nm", text)
                self.assertIn("1.000 ≤ ΔExt450", text)
                for page in document:
                    for word in page.get_text("words"):
                        self.assertTrue(0 <= word[0] < word[2] <= page.rect.width)
                        self.assertTrue(0 <= word[1] < word[3] <= page.rect.height)

    def test_fixed_curve_without_standards_has_a_plot_without_duplicate_metadata(self) -> None:
        result = self.make_result(CalculationMethod.CALIBRATION)
        with TemporaryDirectory() as directory, patch("elisaflow.reporting.calibration.plot_calibration", wraps=plot_calibration) as plot:
            path = PdfReportWriter(result, analyst="operator").write(Path(directory) / "fixed.pdf")
            plot.assert_called_once()
            self.assertFalse(plot.call_args.kwargs["curve_was_fitted"])
            with pymupdf.open(path) as document:
                text = " ".join(" ".join(page.get_text().split()) for page in document)
            self.assertNotIn("Curvemodel:", text)
            self.assertNotIn("Parameters:", text)
            self.assertNotIn("Responsschaal:", text)
            self.assertIn("Vaste kalibratiecurve", text)

    def test_fixed_curve_plot_does_not_require_nominal_concentrations_for_all_standards(self):
        base = self.make_result(CalculationMethod.CALIBRATION)
        samples = samples_from_records(well_table(base.samples).to_dict("records"))
        standards = tuple(
            Sample(f"std-{index}", f"Standard {index}", SampleType.STANDARD,
                   (WellMeasurement(wavelength_nm=450, source_extinction=0.6, well=f"C{index + 1}"),), concentration=concentration)
            for index, concentration in enumerate((None, 1))
        )
        result = AssayAnalyzer(AssayData((*samples, *standards), base.metadata), base.settings).analyze()
        with TemporaryDirectory() as directory:
            report = PdfReportWriter(result, analyst="operator").write(Path(directory) / "fixed.pdf")
            self.assertTrue(report.read_bytes().startswith(b"%PDF"))

    @staticmethod
    def quality_control_row(writer, name):
        return next(row for row in quality_control_table(writer)._cellvalues if row[0].getPlainText() == name)

    @staticmethod
    def patient_names_by_column(page, expected):
        names = [word for word in page.get_text("words") if word[4] in expected]
        positions = sorted({round(word[0], 1) for word in names})
        return [[word[4] for word in sorted(names, key=lambda word: word[1]) if round(word[0], 1) == x] for x in positions]

    @staticmethod
    def make_result(analysis_type: CalculationMethod):
        metadata = assay_metadata()
        extinctions = np.array([0.8, 0.6, 0.5, 0.3])
        if analysis_type == CalculationMethod.CALIBRATION:
            extinctions = FourPLModel().calculate_extinctions([1, 4, 2, 0.5], reference_parameters())
        frame = pd.DataFrame(
            {
                "sample_id": ["blank", "positive", "positive", "cutoff", "extra"],
                "sample_name": ["Blank", "Positive", "Positive", "Cutoff", "Extra"],
                "sample_type": ["Blank", "Control", "Control", "Control", "Control"],
                "source_extinction": [0.075, *(extinctions + 0.075)],
                "well": ["A1", "A2", "A3", "A4", "A5"],
            }
        )
        settings = AssaySettings.from_dict(
            "test.skax",
            {
                "lot_number": "lot",
                **(
                    {
                        "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / cutoff", "references": {"cutoff": "Cutoff"}},
                        "result": {"label": "Ratio"},
                    }
                    if analysis_type == CalculationMethod.FORMULA
                    else {}
                ),
                "controls": {"Positive": {}, "Cutoff": {}, "Missing": {}},
                "blank": {},
                **(
                    {
                        "calculation": {"response_scale": "blank_corrected", "method": "calibration", "fit_model": "4pl", "parameters": {"A": 2.2, "B": -1.3, "C": 2, "D": 0.2}},
                        "result": {"label": "Conc"},
                    }
                    if analysis_type == CalculationMethod.CALIBRATION
                    else {}
                ),
                "validation": {
                    "result": {
                        "controls": {"Positive": {"extinction": {"response_scale": "blank_corrected", "min": 0.1}}, "Missing": {"result": {"min": 0.1}}},
                        "blank": {"extinction": {"response_scale": "raw", "max": 0.1}},
                    }
                },
            },
        )
        return AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), metadata), settings).analyze()

    @classmethod
    def make_fitted_result(cls):
        parameters = reference_parameters()
        nominal = np.repeat([0, 1, 5, 25, 100], 2)
        frame = pd.DataFrame(
            {
                "sample_id": [f"standard-{x}" for x in nominal] + ["blank", "patient", "patient"],
                "sample_name": [f"Standard {x}" for x in nominal] + ["Blank", "Patient 001", "Patient 001"],
                "sample_type": ["Standard"] * len(nominal) + ["Blank", "Unknown", "Unknown"],
                "sample_concentration": [*nominal, np.nan, np.nan, np.nan],
                "source_extinction": [
                    *(FourPLModel().calculate_extinctions(nominal, parameters) + np.tile([-0.002, 0.002], 5) + 0.05),
                    0.05,
                    *FourPLModel().calculate_extinctions([3, 4], parameters) + 0.05,
                ],
                "well": [f"A{index + 1}" for index in range(12)] + ["B1"],
            }
        )
        settings = AssaySettings.from_dict(
            "fitted.skax",
            {"lot_number": "lot", "calculation": {"response_scale": "raw", "method": "calibration", "fit_model": "4pl"}, "result": {"label": "Conc"}},
        )
        return AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), cls.make_result(CalculationMethod.FORMULA).metadata), settings).analyze()

    def test_report_uses_the_entered_full_name_for_the_analyst_and_author(self):
        result = self.make_result(CalculationMethod.FORMULA)
        for analyst in ("Friso van Dijk", "Zo\u00eb de Vries", r"LAB\operator", "  Anne-Marie  "):
            with self.subTest(analyst=analyst), TemporaryDirectory() as directory:
                path = PdfReportWriter(result, analyst=analyst).write(Path(directory) / "report.pdf")
                with pymupdf.open(path) as document:
                    self.assertIn(f"Analist: {analyst.strip()}", document[0].get_text())
                    self.assertEqual(document.metadata["author"], analyst.strip())

    def test_report_requires_a_nonempty_analyst_name(self):
        result = self.make_result(CalculationMethod.FORMULA)
        for analyst in ("", "   ", "Name\nOther"):
            with self.subTest(analyst=analyst), self.assertRaises(ValueError):
                PdfReportWriter(result, analyst=analyst)

    def test_polynomial_reports_include_model_coefficients_and_calibration(self) -> None:
        base = self.make_fitted_result()
        for name, model, params in (
            ("quadratic", QuadraticModel(), QuadraticParameters(0.1, 0.02, 0.0001)),
            ("cubic", CubicModel(), CubicParameters(0.1, 0.02, 0.0001, 0.000001)),
        ):
            with self.subTest(model=name), TemporaryDirectory() as directory:
                frame = well_table(base.samples).copy()
                nominal = frame["sample_concentration"].fillna(3).to_numpy()
                frame["source_extinction"] = model.calculate_extinctions(nominal, params)
                frame.loc[frame["sample_type"] == "Blank", "source_extinction"] = 0.05
                settings = replace(base.settings, calculation=CalibrationSettings(name, response_scale="raw", weighting_exponent=2.0))
                result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), base.metadata), settings).analyze()
                with patch("elisaflow.reporting.calibration.plot_calibration", wraps=plot_calibration) as plot:
                    path = PdfReportWriter(result, analyst="operator").write(Path(directory) / f"{name}.pdf")
                self.assertEqual(plot.call_args.kwargs["model_name"], name)
                self.assertEqual(plot.call_args.args[2], result.calculation_details.curve.parameters)
                with pymupdf.open(path) as document:
                    text = " ".join(" ".join(page.get_text().split()) for page in document)
                    self.assertNotIn("Curvemodel:", text)
                    self.assertNotIn("Parameters:", text)
                    self.assertIn("Kalibratiecurve", text)
                    self.assertIn("Terugvinding%", text)
                    self.assertIn("Patient 001", text)

    def test_five_pl_fitted_calibration_report_contains_plot_and_recoveries(self) -> None:
        base = self.make_fitted_result()
        frame = well_table(base.samples).copy()
        nominal = frame["sample_concentration"].fillna(3).to_numpy()
        frame["source_extinction"] = FivePLModel().calculate_extinctions(nominal, FivePLParameters(2.5, -1.3, 2, 0.1, 0.4)) + 0.05
        frame.loc[frame["sample_type"] == "Blank", "source_extinction"] = 0.05
        settings = replace(base.settings, calculation=CalibrationSettings("5pl", response_scale="raw"))
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), base.metadata), settings).analyze()
        self.assertIsInstance(result.calculation_details.curve.parameters, FivePLParameters)
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(result, analyst="operator")
            with patch("elisaflow.reporting.calibration.plot_calibration", wraps=plot_calibration) as plot:
                path = writer.write(Path(directory) / "five-pl.pdf")
            self.assertEqual(plot.call_args.kwargs["model_name"], "5pl")
            self.assertEqual(plot.call_args.args[2], result.calculation_details.curve.parameters)
            with pymupdf.open(path) as document:
                text = " ".join(" ".join(page.get_text() for page in document).split())
                self.assertIn("Standaard Conc", text)
                self.assertIn("Terugvinding%", text)
                self.assertIn("Patient 001", text)
                self.assertNotIn("Curvemodel:", text)
                self.assertNotIn("Parameters:", text)
                self.assertTrue(any(page.get_images() for page in document))
        table = calibration_table(writer.context, result.samples, 200, scale=writer.measurements[0].calculation_input_scale)
        recoveries = [row[-1].getPlainText() for row in table._cellvalues[1:]]
        self.assertTrue(all(value in ("-", "100.0") for value in recoveries))

    def test_report_precision_is_assay_specific_and_does_not_change_analysis(self) -> None:
        base = self.make_result(CalculationMethod.FORMULA)
        controls = dict(base.settings.validation.result.to_dict()["controls"])
        controls["Positive"] = {"result": {"min": 1.3, "max": 1.5}}
        patients = pd.DataFrame(
            {
                "sample_id": ["patient", "patient"],
                "sample_name": ["Patient 001", "Patient 001"],
                "sample_type": ["Unknown", "Unknown"],
                "sample_dilution": [1, 1],
                "source_extinction": [0.81, 1.03],
                "wavelength_nm": 450,
                "well": ["B1", "B2"],
            }
        )
        data = AssayData(samples_from_records(pd.concat([well_table(base.samples), patients], ignore_index=True).astype(dict.fromkeys(('saturated', 'disabled', 'reference_invalid', 'reference_saturated', 'reference_disabled'), 'boolean')).fillna(dict.fromkeys(('saturated', 'disabled', 'reference_invalid', 'reference_saturated', 'reference_disabled'), False)).to_dict('records')), base.metadata)
        settings = with_validation(base.settings, controls=controls)
        default = AssayAnalyzer(data, settings).analyze()
        report = ReportSettings(result_decimals=4, extinction_decimals=2, percentage_decimals=0)
        precise = AssayAnalyzer(data, replace(settings, report=report)).analyze()
        pd.testing.assert_frame_equal(well_table(default.samples), well_table(precise.samples))
        pd.testing.assert_frame_equal(sample_table(default.samples), sample_table(precise.samples))
        self.assertEqual(default.validation_results, precise.validation_results)
        self.assertEqual(default.calculation_details, precise.calculation_details)
        with TemporaryDirectory() as directory:
            # Interleave two writers to catch accidental changes to module defaults.
            default_writer, precise_writer = PdfReportWriter(default, analyst="operator"), PdfReportWriter(precise, analyst="operator")
            precise_path = precise_writer.write(Path(directory) / "precise.pdf")
            default_path = default_writer.write(Path(directory) / "default.pdf")
            positive = self.quality_control_row(precise_writer, "Positive")
            self.assertEqual(positive[5].getPlainText(), "20")
            self.assertEqual(positive[2].getPlainText(), "1.3000 ≤ Ratio ≤ 1.5000")
            self.assertEqual(positive[-4].getPlainText(), "1.4000")
            default_positive = self.quality_control_row(default_writer, "Positive")
            self.assertEqual(default_positive[5].getPlainText(), "20.2")
            self.assertEqual(default_positive[-4].getPlainText(), "1.4")
            for path, expected in (
                (precise_path, ("Berekening: Ratio (Monster / Cutoff)", "0.88", "1.6900")),
                (default_path, ("Berekening: Ratio (Monster / Cutoff)", "0.875", "1.7")),
            ):
                with pymupdf.open(path) as document:
                    text = "\n".join(page.get_text() for page in document)
                    for value in expected:
                        self.assertIn(value, text)

    def test_custom_cv_threshold_uses_unrounded_values_in_qc_and_patient_tables(self) -> None:
        base = self.make_fitted_result()
        report = ReportSettings(percentage_decimals=0, cv_warning_percent=7.5)
        control_result = self.make_result(CalculationMethod.CALIBRATION)
        controls = sample_table(control_result.samples)
        finite_patient = sample_table(base.samples).loc[sample_table(base.samples)["sample_type"].eq("Unknown")].assign(sample_name="Patient finite")
        samples = pd.concat([controls, sample_table(base.samples), finite_patient], ignore_index=True)
        for cv in (7.5, 7.5001, float("nan")):
            with self.subTest(cv=cv), TemporaryDirectory() as directory:
                samples = samples.copy()
                samples.loc[samples["sample_name"].eq("Positive"), "extinction_cv_percent"] = cv
                # Keep finite CVs so both tables retain their CV columns.
                samples.loc[samples["sample_name"].eq("Extra"), "extinction_cv_percent"] = 1
                samples.loc[samples["sample_type"].eq("Unknown"), "extinction_cv_percent"] = 1
                samples.loc[samples["sample_name"].eq("Patient 001"), "result_cv_percent"] = cv
                checks = tuple(
                    replace(check, cv_percent=cv if np.isfinite(cv) else None) if check.sample_name == "Positive" else check
                    for check in control_result.validation_results
                )
                checks += (replace(next(check for check in checks if check.sample_name == "Positive"), sample_name="Extra", cv_percent=1),)
                result = replace(
                    base,
                    settings=replace(base.settings, report=report),
                    selected=replace(
                        base.selected, samples=samples_with_statistics((samples).to_dict("records"), source_samples=base.samples), validation_results=checks
                    ),
                )
                patient = next(sample for sample in result.samples if sample.name == "Patient 001")
                self.assertIsNotNone(result.patient_extinction(patient))
                if cv == report.cv_warning_percent:
                    # Set equality using the CV of the concrete replicate readings.
                    threshold = max(cv, patient.result.cv_percent)
                    result = replace(result, settings=replace(result.settings, report=replace(report, cv_warning_percent=threshold)))
                writer = PdfReportWriter(result, analyst="operator")
                writer._initialize_report(Path(directory) / "report.pdf")
                qc_cv = self.quality_control_row(writer, "Positive")[5]
                patient_cv = self.patient_table(writer)._cellvalues[2][-1]
                expected_text = "" if np.isnan(cv) else "8"
                expected_color = colors.firebrick if cv > 7.5 else writer._INK
                for cell in (qc_cv, patient_cv):
                    self.assertEqual(cell.getPlainText(), expected_text)
                    self.assertEqual(cell.style.textColor, expected_color)

    def test_calibration_precision_and_custom_recovery_warning_use_unrounded_values(self) -> None:
        base = self.make_fitted_result()
        recoveries = [87.5, 87.4999, 112.5, 112.5001, float("nan"), float("inf")]
        standards = tuple(
            Sample(str(index), "Standard", SampleType.STANDARD, (
                WellMeasurement(
                    well=f"C{index + 1}", wavelength_nm=450, source_extinction=0.56789,
                    curve_extinction=0.56789, back_calculated_concentration=recovery * 1.234 / 100,
                ),
            ), concentration=1.234)
            for index, recovery in enumerate(recoveries)
        )
        for threshold, warned in ((12.5, [False, True, False, True, False, False]), (20, [False] * len(recoveries))):
            with self.subTest(threshold=threshold), TemporaryDirectory() as directory:
                report = ReportSettings(
                    calibration_decimals=3,
                    extinction_decimals=2,
                    percentage_decimals=2,
                    recovery_warning_percent=threshold,
                )
                writer = PdfReportWriter(replace(base, settings=replace(base.settings, report=report)), analyst="operator")
                writer._initialize_report(Path(directory) / "report.pdf")
                rows = calibration_table(writer.context, standards, writer.content_width * 0.45, scale=writer.measurements[0].calculation_input_scale)._cellvalues[1:]
                for index, row in enumerate(rows):
                    expected = ["1.080", "1.080", "1.388", "1.388", "-", "-"][index]
                    self.assertEqual([cell.getPlainText() for cell in row[:3]], ["1.234", "0.57", expected])
                    recovery = row[-1]
                    self.assertEqual(recovery.getPlainText(), ["87.50", "87.50", "112.50", "112.50", "-", "-"][index])
                    self.assertEqual(recovery.style.textColor, colors.firebrick if warned[index] else writer._INK)

    def test_patient_report_shows_dilution_and_keeps_results_separate(self) -> None:
        result = self.make_result(CalculationMethod.CALIBRATION)
        result = replace(result, settings=replace(result.settings, report=ReportSettings(result_decimals=3, calibration_decimals=3)))
        frame = pd.DataFrame(
            {
                "sample_id": ["patient"] * 3,
                "sample_name": ["Patient 001"] * 3,
                "sample_type": ["Unknown"] * 3,
                "sample_dilution": [1, 2.5, 10],
                "source_extinction": FourPLModel().calculate_extinctions([3, 3, 3], result.settings.calculation.parameters),
                "well": ["B1", "B2", "B3"],
            }
        )
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), result.metadata), result.settings).analyze()
        with TemporaryDirectory() as directory:
            path = PdfReportWriter(result, analyst="operator").write(Path(directory) / "report.pdf")
            with pymupdf.open(path) as document:
                text = "\n".join(page.get_text() for page in document)
                self.assertIn("Verd.", text)
                self.assertIn("2.5", text)
                self.assertEqual(text.count("Patient 001"), 3)
                for concentration in ("3.000", "7.500", "30.000"):
                    self.assertIn(concentration, text)

    def test_optional_patient_columns_keep_values_aligned_in_the_pdf(self) -> None:
        data = saturate(dual_data(), 450, "p1", all_replicates=True)
        data = replace(data, samples=tuple(
            replace(sample, name="P405" if sample.sample_id == "p1" else "P450", dilution=2.5 if sample.sample_id == "p1" else 1)
            if sample.sample_type == SampleType.UNKNOWN else sample
            for sample in data.samples
        ))
        config = dual_settings().to_dict()
        config["calculation"]["wavelengths"]["405"] = {
            "fit_model": "4pl", "response_scale": "blank_corrected",
            "parameters": {"A": 1.395, "B": -1.2, "C": 30, "D": 0.020},
        }
        result = AssayAnalyzer(data, AssaySettings.from_dict("test.skax", config)).analyze()
        for visible in (True, False):
            with self.subTest(result_cv=visible), TemporaryDirectory() as directory:
                settings = replace(result.settings, report=replace(result.settings.report, show_result_cv=visible))
                output = PdfReportWriter(replace(result, settings=settings), analyst="operator").write(Path(directory) / "report.pdf")
                with pymupdf.open(output) as document:
                    text = "\n".join(page.get_text() for page in document).replace("∆", "Δ")
                patient_text = text.split("Patiëntresultaten", 1)[1]
                self.assertEqual(patient_text.count("P405"), 1)
                self.assertEqual(patient_text.count("P450"), 1)
                row_405, _, row_450 = patient_text.partition("P405")[2].partition("P450")
                for row, expected in (
                    (row_405, ["E2,F2", "405", "2.5", "ΔExt", "0.543", "0.0", "50.0"]),
                    (row_450, ["G2,H2", "450", "1", "Ext", "0.245", "0.0", "2.0"]),
                ):
                    if visible:
                        expected.append("0.0")
                    self.assertEqual("".join(row.split()), "".join(expected))

    def test_report_displays_bounds_for_controls_and_diluted_patients(self) -> None:
        base = self.make_result(CalculationMethod.CALIBRATION)
        settings = with_validation(
            replace(
                base.settings,
                result=replace(base.settings.result, unit="U/ml"),
                report=ReportSettings(result_decimals=3, calibration_decimals=3),
                calculation=replace(base.settings.calculation, concentration_range=(40, 50)),
                controls={name: SampleSettings.from_dict({}) for name in ("Below", "Above", "Inside", "Missing")},
            ),
            controls={name: {"result": {"min": 30, "max": 88}} for name in ("Below", "Above", "Inside", "Missing")},
        )
        entries = [
            ("Below", "Control", 1, 20),
            ("Below", "Control", 1, 30),
            ("Above", "Control", 1, 60),
            ("Above", "Control", 1, 70),
            ("Inside", "Control", 1, 45),
            ("Patient below", "Unknown", 2, 30),
            ("Patient above", "Unknown", 3, 60),
            ("Patient minimum", "Unknown", 1, 40),
            ("Patient maximum", "Unknown", 1, 50),
            ("Patient missing", "Unknown", 1, np.nan),
        ]
        frame = pd.DataFrame(entries, columns=["sample_name", "sample_type", "sample_dilution", "nominal"])
        frame["sample_id"] = frame["sample_name"]
        frame["source_extinction"] = FourPLModel().calculate_extinctions(frame.pop("nominal"), settings.calculation.parameters)
        frame["well"] = [f"A{index + 1}" for index in range(len(frame))]
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), base.metadata), settings).analyze()
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(result, analyst="operator")
            path = writer.write(Path(directory) / "report.pdf")
            for name, value, status in (
                ("Below", "<40.000", "Niet bepaald"),
                ("Above", ">50.000", "Niet bepaald"),
                ("Inside", "45.000", "OK"),
                ("Missing", "-", "Niet bepaald"),
            ):
                row = self.quality_control_row(writer, name)
                self.assertEqual(row[-3].getPlainText(), value)
                self.assertEqual(row[-1].getPlainText(), status)
            with pymupdf.open(path) as document:
                text = "\n".join(page.get_text() for page in document)
                for value in ("<40.000", ">50.000", "<80.000", ">150.000", "40.000", "50.000", "Conc [U/ml]"):
                    self.assertIn(value, text)
                for page in document:
                    for word in page.get_text("words"):
                        self.assertTrue(0 <= word[0] < word[2] <= page.rect.width)

    def test_control_bounds_outside_acceptance_ranges_are_rejected_in_analysis_and_pdf(self):
        base = self.make_result(CalculationMethod.CALIBRATION)
        settings = with_validation(
            replace(
                base.settings,
                calculation=replace(base.settings.calculation, concentration_range=(0.11, 9)),
                result=replace(base.settings.result, unit="µg/mL"),
                report=ReportSettings(result_decimals=2),
                controls={name: SampleSettings() for name in ("Hoog", "Laag")},
            ),
            controls={"Hoog": {"result": {"min": 2.42, "max": 4.03}}, "Laag": {"result": {"min": 0.44, "max": 0.92}}},
            blank={"extinction": {"response_scale": "raw", "max": 0.2}},
        )
        responses = FourPLModel().calculate_extinctions([12, 0.05], settings.calculation.parameters)
        frame = pd.DataFrame({
            "sample_id": ["blank", "high", "low"],
            "sample_name": ["Blanco", "Hoog", "Laag"],
            "sample_type": ["Blank", "Control", "Control"],
            "well": ["A1", "G1", "H1"],
            "source_extinction": [0.05, *(responses + 0.05)],
        })
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict("records")), base.metadata), settings).analyze()
        self.assertEqual(result.validation_status, ValidationStatus.FAILED)
        self.assertEqual(next(check for check in result.validation_results if check.sample_type == SampleType.BLANK).status, ValidationStatus.PASSED)
        controls = [check for check in result.validation_results if check.sample_type == SampleType.CONTROL]
        self.assertEqual(len(controls), 2)
        self.assertTrue(all(check.status == ValidationStatus.FAILED and check.value is None for check in controls))
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(result, analyst="operator")
            path = writer.write(Path(directory) / "rejected-controls.pdf")
            for name, value in (("Hoog", ">9.00"), ("Laag", "<0.11")):
                row = self.quality_control_row(writer, name)
                self.assertEqual(row[-3].getPlainText(), value)
                self.assertEqual(row[-2].getPlainText(), "µg/mL")
                self.assertEqual(row[-1].getPlainText(), "Afgekeurd")
            with pymupdf.open(path) as document:
                text = " ".join(page.get_text() for page in document)
            self.assertIn(">9.00", text)
            self.assertIn("<0.11", text)
            self.assertEqual(text.count("Afgekeurd"), 2)

    def test_control_result_label_and_unit_are_independent_in_report_criteria(self) -> None:
        base = self.make_result(CalculationMethod.FORMULA)
        for unit in ("index", None):
            with self.subTest(unit=unit), TemporaryDirectory() as directory:
                control_result = {"label": "Controle-index", **({"unit": unit} if unit is not None else {})}
                settings = AssaySettings.from_dict("test.skax", {
                    "lot_number": "lot",
                    "controls": {"Positive": {}},
                    "validation": {"result": {"controls": {
                        "Positive": {"result": {"min": 1.3}}, "Cutoff": {"result": {"min": 0}},
                    }}},
                    "calculation": {
                        "response_scale": "blank_corrected",
                        "method": "formula",
                        "expression": "sample / (cutoff * 0.5)",
                        "references": {"cutoff": "Cutoff"},
                        "control_overrides": {"Positive": {"expression": "sample / cutoff", "response_scale": "blank_corrected"}},
                        "control_results": {"Positive": control_result},
                    },
                    "result": {"label": "Score", "unit": "%"},
                })
                result = AssayAnalyzer(AssayData(samples_from_records(well_table(base.samples).to_dict('records')), base.metadata), settings).analyze()
                writer = PdfReportWriter(result, analyst="operator")
                path = writer.write(Path(directory) / "report.pdf")
                rows = quality_control_table(writer)._cellvalues
                row = next(row for row in rows if row[0].getPlainText() == "Positive")
                self.assertEqual(row[-4].getPlainText(), "1.4")
                self.assertIn("Controle-index (Positive / Cutoff)", row[2].getPlainText())
                self.assertNotIn("Score", row[2].getPlainText())
                self.assertEqual(row[-2].getPlainText(), unit or "")
                self.assertEqual(row[-1].getPlainText(), "OK")
                inherited = next(row for row in rows if row[0].getPlainText() == "Cutoff")
                self.assertIn("Score", inherited[2].getPlainText())
                self.assertEqual(inherited[-2].getPlainText(), "%")
                with pymupdf.open(path) as document:
                    text = " ".join(" ".join(page.get_text().split()) for page in document)
                    self.assertIn("Berekening: Score (Monster / (Cutoff * 0.5))", text)
                    self.assertIn("Positive / Cutoff", text)
                    self.assertIn("0.0 ≤ Score", text)

    def test_reports_pair_criteria_with_values_and_omit_separate_extinction_columns(self) -> None:
        for analysis_type in CalculationMethod:
            with self.subTest(analysis_type=analysis_type), TemporaryDirectory() as directory:
                result = self.make_result(analysis_type)
                writer = PdfReportWriter(result, analyst="operator")
                path = writer.write(Path(directory) / "report.pdf")
                self.assertTrue(path.read_bytes().startswith(b"%PDF-"))
                table = quality_control_table(writer)
                cells = [[cell.getPlainText() for cell in row] for row in table._cellvalues]
                self.assertEqual(cells[0], ["Monster", "Coörd.", "Criterium", "", "Waarde", "", "Eenheid", "Status"])
                self.assertEqual(cells[1], ["", "", "", "", "Gem.", "CV%", "", ""])
                self.assertIn(("SPAN", (4, 0), (5, 0)), table._spanCmds)
                self.assertEqual(next(row[1] for row in cells[2:] if row[0] == "Positive"), "A2, A3")
                rows = {row[0]: row[2:] for row in cells[2:]}
                expected_extinction = (
                    "0.700" if analysis_type == CalculationMethod.FORMULA else f"{sample_table(result.samples).iloc[1]['extinction_mean']:.3f}"
                )
                sample = sample_table(result.samples).iloc[1]
                self.assertEqual(
                    rows["Positive"][2:4],
                    [
                        expected_extinction,
                        f"{sample['extinction_cv_percent']:.{result.settings.report.percentage_decimals}f}",
                    ],
                )
                self.assertEqual(rows["Positive"][-4], expected_extinction)
                self.assertEqual(rows["Cutoff"], ["", "", "", "", "", "-"])
                self.assertEqual(rows["Extra"], ["", "", "", "", "", "-"])
                self.assertEqual(rows["Missing"][1], "")
                self.assertEqual(rows["Missing"][-4], "-")
                self.assertEqual(rows["Missing"][-1], "Niet bepaald")
                self.assertEqual(cells[2], ["Blanco", "A1", "Ext ≤ 0.100", "", "0.075", "", "", "OK"])
                with pymupdf.open(path) as document:
                    for label in ("Waarde", "Eenheid", "Status"):
                        headers = [word for word in document[0].get_text("words") if word[4] == label]
                        self.assertEqual(len(headers), 1)
                        self.assertGreater(headers[0][0], 0)
                    text = " ".join(document[0].get_text().split()).replace("∆", "Δ")
                    self.assertLess(text.index("Blanco"), text.index("Positive"))
                    self.assertIn("ΔExt", text)

    def test_labels_follow_the_actual_extinction_scale_of_each_report_section(self) -> None:
        base = self.make_fitted_result()
        settings = with_validation(
            replace(base.settings, blank=SampleSettings.from_dict({}), standards={25: SampleSettings.from_dict({})}),
            blank={"extinction": {"response_scale": "raw", "max": 0.1}},
            standards={25: {"extinction": {"response_scale": "blank_corrected", "min": 0.1}}},
        )
        for blank_value in (None, 0.0, 0.05):
            with self.subTest(blank=blank_value), TemporaryDirectory() as directory:
                frame = well_table(base.samples).copy()
                blank_mask = frame["sample_type"].eq("Blank")
                if blank_value is None:
                    frame = frame.loc[~blank_mask].copy()
                else:
                    frame.loc[blank_mask, "source_extinction"] = blank_value
                result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), base.metadata), settings).analyze()
                label = "Ext" if blank_value is None else "ΔExt"
                writer = PdfReportWriter(result, analyst="operator")
                with patch("elisaflow.reporting.calibration.plot_calibration", wraps=plot_calibration) as plot:
                    path = writer.write(Path(directory) / "report.pdf")
                self.assertEqual(plot.call_args.kwargs["extinction_label"], "Ext")
                standards = select_samples(result.samples, SampleType.STANDARD)
                calibration = calibration_table(writer.context, standards, writer.content_width * 0.45, scale=writer.measurements[0].calculation_input_scale)
                self.assertEqual(calibration._cellvalues[0][1].getPlainText(), "Ext")
                qc = quality_control_table(writer)._cellvalues
                self.assertEqual(qc[0][2].getPlainText(), "Criterium")
                standard_row = next(row for row in qc if row[0].getPlainText() == "Standaard 25")
                self.assertEqual(standard_row[2].getPlainText(), f"0.100 ≤ {label}")
                self.assertEqual(qc[2][0].getPlainText(), "Blanco")
                self.assertEqual(qc[2][2].getPlainText(), "Ext ≤ 0.100")
                self.assertEqual(qc[2][4].getPlainText(), "-" if blank_value is None else f"{blank_value:.3f}")
                with pymupdf.open(path) as document:
                    text = " ".join(" ".join(page.get_text().split()) for page in document).replace("∆", "Δ")
                    patient_text = text.split("Patiëntresultaten", 1)[1]
                    self.assertIn("Ext", patient_text)
                    self.assertNotIn("ΔExt", patient_text)
                    if blank_value is None:
                        self.assertNotIn("Δ", text)

    def test_patient_extinction_values_and_cvs_follow_the_final_calculation_scale(self):
        base = self.make_fitted_result()
        data = AssayData(samples_from_records(well_table(base.samples).to_dict("records")), base.metadata)
        patient = next(sample for sample in data.samples if sample.sample_type == SampleType.UNKNOWN)
        raw = np.array([well.source_extinction for well in patient.measurements])
        corrected = raw - base.blank_value
        parameters = reference_parameters()
        raw_parameters = replace(parameters, A=parameters.A + base.blank_value, D=parameters.D + base.blank_value)
        for calculation, label, values in (
            (base.settings.calculation, "Ext", raw),
            (CalibrationSettings("4pl", raw_parameters, response_scale="raw"), "Ext", raw),
            (CalibrationSettings("4pl", parameters, response_scale="blank_corrected"), "ΔExt", corrected),
            (FormulaSettings(Formula("sample"), response_scale="blank_corrected"), "ΔExt", corrected),
            (FormulaSettings(Formula("sample"), response_scale="raw"), "Ext", raw),
        ):
            with self.subTest(calculation=calculation), TemporaryDirectory() as directory:
                result = AssayAnalyzer(data, replace(base.settings, calculation=calculation)).analyze()
                writer = PdfReportWriter(result, analyst="operator")
                path = Path(directory) / "report.pdf"
                writer._initialize_report(path)
                table = self.patient_table(writer)
                self.assertEqual(table._cellvalues[0][2].getPlainText(), label)
                self.assertEqual(table._cellvalues[2][2].getPlainText(), f"{values.mean():.3f}")
                self.assertEqual(table._cellvalues[2][3].getPlainText(), f"{100 * values.std(ddof=1) / values.mean():.1f}")
                writer.write(path)
                with pymupdf.open(path) as document:
                    text = " ".join(" ".join(page.get_text().split()) for page in document).replace("∆", "Δ")
                patients = text.split("Patiëntresultaten", 1)[1]
                self.assertIn(label, patients)
                self.assertIn(f"{values.mean():.3f}", patients)
                if label == "Ext":
                    self.assertNotIn("ΔExt", patients)

    def test_patient_extinction_reports_the_preprocessed_input_for_fixed_curves(self):
        base = self.make_fitted_result()
        data = AssayData(samples_from_records(well_table(base.samples).to_dict("records")), base.metadata)
        raw = np.array([well.source_extinction for sample in data.samples if sample.sample_type == SampleType.UNKNOWN for well in sample.measurements])
        parameters = reference_parameters()
        raw_parameters = replace(parameters, A=parameters.A + base.blank_value, D=parameters.D + base.blank_value)
        for calculation, label, inputs in (
            (CalibrationSettings("4pl", replace(raw_parameters, A=2 * raw_parameters.A, D=2 * raw_parameters.D),
                                 response_scale="raw", preprocessing=Formula("sample * 2")), "Ext", raw),
            (CalibrationSettings("4pl", replace(parameters, A=2 * parameters.A, D=2 * parameters.D),
                                 preprocessing=Formula("sample * 2"), response_scale="blank_corrected"), "ΔExt", raw - base.blank_value),
        ):
            with self.subTest(calculation=calculation), TemporaryDirectory() as directory:
                result = AssayAnalyzer(data, replace(base.settings, calculation=calculation)).analyze()
                writer = PdfReportWriter(result, analyst="operator")
                writer._initialize_report(Path(directory) / "report.pdf")
                table = self.patient_table(writer)
                self.assertEqual(table._cellvalues[0][2].getPlainText(), label)
                self.assertEqual(table._cellvalues[2][2].getPlainText(), f"{2 * inputs.mean():.3f}")
                with patch("elisaflow.reporting.calibration.plot_calibration", wraps=plot_calibration) as plot:
                    writer.write(Path(directory) / "report.pdf")
                self.assertEqual(plot.call_args.kwargs["preprocessing"], calculation.preprocessing)
                with pymupdf.open(Path(directory) / "report.pdf") as document:
                    text = " ".join(page.get_text() for page in document)
                self.assertNotIn("(curve-invoer)", text)

    def test_patient_extinction_uses_each_selected_wavelengths_scale(self):
        data = saturate(dual_data(), 450, "p1", all_replicates=True)
        for corrected_wavelength, top, baseline, blank in ((405, 1.4, 0.025, 0.005), (450, 4.5, 0.08, 0.01)):
            config = dual_settings().to_dict()
            config["calculation"]["wavelengths"][str(corrected_wavelength)] = {
                "fit_model": "4pl", "response_scale": "blank_corrected",
                "parameters": {"A": top - blank, "B": -1.2, "C": 30, "D": baseline - blank},
            }
            result = AssayAnalyzer(data, AssaySettings.from_dict("test.skax", config)).analyze()
            patients = [sample for sample in result.samples if sample.sample_type == SampleType.UNKNOWN]
            self.assertEqual([sample.selected_wavelength_nm for sample in patients], [405, 450])
            with self.subTest(corrected_wavelength=corrected_wavelength), TemporaryDirectory() as directory:
                writer = PdfReportWriter(result, analyst="operator")
                path = Path(directory) / "report.pdf"
                writer._initialize_report(path)
                table = self.patient_table(writer)
                self.assertEqual(table._cellvalues[0][3].getPlainText(), "Extinctie")
                for row, patient in zip(table._cellvalues[2:], patients, strict=True):
                    corrected = patient.selected_wavelength_nm == corrected_wavelength
                    values = np.array([well.source_extinction for well in patient.retained]) - (blank if corrected else 0)
                    label = "ΔExt" if corrected else "Ext"
                    self.assertEqual(row[2].getPlainText(), str(patient.selected_wavelength_nm))
                    self.assertEqual(row[3].getPlainText(), f"{label} {values.mean():.3f}")
                writer.write(path)
                with pymupdf.open(path) as document:
                    text = " ".join(" ".join(page.get_text().split()) for page in document)
                    patient_text = text.split("Patiëntresultaten", 1)[1]
                    self.assertIn("405", patient_text)
                    self.assertIn("450", patient_text)
                    for page in document:
                        for word in page.get_text("words"):
                            self.assertTrue(0 <= word[0] < word[2] <= page.rect.width)

    def test_patient_extinction_uses_reference_corrected_and_converted_inputs(self):
        from test_reference_wavelength import reference_case
        from test_saturation_conversion import INTERCEPT, SLOPE, conversion_settings

        reference_data, reference_settings = reference_case(switching=False)
        converted_data = saturate(dual_data(), 450, "p1", all_replicates=True)
        frame = well_table(converted_data.samples)
        frame.loc[frame["sample_id"].eq("p1") & frame["wavelength_nm"].eq(450), "source_extinction"] = 99
        converted_data = replace(converted_data, samples=samples_from_records(frame.to_dict("records")))
        for data, settings in ((reference_data, reference_settings), (converted_data, conversion_settings())):
            result = AssayAnalyzer(data, settings).analyze()
            patient = next(sample for sample in result.samples if sample.sample_id == "p1")
            source = well_table(data.samples)
            source = source.loc[source["sample_id"].eq("p1")].set_index(["wavelength_nm", "well"])["source_extinction"]
            if settings.measurement.reference_wavelength_nm is not None:
                inputs = source.loc[450] - source.loc[620]
            else:
                self.assertEqual(patient.selected_wavelength_nm, 450)
                inputs = INTERCEPT + SLOPE * source.loc[405]
            with self.subTest(reference=settings.measurement.reference_wavelength_nm), TemporaryDirectory() as directory:
                writer = PdfReportWriter(result, analyst="operator")
                writer._initialize_report(Path(directory) / "report.pdf")
                table = self.patient_table(writer)
                column = 3 if writer.context.multiple_measurements else 2
                self.assertEqual(table._cellvalues[0][column].getPlainText(), "Ext")
                self.assertEqual(table._cellvalues[2][column].getPlainText(), f"{inputs.mean():.3f}")

    def test_patient_without_a_selected_result_does_not_display_an_arbitrary_extinction(self):
        data = saturate(saturate(dual_data(), 450, "p1", all_replicates=True), 405, "p1", all_replicates=True)
        result = AssayAnalyzer(data, dual_settings()).analyze()
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(result, analyst="operator")
            writer._initialize_report(Path(directory) / "report.pdf")
            table = self.patient_table(writer)
            self.assertEqual(table._cellvalues[2][3].getPlainText(), "-")
            self.assertEqual(table._cellvalues[2][4].getPlainText(), "")

    def test_raw_calibration_labels_and_exclusion_values_keep_extinction_qc_separate(self) -> None:
        base = self.make_fitted_result()
        frame = well_table(base.samples).copy()
        outlier = frame.loc[frame["sample_concentration"] == 0].iloc[[0]].copy()
        outlier["source_extinction"] = 0.4
        outlier["well"] = "C1"
        frame = pd.concat([frame, outlier], ignore_index=True)
        settings = with_validation(
            replace(
                base.settings,
                calculation=replace(base.settings.calculation, response_scale="raw", weighting_exponent=2.0),
                blank=SampleSettings.from_dict({}),
                standards={0: SampleSettings(replicate_filter=ReplicateFilterSettings(25, 1, response_scale="blank_corrected")), 25: SampleSettings.from_dict({})},
            ),
            blank={"extinction": {"response_scale": "raw", "max": 0.1}},
            standards={25: {"extinction": {"response_scale": "blank_corrected", "min": 0.1}}},
        )
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), base.metadata), settings).analyze()
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(result, analyst="operator")
            with patch("elisaflow.reporting.calibration.plot_calibration", wraps=plot_calibration) as plot:
                path = writer.write(Path(directory) / "raw-calibration.pdf")
            self.assertEqual(plot.call_args.kwargs["extinction_label"], "Ext")
            standards = select_samples(result.samples, SampleType.STANDARD)
            calibration = calibration_table(writer.context, standards, writer.content_width * 0.45, scale=writer.measurements[0].calculation_input_scale)
            self.assertEqual(calibration._cellvalues[0][1].getPlainText(), "Ext")
            qc = quality_control_table(writer)._cellvalues
            standard_row = next(row for row in qc if row[0].getPlainText() == "Standaard 25")
            self.assertIn("ΔExt", standard_row[2].getPlainText())
            self.assertEqual(qc[2][0].getPlainText(), "Blanco")
            self.assertEqual(qc[2][2].getPlainText(), "Ext ≤ 0.100")
            writer._write_replicate_exclusions(writer.measurements[0])
            exclusions = writer.story[-1]._cellvalues
            self.assertEqual(exclusions[0][2], "ΔExt")
            self.assertEqual(len(exclusions), 2)
            self.assertEqual(exclusions[1][2], "0.350")
            with pymupdf.open(path) as document:
                text = " ".join(" ".join(page.get_text().split()) for page in document)
                self.assertIn("Kalibratiecurve", text)
                self.assertIn("Uitgesloten replicaten", text)

    def test_blank_quality_control_can_be_the_only_qc_table(self) -> None:
        result = self.make_fitted_result()
        settings = with_validation(replace(result.settings, blank=SampleSettings.from_dict({})), blank={"extinction": {"response_scale": "raw", "max": 0.1}})
        result = AssayAnalyzer(AssayData(samples_from_records(well_table(result.samples).to_dict('records')), result.metadata), settings).analyze()
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(result, analyst="operator")
            path = writer.write(Path(directory) / "report.pdf")
            table = quality_control_table(writer)
            self.assertIsNotNone(table)
            self.assertEqual([row[0].getPlainText() for row in table._cellvalues[table.repeatRows:]], ["Blanco"])
            with pymupdf.open(path) as document:
                text = " ".join(document[0].get_text().split())
                self.assertIn("Kwaliteitscontrole", text)
                self.assertIn("Blanco", text)

    def test_flagged_readings_are_marked_and_only_usable_patient_results_are_calculated(self) -> None:
        base = self.make_result(CalculationMethod.FORMULA)
        frame = well_table(base.samples).assign(saturated=False, disabled=False)
        patients = pd.DataFrame(
            {
                "sample_id": ["paired", "paired", "disabled", "saturated", "both", "healthy"],
                "sample_name": ["Paired", "Paired", "Disabled", "Unavailable", "Both", "Healthy"],
                "sample_type": "Unknown",
                "sample_dilution": 1,
                "well": [f"B{index + 1}" for index in range(6)],
                "source_extinction": [1.075, 6.1234, 2.075, 7.2345, 8.3456, 3.075],
                "wavelength_nm": 450,
                "saturated": [False, True, False, True, True, False],
                "disabled": [False, False, True, False, True, False],
            }
        )
        result = AssayAnalyzer(AssayData(samples_from_records(pd.concat([frame, patients], ignore_index=True).astype(dict.fromkeys(('saturated', 'disabled', 'reference_invalid', 'reference_saturated', 'reference_disabled'), 'boolean')).fillna(dict.fromkeys(('saturated', 'disabled', 'reference_invalid', 'reference_saturated', 'reference_disabled'), False)).to_dict('records')), base.metadata), base.settings).analyze()
        with TemporaryDirectory() as directory:
            path = PdfReportWriter(result, analyst="operator").write(Path(directory) / "report.pdf")
            with pymupdf.open(path) as document:
                text = " ".join(" ".join(page.get_text().split()) for page in document)
                self.assertIn("SAT", text)
                self.assertIn("UIT", text)
                self.assertIn("Meting verzadigd in SkanIt", text)
                self.assertIn("Well uitgeschakeld in SkanIt.", text)
                self.assertNotIn("6.123", text)
                self.assertNotIn("7.234", text)
                self.assertNotIn("8.346", text)
                patient_text = text.split("Patiëntresultaten", 1)[1]
                self.assertIn("Paired", patient_text)
                self.assertIn("Unavailable", patient_text)
                self.assertIn("Healthy", patient_text)
                self.assertNotIn("Disabled", patient_text)
                self.assertNotIn("Both", patient_text)
        samples = sample_table(result.samples).set_index("sample_id")
        self.assertAlmostEqual(samples.loc["paired", "result_mean"], 2)
        self.assertAlmostEqual(samples.loc["healthy", "result_mean"], 6)
        self.assertTrue(np.isnan(samples.loc["saturated", "result_mean"]))

    def test_multiple_criteria_keep_their_values_and_statuses_in_separate_rows(self) -> None:
        for analysis_type, minimum, expected in (
            (CalculationMethod.FORMULA, 1.5, "1.400"),
            (CalculationMethod.CALIBRATION, 3, "2.500"),
        ):
            with self.subTest(analysis_type=analysis_type), TemporaryDirectory() as directory:
                base = self.make_result(analysis_type)
                controls = {
                    "Positive": {
                        "extinction": {"response_scale": "blank_corrected", "min": 0.1, "max": 2},
                        "result": {"min": minimum},
                    }
                }
                settings = with_validation(
                    replace(base.settings, report=ReportSettings(result_decimals=3, calibration_decimals=3)), controls=controls
                )
                result = AssayAnalyzer(AssayData(samples_from_records(well_table(base.samples).to_dict('records')), base.metadata), settings).analyze()
                writer = PdfReportWriter(result, analyst="operator")
                path = writer.write(Path(directory) / "report.pdf")
                table = quality_control_table(writer)
                rows = table._cellvalues
                index = next(i for i, row in enumerate(rows) if row[0].getPlainText() == "Positive")
                extinction_row, outcome_row = rows[index : index + 2]
                self.assertEqual(extinction_row[2].getPlainText(), "0.100 ≤ ΔExt ≤ 2.000")
                positive = sample_table(result.samples).loc[sample_table(result.samples)["sample_name"].eq("Positive")].iloc[0]
                self.assertEqual(extinction_row[-4].getPlainText(), f"{positive['extinction_mean']:.3f}")
                self.assertEqual(extinction_row[-1].getPlainText(), "OK")
                self.assertIn("Ratio" if analysis_type == CalculationMethod.FORMULA else "Conc", outcome_row[2].getPlainText())
                self.assertEqual(outcome_row[-4].getPlainText(), expected)
                self.assertEqual(extinction_row[-3].getPlainText(), f"{positive['extinction_cv_percent']:.1f}")
                self.assertEqual(outcome_row[-3].getPlainText(), f"{positive['result_cv_percent']:.1f}")
                self.assertEqual(outcome_row[-1].getPlainText(), "Afgekeurd")
                self.assertIn(("SPAN", (0, index), (0, index + 1)), table._spanCmds)
                with pymupdf.open(path) as document:
                    words = document[0].get_text("words")
                    value_word = next(word for word in words if word[4] == expected)
                    status_word = next(word for word in words if word[4] == "Afgekeurd")
                    self.assertAlmostEqual(value_word[1], status_word[1], places=2)

    def test_dengue_quality_controls_follow_the_skax_order_in_the_table_and_pdf(self) -> None:
        root = Path(__file__).resolve().parents[1]
        protocol = root / "protocols" / "Dengue_NS1_as.skax"
        settings = replace(AssayRegistry(root / "configs").get(protocol), measurement=MeasurementSettings())
        wells = SkaxReader.from_file(protocol).wells
        frame = wells.loc[wells["sample_type"].eq("Control")].copy()
        frame["source_extinction"] = frame["sample_name"].map(
            {
                "Positief": 0.6,
                "Grens": 0.2,
                "Negatief": 0.1,
            }
        )
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.make_result(CalculationMethod.FORMULA).metadata), settings).analyze()
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(result, analyst="operator")
            path = writer.write(Path(directory) / "report.pdf")
            table = quality_control_table(writer)
            rows = table._cellvalues[table.repeatRows :]
            self.assertEqual([row[0].getPlainText() for row in rows], ["Positief", "Grens", "Negatief"])
            with pymupdf.open(path) as document:
                words = document[0].get_text("words")
                controls = sorted(
                    [
                        word
                        for word in words
                        if word[0] >= document[0].search_for("Kwaliteitscontrole")[0].x0 and word[4] in {"Positief", "Grens", "Negatief"}
                    ],
                    key=lambda word: word[1],
                )
                self.assertEqual([word[4] for word in controls], ["Positief", "Grens", "Negatief"])

    def test_quality_control_order_runs_down_numeric_columns_and_keeps_missing_targets(self) -> None:
        frame = pd.DataFrame(
            {
                "sample_id": ["standard-1", "positive", "extra", "cutoff", "standard-10", "positive", "blank"],
                "sample_name": ["Cal 1", "Positive", "Extra", "Cutoff", "Cal 10", "Positive", "Blank"],
                "sample_type": ["Standard", "Control", "Control", "Control", "Standard", "Control", "Blank"],
                "sample_concentration": [1, np.nan, np.nan, np.nan, 10, np.nan, np.nan],
                "source_extinction": [0.8, 0.9, 0.4, 0.5, 0.7, 0.8, 0.05],
                "well": ["A10", "B10", "B2", "H1", "A2", "C2", "A1"],
                "well_x": [9, 9, 1, 0, 1, 1, 0],
                "well_y": [0, 1, 1, 7, 0, 2, 0],
                "plate_index": [0] * 7,
            }
        )
        settings = AssaySettings.from_dict(
            "test.skax",
            {
                "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / cutoff", "references": {"cutoff": "Cutoff"}},
                "lot_number": "lot",
                "controls": {"Positive": {}, "Cutoff": {}, "Missing": {}},
                "blank": {},
                "validation": {
                    "result": {
                        "controls": {"Positive": {"result": {"min": 1}}, "Cutoff": {"extinction": {"response_scale": "blank_corrected", "min": 0.1}}, "Missing": {"result": {"min": 1}}},
                        "standards": {str(nominal): {"extinction": {"response_scale": "blank_corrected", "min": 0}} for nominal in (1, 10, 100)},
                        "blank": {"extinction": {"response_scale": "raw", "max": 0.1}},
                    }
                },
            },
        )
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), self.make_result(CalculationMethod.FORMULA).metadata), settings).analyze()
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(result, analyst="operator")
            writer._initialize_report(Path(directory) / "report.pdf")
            rows = quality_control_table(writer)._cellvalues[2:]
            self.assertEqual(
                [row[0].getPlainText() for row in rows],
                [
                    "Blanco",
                    "Cutoff",
                    "Standaard 10",
                    "Extra",
                    "Positive",
                    "Standaard 1",
                    "Missing",
                    "Standaard 100",
                ],
            )
            self.assertEqual(rows[-2][-1].getPlainText(), "Niet bepaald")
            self.assertEqual(rows[-1][-1].getPlainText(), "Niet bepaald")
            coordinates = {row[0].getPlainText(): row[1].getPlainText() for row in rows}
            self.assertEqual(coordinates["Positive"], "C2, B10")
            self.assertEqual(coordinates["Standaard 10"], "A2")
            self.assertEqual(coordinates["Missing"], "")
            self.assertEqual(coordinates["Standaard 100"], "")

    def test_controls_without_acceptance_rules_keep_names_and_coordinates(self) -> None:
        result = self.make_result(CalculationMethod.FORMULA)
        result = replace(result, selected=replace(result.selected, validation_results=()))
        self.assertEqual(result.validation_status, ValidationStatus.NOT_EVALUATED)
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(result, analyst="operator")
            writer.write(Path(directory) / "report.pdf")
            table = quality_control_table(writer)
            cells = [[cell.getPlainText() for cell in row] for row in table._cellvalues[table.repeatRows:]]
            self.assertEqual(
                cells,
                [
                    ["Positive", "A2, A3", "", "", "", "", "-"],
                    ["Cutoff", "A4", "", "", "", "", "-"],
                    ["Extra", "A5", "", "", "", "", "-"],
                ],
            )

    def test_unavailable_criterion_values_display_a_dash(self) -> None:
        base = self.make_result(CalculationMethod.FORMULA)
        for value in (None, np.nan, np.inf):
            with self.subTest(value=value), TemporaryDirectory() as directory:
                checks = tuple(
                    replace(check, value=value, status=ValidationStatus.NOT_EVALUATED) if check.sample_name == "Positive" else check
                    for check in base.validation_results
                )
                writer = PdfReportWriter(
                    replace(
                        base,
                        selected=replace(base.selected, validation_results=checks),
                    )
                , analyst="operator")
                writer.write(Path(directory) / "report.pdf")
                positive = self.quality_control_row(writer, "Positive")
                self.assertEqual(positive[2].getPlainText(), "0.100 ≤ ΔExt")
                self.assertEqual(positive[-3].getPlainText(), "-")
                self.assertEqual(positive[-1].getPlainText(), "Niet bepaald")

    def test_report_renders_stored_quality_control_without_recalculating(self) -> None:
        result = self.make_result(CalculationMethod.FORMULA)
        samples = sample_table(result.samples)
        samples["extinction_mean"] = 999
        samples["extinction_cv_percent"] = 999
        wells = well_table(result.samples)
        wells["extinction"] = 999
        with (
            TemporaryDirectory() as directory,
            patch("elisaflow.quality_control.summary.summarize_quality_control", side_effect=AssertionError("QC was recalculated")),
        ):
            writer = PdfReportWriter(result, analyst="operator")
            writer.write(Path(directory) / "report.pdf")
            positive = self.quality_control_row(writer, "Positive")
            self.assertEqual(positive[4].getPlainText(), "0.700")
            self.assertEqual(positive[5].getPlainText(), "20.2")

    def test_report_displays_instrument_name_or_missing_marker(self) -> None:
        result = self.make_result(CalculationMethod.FORMULA)
        for name, expected in ((None, "-"), ("Multiskan SkyHigh SN 1234", "Multiskan SkyHigh SN 1234")):
            with self.subTest(name=name), TemporaryDirectory() as directory:
                named = replace(result, metadata=replace(result.metadata, instrument_name=name))
                path = PdfReportWriter(named, analyst="operator").write(Path(directory) / "report.pdf")
                with pymupdf.open(path) as document:
                    text = " ".join(document[0].get_text().split())
                self.assertIn(f"Instrument: {expected}", text)

    def test_quality_control_extinction_cv_warning_uses_unrounded_values(self) -> None:
        for analysis_type in CalculationMethod:
            for extinction_cv in (20, 20.0001, 21, 0, 19.999):
                with self.subTest(analysis_type=analysis_type, cv=extinction_cv):
                    result = self.make_result(analysis_type)
                    checks = tuple(
                        replace(check, cv_percent=extinction_cv) if check.sample_name == "Positive" else check
                        for check in result.validation_results
                    )
                    result = replace(
                        result,
                        selected=replace(result.selected, validation_results=checks),
                    )
                    with TemporaryDirectory() as directory:
                        writer = PdfReportWriter(result, analyst="operator")
                        path = writer.write(Path(directory) / "report.pdf")
                        positive = self.quality_control_row(writer, "Positive")
                        cell = positive[5]
                        self.assertEqual(cell.getPlainText(), f"{extinction_cv:.{result.settings.report.percentage_decimals}f}")
                        self.assertEqual(cell.style.textColor, colors.firebrick if extinction_cv > 20 else writer._INK)
                        with pymupdf.open(path) as document:
                            spans = [
                                span
                                for page in document
                                for block in page.get_text("dict")["blocks"]
                                if "lines" in block
                                for line in block["lines"]
                                for span in line["spans"]
                            ]
                            expected_color = 0xB22222 if extinction_cv > 20 else 0x263442
                            self.assertTrue(
                                any(
                                    span["text"] == f"{extinction_cv:.{result.settings.report.percentage_decimals}f}"
                                    and span["color"] == expected_color
                                    for span in spans
                                )
                            )

    def test_single_measurements_and_unavailable_cvs_omit_mean_and_cv_subheaders(self) -> None:
        for analysis_type in CalculationMethod:
            result = self.make_result(analysis_type)
            single = AssayAnalyzer(AssayData(samples_from_records(well_table(result.samples).drop(index=2).to_dict('records')), result.metadata), result.settings).analyze()
            result = replace(
                result,
                selected=replace(
                    result.selected,
                    validation_results=tuple(replace(check, cv_percent=None) for check in result.selected.validation_results),
                ),
            )
            for unavailable in (single, result):
                with self.subTest(analysis_type=analysis_type), TemporaryDirectory() as directory:
                    writer = PdfReportWriter(unavailable, analyst="operator")
                    writer.write(Path(directory) / "report.pdf")
                    table = quality_control_table(writer)
                    headers = [cell.getPlainText() for cell in table._cellvalues[0]]
                    self.assertEqual(
                        headers,
                        ["Monster", "Coörd.", "Criterium", "", "Waarde", "Eenheid", "Status"],
                    )
                    self.assertEqual(table.repeatRows, 1)
                    self.assertEqual(table._cellvalues[1][0].getPlainText(), "Blanco")
                    positive = next(row for row in table._cellvalues if row[0].getPlainText() == "Positive")
                    self.assertEqual(positive[4].getPlainText(), f"{sample_table(unavailable.samples).iloc[1]['extinction_mean']:.3f}")

    def test_control_cv_combines_measurements_with_different_sample_ids(self) -> None:
        result = self.make_result(CalculationMethod.FORMULA)
        frame = well_table(result.samples).copy()
        frame.loc[2, "sample_id"] = "another-positive"
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), result.metadata), result.settings).analyze()
        controls = sample_table(result.samples).loc[sample_table(result.samples)["sample_name"].eq("Positive")]
        self.assertEqual(controls["count"].tolist(), [1, 1])
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(result, analyst="operator")
            writer.write(Path(directory) / "report.pdf")
            positive = self.quality_control_row(writer, "Positive")
            cv = f"{20.203050891044217:.{result.settings.report.percentage_decimals}f}"
            self.assertEqual(positive[1].getPlainText(), "A2, A3")
            self.assertEqual(positive[5].getPlainText(), cv)

    def test_blank_cv_uses_raw_extinctions(self) -> None:
        result = self.make_result(CalculationMethod.FORMULA)
        frame = well_table(result.samples).copy()
        blank = frame.iloc[[0]].copy()
        blank["source_extinction"] = 0.15
        blank["well"] = "B1"
        frame = pd.concat([frame, blank], ignore_index=True)
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), result.metadata), result.settings).analyze()
        self.assertAlmostEqual(sample_table(result.samples).iloc[0]["extinction_mean"], 0)
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(result, analyst="operator")
            path = writer.write(Path(directory) / "report.pdf")
            row = next(row for row in quality_control_table(writer)._cellvalues if row[0].getPlainText() == "Blanco")
            self.assertEqual(row[1].getPlainText(), "A1, B1")
            self.assertEqual(row[4].getPlainText(), "0.112")
            cv = f"{47.14045207910317:.{result.settings.report.percentage_decimals}f}"
            self.assertEqual(row[5].getPlainText(), cv)
            self.assertEqual(row[5].style.textColor, colors.firebrick)
            self.assertEqual(row[-4].getPlainText(), "0.112")
            self.assertEqual(row[-3].getPlainText(), cv)
            with pymupdf.open(path) as document:
                self.assertIn(cv, document[0].get_text())

    def test_standard_cv_matches_nominal_concentration_and_sample_type(self) -> None:
        result = self.make_fitted_result()
        frame = well_table(result.samples).copy()
        settings = with_validation(result.settings, standards={"25": {"extinction": {"response_scale": "blank_corrected", "min": 0.1}}})
        result = AssayAnalyzer(AssayData(samples_from_records(frame.to_dict('records')), result.metadata), settings).analyze()
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(result, analyst="operator")
            writer._initialize_report(Path(directory) / "report.pdf")
            table = quality_control_table(writer)
            standard = sample_table(result.samples).loc[sample_table(result.samples)["sample_concentration"].eq(25)].iloc[0]
            row = next(row for row in table._cellvalues if row[0].getPlainText() == "Standaard 25")
            self.assertEqual(row[1].getPlainText(), "A7, A8")
            self.assertEqual(row[5].getPlainText(), f"{standard['extinction_cv_percent']:.{result.settings.report.percentage_decimals}f}")

    def test_calibration_table_uses_unrestricted_sample_means_and_zero_has_no_recovery(self) -> None:
        result = self.make_fitted_result()
        result = replace(result, settings=replace(result.settings, report=ReportSettings(result_decimals=3, calibration_decimals=3)))
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(result, analyst="operator")
            with patch("elisaflow.calibration.curves.fit_curve", side_effect=AssertionError("Refitting")):
                path = writer.write(Path(directory) / "report.pdf")
            self.assertTrue(path.read_bytes().startswith(b"%PDF-"))
            standards = select_samples(result.samples, SampleType.STANDARD)
            table = calibration_table(writer.context, standards, writer.content_width * 0.45, scale=writer.measurements[0].calculation_input_scale)
            cells = [[cell.getPlainText() for cell in row] for row in table._cellvalues]
            self.assertEqual(cells[0], ["Standaard Conc", "Ext", "Berekende Conc", "Terugvinding%"])
            self.assertEqual(len(cells), 6)
            self.assertEqual(cells[1][0], "0.000")
            self.assertEqual(cells[1][-1], "-")
            low = standards[1]
            self.assertTrue(np.isnan(low.result.mean))
            self.assertEqual(cells[2][2], f"{low.back_calculated_concentration.mean:.3f}")
            self.assertEqual(cells[2][3], f"{low.recovery_percent:.{result.settings.report.percentage_decimals}f}")

    def test_formula_omits_calibration_section(self) -> None:
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(self.make_result(CalculationMethod.FORMULA), analyst="operator")
            writer._initialize_report(Path(directory) / "report.pdf")
            self.assertEqual(calibration_section(writer.context, writer.measurements[0]), [])

    def test_concentration_units_reach_column_headers_quality_control_and_the_graph(self) -> None:
        base = self.make_fitted_result()
        for unit in ("U/ml", "U<ref>/ml", "%"):
            with self.subTest(unit=unit), TemporaryDirectory() as directory:
                settings = with_validation(
                    replace(
                        base.settings,
                        result=replace(base.settings.result, unit=unit),
                        report=ReportSettings(result_decimals=3, calibration_decimals=3),
                        standards={25: SampleSettings.from_dict({})},
                    ),
                    standards={25: {"extinction": {"response_scale": "blank_corrected", "min": 0.1}, "result": {"min": 20, "max": 30}}},
                )
                result = AssayAnalyzer(AssayData(samples_from_records(well_table(base.samples).to_dict('records')), base.metadata), settings).analyze()
                writer = PdfReportWriter(result, analyst="operator")
                with patch("elisaflow.reporting.calibration.plot_calibration", wraps=plot_calibration) as plot:
                    path = writer.write(Path(directory) / "report.pdf")
                self.assertEqual(plot.call_args.kwargs["concentration_unit"], unit)
                standards = select_samples(result.samples, SampleType.STANDARD)
                table = calibration_table(writer.context, standards, writer.content_width * 0.45, scale=writer.measurements[0].calculation_input_scale)
                self.assertEqual(
                    [cell.getPlainText() for cell in table._cellvalues[0]],
                    [
                        f"Standaard Conc [{unit}]",
                        "Ext",
                        f"Berekende Conc [{unit}]",
                        "Terugvinding%",
                    ],
                )
                qc = quality_control_table(writer)._cellvalues
                result_row = next(row for row in qc if "Conc" in row[2].getPlainText())
                self.assertEqual(result_row[2].getPlainText(), "20.000 ≤ Conc ≤ 30.000")
                self.assertEqual(result_row[-2].getPlainText(), unit)
                extinction_row = next(row for row in qc if "Ext" in row[2].getPlainText())
                self.assertEqual(extinction_row[-2].getPlainText(), "")
                with pymupdf.open(path) as document:
                    text = " ".join(" ".join(page.get_text().split()) for page in document)
                    self.assertIn(f"Conc [{unit}]", text)
                    self.assertIn("Eenheid", text)
                    self.assertIn(f"Standaard Conc [{unit}]", text)
                    self.assertIn(f"Berekende Conc [{unit}]", text)
                    for page in document:
                        for word in page.get_text("words"):
                            self.assertTrue(0 <= word[0] < word[2] <= page.rect.width)
                            self.assertTrue(0 <= word[1] < word[3] <= page.rect.height)

    def test_large_fitted_report_continues_on_additional_pages(self) -> None:
        result = self.make_fitted_result()
        standards = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].ne("Unknown")]
        patients = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")]
        patients = pd.concat([patients] * 160, ignore_index=True)
        patients["sample_name"] = [f"Patient {index:03d}" for index in range(160)]
        result = replace(
            result,
            selected=replace(result.selected, samples=samples_with_statistics(
                pd.concat([standards, patients], ignore_index=True).to_dict("records"), source_samples=result.samples,
            )),
        )
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(result, analyst="operator")
            pages = []
            draw_footer = writer._draw_footer

            def record_page(canvas, document):
                pages.append(document.page)
                draw_footer(canvas, document)

            with patch.object(writer, "_draw_footer", side_effect=record_page):
                writer.write(Path(directory) / "report.pdf")
            self.assertGreater(len(pages), 1)
            self.assertEqual(pages, list(range(1, len(pages) + 1)))

    def test_patient_results_read_down_columns_for_uneven_counts(self) -> None:
        result = self.make_fitted_result()
        others = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].ne("Unknown")]
        patient = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")]
        for count in (0, 1, 2, 4, 8, 13):
            with self.subTest(count=count), TemporaryDirectory() as directory:
                patients = pd.concat([patient] * count, ignore_index=True) if count else patient.iloc[:0].copy()
                expected = [f"Patient{index:03d}" for index in range(count)]
                patients["sample_name"] = expected
                samples = pd.concat([others, patients], ignore_index=True)
                writer = PdfReportWriter(
                    replace(
                        result,
                        selected=replace(result.selected, samples=samples_with_statistics(samples.to_dict("records"), source_samples=result.samples)),
                    )
                , analyst="operator")
                path = writer.write(Path(directory) / "report.pdf")
                with pymupdf.open(path) as document:
                    observed = []
                    for page in document:
                        groups = self.patient_names_by_column(page, expected)
                        observed.extend(name for group in groups for name in group)
                        if count:
                            self.assertEqual(len(groups), 1)
                    self.assertEqual(observed, expected)
                    if not count:
                        self.assertNotIn("Patiëntresultaten", document[0].get_text())

    def test_dual_wavelength_patients_fill_columns_in_reading_order(self) -> None:
        result = AssayAnalyzer(dual_data(), dual_settings()).analyze()
        others = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].ne("Unknown")]
        patient = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")].iloc[:1]
        patients = pd.concat([patient] * 88, ignore_index=True)
        expected = [f"P{index:03d}" for index in range(len(patients))]
        patients["sample_name"] = expected
        patients[["extinction_cv_percent", "result_cv_percent"]] = np.nan
        result = replace(result, selected=replace(result.selected, samples=samples_with_statistics(
            pd.concat([others, patients], ignore_index=True).to_dict("records"), source_samples=result.samples,
        )))
        with TemporaryDirectory() as directory:
            writer = PdfReportWriter(result, analyst="operator")
            path = writer.write(Path(directory) / "report.pdf")
            with pymupdf.open(path) as document:
                observed = []
                most_columns = 0
                for page in document:
                    groups = self.patient_names_by_column(page, expected)
                    if groups:
                        self.assertLessEqual(len(groups), 3)
                        most_columns = max(most_columns, len(groups))
                        positions = {word[4]: word for word in page.get_text("words") if word[4] in expected}
                        page_bottom = page.rect.height - writer.document.bottomMargin - writer._FRAME_PADDING
                        for group in groups[:-1]:
                            row_spacing = positions[group[-1]][1] - positions[group[-2]][1]
                            self.assertLess(page_bottom - positions[group[-1]][3], row_spacing + 6)
                        observed.extend(name for group in groups for name in group)
                self.assertGreater(most_columns, 1)
                self.assertEqual(observed, expected)

    def test_patient_column_is_filled_to_the_page_bottom_before_moving_right(self) -> None:
        result = self.make_fitted_result()
        others = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].ne("Unknown")]
        patient = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")]
        expected = [f"Patient{index:03d}" for index in range(160)]
        patients = pd.concat([patient] * len(expected), ignore_index=True)
        patients["sample_name"] = expected
        samples = pd.concat([others, patients], ignore_index=True)
        writer = PdfReportWriter(
            replace(
                result,
                selected=replace(result.selected, samples=samples_with_statistics(samples.to_dict("records"), source_samples=result.samples)),
            )
        , analyst="operator")
        with TemporaryDirectory() as directory:
            path = writer.write(Path(directory) / "report.pdf")
            with pymupdf.open(path) as document:
                filled_columns = 0
                for page in document:
                    groups = self.patient_names_by_column(page, expected)
                    positions = {word[4]: word for word in page.get_text("words") if word[4] in expected}
                    page_bottom = page.rect.height - writer.document.bottomMargin - writer._FRAME_PADDING
                    for group in groups[:-1]:
                        self.assertGreater(len(group), 1)
                        row_spacing = positions[group[-1]][1] - positions[group[-2]][1]
                        self.assertLess(page_bottom - positions[group[-1]][3], row_spacing + 6)
                        filled_columns += 1
                self.assertGreater(filled_columns, 2)

    def test_patient_columns_continue_in_sequence_across_page_splits(self) -> None:
        result = self.make_fitted_result()
        others = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].ne("Unknown")]
        patient = sample_table(result.samples).loc[sample_table(result.samples)["sample_type"].eq("Unknown")]
        patients = pd.concat([patient] * 160, ignore_index=True)
        expected = [f"Patient{index:03d}" for index in range(160)]
        patients["sample_name"] = expected
        patients.loc[20, "sample_name"] += " with a long name and identifier ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        patients.loc[90, "sample_name"] += " with another long name and identifier ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        with TemporaryDirectory() as directory:
            samples = pd.concat([others, patients], ignore_index=True)
            writer = PdfReportWriter(
                replace(
                    result,
                    selected=replace(result.selected, samples=samples_with_statistics(samples.to_dict("records"), source_samples=result.samples)),
                )
            , analyst="operator")
            path = writer.write(Path(directory) / "report.pdf")
            with pymupdf.open(path) as document:
                observed = []
                for page in document:
                    groups = self.patient_names_by_column(page, expected)
                    observed.extend(name for group in groups for name in group)
                    for word in page.get_text("words"):
                        self.assertTrue(0 <= word[0] < word[2] <= page.rect.width)
                        self.assertTrue(0 <= word[1] < word[3] <= page.rect.height)
                self.assertGreater(len(document), 1)
                self.assertIn("Patient000", document[0].get_text())
                self.assertEqual(observed, expected)


if __name__ == "__main__":
    unittest.main()
