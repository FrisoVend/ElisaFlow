from __future__ import annotations

import json
import os
import stat
import unittest
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from fixtures import add_measurement_steps

from elisaflow.assays.access_control import AccessControl
from elisaflow.assays.management import AssayManagement
from elisaflow.assays.protocol_samples import ProtocolSamples
from elisaflow.configuration.assay_settings import AssaySettings, ReplicateFilterSettings, ReportSettings, ResultSettings

try:
    from PySide6.QtCore import QPoint, Qt
    from PySide6.QtGui import QColor, QPalette
    from PySide6.QtTest import QTest
    from PySide6.QtWidgets import QAbstractSpinBox, QApplication, QDialogButtonBox, QLabel, QLineEdit, QPushButton

    from elisaflow.ui.assay_editor import AssayEditor
    from elisaflow.ui.assay_settings import AssaySettingsWidget
    from elisaflow.ui.editors.formula import FormulaDialog, FormulaEditor
except ModuleNotFoundError as error:
    if error.name != "PySide6":
        raise
    QApplication = None


def write_protocol(path: Path, samples: list[tuple[str, str | None, object]], wavelengths=(405,)) -> None:
    """Create a real SKAX archive so tests exercise the production reader."""
    session = ET.Element("Session")
    add_measurement_steps(session, wavelengths)
    plates = ET.SubElement(ET.SubElement(session, "PlateLayout"), "Plates")
    plate = ET.SubElement(plates, "Plate")
    ET.SubElement(plate, "Id").text = "plate"
    ET.SubElement(plate, "Name").text = "plate"
    wells = ET.SubElement(plate, "Wells")
    for index, (sample_type, name, concentration) in enumerate(samples):
        well = ET.SubElement(wells, "Well")
        ET.SubElement(well, "WellLocationX").text = str(index // 8)
        ET.SubElement(well, "WellLocationY").text = str(index % 8)
        sample = ET.SubElement(
            ET.SubElement(ET.SubElement(well, "Replicate"), "Contents"),
            "Sample",
        )
        ET.SubElement(sample, "Id").text = f"sample-{index}"
        ET.SubElement(sample, "Name").text = name
        ET.SubElement(sample, "SampleType").text = sample_type
        if concentration is not None:
            ET.SubElement(sample, "Concentration").text = str(concentration)
    with ZipFile(path, "w") as archive:
        archive.writestr("Header.xml", "<Header />")
        archive.writestr("Session.xml", ET.tostring(session))


class TestAccess(AccessControl):
    def __init__(self, editable: bool = True) -> None:
        self.editable = editable

    def can_edit_assays(self) -> bool:
        return self.editable

    def has_password(self) -> bool:
        return True

    @contextmanager
    def edit_assays(self, entered_password: str):
        self.require_assay_editor()
        yield


class ProtocolFixture:
    def setUp(self) -> None:
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.path = self.directory / "test.skax"
        self.samples = [
            ("Blank", "Blanco", 0),
            ("Control", "Grens", 0),
            ("Control", "Positief", 0),
            ("Control", "Positief", 0),
            ("Standard", "Standaard1", "2.5000"),
            ("Standard", "Standaard1", "2.5"),
            ("Standard", "Standaard2", 0),
            ("Unknown", "Patient", 10),
        ]
        write_protocol(self.path, self.samples)
        self.config = {
            "lot_number": "kit lot",
            "controls": {
                "Positief": {"lot_number": "control lot"}
            },
            "standards": {"2.5000": {}},
            "validation": {"result": {"controls": {"Positief": {"result": {"min": 1.2}}}, "standards": {"2.5000": {"extinction": {"response_scale": "blank_corrected", "min": 0.1}}}}},
            "calculation": {
                "response_scale": "blank_corrected",
                "method": "formula",
                "expression": "sample / grens",
                "references": {"grens": "Grens"},
                "control_overrides": {"Positief": {"expression": "sample / grens", "response_scale": "blank_corrected"}},
            },
            "result": {"label": "Ratio"},
        }
        self.configs = self.directory / "configs"
        self.configs.mkdir()
        self.source = self.configs / "library.json"
        self.source.write_text(json.dumps({"test.skax": self.config}), encoding="utf-8")
        self.access = TestAccess()
        self.assays = AssayManagement(self.configs, self.access, protocols_dir=self.directory)
        self.settings = self.assays.load_assay("test.skax")

    def incoming_protocol(self, samples: list | None = None, name: str = "new.skax", wavelengths=(405,)) -> Path:
        incoming = self.directory / "incoming"
        incoming.mkdir(exist_ok=True)
        path = incoming / name
        write_protocol(path, self.samples if samples is None else samples, wavelengths)
        return path


class ProtocolSamplesTests(ProtocolFixture, unittest.TestCase):
    def test_protocol_deduplicates_replicates_and_groups_standards_by_concentration(self) -> None:
        write_protocol(self.path, self.samples + [("Standard", "Another standard", 2.5)])
        samples = self.assays.load_protocol_samples("test.skax")
        self.assertEqual(samples.controls, ("Grens", "Positief"))
        self.assertEqual(samples.standards, {2.5: ("Standaard1", "Another standard"), 0: ("Standaard2",)})
        samples.validate_settings(self.settings)

    def test_wavelengths_use_filter_metadata_in_protocol_order(self) -> None:
        write_protocol(self.path, self.samples, (620, 450, 405))
        protocol = self.assays.load_protocol_samples("test.skax")
        self.assertEqual(protocol.wavelengths_nm, (620, 450, 405))
        with self.assertRaisesRegex(ValueError, 'precies 1 analysegolflengte'):
            protocol.validate_measurements(self.settings)

    def test_missing_or_duplicate_skax_wavelengths_have_no_config_fallback(self) -> None:
        for wavelengths in ((), (405, 405)):
            with self.subTest(wavelengths=wavelengths):
                write_protocol(self.path, self.samples, wavelengths)
                with self.assertRaisesRegex(ValueError, 'verschillende ruwe absorbantiegolflengten'):
                    self.assays.load_protocol_samples("test.skax")

    def test_missing_and_unreadable_protocols_fail_clearly(self) -> None:
        with self.assertRaises(FileNotFoundError):
            ProtocolSamples.from_file(self.directory / "missing.skax")
        self.path.write_bytes(b"not a zip archive")
        with self.assertRaisesRegex(ValueError, 'kon niet worden gelezen'):
            self.assays.load_protocol_samples("test.skax")

    def test_invalid_sample_identities_are_rejected(self) -> None:
        for sample_type in ("Control", "Standard"):
            for name in (None, "", " "):
                with self.subTest(type=sample_type, name=name):
                    write_protocol(self.path, [(sample_type, name, 1)])
                    with self.assertRaisesRegex(ValueError, 'moeten niet-lege teksten zijn'):
                        ProtocolSamples.from_file(self.path)
        for concentration in (None, "bad", -1, "NaN", "inf"):
            with self.subTest(concentration=concentration):
                write_protocol(self.path, [("Standard", "Standaard", concentration)])
                with self.assertRaisesRegex(ValueError, 'eindige, niet-negatieve concentratie'):
                    ProtocolSamples.from_file(self.path)

    def test_save_rejects_sample_settings_absent_from_protocol_without_writing(self) -> None:
        original = self.source.read_bytes()
        for change in ({"controls": {"Invented": {}}}, {"standards": {"3": {}}}):
            with self.subTest(change=change):
                draft = AssaySettings.from_dict("test.skax", self.config | change)
                with self.assertRaisesRegex(ValueError, 'ontbreken in het SKAX-protocol'):
                    self.assays.save_assay(draft, expected_settings=self.settings)
                self.assertEqual(self.source.read_bytes(), original)
        self.assertFalse((self.configs / ".history").exists())

    def test_control_results_require_exact_protocol_controls(self) -> None:
        config = self.config | {"calculation": self.config["calculation"] | {"control_results": {"Invented": {"label": "Controle-index", "unit": "index"}}}}
        with self.assertRaisesRegex(ValueError, "Ingestelde controles ontbreken.*Invented"):
            self.assays.validate_assay("test.skax", config)

    def test_all_control_reference_locations_are_validated(self) -> None:
        cases = [
            {
                "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / invented", "references": {"invented": "Invented"}},
                "result": {"label": "Ratio"},
            },
            {
                "calculation": {
                    "response_scale": "blank_corrected",
                    "method": "formula",
                    "expression": "sample / (grens - invented)",
                    "references": {"grens": "Grens", "invented": "Invented"},
                },
                "result": {"label": "Ratio"},
            },
            {
                "calculation": {
                    "response_scale": "blank_corrected",
                    "method": "formula",
                    "expression": "sample",
                    "references": {"ref": "Invented"},
                    "control_overrides": {"Positief": {"expression": "sample / ref", "response_scale": "blank_corrected"}},
                }
            },
            {
                "controls": {},
                "calculation": {
                    "response_scale": "raw",
                    "method": "calibration",
                    "fit_model": "4pl",
                    "preprocessing": {"expression": "(sample - invented)", "references": {"invented": "Invented"}},
                    "parameters": {"A": 2, "B": -1, "C": 1, "D": 0},
                },
                "result": {"label": "Conc"},
            },
        ]
        for change in cases:
            with (
                self.subTest(change=change),
                self.assertRaisesRegex(ValueError, 'verwijst naar een controle die ontbreekt'),
            ):
                self.assays.validate_assay("test.skax", self.config | change)
        self.assays.validate_assay("test.skax", self.config)

    def test_each_wavelength_preprocessing_checks_its_own_control_references(self) -> None:
        write_protocol(self.path, self.samples, (405, 450))
        curves = {
            str(wavelength): {
                "response_scale": "raw",
                "fit_model": "4pl",
                "preprocessing": {"expression": "sample / ref", "references": {"ref": control}},
                "parameters": {"A": 2, "B": -1, "C": 1, "D": 0},
            }
            for wavelength, control in ((405, "Grens"), (450, "Positief"))
        }
        config = self.config | {
            "measurement": {"mode": "switch_on_concentration", "cutoff_concentration": 100, "preferred_wavelength_nm": 450},
            "calculation": {"method": "calibration", "wavelengths": curves},
            "validation": {},
        }
        self.assays.validate_assay("test.skax", config)
        for wavelength in (405, 450):
            formula = curves[str(wavelength)]["preprocessing"]
            control = formula["references"]["ref"]
            formula["references"]["ref"] = "Invented"
            with self.subTest(wavelength=wavelength), self.assertRaisesRegex(ValueError, f"Kalibratievoorbewerking bij {wavelength} nm.*Invented"):
                self.assays.validate_assay("test.skax", config)
            formula["references"]["ref"] = control

    def test_save_rechecks_protocol_after_it_changes(self) -> None:
        original = self.source.read_bytes()
        write_protocol(self.path, [sample for sample in self.samples if sample[1] != "Positief"])
        with self.assertRaisesRegex(ValueError, 'Ingestelde controles ontbreken'):
            self.assays.save_assay(replace(self.settings, lot_number="new lot"), expected_settings=self.settings)
        self.assertEqual(self.source.read_bytes(), original)


class NewAssayTests(ProtocolFixture, unittest.TestCase):
    def test_prepare_is_read_only_and_save_imports_reloadable_assay(self) -> None:
        path = self.incoming_protocol()
        original_library = self.source.read_bytes()
        original_protocol = path.read_bytes()
        samples = self.assays.prepare_new_assay(path)
        self.assertEqual(samples.controls, ("Grens", "Positief"))
        self.assertFalse((self.directory / "new.skax").exists())
        self.assertEqual(list(self.configs.iterdir()), [self.source])
        settings = AssaySettings.from_dict(path.name, self.config)
        self.assays.create_assay(settings, protocol_path=path)
        self.assertEqual(self.assays.load_assay(path.name), settings)
        self.assertEqual(self.source.read_bytes(), original_library)
        self.assertEqual(path.read_bytes(), original_protocol)
        self.assertEqual((self.directory / path.name).read_bytes(), original_protocol)
        saved = json.loads((self.configs / "new.json").read_text(encoding="utf-8"))
        self.assertEqual(saved, {path.name: settings.to_dict()})
        self.assertFalse((self.configs / "new.json").stat().st_mode & stat.S_IWUSR)
        self.assertTrue((self.directory / path.name).stat().st_mode & stat.S_IWUSR)
        if os.name == "nt":
            self.assertEqual(list(self.configs.glob("*.json.lock")), [])
        reopened = AssayManagement(self.configs, self.access, protocols_dir=self.directory)
        self.assertIn(path.name, reopened.list_assays())
        self.assertEqual(reopened.load_protocol_samples(path.name), samples)

    def test_failed_configuration_protection_cleans_up_import_and_staged_file(self) -> None:
        path = self.incoming_protocol()
        settings = AssaySettings.from_dict(path.name, self.config)
        with patch("elisaflow.assays.management.make_read_only", side_effect=PermissionError("protection error")):
            with self.assertRaisesRegex(PermissionError, "protection error"):
                self.assays.create_assay(settings, protocol_path=path)
        self.assertFalse((self.configs / "new.json").exists())
        self.assertFalse((self.directory / path.name).exists())
        self.assertEqual(list(self.configs.glob("*.tmp")), [])

    def test_edit_permission_is_required_to_prepare_and_to_create(self) -> None:
        path = self.incoming_protocol()
        self.access.editable = False
        settings = AssaySettings.from_dict(path.name, self.config)
        with self.assertRaises(PermissionError):
            self.assays.prepare_new_assay(path)
        with self.assertRaises(PermissionError):
            self.assays.create_assay(settings, protocol_path=path)
        self.assertEqual(list(self.configs.iterdir()), [self.source])
        self.assertFalse((self.directory / path.name).exists())

    def test_registered_names_and_configuration_collisions_are_preserved(self) -> None:
        path = self.incoming_protocol(name="test.skax")
        with self.assertRaisesRegex(ValueError, 'al een assay geregistreerd'):
            self.assays.prepare_new_assay(path)
        path = self.incoming_protocol()
        collision = self.configs / "new.json"
        contents = b'{"unrelated.skax": {"lot_number": "lot", "calculation": {"method": "formula", "expression": "sample", "response_scale": "raw"}}}'
        collision.write_bytes(contents)
        with self.assertRaisesRegex(FileExistsError, 'bestaat al een assayconfiguratie'):
            self.assays.prepare_new_assay(path)
        self.assertEqual(collision.read_bytes(), contents)

    def test_existing_protocol_is_reused_only_when_identical(self) -> None:
        path = self.incoming_protocol()
        destination = self.directory / path.name
        destination.write_bytes(b"a different protocol")
        with self.assertRaisesRegex(FileExistsError, 'bestaat al een ander protocol'):
            self.assays.prepare_new_assay(path)
        self.assertEqual(destination.read_bytes(), b"a different protocol")
        destination.write_bytes(path.read_bytes())
        settings = AssaySettings.from_dict(path.name, self.config)
        self.assays.create_assay(settings, protocol_path=destination)
        self.assertEqual(self.assays.load_assay(path.name), settings)
        self.assertEqual(destination.read_bytes(), path.read_bytes())

    def test_invalid_protocols_and_settings_do_not_create_an_assay(self) -> None:
        path = self.incoming_protocol([("Unknown", "Patient", 1)])
        with self.assertRaisesRegex(ValueError, 'controles of standaarden bevatten'):
            self.assays.prepare_new_assay(path)
        path.write_bytes(b"invalid archive")
        with self.assertRaisesRegex(ValueError, 'kon niet worden gelezen'):
            self.assays.prepare_new_assay(path)
        path = self.incoming_protocol()
        settings = AssaySettings.from_dict(
            path.name,
            self.config
            | {
                "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / invented", "references": {"invented": "Invented"}},
                "result": {"label": "Ratio"},
            },
        )
        with self.assertRaisesRegex(ValueError, 'verwijst naar een controle die ontbreekt'):
            self.assays.create_assay(settings, protocol_path=path)
        self.assertEqual(list(self.configs.iterdir()), [self.source])
        self.assertFalse((self.directory / path.name).exists())

    def test_failed_configuration_publish_removes_only_the_new_import(self) -> None:
        from elisaflow.assays.management import _publish_new_file

        path = self.incoming_protocol()
        settings = AssaySettings.from_dict(path.name, self.config)
        for existing in (False, True):
            with self.subTest(existing=existing):
                destination = self.directory / path.name
                if existing:
                    destination.write_bytes(path.read_bytes())

                def fail_config(target: Path, contents: bytes, *, read_only: bool = False) -> None:
                    if target.suffix == ".json":
                        raise OSError("disk error")
                    _publish_new_file(target, contents, read_only=read_only)

                with patch("elisaflow.assays.management._publish_new_file", side_effect=fail_config):
                    with self.assertRaisesRegex(OSError, "disk error"):
                        self.assays.create_assay(settings, protocol_path=path)
                self.assertEqual(destination.exists(), existing)
                self.assertFalse((self.configs / "new.json").exists())
                if os.name == "nt":
                    self.assertEqual(list(self.configs.glob("*.json.lock")), [])
                self.assertEqual(list(self.directory.glob("*.tmp")), [])
        self.assays.create_assay(settings, protocol_path=path)
        self.assertEqual(self.assays.load_assay(path.name), settings)

    def test_concurrent_configuration_creation_is_never_overwritten(self) -> None:
        from elisaflow.assays.management import _publish_new_file

        path = self.incoming_protocol()
        settings = AssaySettings.from_dict(path.name, self.config)
        collision = self.configs / "new.json"
        external = b'{"external.skax": {"lot_number": "lot", "calculation": {"method": "formula", "expression": "sample", "response_scale": "raw"}}}'

        def create_concurrently(target: Path, contents: bytes, *, read_only: bool = False) -> None:
            if target.suffix == ".json":
                collision.write_bytes(external)
            _publish_new_file(target, contents, read_only=read_only)

        with patch("elisaflow.assays.management._publish_new_file", side_effect=create_concurrently):
            with self.assertRaises(FileExistsError):
                self.assays.create_assay(settings, protocol_path=path)
        self.assertEqual(collision.read_bytes(), external)
        self.assertFalse((self.directory / path.name).exists())
        self.assertEqual(list(self.configs.glob("*.tmp")), [])

    def test_changed_source_protocol_is_not_imported_with_stale_validation(self) -> None:
        path = self.incoming_protocol()
        settings = AssaySettings.from_dict(path.name, self.config)
        prepare = self.assays.prepare_new_assay

        def change_after_read(source: Path) -> ProtocolSamples:
            samples = prepare(source)
            write_protocol(source, [("Control", "Grens", 0)])
            return samples

        with patch.object(self.assays, "prepare_new_assay", side_effect=change_after_read):
            with self.assertRaisesRegex(ValueError, 'geselecteerde protocol is gewijzigd'):
                self.assays.create_assay(settings, protocol_path=path)
        self.assertFalse((self.configs / "new.json").exists())
        self.assertFalse((self.directory / path.name).exists())


@unittest.skipIf(QApplication is None, "Install the optional PySide6 UI dependency")
class AssayEditorTests(ProtocolFixture, unittest.TestCase):
    def setUp(self) -> None:
        super().setUp()
        self.enterContext(patch("elisaflow.ui.assay_settings.QInputDialog.getText", return_value=("editor password", True)))

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        environment = patch.dict(os.environ, {"QT_QPA_PLATFORM": "offscreen"})
        environment.start()
        cls.addClassCleanup(environment.stop)
        cls.application = QApplication.instance() or QApplication([])

    def open_editor(self, settings: AssaySettings | None = None) -> AssayEditor:
        editor = AssayEditor(settings or self.settings, self.assays, self.access)
        self.addCleanup(editor.deleteLater)
        return editor

    def edit_formula(self, editor: FormulaEditor, expression: str, control: str | None = None) -> None:
        def save(dialog):
            dialog.input.expression.setText(expression)
            if control is not None:
                field = dialog.input
                field.operand.setCurrentIndex(field.operand.findData(control))
                field.insert_button.click()
            dialog._save()
            self.assertEqual(dialog.result(), dialog.DialogCode.Accepted)
            return dialog.result()

        with patch.object(FormulaDialog, "exec", save):
            editor._edit()

    def open_switching_editor(self, reference=None, *, fixed=False) -> AssayEditor:
        config = self.settings.to_dict()
        write_protocol(self.path, self.samples, (405, 450) + ((reference,) if reference is not None else ()))
        config["measurement"] = {"mode": "switch_on_concentration", "cutoff_concentration": 1000.0, "preferred_wavelength_nm": 450}
        if reference is not None:
            config["measurement"]["reference_wavelength_nm"] = reference
        config["calculation"] = {"method": "calibration", "wavelengths": {
            "405": {"response_scale": "raw", "fit_model": "quadratic", "weighting_exponent": 1.5},
            "450": {"response_scale": "raw", "fit_model": "4pl", "weighting_exponent": 2.0},
        }}
        if fixed:
            config["calculation"]["concentration_range"] = [1, 300]
            for wavelength, parameters in (
                ("405", {"A": 0.03, "B": 0.003, "C": 0.000005}),
                ("450", {"A": 4.5, "B": -1.2, "C": 30, "D": 0.08}),
            ):
                curve = config["calculation"]["wavelengths"][wavelength]
                curve.pop("weighting_exponent")
                curve["parameters"] = parameters
        config["validation"] = {"wavelengths": {
            "405": {"blank": {"extinction": {"response_scale": "raw", "max": 0.15}}},
            "450": {"blank": {"extinction": {"response_scale": "raw", "max": 0.2}}},
        }}
        self.source.chmod(self.source.stat().st_mode | stat.S_IWUSR)
        self.source.write_text(json.dumps({"test.skax": config}), encoding="utf-8")
        self.assays = AssayManagement(self.configs, self.access, protocols_dir=self.directory)
        return self.open_editor(self.assays.load_assay("test.skax"))

    def test_sample_categories_and_calculation_share_the_same_control_choices(self) -> None:
        editor = self.open_editor()
        self.assertEqual(
            [editor.navigation.item(index).text() for index in range(editor.navigation.count())],
            ["Meting", "Monsterinstellingen", "Berekening", "QC-criteria", "Rapport"],
        )
        for index, (title, page) in enumerate((
            ("Meting", editor.measurement),
            ("Monsterinstellingen", editor.sample_tabs),
            ("Berekening", editor._calculation_stack),
            ("QC-criteria", editor.validation),
            ("Rapport", editor.report),
        )):
            editor.navigation.setCurrentRow(index)
            self.assertEqual(editor.page_title.text(), title)
            current = editor.pages.currentWidget()
            self.assertTrue(current is page or current.isAncestorOf(page))
        editor.navigation.setCurrentRow(1)
        self.assertIs(editor.pages.currentWidget(), editor.sample_tabs)
        self.assertEqual(
            [editor.sample_tabs.tabText(index) for index in range(editor.sample_tabs.count())],
            ["Patiënten", "Controles", "Blanco", "Standaarden"],
        )
        self.assertIs(editor.sample_tabs.widget(0).widget(), editor.patients)
        self.assertEqual(editor.patients.findChildren(QLineEdit), list(editor.patients.replicate_filter.fields.values()))
        editor.navigation.setCurrentRow(2)
        self.assertIs(editor.pages.currentWidget().widget(), editor._calculation_stack)
        calculation = editor.calculation
        self.assertEqual([calculation.tabs.tabText(index) for index in range(calculation.tabs.count())], ["Patiënten", "Controles"])
        self.assertTrue(calculation.tabs.widget(0).isAncestorOf(calculation.formula))
        self.assertTrue(calculation.tabs.widget(0).isAncestorOf(calculation.result_label))
        self.assertTrue(calculation.tabs.widget(0).isAncestorOf(calculation.result_unit))
        self.assertEqual(tuple(calculation.control_tabs.tabText(i) for i in range(calculation.control_tabs.count())), ("Grens", "Positief"))
        self.assertEqual(tuple(calculation.overrides), tuple(editor.controls.editors))
        for name in ("Grens", "Positief"):
            calculation.control_tabs.setCurrentWidget(calculation.overrides[name])
            self.assertIs(calculation.control_tabs.currentWidget(), calculation.overrides[name])

    def test_read_only_calculation_tabs_allow_browsing_without_enabling_changes(self) -> None:
        original = self.source.read_bytes()
        for read_only, editable in ((True, True), (False, False)):
            with self.subTest(read_only=read_only, editable=editable):
                self.access.editable = editable
                editor = AssayEditor(self.settings, self.assays, self.access, read_only=read_only)
                self.addCleanup(editor.deleteLater)
                editor.show()
                editor.navigation.setCurrentRow(2)
                self.application.processEvents()
                calculation = editor.calculation
                QTest.mouseClick(calculation.tabs.tabBar(), Qt.MouseButton.LeftButton, pos=calculation.tabs.tabBar().tabRect(1).center())
                self.assertEqual(calculation.tabs.currentIndex(), 1)
                for index, name in enumerate(("Grens", "Positief")):
                    selector = calculation.control_tabs
                    QTest.mouseClick(selector.tabBar(), Qt.MouseButton.LeftButton, pos=selector.tabBar().tabRect(index).center())
                    self.assertEqual(selector.tabText(selector.currentIndex()), name)
                    control = calculation.overrides[name]
                    self.assertIs(selector.currentWidget(), control)
                    self.assertFalse(control.use_sample_result.isEnabled())
                    self.assertFalse(control.result_label.isEnabled())
                    self.assertFalse(control.result_unit.isEnabled())
                    self.assertFalse(control.formula.edit_button.isEnabled())
                    self.assertFalse(control.formula.clear_button.isEnabled())
                    checked = control.use_sample_result.isChecked()
                    QTest.mouseClick(control.use_sample_result, Qt.MouseButton.LeftButton)
                    self.assertEqual(control.use_sample_result.isChecked(), checked)
                QTest.mouseClick(calculation.tabs.tabBar(), Qt.MouseButton.LeftButton, pos=calculation.tabs.tabBar().tabRect(0).center())
                self.assertEqual(calculation.tabs.currentIndex(), 0)
                self.assertFalse(calculation.result_label.isEnabled())
                self.assertFalse(calculation.result_unit.isEnabled())
                self.assertFalse(calculation.formula.edit_button.isEnabled())
                self.assertEqual(self.source.read_bytes(), original)
                editor.close()

    def test_skax_wavelengths_cannot_be_edited_with_or_without_edit_permission(self) -> None:
        for editable in (True, False):
            with self.subTest(editable=editable):
                self.access.editable = editable
                editor = self.open_switching_editor()
                self.assertIsInstance(editor.measurement.wavelengths_nm, QLabel)
                self.assertEqual(editor.measurement.wavelengths_nm.text(), "405, 450")
                parameter_fields = [editor.measurement.cutoff_concentration]
                for saturation in editor.measurement.saturation.values():
                    parameter_fields.extend((saturation.slope, saturation.intercept))
                    self.assertFalse(saturation.source_wavelength_nm.isEditable())
                self.assertCountEqual(editor.measurement.findChildren(QLineEdit), parameter_fields)
                self.assertFalse(editor.measurement.preferred_wavelength_nm.isEditable())
                self.assertFalse(editor.measurement.reference_wavelength_nm.isEditable())
                self.assertEqual(editor.measurement.preferred_wavelength_nm.count(), 2)
                if editable:
                    fits = editor.calculation.to_config()["calculation"]["wavelengths"]
                    with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
                        editor._save()
                    warning.assert_not_called()
                    saved = self.assays.load_assay("test.skax")
                    self.assertNotIn("wavelengths_nm", saved.to_dict()["measurement"])
                    reopened = self.open_editor(saved)
                    self.assertEqual(reopened.measurement.wavelengths_nm.text(), "405, 450")
                    self.assertEqual(reopened.calculation.to_config()["calculation"]["wavelengths"], fits)

    def test_save_rechecks_wavelengths_when_the_skax_protocol_changes(self) -> None:
        editor = self.open_switching_editor()
        original = self.source.read_bytes()
        write_protocol(self.path, self.samples, (405, 620))
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            editor._save()
        warning.assert_called_once()
        self.assertEqual(self.source.read_bytes(), original)

    def test_rounding_fields_can_be_typed_with_the_theme_applied(self) -> None:
        editor = self.open_switching_editor()
        editor.navigation.setCurrentRow(4)
        editor.show()
        self.application.processEvents()
        spinners = editor.findChildren(QAbstractSpinBox)
        self.assertTrue(spinners)
        for field in editor.report.decimals.values():
            self.assertGreater(field.lineEdit().width(), 20)
            for value in (0, 6, 2):
                field.setFocus()
                field.selectAll()
                QTest.keyClicks(field, str(value))
                QTest.keyClick(field, Qt.Key.Key_Tab)
                self.assertEqual(field.value(), value)
        self.assertEqual(editor.report.to_config()["result_decimals"], 2)
        self.assertFalse(editor.measurement.preferred_wavelength_nm.isEditable())
        self.assertEqual(editor.measurement.preferred_wavelength_nm.count(), 2)
        self.assertEqual(editor.validation.findChildren(QPushButton), [])
        self.assertEqual(editor.calculation.wavelength_tabs.count(), 2)
        self.assertFalse(editor.calculation.wavelength_tabs.tabBar().isMovable())
        self.assertEqual(set(editor.calculation.wavelength_curves), {405, 450})

    def test_reference_choices_use_skax_and_preserve_settings_by_wavelength(self) -> None:
        editor = self.open_switching_editor(reference=620)
        fits = editor.calculation.to_config()["calculation"]["wavelengths"]
        qc = editor.validation.editors[405]
        reference = editor.measurement.reference_wavelength_nm
        self.assertEqual(editor.measurement.wavelengths_nm.text(), "405, 450, 620")
        self.assertEqual(reference.currentData(), 620)
        self.assertEqual(editor.measurement.preferred_wavelength_nm.count(), 2)
        self.assertEqual([editor.validation.tabs.tabText(i) for i in range(editor.validation.tabs.count())],
                         ["Resultaat", "Meting 405 − 620 nm", "Meting 450 − 620 nm"])
        self.assertEqual(editor.calculation.wavelength_tabs.tabText(0), "405 − 620 nm")
        reference.setCurrentIndex(reference.findData(405))
        self.assertEqual(set(editor.calculation.to_config()["calculation"]["wavelengths"]), {"450", "620"})
        reference.setCurrentIndex(reference.findData(620))
        self.assertEqual(editor.calculation.to_config()["calculation"]["wavelengths"], fits)
        self.assertIs(editor.validation.editors[405], qc)
        original = self.source.read_bytes()
        reference.setCurrentIndex(reference.findData(None))
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            editor._save()
        warning.assert_called_once()
        self.assertEqual(self.source.read_bytes(), original)
        reference.setCurrentIndex(reference.findData(620))
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            editor._save()
        warning.assert_not_called()
        saved = self.assays.load_assay("test.skax")
        self.assertEqual(saved.measurement.reference_wavelength_nm, 620)
        self.assertNotIn("wavelengths_nm", saved.to_dict()["measurement"])
        reopened = self.open_editor(saved)
        self.assertEqual(reopened.measurement.reference_wavelength_nm.currentData(), 620)
        self.assertEqual(reopened.calculation.to_config()["calculation"]["wavelengths"], fits)

    def test_wavelength_switching_and_scoped_validation_save_and_reopen(self) -> None:
        write_protocol(self.path, self.samples, (405, 450))
        config = {
            "lot_number": "kit lot",
            "controls": {"Positief": {}},
            "validation": {
                "result": {"controls": {"Positief": {"result": {"min": 1, "max": 2}}}},
                "wavelengths": {
                    "450": {"blank": {"extinction": {"response_scale": "raw", "max": 0.2}}, "standards": {"2.5": {"extinction": {"response_scale": "blank_corrected", "min": 0.1}}}},
                    "405": {"blank": {"extinction": {"response_scale": "raw", "max": 0.1}}},
                },
            },
            "measurement": {"mode": "switch_on_concentration", "cutoff_concentration": 1000.0, "preferred_wavelength_nm": 450},
            "calculation": {
                "method": "calibration",
                "wavelengths": {
                    "405": {"response_scale": "raw", "fit_model": "quadratic", "weighting_exponent": 1.0},
                    "450": {"response_scale": "raw", "fit_model": "4pl", "weighting_exponent": 2.0},
                },
            },
            "result": {"label": "Conc"},
        }
        self.source.chmod(self.source.stat().st_mode | stat.S_IWUSR)
        self.source.write_text(json.dumps({"test.skax": config}), encoding="utf-8")
        self.assays = AssayManagement(self.configs, self.access, protocols_dir=self.directory)
        settings = self.assays.load_assay("test.skax")
        editor = self.open_editor(settings)
        self.assertEqual(editor.measurement.preferred_wavelength_nm.currentData(), 450)
        self.assertTrue(editor.calculation._single_curve_stack.isHidden())
        self.assertFalse(editor.calculation.wavelength_tabs.isHidden())
        self.assertTrue(editor.calculation.wavelength_curves[405].fixed_parameters.isEnabled())
        self.assertEqual(editor.calculation.wavelength_curves[405].fit_model.currentText(), "Kwadratisch")
        editor.calculation.wavelength_curves[405].fit_model.setCurrentText("Kubisch")
        editor.calculation.wavelength_curves[405].weighting_exponent.setText("1.5")
        self.assertEqual(editor.validation.tabs.count(), 3)
        editor.validation.editors[405].fields[("blank", None, "extinction")][1].setText("0.15")
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            editor._save()
        warning.assert_not_called()
        saved = self.assays.load_assay("test.skax")
        self.assertEqual(saved.measurement.preferred_wavelength_nm, 450)
        self.assertEqual(saved.calculation.wavelengths[405].fit_model, "cubic")
        self.assertEqual(saved.calculation.wavelengths[405].weighting_exponent, 1.5)
        self.assertEqual(saved.calculation.wavelengths[450].fit_model, "4pl")
        self.assertEqual(saved.calculation.wavelengths[450].weighting_exponent, 2.0)
        self.assertEqual(saved.validation.wavelengths[405].to_dict()["blank"], {"extinction": {"response_scale": "raw", "max": 0.15}})
        reopened = self.open_editor(saved)
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            reopened._save()
        warning.assert_not_called()
        self.assertEqual(self.assays.load_assay("test.skax"), saved)

    def test_switching_cannot_invent_another_wavelength(self) -> None:
        config = self.settings.to_dict()
        config["calculation"] = {"response_scale": "raw", "method": "calibration", "fit_model": "4pl"}
        self.source.chmod(self.source.stat().st_mode | stat.S_IWUSR)
        self.source.write_text(json.dumps({"test.skax": config}), encoding="utf-8")
        editor = self.open_editor(self.assays.load_assay("test.skax"))
        original = self.source.read_bytes()
        editor.measurement.cutoff_concentration.setText("1000")
        editor.measurement.mode.setCurrentIndex(editor.measurement.mode.findData("switch_on_concentration"))
        self.assertEqual(editor.measurement.wavelengths_nm.text(), "405")
        self.assertEqual(editor.measurement.preferred_wavelength_nm.count(), 1)
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            editor._save()
        warning.assert_called_once()
        self.assertEqual(self.source.read_bytes(), original)

    def test_wavelengths_can_mix_fixed_and_fitted_curves_and_reopen_without_changes(self) -> None:
        editor = self.open_switching_editor()
        editor.navigation.setCurrentRow(2)
        editor.show()
        self.application.processEvents()
        curve = editor.calculation.wavelength_curves[405]
        QTest.mouseClick(curve.fixed_parameters, Qt.MouseButton.LeftButton, pos=QPoint(8, curve.fixed_parameters.height() // 2))
        self.assertTrue(curve.fixed_parameters.isChecked())
        self.assertFalse(curve.response_scale.isHidden())
        self.assertTrue(curve.weighting_exponent.isHidden())
        curve.response_scale.setCurrentIndex(curve.response_scale.findData("raw"))
        parameters = {"A": 0.03, "B": 0.003, "C": 0.000005}
        for name, value in parameters.items():
            curve.parameters[name].setText(str(value))
        editor.calculation.concentration_minimum.setText("1")
        editor.calculation.concentration_maximum.setText("300")
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            editor._save()
        warning.assert_not_called()
        saved = self.assays.load_assay("test.skax")
        fixed = saved.calculation_for_wavelength(405)
        self.assertEqual(fixed.to_dict()["parameters"], parameters)
        self.assertEqual(fixed.response_scale, "raw")
        self.assertEqual(fixed.concentration_range, (1, 300))
        self.assertIsNone(saved.calculation_for_wavelength(450).parameters)
        self.assertEqual(saved.calculation_for_wavelength(450).weighting_exponent, 2)
        reopened = self.open_editor(saved)
        self.assertTrue(reopened.calculation.wavelength_curves[405].fixed_parameters.isChecked())
        self.assertFalse(reopened.calculation.wavelength_curves[450].fixed_parameters.isChecked())
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            reopened._save()
        warning.assert_not_called()
        self.assertEqual(self.assays.load_assay("test.skax"), saved)

        second = reopened.calculation.wavelength_curves[450]
        second.fixed_parameters.setChecked(True)
        for name, value in {"A": 4.5, "B": -1.2, "C": 30, "D": 0.08}.items():
            second.parameters[name].setText(str(value))
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            reopened._save()
        warning.assert_not_called()
        both_fixed = self.assays.load_assay("test.skax")
        self.assertIsNotNone(both_fixed.calculation_for_wavelength(450).parameters)
        self.assertEqual(both_fixed.calculation_for_wavelength(450).response_scale, "raw")
        self.assertEqual(both_fixed.calculation_for_wavelength(405), fixed)

    def test_fixed_curve_drafts_survive_reference_and_measurement_mode_changes(self) -> None:
        editor = self.open_switching_editor(reference=620)
        curve = editor.calculation.wavelength_curves[405]
        curve.fixed_parameters.setChecked(True)
        curve.fit_model.setCurrentText("5pl")
        curve.parameters["E"].setText("0.4")
        curve.parameters["A"].setText("unfinished")
        curve.response_scale.setCurrentIndex(curve.response_scale.findData("raw"))
        reference = editor.measurement.reference_wavelength_nm
        reference.setCurrentIndex(reference.findData(405))
        reference.setCurrentIndex(reference.findData(620))
        editor.measurement.mode.setCurrentIndex(0)
        editor.measurement.mode.setCurrentIndex(1)
        self.assertIs(editor.calculation.wavelength_curves[405], curve)
        self.assertTrue(curve.fixed_parameters.isChecked())
        self.assertEqual(curve.parameters["A"].text(), "unfinished")
        self.assertEqual(curve.parameters["E"].text(), "0.4")
        self.assertEqual(curve.response_scale.currentData(), "raw")
        curve.fixed_parameters.setChecked(False)
        self.assertEqual(editor.calculation.to_config()["calculation"]["wavelengths"]["405"]["weighting_exponent"], 1.5)
        curve.fixed_parameters.setChecked(True)
        original = self.source.read_bytes()
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            editor._save()
        warning.assert_called_once()
        self.assertIn("405 nm", warning.call_args.args[2])
        self.assertEqual(self.source.read_bytes(), original)
        self.assertIs(editor.calculation.wavelength_tabs.currentWidget(), curve)

    def test_read_only_wavelength_choices_can_be_browsed_with_curve_controls_locked(self) -> None:
        editable = self.open_switching_editor()
        readonly = AssayEditor(editable._settings, self.assays, self.access, read_only=True)
        self.addCleanup(readonly.deleteLater)
        readonly.navigation.setCurrentRow(2)
        readonly.show()
        self.application.processEvents()
        selector = readonly.calculation.wavelength_tabs
        for index, wavelength in enumerate((405, 450)):
            QTest.mouseClick(selector.tabBar(), Qt.MouseButton.LeftButton, pos=selector.tabBar().tabRect(index).center())
            self.assertIs(selector.currentWidget(), readonly.calculation.wavelength_curves[wavelength])
            self.assertFalse(selector.currentWidget().fixed_parameters.isEnabled())
            self.assertFalse(selector.currentWidget().fit_model.isEnabled())

    def test_new_wavelengths_preserve_the_known_fixed_draft_and_start_the_other_with_fitting(self) -> None:
        editor = AssayEditor(None, self.assays, self.access, protocol_path=self.incoming_protocol(wavelengths=(405, 450)))
        self.addCleanup(editor.deleteLater)
        reference = editor.measurement.reference_wavelength_nm
        reference.setCurrentIndex(reference.findData(450))
        editor.calculation.curve.fixed_parameters.setChecked(True)
        editor.calculation.curve.parameters["A"].setText("unfinished")
        editor.measurement.cutoff_concentration.setText("100")
        reference.setCurrentIndex(reference.findData(None))
        editor.measurement.mode.setCurrentIndex(1)
        curves = editor.calculation.wavelength_curves
        self.assertTrue(curves[405].fixed_parameters.isChecked())
        self.assertEqual(curves[405].parameters["A"].text(), "unfinished")
        self.assertFalse(curves[450].fixed_parameters.isChecked())
        self.assertEqual(curves[450].parameters["A"].text(), "")
        curves[405].parameters["A"].setText("1.4")
        self.assertEqual(curves[450].parameters["A"].text(), "")
        self.assertIs(editor.calculation.curve, curves[405])
        self.assertEqual(editor.calculation.curve.parameters["A"].text(), "1.4")

    def test_shared_calibration_settings_stay_above_the_wavelengths_when_curves_change(self) -> None:
        editor = self.open_switching_editor()
        editor.navigation.setCurrentRow(2)
        editor.show()
        QTest.qWait(30)
        minimum = editor.calculation.concentration_minimum
        shared = editor.calculation.shared_settings
        tabs = editor.calculation.wavelength_tabs
        self.assertTrue(shared.isAncestorOf(minimum))
        for curve in editor.calculation.wavelength_curves.values():
            self.assertTrue(curve.isAncestorOf(curve.preprocessing))
            self.assertTrue(tabs.isAncestorOf(curve.preprocessing))
            self.assertFalse(shared.isAncestorOf(curve.preprocessing))
        original_y = minimum.mapTo(editor.calculation, minimum.rect().topLeft()).y()
        self.assertLess(original_y, tabs.mapTo(editor.calculation, tabs.rect().topLeft()).y())
        editor.calculation.wavelength_curves[405].fixed_parameters.setChecked(True)
        QTest.qWait(30)
        expanded_y = minimum.mapTo(editor.calculation, minimum.rect().topLeft()).y()
        self.assertEqual(expanded_y, original_y)
        editor.calculation.wavelength_tabs.setCurrentWidget(editor.calculation.wavelength_curves[450])
        QTest.qWait(30)
        fitted_y = minimum.mapTo(editor.calculation, minimum.rect().topLeft()).y()
        self.assertEqual(fitted_y, original_y)

    def test_free_fitting_removes_preprocessing_without_restoring_it(self) -> None:
        editor = self.open_switching_editor(fixed=True)
        for curve in editor.calculation.wavelength_curves.values():
            self.assertFalse(curve.preprocessing.isHidden())
            self.edit_formula(curve.preprocessing, "Monster * 2")
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            editor._save()
        warning.assert_not_called()
        saved = self.assays.load_assay("test.skax")
        editor = self.open_editor(saved)
        curve = editor.calculation.wavelength_curves[405]
        curve.fixed_parameters.setChecked(False)
        self.assertTrue(curve.preprocessing.isHidden())
        self.assertTrue(curve._form.labelForField(curve.preprocessing).isHidden())
        self.assertIsNone(curve.preprocessing.to_config())
        self.assertIsNone(curve.draft()["preprocessing"])
        self.assertNotIn("preprocessing", curve.to_config())
        curve.fixed_parameters.setChecked(True)
        self.assertFalse(curve.preprocessing.isHidden())
        self.assertIsNone(curve.preprocessing.to_config())
        curve.fixed_parameters.setChecked(False)
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            editor._save()
        warning.assert_not_called()
        updated = self.assays.load_assay("test.skax")
        self.assertIsNone(updated.calculation_for_wavelength(405).preprocessing)
        self.assertIsNone(updated.calculation_for_wavelength(405).parameters)
        self.assertEqual(updated.calculation_for_wavelength(450).preprocessing, saved.calculation_for_wavelength(450).preprocessing)
        reopened = self.open_editor(updated)
        curve = reopened.calculation.wavelength_curves[405]
        self.assertTrue(curve.preprocessing.isHidden())
        curve.fixed_parameters.setChecked(True)
        self.assertIsNone(curve.preprocessing.to_config())

    def test_wavelength_preprocessing_round_trips_and_can_be_cleared_independently(self) -> None:
        editor = self.open_switching_editor(fixed=True)
        curves = editor.calculation.wavelength_curves
        self.edit_formula(curves[405].preprocessing, "Monster * 2 / ", "Grens")
        self.edit_formula(curves[450].preprocessing, "Monster / ", "Positief")
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        saved = self.assays.load_assay("test.skax")
        self.assertNotIn("preprocessing", saved.to_dict()["calculation"])
        self.assertEqual(saved.calculation_for_wavelength(405).preprocessing.expression, "sample * 2 / grens")
        self.assertEqual(saved.calculation_for_wavelength(450).preprocessing.references, {"positief": "Positief"})

        reopened = self.open_editor(saved)
        curves = reopened.calculation.wavelength_curves
        self.assertEqual(curves[405].preprocessing.to_config(), saved.calculation_for_wavelength(405).preprocessing.to_dict())
        self.assertEqual(curves[450].preprocessing.to_config(), saved.calculation_for_wavelength(450).preprocessing.to_dict())
        curves[405].preprocessing.clear()
        reopened._save()
        self.assertEqual(reopened.result(), reopened.DialogCode.Accepted)
        updated = self.assays.load_assay("test.skax")
        self.assertIsNone(updated.calculation_for_wavelength(405).preprocessing)
        self.assertEqual(updated.calculation_for_wavelength(450).preprocessing, saved.calculation_for_wavelength(450).preprocessing)

    def test_preprocessing_drafts_follow_their_wavelength_through_mode_and_reference_changes(self) -> None:
        editor = AssayEditor(None, self.assays, self.access, protocol_path=self.incoming_protocol(wavelengths=(405, 450)))
        self.addCleanup(editor.deleteLater)
        reference = editor.measurement.reference_wavelength_nm
        reference.setCurrentIndex(reference.findData(450))
        editor.calculation.curve.fixed_parameters.setChecked(True)
        for name, value in {"A": 2, "B": -1, "C": 1, "D": 0}.items():
            editor.calculation.curve.parameters[name].setText(str(value))
        self.edit_formula(editor.calculation.curve.preprocessing, "Monster / ", "Grens")
        known = editor.calculation.curve.preprocessing.to_config()
        editor.measurement.cutoff_concentration.setText("100")
        reference.setCurrentIndex(reference.findData(None))
        editor.measurement.mode.setCurrentIndex(1)
        curves = editor.calculation.wavelength_curves
        self.assertEqual(curves[405].preprocessing.to_config(), known)
        self.assertIsNone(curves[450].preprocessing.to_config())
        curves[450].fixed_parameters.setChecked(True)
        for name, value in {"A": 2, "B": -1, "C": 1, "D": 0}.items():
            curves[450].parameters[name].setText(str(value))
        self.edit_formula(curves[450].preprocessing, "Monster / ", "Positief")
        other = curves[450].preprocessing.to_config()
        editor.measurement.mode.setCurrentIndex(0)
        reference.setCurrentIndex(reference.findData(450))
        self.assertEqual(editor.calculation.curve.preprocessing.to_config(), known)
        reference.setCurrentIndex(reference.findData(None))
        editor.measurement.mode.setCurrentIndex(1)
        serialized = editor.calculation.to_config()["calculation"]
        self.assertNotIn("preprocessing", serialized)
        self.assertEqual(serialized["wavelengths"]["405"]["preprocessing"], known)
        self.assertEqual(serialized["wavelengths"]["450"]["preprocessing"], other)

    def test_single_wavelength_selects_its_own_curve_and_keeps_edits_when_switching_back(self) -> None:
        editor = self.open_switching_editor(fixed=True)
        curves = editor.calculation.wavelength_curves
        self.edit_formula(curves[405].preprocessing, "Monster / ", "Grens")
        self.edit_formula(curves[450].preprocessing, "Monster / ", "Positief")
        other = curves[450].preprocessing.to_config()
        reference = editor.measurement.reference_wavelength_nm
        editor.measurement.mode.setCurrentIndex(0)
        reference.setCurrentIndex(reference.findData(450))
        self.assertIs(editor.calculation.curve, curves[405])
        single = editor.calculation.to_config()["calculation"]
        self.assertEqual(single["fit_model"], "quadratic")
        self.assertEqual(single["preprocessing"], curves[405].preprocessing.to_config())
        self.edit_formula(editor.calculation.curve.preprocessing, "Monster * 3 / ", "Grens")
        changed = editor.calculation.curve.preprocessing.to_config()
        reference.setCurrentIndex(reference.findData(405))
        self.assertIs(editor.calculation.curve, curves[450])
        self.assertEqual(editor.calculation.to_config()["calculation"]["preprocessing"], other)
        reference.setCurrentIndex(reference.findData(None))
        editor.measurement.mode.setCurrentIndex(1)
        serialized = editor.calculation.to_config()["calculation"]["wavelengths"]
        self.assertEqual(serialized["405"]["preprocessing"], changed)
        self.assertEqual(serialized["450"]["preprocessing"], other)
        self.assertIs(editor.calculation.wavelength_tabs.widget(0), curves[405])
        self.assertIs(editor.calculation.wavelength_tabs.widget(1), curves[450])

    def test_wavelength_fitting_drafts_survive_measurement_mode_changes(self) -> None:
        editor = AssayEditor(None, self.assays, self.access, protocol_path=self.incoming_protocol(wavelengths=(405, 450)))
        self.addCleanup(editor.deleteLater)
        editor.calculation.curve.fit_model.setCurrentText("Kwadratisch")
        editor.calculation.curve.weighting_exponent.setText("0.75")
        mode = editor.measurement.mode
        editor.measurement.cutoff_concentration.setText("1000")
        mode.setCurrentIndex(mode.findData("switch_on_concentration"))
        fits = editor.calculation.to_config()["calculation"]["wavelengths"]
        self.assertEqual(fits, {str(w): {"response_scale": "raw", "fit_model": "quadratic", "weighting_exponent": 0.75} for w in (405, 450)})
        editor.calculation.wavelength_curves[405].fit_model.setCurrentText("Kubisch")
        editor.calculation.wavelength_curves[405].weighting_exponent.setText("1.5")
        mode.setCurrentIndex(mode.findData("single"))
        single = editor.calculation.to_config()["calculation"]
        self.assertNotIn("wavelengths", single)
        self.assertEqual(single["fit_model"], "quadratic")
        self.assertEqual(single["weighting_exponent"], 0.75)
        editor.measurement.cutoff_concentration.setText("1000")
        mode.setCurrentIndex(mode.findData("switch_on_concentration"))
        fits["405"] = {"response_scale": "raw", "fit_model": "cubic", "weighting_exponent": 1.5}
        self.assertEqual(editor.calculation.to_config()["calculation"]["wavelengths"], fits)

    def test_analysis_type_change_refreshes_measurement_and_validation_contexts(self) -> None:
        editor = AssayEditor(None, self.assays, self.access, protocol_path=self.incoming_protocol())
        self.addCleanup(editor.deleteLater)
        editor.measurement.cutoff_concentration.setText("1000")
        editor.measurement.mode.setCurrentIndex(editor.measurement.mode.findData("switch_on_concentration"))
        blank_limit = editor.validation.editors[None].fields[("blank", None, "extinction")][1]
        self.assertFalse(blank_limit.isEnabled())
        editor.calculation_method.setCurrentText("Berekening")
        self.assertEqual(editor.measurement.to_config()["mode"], "single")
        self.assertTrue(editor.validation.editors[None].fields[("blank", None, "extinction")][1].isEnabled())
        editor.calculation_method.setCurrentText("Standaardcurve")
        self.assertEqual(editor.measurement.selection().mode, "switch_on_concentration")
        self.assertFalse(editor.validation.editors[None].fields[("blank", None, "extinction")][1].isEnabled())

    def test_ui_shows_every_protocol_sample_and_locks_identities(self) -> None:
        editor = self.open_editor()
        self.assertEqual(tuple(editor.controls.tabs.tabText(i) for i in range(editor.controls.tabs.count())), ("Grens", "Positief"))
        self.assertEqual(len(editor.standards.editors), 2)
        standard = editor.standards.editors[2.5]
        self.assertEqual(editor.standards.tabs.tabText(0), "Standaard1")
        self.assertEqual(float(standard.concentration.text()), 2.5)
        for field in (standard.concentration,):
            before = field.text()
            field.selectAll()
            QTest.keyClicks(field, "Invented")
            self.assertEqual(field.text(), before)
        buttons = {button.text() for button in editor.findChildren(QPushButton)}
        self.assertTrue(buttons.isdisjoint({"Add control", "Remove control", "Add standard", "Remove standard"}))
        self.assertFalse(editor.controls.editors["Positief"].lot_number.isReadOnly())

    def test_validation_shows_duplo_samples_once_with_independent_metric_limits(self) -> None:
        editor = self.open_editor()
        selected = editor.validation.editors[None]
        self.assertEqual(
            [selected.item(row, 0).text() for row in range(selected.rowCount())],
            ["Blanco", "Grens", "Positief", "2.5 (Standaard1)", "0 (Standaard2)"],
        )
        selected.fields[("controls", "Positief", "extinction")][1].setText("2.0")
        selected.fields[("controls", "Positief", "result")][0].setText("1.3")
        for fields in (selected.fields[("controls", "Positief", metric)] for metric in ("extinction", "result")):
            for field in fields:
                self.assertIn(field, [selected.cellWidget(2, column) for column in range(selected.columnCount())])
        self.assertEqual(editor.validation.tabs.count(), 1)
        self.assertTrue(editor.validation.tabs.tabBar().isHidden())
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            editor._save()
        warning.assert_not_called()
        reopened = self.open_editor(self.assays.load_assay("test.skax"))
        self.assertFalse(self.assays.load_assay("test.skax").validation.wavelengths)
        self.assertEqual(reopened.validation.editors[None].to_config(), selected.to_config())

    def test_control_pages_remain_readable_with_a_dark_system_palette(self) -> None:
        original = self.application.palette()
        try:
            for editable in (True, False):
                with self.subTest(editable=editable):
                    dark = QPalette(original)
                    for role in (
                        QPalette.ColorRole.Window,
                        QPalette.ColorRole.Base,
                        QPalette.ColorRole.Button,
                    ):
                        dark.setColor(role, QColor("#20242A"))
                    dark.setColor(QPalette.ColorRole.WindowText, QColor("#F5F5F5"))
                    dark.setColor(QPalette.ColorRole.Text, QColor("#F5F5F5"))
                    self.application.setPalette(dark)
                    self.access.editable = editable
                    editor = self.open_editor()
                    editor.show()
                    self.application.processEvents()
                    # Check the content surface below the attached tab rows;
                    # the top-left pixel now belongs to the selected subtab.
                    for surface in (editor.controls, editor.controls.editors["Grens"]):
                        background = surface.grab().toImage().pixelColor(12, surface.height() - 12)
                        self.assertGreaterEqual(background.lightness(), 240)
                    field = editor.standards.editors[2.5].concentration
                    for group in (QPalette.ColorGroup.Active, QPalette.ColorGroup.Inactive):
                        foreground = field.palette().color(group, QPalette.ColorRole.Text)
                        self.assertLess(foreground.lightness(), 140)
                    editor.close()
        finally:
            self.access.editable = True
            self.application.setPalette(original)

    def test_dual_result_limits_are_exclusive_to_the_result_tab(self) -> None:
        editor = self.open_switching_editor()
        editor.calculation.result_label.setText("Conc")
        editor.calculation.result_unit.setText("U/ml")
        validation = editor.validation
        self.assertEqual([validation.tabs.tabText(i) for i in range(validation.tabs.count())],
                         ["Resultaat", "Meting 405 nm", "Meting 450 nm"])
        result = validation.editors[None]
        self.assertTrue(result.isColumnHidden(1))
        self.assertTrue(result.isColumnHidden(2))
        self.assertTrue(result.isRowHidden(0))
        self.assertTrue(result.isColumnHidden(3))
        self.assertEqual(result.horizontalHeaderItem(4).text(), "Resultaat\nMinimum")
        self.assertEqual(result.horizontalHeaderItem(5).text(), "Resultaat\nMaximum")
        result.fields[("controls", "Positief", "result")][0].setText("1.25")
        for wavelength in (405, 450):
            measurement = validation.editors[wavelength]
            self.assertEqual(measurement.columnCount(), 4)
            self.assertTrue(all(metric == "extinction" for _, _, metric in measurement.fields))
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            editor._save()
        warning.assert_not_called()
        saved = self.assays.load_assay("test.skax").validation.to_dict()
        self.assertEqual(saved["result"]["controls"]["Positief"], {"result": {"min": 1.25}})
        self.assertNotIn("selected", saved)
        self.assertTrue(all("result" not in scope.get("controls", {}).get("Positief", {}) for scope in saved["wavelengths"].values()))

    def test_single_reference_measurement_has_one_criteria_table(self) -> None:
        write_protocol(self.path, self.samples, (450, 620))
        config = self.settings.to_dict()
        config["measurement"] = {"mode": "single", "reference_wavelength_nm": 620}
        settings = AssaySettings.from_dict("test.skax", config)
        editor = self.open_editor(settings)
        self.assertEqual(editor.validation.tabs.count(), 1)
        self.assertTrue(editor.validation.tabs.tabBar().isHidden())
        table = editor.validation.editors[None]
        self.assertIn("450 − 620 nm", table.horizontalHeaderItem(1).text())
        self.assertFalse(table.isColumnHidden(1))
        self.assertFalse(table.isColumnHidden(3))
        self.assertNotIn("wavelengths", editor.validation.to_config())

    def test_qc_drafts_survive_switching_between_single_and_dual_modes(self) -> None:
        editor = self.open_switching_editor()
        mode = editor.measurement.mode
        mode.setCurrentIndex(mode.findData("single"))
        single = editor.validation.editors[None]
        blank = ("blank", None, "extinction")
        self.assertEqual(single.fields[blank][1].text(), "0.2")
        single.fields[blank][1].setText("0.25")
        editor.measurement.cutoff_concentration.setText("1000")
        mode.setCurrentIndex(mode.findData("switch_on_concentration"))
        self.assertEqual(editor.validation.editors[450].fields[blank][1].text(), "0.25")
        self.assertEqual(editor.validation.editors[405].fields[blank][1].text(), "0.15")
        self.assertNotIn("blank", editor.validation.to_config().get("result", {}))
        mode.setCurrentIndex(mode.findData("single"))
        self.assertEqual(single.fields[blank][1].text(), "0.25")
        self.assertEqual(editor.validation.tabs.count(), 1)
        single.fields[blank][1].clear()
        editor.measurement.cutoff_concentration.setText("1000")
        mode.setCurrentIndex(mode.findData("switch_on_concentration"))
        self.assertNotIn("450", editor.validation.to_config()["wavelengths"])
        self.assertEqual(editor.validation.editors[405].fields[blank][1].text(), "0.15")

    def test_invalid_criteria_are_highlighted_and_block_save_until_corrected(self) -> None:
        editor = self.open_editor()
        table = editor.validation.editors[None]
        minimum, maximum = table.fields[("controls", "Positief", "result")]
        original = self.source.read_bytes()
        minimum.setText("2,5")
        maximum.setText("2")
        self.assertTrue(minimum.property("invalid"))
        self.assertIn("minimum", editor.validation.error.text())
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            editor._save()
        warning.assert_called_once()
        self.assertEqual(self.source.read_bytes(), original)
        maximum.setText("3")
        self.assertFalse(minimum.property("invalid"))
        self.assertTrue(editor.validation.error.isHidden())
        self.assertEqual(table.to_config()["controls"]["Positief"]["result"], {"min": 2.5, "max": 3})
        for invalid in ("text", "nan", "inf"):
            with self.subTest(value=invalid):
                maximum.setText(invalid)
                self.assertTrue(maximum.property("invalid"))
                with self.assertRaises(ValueError):
                    table.to_config()

    def test_formula_report_hides_calibration_options_and_restores_their_drafts(self) -> None:
        editor = AssayEditor(None, self.assays, self.access, protocol_path=self.incoming_protocol())
        self.addCleanup(editor.deleteLater)
        row = editor.report.decimals["calibration_decimals"]
        self.assertFalse(row.isHidden())
        editor.report.decimals["calibration_decimals"].setValue(4)
        editor.calculation_method.setCurrentText("Berekening")
        self.assertTrue(row.isHidden())
        self.assertTrue(editor.report.recovery_warning_percent.isHidden())
        editor.calculation_method.setCurrentText("Standaardcurve")
        self.assertFalse(row.isHidden())
        self.assertEqual(editor.report.decimals["calibration_decimals"].value(), 4)

    def test_unchanged_editor_preserves_exact_source_and_sparse_settings(self) -> None:
        original = self.source.read_bytes()
        editor = self.open_editor()
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        self.assertEqual(self.source.read_bytes(), original)
        self.assertFalse((self.configs / ".history").exists())

    def test_editing_qc_saves_under_skax_names_and_nominal_concentrations(self) -> None:
        editor = self.open_editor()
        editor.controls.editors["Grens"].lot_number.setText("new control lot")
        standard = editor.standards.editors[0.0]
        standard.lot_number.setText("standard lot")
        editor.validation.editors[None].fields[("standards", "0.0", "extinction")][1].setText("0.5")
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        saved = self.assays.load_assay("test.skax")
        self.assertEqual(saved.controls["Grens"].lot_number, "new control lot")
        self.assertEqual(saved.controls["Positief"], self.settings.controls["Positief"])
        self.assertEqual(saved.standards[0].lot_number, "standard lot")
        self.assertEqual(saved.standards[2.5], self.settings.standards[2.5])
        self.assertEqual(set(saved.controls), {"Grens", "Positief"})
        self.assertEqual(set(saved.standards), {0, 2.5})

    def test_replicate_filter_editor_saves_reopens_and_disables_rules(self) -> None:
        editor = self.open_editor()
        patient = editor.patients.replicate_filter
        self.assertFalse(patient.enabled.isChecked())
        patient.enabled.setChecked(True)
        patient.fields["max_deviation_percent"].setText("10")
        patient.response_scale.setCurrentIndex(patient.response_scale.findData("raw"))
        control = editor.controls.editors["Grens"].replicate_filter
        self.assertFalse(control.enabled.isChecked())
        control.enabled.setChecked(True)
        for name, value in {"max_deviation_percent": 15, "max_exclusions": 2}.items():
            control.fields[name].setText(str(value))
        standard = editor.standards.editors[2.5].replicate_filter
        standard.enabled.setChecked(True)
        standard.fields["max_deviation_percent"].setText("20")
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        saved = self.assays.load_assay("test.skax")
        self.assertEqual(saved.patients.replicate_filter, ReplicateFilterSettings(10, response_scale="raw"))
        self.assertEqual(saved.controls["Grens"].replicate_filter, ReplicateFilterSettings(15, 2, response_scale="blank_corrected"))
        self.assertEqual(saved.standards[2.5].replicate_filter, ReplicateFilterSettings(20, response_scale="blank_corrected"))
        reopened = self.open_editor(saved)
        self.assertEqual(reopened.patients.replicate_filter.to_config(), saved.patients.replicate_filter.to_dict())
        self.assertEqual(reopened.controls.editors["Grens"].replicate_filter.to_config(), saved.controls["Grens"].replicate_filter.to_dict())
        # Saving unrelated settings must preserve all filtering rules.
        reopened.lot_number.setText("updated lot")
        reopened._save()
        updated = self.assays.load_assay("test.skax")
        self.assertEqual(updated.patients, saved.patients)
        self.assertEqual(updated.controls, saved.controls)
        self.assertEqual(updated.standards, saved.standards)
        disabled = self.open_editor(updated)
        disabled.patients.replicate_filter.enabled.setChecked(False)
        disabled.controls.editors["Grens"].replicate_filter.enabled.setChecked(False)
        disabled._save()
        self.assertIsNone(self.assays.load_assay("test.skax").patients.replicate_filter)
        self.assertNotIn("patients", self.assays.load_assay("test.skax").to_dict())
        self.assertIsNone(self.assays.load_assay("test.skax").controls["Grens"].replicate_filter)

    def test_report_tab_loads_defaults_and_saves_all_options(self) -> None:
        editor = self.open_editor()
        self.assertEqual(editor.report.to_config(), ReportSettings().to_dict())
        options = ReportSettings(
            result_decimals=2,
            extinction_decimals=4,
            calibration_decimals=3,
            percentage_decimals=0,
            cv_warning_percent=12.5,
            recovery_warning_percent=15,
        )
        for name, field in editor.report.decimals.items():
            field.setValue(getattr(options, name))
        editor.report.cv_warning_percent.setText("12.5")
        editor.report.recovery_warning_percent.setText("15")
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        saved = self.assays.load_assay("test.skax")
        self.assertEqual(saved.report, options)
        self.assertEqual(saved.controls["Positief"], self.settings.controls["Positief"])
        self.assertEqual(json.loads(self.source.read_text())["test.skax"]["report"], options.to_dict())
        reopened = self.open_editor(saved)
        self.assertEqual(reopened.report.to_config(), options.to_dict())

    def test_invalid_report_thresholds_do_not_save(self) -> None:
        editor = self.open_editor()
        original = self.source.read_bytes()
        for field in (editor.report.cv_warning_percent, editor.report.recovery_warning_percent):
            for value in ("", "nan", "inf", "-1", "invalid"):
                with self.subTest(value=value), patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
                    field.setText(value)
                    editor._save()
                    warning.assert_called_once()
                    self.assertEqual(self.source.read_bytes(), original)
                    self.assertNotEqual(editor.result(), editor.DialogCode.Accepted)
            field.setText("20")

    def test_missing_reference_requires_explicit_correction(self):
        self.config["calculation"]["references"]["grens"] = "Missing"
        self.source.chmod(self.source.stat().st_mode | stat.S_IWUSR)
        self.source.write_text(json.dumps({"test.skax": self.config}), encoding="utf-8")
        self.settings = self.assays.load_assay("test.skax")
        editor = self.open_editor()
        original = self.source.read_bytes()
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            editor._save()
        warning.assert_called_once()
        self.assertEqual(self.source.read_bytes(), original)
        self.edit_formula(editor.calculation.formula, "Patient / ", "Grens")
        self.edit_formula(editor.calculation.overrides["Positief"].formula, "Positief / ", "Grens")
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        self.assertEqual(set(self.assays.load_assay("test.skax").calculation.formula.references.values()), {"Grens"})

    def test_formula_overrides_use_protocol_controls(self):
        editor = self.open_editor()
        self.assertEqual(editor.calculation.formula.summary.text(), "Patient / Grens")
        self.assertEqual(editor.calculation.overrides["Positief"].formula.summary.text(), "Positief / Grens")
        self.assertEqual(editor.calculation.formula.to_config()["expression"], "sample / grens")
        self.assertEqual(editor.calculation.overrides["Positief"].formula.to_config()["expression"], "sample / grens")
        self.assertEqual(editor.calculation.to_config()["calculation"]["control_overrides"]["Positief"]["expression"], "sample / grens")
        control = editor.calculation.overrides["Positief"]
        self.assertFalse(control.use_sample_result.isChecked())
        self.assertTrue(control.formula.edit_button.isEnabled())
        self.assertEqual(control.result_label.text(), "Ratio")

    def test_unused_formula_references_are_removed_when_saving(self):
        preprocessing = FormulaEditor(
            {"expression": "sample / ref", "references": {"unused": "Positief", "ref": "Grens"}}, ("Grens", "Positief"), sample_name="Monster",
        )
        self.addCleanup(preprocessing.deleteLater)
        self.assertEqual(preprocessing.to_config(), {"expression": "sample / ref", "references": {"ref": "Grens"}})
        self.edit_formula(preprocessing, "Monster")
        self.assertEqual(preprocessing.to_config(), {"expression": "sample"})

    def test_control_overrides_insert_controls_and_keep_only_used_references(self):
        editor = self.open_editor()
        self.edit_formula(editor.calculation.formula, "Patient")
        self.edit_formula(editor.calculation.overrides["Positief"].formula, "Positief / ", "Grens")
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        saved = self.assays.load_assay("test.skax")
        self.assertEqual(saved.calculation.formula.expression, "sample")
        self.assertEqual(saved.calculation.formula.references, {"grens": "Grens"})
        self.assertEqual(saved.calculation.control_overrides["Positief"].expression, "sample / grens")
        self.assertEqual(saved.calculation.overrides["Positief"].evaluate(4, {"Grens": 2}), 2)
        reopened = self.open_editor(saved)
        reopened.calculation.overrides["Positief"].formula.clear()
        config = reopened.calculation.to_config()["calculation"]
        self.assertNotIn("references", config)
        self.assertEqual(config["control_overrides"]["Positief"], {"response_scale": "blank_corrected"})

    def test_control_calculation_checkbox_applies_formula_label_and_unit_together(self):
        editor = self.open_editor()
        control = editor.calculation.overrides["Grens"]
        self.assertTrue(control.use_sample_result.isChecked())
        self.assertEqual(control.formula.summary.text(), "Patientformule")
        self.assertFalse(control.formula.edit_button.isEnabled())
        self.assertFalse(control.formula.clear_button.isEnabled())
        control.use_sample_result.setChecked(False)
        self.assertTrue(control.result_label.isEnabled())
        self.assertTrue(control.result_unit.isEnabled())
        self.assertTrue(control.formula.edit_button.isEnabled())
        self.edit_formula(control.formula, "Grens * 2")
        control.result_label.setText("Controle-index")
        control.result_unit.setText("index")
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        saved = self.assays.load_assay("test.skax")
        self.assertEqual(saved.calculation.control_overrides["Grens"].expression, "sample * 2")
        self.assertEqual(saved.calculation.control_results["Grens"], ResultSettings("Controle-index", "index"))

        reopened = self.open_editor(saved)
        control = reopened.calculation.overrides["Grens"]
        self.assertFalse(control.use_sample_result.isChecked())
        self.assertEqual(control.formula.summary.text(), "Grens * 2")
        control.use_sample_result.setChecked(True)
        self.assertEqual(control.formula.summary.text(), "Patientformule")
        self.assertFalse(control.formula.edit_button.isEnabled())
        self.assertFalse(control.result_label.isEnabled())
        self.assertFalse(control.result_unit.isEnabled())
        reopened.calculation.result_label.setText("Patiëntscore")
        reopened.calculation.result_unit.setText("%")
        self.assertEqual(control.result_label.text(), "Patiëntscore")
        self.assertEqual(control.result_unit.text(), "%")
        reopened._save()
        self.assertEqual(reopened.result(), reopened.DialogCode.Accepted)
        saved = self.assays.load_assay("test.skax")
        self.assertNotIn("Grens", saved.calculation.control_overrides)
        self.assertNotIn("Grens", saved.calculation.control_results)
        self.assertTrue(self.open_editor(saved).calculation.overrides["Grens"].use_sample_result.isChecked())

    def test_unchecking_patient_calculation_keeps_identical_control_metadata_independent(self):
        editor = self.open_editor()
        editor.calculation.overrides["Grens"].use_sample_result.setChecked(False)
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        reopened = self.open_editor(self.assays.load_assay("test.skax"))
        control = reopened.calculation.overrides["Grens"]
        self.assertFalse(control.use_sample_result.isChecked())
        self.assertTrue(control.formula.edit_button.isEnabled())
        self.assertEqual(control.formula.summary.text(), "Patientformule")
        reopened.calculation.result_label.setText("Nieuwe score")
        reopened.calculation.result_unit.setText("%")
        self.assertEqual(control.result_label.text(), "Ratio")
        self.assertEqual(control.result_unit.text(), "")

    def test_control_results_save_reopen_and_stay_independent_of_sample_results(self):
        editor = self.open_editor()
        editor.calculation.result_label.setText("Score")
        editor.calculation.result_unit.setText("%")
        card = editor.calculation.overrides["Positief"]
        card.use_sample_result.setChecked(False)
        card.result_label.setText("Controle-index")
        card.result_unit.setText("index")
        linked = editor.calculation.overrides["Grens"]
        self.assertTrue(linked.use_sample_result.isChecked())
        self.assertEqual(linked.result_label.text(), "Score")
        self.assertEqual(linked.result_unit.text(), "%")
        self.assertFalse(linked.result_label.isEnabled())
        self.assertFalse(linked.result_unit.isEnabled())
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        saved = self.assays.load_assay("test.skax")
        self.assertEqual(saved.calculation.control_results, {"Positief": ResultSettings("Controle-index", "index")})
        self.assertEqual(saved.result.label, "Score")
        self.assertEqual(saved.result.unit, "%")
        reopened = self.open_editor(saved)
        custom = reopened.calculation.overrides["Positief"]
        self.assertFalse(custom.use_sample_result.isChecked())
        self.assertEqual(custom.result_label.text(), "Controle-index")
        self.assertEqual(custom.result_unit.text(), "index")
        self.assertTrue(custom.result_label.isEnabled())
        self.assertTrue(custom.result_unit.isEnabled())
        reopened.calculation.result_label.setText("Monsterscore")
        reopened.calculation.result_unit.setText("U/ml")
        self.assertEqual(custom.result_label.text(), "Controle-index")
        self.assertEqual(custom.result_unit.text(), "index")
        self.assertEqual(reopened.calculation.overrides["Grens"].result_label.text(), "Monsterscore")
        self.assertEqual(reopened.calculation.overrides["Grens"].result_unit.text(), "U/ml")
        self.assertEqual(reopened.validation.editors[None].horizontalHeaderItem(4).text(), "Resultaat\nMinimum")
        custom.use_sample_result.setChecked(True)
        self.assertEqual(custom.result_label.text(), "Monsterscore")
        self.assertEqual(custom.result_unit.text(), "U/ml")
        self.assertNotIn("control_results", reopened.calculation.to_config()["calculation"])

    def test_dimensionless_control_result_is_preserved_and_empty_label_blocks_save(self):
        editor = self.open_editor()
        editor.calculation.result_label.setText("Score")
        editor.calculation.result_unit.setText("%")
        card = editor.calculation.overrides["Positief"]
        card.use_sample_result.setChecked(False)
        card.result_label.setText("Controle-index")
        card.result_unit.clear()
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        saved = self.assays.load_assay("test.skax")
        self.assertEqual(saved.calculation.control_results["Positief"], ResultSettings("Controle-index"))
        self.assertEqual(saved.to_dict()["calculation"]["control_results"]["Positief"], {"label": "Controle-index"})
        reopened = self.open_editor(saved)
        custom = reopened.calculation.overrides["Positief"]
        self.assertEqual(custom.result_unit.text(), "")
        self.assertFalse(custom.use_sample_result.isChecked())
        custom.result_label.clear()
        original = self.source.read_bytes()
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            reopened._save()
        warning.assert_called_once()
        self.assertEqual(self.source.read_bytes(), original)

    def test_concentration_editor_round_trips_and_control_free_preprocessing_uses_numbers(
        self,
    ) -> None:
        write_protocol(self.path, [sample for sample in self.samples if sample[0] != "Control"])
        self.config = {
            "lot_number": "kit lot",
            "calculation": {"response_scale": "raw", "method": "calibration", "fit_model": "4pl",
                            "parameters": {"A": 2.5, "B": -1.3, "C": 2, "D": 0.1}, "preprocessing": {"expression": "(sample * 2)"}},
            "result": {"label": "Conc"},
        }
        self.source.chmod(self.source.stat().st_mode | stat.S_IWUSR)
        self.source.write_text(json.dumps({"test.skax": self.config}), encoding="utf-8")
        self.settings = self.assays.load_assay("test.skax")
        editor = self.open_editor()
        self.assertEqual(len(editor.controls.editors), 0)
        self.assertEqual(len(editor.standards.editors), 2)
        def inspect(dialog):
            self.assertEqual(dialog.input.operand.count(), 1)
            self.assertEqual(dialog.input.operand.itemText(0), "Monster (meetwaarde)")
            self.assertFalse(dialog.input.insert_button.isEnabled())
            return dialog.DialogCode.Rejected

        with patch.object(FormulaDialog, "exec", inspect):
            editor.calculation.curve.preprocessing._edit()
        original = self.source.read_bytes()
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        self.assertEqual(self.source.read_bytes(), original)

        self.assertFalse(editor.calculation.curve.response_scale.isHidden())
        self.assertEqual(self.assays.load_assay("test.skax").calculation.response_scale, "raw")
        editor.calculation.curve.fixed_parameters.setChecked(True)
        self.assertFalse(editor.calculation.curve.response_scale.isHidden())
        self.assertEqual(editor.calculation.curve.response_scale.currentData(), "raw")
        editor.calculation.curve.response_scale.setCurrentIndex(editor.calculation.curve.response_scale.findData("blank_corrected"))
        self.assertEqual(editor.calculation.curve.response_scale.currentData(), "blank_corrected")
        for name, value in {"A": 2.5, "B": -1.3, "C": 2, "D": 0.1}.items():
            editor.calculation.curve.parameters[name].setText(str(value))
        editor._save()
        saved = self.assays.load_assay("test.skax")
        self.assertEqual(saved.calculation.response_scale, "blank_corrected")
        reopened = self.open_editor(saved)
        self.assertFalse(reopened.calculation.curve.response_scale.isHidden())
        reopened.calculation.curve.response_scale.setCurrentIndex(reopened.calculation.curve.response_scale.findData("raw"))
        reopened._save()
        self.assertEqual(self.assays.load_assay("test.skax").calculation.response_scale, "raw")
        fitting = self.open_editor(self.assays.load_assay("test.skax"))
        fitting.calculation.curve.fixed_parameters.setChecked(False)
        self.assertFalse(fitting.calculation.curve.response_scale.isHidden())
        fitting._save()
        fitted = self.assays.load_assay("test.skax").calculation
        self.assertEqual(fitted.response_scale, "raw")
        self.assertEqual(fitted.to_dict()["response_scale"], "raw")

    def test_polynomial_editor_saves_reopens_and_switches_coefficient_fields(self) -> None:
        path = self.incoming_protocol()
        editor = AssayEditor(None, self.assays, self.access, protocol_path=path)
        self.addCleanup(editor.deleteLater)
        editor.lot_number.setText("kit lot")
        editor.calculation.curve.fit_model.setCurrentText("Kwadratisch")
        self.assertTrue(editor.calculation.curve.parameters["D"].isHidden())
        self.assertTrue(editor.calculation.curve.parameters["E"].isHidden())
        editor.calculation.curve.fixed_parameters.setChecked(True)
        for name, value in {"A": 0.1, "B": 0.3, "C": 0.02}.items():
            editor.calculation.curve.parameters[name].setText(str(value))
        editor.calculation.concentration_minimum.setText("0")
        editor.calculation.concentration_maximum.setText("10")
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        saved = self.assays.load_assay(path.name)
        self.assertEqual(saved.calculation.fit_model, "quadratic")
        self.assertEqual(saved.calculation.to_dict()["parameters"], {"A": 0.1, "B": 0.3, "C": 0.02})
        reopened = self.open_editor(saved)
        reopened.calculation.curve.fit_model.setCurrentText("Kubisch")
        self.assertFalse(reopened.calculation.curve.parameters["D"].isHidden())
        self.assertTrue(reopened.calculation.curve.parameters["E"].isHidden())
        reopened.calculation.curve.parameters["D"].setText("0.003")
        reopened._save()
        saved = self.assays.load_assay(path.name)
        self.assertEqual(saved.calculation.parameters.D, 0.003)
        fitted = self.open_editor(saved)
        fitted.calculation.curve.fixed_parameters.setChecked(False)
        fitted.calculation.curve.weighting_exponent.setText("0.75")
        fitted._save()
        saved = self.assays.load_assay(path.name)
        self.assertEqual(saved.calculation.fit_model, "cubic")
        self.assertEqual(saved.calculation.response_scale, "raw")
        self.assertIsNone(saved.calculation.parameters)
        self.assertEqual(saved.calculation.weighting_exponent, 0.75)

    def test_power_law_editor_saves_reopens_and_switches_between_fixed_and_fitted(self) -> None:
        path = self.incoming_protocol()
        editor = AssayEditor(None, self.assays, self.access, protocol_path=path)
        self.addCleanup(editor.deleteLater)
        editor.lot_number.setText("kit lot")
        curve = editor.calculation.curve
        curve.fit_model.setCurrentText("Machtsfunctie met offset")
        self.assertEqual(curve.fit_model.currentData(), "power_law")
        curve.fixed_parameters.setChecked(True)
        self.assertEqual(curve.equation.text(), "y = C + A·x^B")
        self.assertFalse(curve.equation.isHidden())
        for name, value in {"A": 0.3, "B": 1.7, "C": 0.12}.items():
            self.assertFalse(curve.parameters[name].isHidden())
            curve.parameters[name].setText(str(value))
        for name in ("D", "E"):
            self.assertTrue(curve.parameters[name].isHidden())
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        saved = self.assays.load_assay(path.name)
        self.assertEqual(saved.calculation.fit_model, "power_law")
        self.assertEqual(saved.calculation.to_dict()["parameters"], {"A": 0.3, "B": 1.7, "C": 0.12})
        reopened = self.open_editor(saved)
        self.assertEqual(reopened.calculation.curve.fit_model.currentData(), "power_law")
        reopened.calculation.curve.fixed_parameters.setChecked(False)
        reopened.calculation.curve.weighting_exponent.setText("0.75")
        reopened._save()
        saved = self.assays.load_assay(path.name)
        self.assertEqual(saved.calculation.fit_model, "power_law")
        self.assertIsNone(saved.calculation.parameters)
        self.assertEqual(saved.calculation.response_scale, "raw")
        self.assertEqual(saved.calculation.weighting_exponent, 0.75)

    def test_viewer_permissions_still_disable_editing_and_saving(self) -> None:
        self.access.editable = False
        editor = self.open_editor()
        self.assertTrue(editor.lot_number.isReadOnly())
        self.assertFalse(editor.controls.editors["Positief"].lot_number.isEnabled())
        self.assertFalse(editor.report.isEnabled())
        self.assertTrue(all(not field.isEnabled() for field in editor.report.decimals.values()))
        self.assertFalse(editor.report.cv_warning_percent.isEnabled())
        self.assertFalse(editor.report.recovery_warning_percent.isEnabled())
        buttons = editor.findChild(QDialogButtonBox)
        self.assertFalse(buttons.button(QDialogButtonBox.StandardButton.Save).isEnabled())

    def test_model_switching_preserves_five_pl_draft_and_saves_matching_parameters(self) -> None:
        path = self.incoming_protocol()
        editor = AssayEditor(None, self.assays, self.access, protocol_path=path)
        self.addCleanup(editor.deleteLater)
        editor.lot_number.setText("kit lot")
        self.assertFalse(editor.calculation.curve.response_scale.isHidden())
        self.assertTrue(editor.calculation.curve.parameters["E"].isHidden())
        editor.calculation.curve.fit_model.setCurrentText("5pl")
        self.assertTrue(editor.calculation.curve.parameters["E"].isHidden())
        self.assertEqual(float(editor.calculation.curve.parameters["E"].text()), 1)
        editor.calculation.curve.fixed_parameters.setChecked(True)
        self.assertFalse(editor.calculation.curve.parameters["E"].isHidden())
        self.assertFalse(editor.calculation.curve.response_scale.isHidden())
        editor.calculation.curve.response_scale.setCurrentIndex(editor.calculation.curve.response_scale.findData("blank_corrected"))
        self.assertEqual(editor.calculation.curve.response_scale.currentData(), "blank_corrected")
        parameters = {"A": 2.5, "B": -1.3, "C": 2, "D": 0.1, "E": 0.4}
        for name, value in parameters.items():
            editor.calculation.curve.parameters[name].setText(str(value))
        editor.calculation.curve.fit_model.setCurrentText("4pl")
        self.assertTrue(editor.calculation.curve.parameters["E"].isHidden())
        editor.calculation.curve.fit_model.setCurrentText("5pl")
        self.assertEqual(editor.calculation.curve.parameters["E"].text(), "0.4")
        editor.calculation_method.setCurrentText("Berekening")
        editor.calculation_method.setCurrentText("Standaardcurve")
        self.assertEqual(editor.calculation.curve.response_scale.currentData(), "blank_corrected")
        self.assertEqual(editor.calculation.curve.fit_model.currentText(), "5pl")
        self.assertEqual(editor.calculation.curve.parameters["E"].text(), "0.4")
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        saved = self.assays.load_assay(path.name)
        self.assertEqual(saved.calculation.to_dict()["parameters"], parameters)
        self.assertEqual(saved.calculation.response_scale, "blank_corrected")
        reopened = self.open_editor(saved)
        self.assertFalse(reopened.calculation.curve.parameters["E"].isHidden())
        self.assertEqual(reopened.calculation.curve.parameters["E"].text(), "0.4")
        reopened.calculation.curve.fit_model.setCurrentText("4pl")
        reopened._save()
        saved = self.assays.load_assay(path.name)
        self.assertEqual(saved.calculation.fit_model, "4pl")
        self.assertEqual(set(saved.calculation.to_dict()["parameters"]), {"A", "B", "C", "D"})
        fitted = self.open_editor(saved)
        fitted.calculation.curve.fit_model.setCurrentText("5pl")
        fitted.calculation.curve.fixed_parameters.setChecked(False)
        fitted.calculation.curve.weighting_exponent.setText("1.25")
        fitted._save()
        saved = self.assays.load_assay(path.name)
        self.assertEqual(saved.calculation.fit_model, "5pl")
        self.assertIsNone(saved.calculation.parameters)
        self.assertEqual(saved.calculation.weighting_exponent, 1.25)

    def test_new_concentration_assay_requires_a_lot_and_saves_protocol_samples(self) -> None:
        path = self.incoming_protocol()
        editor = AssayEditor(None, self.assays, self.access, protocol_path=path)
        self.addCleanup(editor.deleteLater)
        self.assertEqual(editor.lot_number.text(), "")
        self.assertEqual(editor.report.to_config(), ReportSettings().to_dict())
        self.assertEqual(len(editor.standards.editors), 2)
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            editor._save()
        warning.assert_called_once()
        self.assertNotIn(path.name, self.assays.list_assays())
        self.assertFalse((self.directory / path.name).exists())
        editor.lot_number.setText("new kit lot")
        editor.calculation.result_unit.setText("U/mL")
        editor.report.decimals["calibration_decimals"].setValue(2)
        editor.standards.editors[2.5].lot_number.setText("standard lot")
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        saved = self.assays.load_assay(path.name)
        self.assertEqual(saved.calculation.method.value, "calibration")
        self.assertEqual(saved.lot_number, "new kit lot")
        self.assertEqual(saved.result.unit, "U/mL")
        self.assertEqual(saved.report.calibration_decimals, 2)
        self.assertEqual(saved.standards[2.5].lot_number, "standard lot")
        self.assertEqual(saved.calculation.fit_model, "4pl")

    def test_new_formula_assay_rechecks_permissions_on_save(self):
        path = self.incoming_protocol([sample for sample in self.samples if sample[0] != "Standard"])
        editor = AssayEditor(None, self.assays, self.access, protocol_path=path)
        self.addCleanup(editor.deleteLater)
        editor.lot_number.setText("lot")
        self.edit_formula(editor.calculation.formula, "Patient / ", "Grens")
        self.access.editable = False
        with patch("elisaflow.ui.assay_editor.QMessageBox.warning") as warning:
            editor._save()
        warning.assert_called_once()
        self.assertFalse((self.directory / path.name).exists())
        self.access.editable = True
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        saved = self.assays.load_assay(path.name)
        self.assertEqual(saved.calculation.method.value, "formula")
        self.assertEqual(saved.calculation.formula.references["grens"], "Grens")

    def test_calculation_switch_preserves_drafts_and_supports_fixed_curves_without_standards(self):
        path = self.incoming_protocol([sample for sample in self.samples if sample[0] != "Standard"])
        editor = AssayEditor(None, self.assays, self.access, protocol_path=path)
        self.addCleanup(editor.deleteLater)
        editor.lot_number.setText("kit lot")
        positive = editor.controls.editors["Positief"]
        positive.lot_number.setText("control lot")
        self.edit_formula(editor.calculation.formula, "Patient * 0.43")
        self.edit_formula(editor.calculation.overrides["Positief"].formula, "Positief")
        editor.calculation_method.setCurrentText("Standaardcurve")
        editor.calculation.curve.fixed_parameters.setChecked(True)
        for name, value in {"A": 0.1, "B": 1, "C": 90, "D": 2.1}.items():
            editor.calculation.curve.parameters[name].setText(str(value))
        editor.calculation.result_unit.setText("U/ml")
        self.edit_formula(editor.calculation.curve.preprocessing, "Monster / ", "Grens")
        self.assertEqual(editor.calculation.curve.preprocessing.summary.text(), "Monster / Grens")
        editor.calculation_method.setCurrentText("Berekening")
        self.assertEqual(editor.calculation.formula.to_config()["expression"], "sample * 0.43")
        self.assertEqual(editor.calculation.overrides["Positief"].formula.to_config()["expression"], "sample")
        editor.calculation_method.setCurrentText("Standaardcurve")
        self.assertTrue(editor.calculation.curve.fixed_parameters.isChecked())
        self.assertEqual(editor.calculation.curve.parameters["C"].text(), "90")
        self.assertEqual(editor.calculation.curve.preprocessing.to_config()["references"], {"grens": "Grens"})
        editor._save()
        self.assertEqual(editor.result(), editor.DialogCode.Accepted)
        saved = self.assays.load_assay(path.name)
        self.assertEqual(saved.calculation.parameters.C, 90)
        self.assertEqual(saved.result.unit, "U/ml")
        self.assertEqual(saved.controls["Positief"].lot_number, "control lot")
        self.assertFalse(saved.standards)

    def test_management_actions_allow_password_prompts_and_require_selection(self) -> None:
        for editable in (True, False):
            with self.subTest(editable=editable):
                self.access.editable = editable
                window = AssaySettingsWidget(self.assays, self.access)
                self.addCleanup(window.deleteLater)
                self.assertEqual(window.view_button.text(), "Bekijk assay")
                self.assertEqual(window.edit_button.text(), "Verander assay")
                self.assertEqual(window.new_button.text(), "Nieuwe assay")
                self.assertFalse(window.view_button.isEnabled())
                self.assertFalse(window.edit_button.isEnabled())
                self.assertTrue(window.new_button.isEnabled())
                window.assay_list.setCurrentRow(0)
                self.assertTrue(window.view_button.isEnabled())
                self.assertTrue(window.edit_button.isEnabled())
                window.assay_list.setCurrentRow(-1)
                self.assertFalse(window.view_button.isEnabled())
                self.assertFalse(window.edit_button.isEnabled())

    def test_view_button_and_double_click_never_allow_changes_even_for_editors(self) -> None:
        original = self.source.read_bytes()
        viewed = []

        def inspect_editor(editor: AssayEditor) -> int:
            self.addCleanup(editor.deleteLater)
            self.assertEqual(editor.windowTitle(), "Bekijk assay: test")
            self.assertTrue(editor.lot_number.isReadOnly())
            self.assertFalse(editor.controls.editors["Positief"].isEnabled())
            self.assertFalse(editor.blank.isEnabled())
            self.assertFalse(editor.patients.replicate_filter.enabled.isEnabled())
            self.assertFalse(editor.standards.editors[2.5].isEnabled())
            self.assertFalse(editor.calculation.formula.isEnabled())
            self.assertFalse(editor.calculation.formula.edit_button.isEnabled())
            self.assertFalse(editor.report.isEnabled())
            buttons = editor.findChild(QDialogButtonBox)
            self.assertIsNone(buttons.button(QDialogButtonBox.StandardButton.Save))
            self.assertTrue(buttons.button(QDialogButtonBox.StandardButton.Close).isEnabled())
            # The view mode must also block writes if the save handler is invoked directly.
            editor.lot_number.setText("accidental change")
            editor._save()
            self.assertNotEqual(editor.result(), editor.DialogCode.Accepted)
            viewed.append(editor)
            editor.reject()
            return editor.result()

        with (
            patch.object(AssayEditor, "exec", inspect_editor),
            patch.object(self.assays, "save_assay", wraps=self.assays.save_assay) as save,
        ):
            for editable in (True, False):
                with self.subTest(editable=editable):
                    self.access.editable = editable
                    window = AssaySettingsWidget(self.assays, self.access)
                    self.addCleanup(window.deleteLater)
                    window.assay_list.setCurrentRow(0)
                    window.view_button.click()
                    window.assay_list.itemDoubleClicked.emit(window.assay_list.currentItem())
                    self.assertEqual(self.source.read_bytes(), original)
        self.assertEqual(len(viewed), 4)
        save.assert_not_called()
        self.assertFalse((self.configs / ".history").exists())

    def test_change_button_opens_editable_assay_and_can_save(self) -> None:
        window = AssaySettingsWidget(self.assays, self.access)
        self.addCleanup(window.deleteLater)
        window.assay_list.setCurrentRow(0)

        def change_editor(editor: AssayEditor) -> int:
            self.addCleanup(editor.deleteLater)
            self.assertEqual(editor.windowTitle(), "Verander assay: test")
            self.assertFalse(editor.lot_number.isReadOnly())
            self.assertTrue(editor.controls.editors["Positief"].isEnabled())
            buttons = editor.findChild(QDialogButtonBox)
            self.assertTrue(buttons.button(QDialogButtonBox.StandardButton.Save).isEnabled())
            editor.lot_number.setText("updated kit lot")
            editor._save()
            self.assertEqual(editor.result(), editor.DialogCode.Accepted)
            return editor.result()

        with (
            patch.object(AssayEditor, "exec", change_editor),
            patch("elisaflow.ui.assay_settings.QMessageBox.warning") as warning,
        ):
            window.edit_button.click()
        warning.assert_not_called()
        saved = self.assays.load_assay("test.skax")
        self.assertEqual(saved.lot_number, "updated kit lot")
        self.assertEqual(saved.controls, self.settings.controls)

    def test_change_button_rechecks_permission_before_opening_editor(self) -> None:
        window = AssaySettingsWidget(self.assays, self.access)
        self.addCleanup(window.deleteLater)
        window.assay_list.setCurrentRow(0)
        self.assertTrue(window.edit_button.isEnabled())
        original = self.source.read_bytes()
        self.access.editable = False
        with (
            patch("elisaflow.ui.assay_settings.AssayEditor") as editor,
            patch("elisaflow.ui.assay_settings.QMessageBox.warning") as warning,
        ):
            window.edit_button.click()
        editor.assert_not_called()
        warning.assert_called_once()
        self.assertEqual(self.source.read_bytes(), original)

    def test_new_assay_requires_edit_permission_before_showing_picker(self) -> None:
        self.access.editable = False
        window = AssaySettingsWidget(self.assays, self.access)
        self.addCleanup(window.deleteLater)
        self.assertTrue(window.new_button.isEnabled())
        with (
            patch("elisaflow.ui.assay_settings.QFileDialog.getOpenFileName") as picker,
            patch("elisaflow.ui.assay_settings.QMessageBox.warning") as warning,
        ):
            window._new_assay()
        picker.assert_not_called()
        warning.assert_called_once()
        self.assertEqual(list(self.configs.iterdir()), [self.source])

    def test_new_assay_picker_cancellation_does_nothing(self) -> None:
        window = AssaySettingsWidget(self.assays, self.access)
        self.addCleanup(window.deleteLater)
        with (
            patch("elisaflow.ui.assay_settings.QFileDialog.getOpenFileName", return_value=("", "")) as picker,
            patch("elisaflow.ui.assay_settings.AssayEditor") as editor,
        ):
            window._new_assay()
        picker.assert_called_once()
        self.assertIn("*.skax", picker.call_args.args[3])
        editor.assert_not_called()
        self.assertEqual(list(self.configs.iterdir()), [self.source])

    def test_new_assay_editor_cancellation_does_not_import_or_register(self) -> None:
        path = self.incoming_protocol()
        window = AssaySettingsWidget(self.assays, self.access)
        self.addCleanup(window.deleteLater)

        def cancel_editor(editor: AssayEditor) -> int:
            self.addCleanup(editor.deleteLater)
            editor.reject()
            return editor.result()

        with (
            patch("elisaflow.ui.assay_settings.QFileDialog.getOpenFileName", return_value=(str(path), "")),
            patch.object(AssayEditor, "exec", cancel_editor),
        ):
            window._new_assay()
        self.assertNotIn(path.name, self.assays.list_assays())
        self.assertFalse((self.directory / path.name).exists())
        self.assertEqual(list(self.configs.iterdir()), [self.source])

    def test_new_assay_picker_save_updates_and_selects_list_entry(self) -> None:
        path = self.incoming_protocol(name="new_assay.skax")
        window = AssaySettingsWidget(self.assays, self.access)
        self.addCleanup(window.deleteLater)
        self.assertTrue(window.new_button.isEnabled())

        def save_editor(editor: AssayEditor) -> int:
            self.addCleanup(editor.deleteLater)
            editor.lot_number.setText("new kit lot")
            editor._save()
            return editor.result()

        with (
            patch("elisaflow.ui.assay_settings.QFileDialog.getOpenFileName", return_value=(str(path), "")),
            patch.object(AssayEditor, "exec", save_editor),
            patch("elisaflow.ui.assay_settings.QMessageBox.warning") as warning,
        ):
            window.new_button.click()
        warning.assert_not_called()
        self.assertEqual(window.assay_list.currentItem().text(), "new assay")
        self.assertEqual(window.assay_list.count(), 2)
        self.assertTrue(window.view_button.isEnabled())
        self.assertTrue(window.edit_button.isEnabled())
        self.assertEqual(self.assays.load_assay(path.name).lot_number, "new kit lot")

    def test_new_assay_invalid_selection_reports_error_and_preserves_library(self) -> None:
        path = self.incoming_protocol(name="test.skax")
        window = AssaySettingsWidget(self.assays, self.access)
        self.addCleanup(window.deleteLater)
        original = self.source.read_bytes()
        with (
            patch("elisaflow.ui.assay_settings.QFileDialog.getOpenFileName", return_value=(str(path), "")),
            patch("elisaflow.ui.assay_settings.QMessageBox.warning") as warning,
        ):
            window._new_assay()
        warning.assert_called_once()
        self.assertEqual(self.source.read_bytes(), original)
        self.assertEqual(window.assay_list.count(), 1)

    def test_open_errors_are_reported_without_changing_settings(self) -> None:
        window = AssaySettingsWidget(self.assays, self.access)
        self.addCleanup(window.deleteLater)
        window.assay_list.setCurrentRow(0)
        original = self.source.read_bytes()
        for samples in ([sample for sample in self.samples if sample[1] != "Positief"], None):
            with self.subTest(samples=samples):
                if samples is None:
                    self.path.write_bytes(b"invalid archive")
                else:
                    write_protocol(self.path, samples)
                with patch("elisaflow.ui.assay_settings.QMessageBox.warning") as warning:
                    window._open_assay()
                warning.assert_called_once()
                self.assertEqual(self.source.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
