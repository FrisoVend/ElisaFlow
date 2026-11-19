import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from skax_fixtures import write_members, write_skax

from elisaflow.assays.protocol_validation import validate_protocol_measurement
from elisaflow.assays.registry import AssayRegistry
from elisaflow.configuration.application import ApplicationSettings
from elisaflow.configuration.assay_settings import MeasurementSettings
from elisaflow.readers.skax_document import SkaxDocument
from elisaflow.runs.workflow import run_workflow

ROOT = Path(__file__).resolve().parents[1]


class ProtocolValidationTests(unittest.TestCase):
    def setUp(self):
        self.registry = AssayRegistry(ROOT / "configs")
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)

    def protocol(self, name="test_ZnT8_as.skax", change=lambda session: None):
        source = ROOT / "protocols" / name
        with ZipFile(source) as archive:
            session = ET.fromstring(archive.read("Session.xml"))
            header = archive.read("Header.xml")
        change(session)
        target = self.directory / name
        write_members(target, {"Session.xml": ET.tostring(session), "Header.xml": header})
        return target

    def test_all_bundled_protocols_match_their_current_configs(self):
        for name in self.registry.protocol_names():
            with self.subTest(protocol=name):
                validate_protocol_measurement(ROOT / "protocols" / name, self.registry.get(name).measurement)

    def test_bundled_templates_disable_automatic_report_export(self):
        for path in (ROOT / "protocols").glob("*.skax"):
            with self.subTest(protocol=path.name):
                document = SkaxDocument.from_file(path)
                document.validate_hashes()
                inner = ET.fromstring(document.session.find("{*}XmlData").text)
                self.assertEqual(inner.find(".//{*}AutoExportSettings/{*}IsSaveFileChecked").text, "false")

    def test_live_templates_require_no_measurements_and_no_export_settings(self):
        template = write_skax(self.directory / "session.skax", completed=False)
        validate_protocol_measurement(template, MeasurementSettings())
        write_skax(template)
        with self.assertRaisesRegex(ValueError, "meetgegevens"):
            validate_protocol_measurement(template, MeasurementSettings())

    def test_switching_uses_filter_metadata_and_ignores_step_names(self):
        def rename(session):
            for step in session.findall("{*}Protocol/{*}ProtocolSteps/*"):
                step.find("{*}Name").text = "Uninformative label"

        validate_protocol_measurement(self.protocol(change=rename), MeasurementSettings("switch_on_concentration", 450, cutoff_concentration=1000))

    def test_rejects_a_mislabeled_duplicate_filter_and_wrong_hardware_position(self):
        def duplicate(session):
            node = session.find(".//{*}PhotometricFilter[{*}Wavelength='450']")
            node.find("{*}Wavelength").text = "405"
            node.find("{*}FilterPosition").text = "1"

        def wrong_position(session):
            session.find(".//{*}PhotometricFilter/{*}FilterPosition").text = "3"

        for change, message in ((duplicate, "verschillende ruwe"), (wrong_position, "filterpositie")):
            with self.subTest(change=change.__name__), self.assertRaisesRegex(ValueError, message):
                validate_protocol_measurement(self.protocol(change=change), MeasurementSettings("switch_on_concentration", 450, cutoff_concentration=1000))

    def test_filter_mismatch_stops_workflow_before_creating_files_or_launching_skanit(self):
        configs = self.directory / "configs"
        configs.mkdir()
        settings = self.registry.get("test_ZnT8_as.skax").to_dict()
        settings["measurement"] = {"mode": "single"}
        settings["calculation"].update(settings["calculation"].pop("wavelengths")["450"])
        settings["validation"].pop("wavelengths")
        (configs / "test.json").write_text(json.dumps({"test_ZnT8_as.skax": settings}), encoding="utf-8")
        config = ApplicationSettings(
            protocols_dir=ROOT / "protocols",
            configs_dir=configs,
            results_dir=self.directory / "results",
            password_file=self.directory / "editor.hash",
        )
        with patch("elisaflow.runs.workflow.open_protocol") as launch, self.assertRaisesRegex(ValueError, 'precies 1'):
            run_workflow(config, protocol_template=ROOT / "protocols/test_ZnT8_as.skax", open_report=False, analyst="operator")
        launch.assert_not_called()
        self.assertFalse(config.results_dir.exists())
