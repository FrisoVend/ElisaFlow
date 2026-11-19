import base64
import json
import stat
import unittest
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from elisaflow.assays.registry import AssayRegistry
from elisaflow.calibration.curves import create_calibration_curve
from elisaflow.configuration.assay_settings import ValidationScope
from elisaflow.xml_utils import xml_attribute


class RegistryTests(unittest.TestCase):
    def test_library_protection_leaves_other_application_files_writable(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            configs = directory / "configs"
            configs.mkdir()
            source = configs / "assay.json"
            source.write_text(json.dumps({"test.skax": {
                "lot_number": "lot",
                "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample"},
            }}), encoding="utf-8")
            other_files = (directory / "app.toml", directory / "settings.json", configs / "readme.txt")
            for path in other_files:
                path.write_text("unchanged", encoding="utf-8")
            registry = AssayRegistry(configs)
            self.assertEqual(registry.get("test.skax").lot_number, "lot")
            self.assertFalse(source.stat().st_mode & stat.S_IWUSR)
            for path in other_files:
                self.assertTrue(path.stat().st_mode & stat.S_IWUSR)
                self.assertEqual(path.read_text(encoding="utf-8"), "unchanged")

    def test_switching_assays_match_actual_protocol_filters(self) -> None:
        project = Path(__file__).resolve().parents[1]
        registry = AssayRegistry(project / "configs")
        names = [name for name in registry.protocol_names() if registry.get(name).measurement.switching]
        expected_cutoffs = {"test_GAD_as.skax": 120, "test_IA-2_as.skax": 120, "test_ZnT8_as.skax": 100}
        self.assertTrue(names)
        for name in names:
            with self.subTest(assay=name):
                settings = registry.get(name)
                self.assertTrue(settings.measurement.switching)
                self.assertEqual(settings.measurement.preferred_wavelength_nm, 450)
                self.assertEqual(settings.measurement.cutoff_concentration, expected_cutoffs[name])
                self.assertTrue(all(curve.fit_model == "4pl" for curve in settings.calculation.wavelengths.values()))
                self.assertLessEqual(set(settings.validation.wavelengths), {405, 450})
                self.assertFalse(settings.validation.result.blank)
                self.assertTrue(settings.validation.wavelengths[450].blank)
                self.assertEqual(settings.validation.wavelengths.get(405, ValidationScope()).to_dict(), {})
                with ZipFile(project / "protocols" / name) as archive:
                    raw = archive.read("Session.xml")
                    self.assertEqual(
                        ET.fromstring(archive.read("Header.xml")).findtext("SessionHash"), base64.b64encode(sha256(raw).digest()).decode("ascii")
                    )
                session = ET.fromstring(raw)
                references = {xml_attribute(node, "Id"): node for node in session.iter() if xml_attribute(node, "Id") is not None}
                hardware = {int(node.findtext("{*}Wavelength")): node for node in session.findall(".//{*}Filter")}
                steps = [
                    step for step in session.find("{*}Protocol/{*}ProtocolSteps") if step.findtext("{*}Type") == "PhotometricFilterMeasurementStep"
                ]
                self.assertEqual(len(steps), 2)
                self.assertEqual(len(session.find("{*}Protocol/{*}ResultSteps")), 0)
                wavelengths = []
                for step in steps:
                    filters = step.findall("{*}Filters/{*}PhotometricFilter")
                    self.assertEqual(len(filters), 1)
                    measurement = filters[0]
                    wavelength = int(measurement.findtext("{*}Wavelength"))
                    wavelengths.append(wavelength)
                    self.assertEqual(step.findtext("{*}Name"), f"{wavelength}nm")
                    self.assertEqual(measurement.findtext("{*}FilterPosition"), hardware[wavelength].findtext("{*}Position"))
                    filter_name = measurement.find("{*}Name")
                    self.assertEqual(references[xml_attribute(filter_name, "Ref")].text, hardware[wavelength].findtext("{*}Name"))
                self.assertEqual(set(wavelengths), {405, 450})


    def test_every_protocol_has_registered_settings(self) -> None:
        project = Path(__file__).resolve().parents[1]
        registry = AssayRegistry(project / "configs")
        for protocol in (project / "protocols").glob("*.skax"):
            with self.subTest(protocol=protocol.name):
                self.assertEqual(registry.get(protocol).protocol_name, protocol.name)

    def test_every_bundled_configuration_loads(self) -> None:
        directory = Path(__file__).resolve().parents[1] / "configs"
        registry = AssayRegistry(directory)
        for path in directory.glob("*.json"):
            for protocol_name, config in json.loads(path.read_text(encoding="utf-8-sig")).items():
                with self.subTest(protocol=protocol_name):
                    settings = registry.get(Path("protocols") / protocol_name)
                    self.assertEqual(settings.protocol_name, protocol_name)
                    self.assertEqual(settings.result.unit, config.get("result", {}).get("unit"))
                    if settings.calculation.method.value == "calibration" and not settings.measurement.switching and settings.calculation.parameters:
                        curve = create_calibration_curve(settings.calculation)
                        self.assertAlmostEqual(curve.calculate_concentration(curve.calculate_extinction(1.0)), 1.0)

    def test_invalid_settings_include_file_and_protocol_context(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            path = directory / "bad.json"
            for invalid in (
                {"response_scale": "raw", "fit_model": []},
                {"response_scale": "raw", "fit_model": "4pl", "weighting_exponent": {}},
                {"response_scale": "blank_corrected", "fit_model": "4pl", "parameters": {}},
            ):
                if path.exists():
                    path.chmod(path.stat().st_mode | stat.S_IWUSR)
                path.write_text(
                    json.dumps(
                        {
                            "test.skax": {
                                "lot_number": "lot",
                                "calculation": {"method": "calibration", **invalid},
                            }
                        }
                    ),
                    encoding="utf-8",
                )
                with self.subTest(settings=invalid), self.assertRaises(ValueError) as caught:
                    AssayRegistry(directory)
                self.assertIn("bad.json", str(caught.exception))
                self.assertIn("test.skax", str(caught.exception))

    def test_duplicate_protocols_are_rejected_and_bom_is_supported(self) -> None:
        config = {
            "test.skax": {
                "lot_number": "lot",
                "calculation": {"response_scale": "blank_corrected", "method": "formula", "expression": "sample / c", "references": {"c": "c"}},
                "result": {"label": "Ratio"},
            }
        }
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "first.json").write_text(json.dumps(config), encoding="utf-8-sig")
            self.assertEqual(AssayRegistry(directory).get("test.skax").calculation.formula.references["c"], "c")
            (directory / "second.json").write_text(json.dumps(config), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, 'Dubbel protocol'):
                AssayRegistry(directory)


if __name__ == "__main__":
    unittest.main()
