import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from fixtures import add_measurement_steps
from skax_fixtures import write_skax

from elisaflow.readers.skax_reader import SkaxReader


class SkaxReaderTests(unittest.TestCase):
    def test_bundled_template_supplies_layout_and_software_metadata_without_measurements(self) -> None:
        protocol = Path(__file__).resolve().parents[1] / "protocols" / "Borellia_C10_as.skax"
        reader = SkaxReader.from_file(protocol)
        self.assertIn("Kalibratie", reader.wells["sample_name"].tolist())
        self.assertEqual(reader.metadata["session"]["software"]["version"], "8.0.0.12")
        self.assertEqual(reader.path, protocol.resolve())
        self.assertFalse(reader.completed)
        self.assertTrue(reader.measurements.empty)

    def test_reading_replaced_session_supplies_independent_results(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "test.skax"
            write_skax(path)
            completed = SkaxReader.from_file(path)
            write_skax(path, completed=False)
            template = SkaxReader.from_file(path)
        self.assertTrue(completed.completed)
        self.assertEqual(len(completed.measurements), 3)
        self.assertFalse(template.completed)
        self.assertTrue(template.measurements.empty)
        self.assertEqual(template.metadata["session"]["execution"], {"start": None, "end": None})

    def test_dilution_preserves_numeric_factors_and_defaults_only_when_unspecified(self) -> None:
        for value, expected in ((None, 1), ("", 1), ("0", 0), ("1", 1), ("2.5", 2.5), ("bad", None)):
            session = ET.fromstring("""
                <Session><PlateLayout><Plates><Plate><Id>plate</Id><Name>plate</Name>
                <Wells><Well><WellLocationX>0</WellLocationX><WellLocationY>0</WellLocationY>
                <Replicate><Contents><Sample><Id>patient</Id><Name>Patient</Name>
                <SampleType>Unknown</SampleType></Sample></Contents></Replicate>
                </Well></Wells></Plate></Plates></PlateLayout></Session>
            """)
            if value is not None:
                ET.SubElement(session.find(".//Sample"), "Dilution").text = value
            add_measurement_steps(session)
            with self.subTest(dilution=value), TemporaryDirectory() as directory:
                path = Path(directory) / "test.skax"
                with ZipFile(path, "w") as archive:
                    archive.writestr("Header.xml", "<Header />")
                    archive.writestr("Session.xml", ET.tostring(session))
                if expected is None:
                    with self.assertRaisesRegex(ValueError, 'Verdunningsfactoren moeten numeriek zijn'):
                        SkaxReader.from_file(path)
                else:
                    reader = SkaxReader.from_file(path)
                    self.assertEqual(reader.wells.iloc[0]["sample_dilution"], expected)

    def test_wells_read_down_columns_in_numeric_order_on_each_plate(self) -> None:
        header = ET.Element("Header")
        indexes = ET.SubElement(header, "Plates")
        session = ET.Element("Session")
        add_measurement_steps(session)
        plates = ET.SubElement(ET.SubElement(session, "PlateLayout"), "Plates")
        for plate_index in (1, 0):
            plate_id = f"plate-{plate_index}"
            pair = ET.SubElement(indexes, "Pair")
            ET.SubElement(pair, "Key").text = str(plate_index)
            ET.SubElement(pair, "Value").text = plate_id
            plate = ET.SubElement(plates, "Plate")
            ET.SubElement(plate, "Id").text = plate_id
            ET.SubElement(plate, "Name").text = plate_id
            wells = ET.SubElement(plate, "Wells")
            for row, column in ((0, 9), (1, 1), (0, 1), (7, 0), (1, 0), (0, 0)):
                well = ET.SubElement(wells, "Well")
                ET.SubElement(well, "WellLocationX").text = str(column)
                ET.SubElement(well, "WellLocationY").text = str(row)
                sample = ET.SubElement(ET.SubElement(ET.SubElement(well, "Replicate"), "Contents"), "Sample")
                ET.SubElement(sample, "Id").text = f"sample-{plate_index}-{row}-{column}"
                ET.SubElement(sample, "Name").text = f"Patient {row}-{column}"
                ET.SubElement(sample, "SampleType").text = "Unknown"
        with TemporaryDirectory() as directory:
            path = Path(directory) / "test.skax"
            with ZipFile(path, "w") as archive:
                archive.writestr("Header.xml", ET.tostring(header))
                archive.writestr("Session.xml", ET.tostring(session))
            reader = SkaxReader.from_file(path)
        self.assertEqual(reader.wells["plate_index"].tolist(), [0] * 6 + [1] * 6)
        self.assertEqual(reader.wells["well"].tolist(), ["A1", "B1", "H1", "A2", "B2", "A10"] * 2)


if __name__ == "__main__":
    unittest.main()
