"""Measurement controls require an explicit cutoff and fixed optical coefficients."""

import unittest

try:
    from PySide6.QtWidgets import QApplication

    from elisaflow.ui.editors.measurement import MeasurementEditor
except ImportError:
    QApplication = None

from elisaflow.configuration.assay_settings import MeasurementSettings, SaturationSettings


@unittest.skipIf(QApplication is None, "PySide6 is not installed")
class MeasurementEditorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def editor(self, settings=None, *, editable=True, wavelengths_nm=(405, 450)):
        editor = MeasurementEditor(settings or MeasurementSettings(), editable, wavelengths_nm=wavelengths_nm)
        self.addCleanup(editor.deleteLater)
        return editor

    def test_switching_requires_user_entered_cutoff(self):
        editor = self.editor()
        editor.mode.setCurrentIndex(editor.mode.findData("switch_on_concentration"))
        with self.assertRaisesRegex(ValueError, "Omschakelconcentratie"):
            editor.to_config()
        editor.cutoff_concentration.setText("100")
        self.assertEqual(editor.to_config()["cutoff_concentration"], 100)
        editor.mode.setCurrentIndex(editor.mode.findData("single"))
        self.assertNotIn("cutoff_concentration", editor.selection().to_dict())
        editor.mode.setCurrentIndex(editor.mode.findData("switch_on_concentration"))
        self.assertEqual(editor.to_config()["cutoff_concentration"], 100)

    def test_conversion_requires_source_slope_and_intercept_and_preserves_drafts(self):
        editor = self.editor(MeasurementSettings("switch_on_concentration", 450, cutoff_concentration=100))
        conversion = editor.saturation[450]
        self.assertEqual(conversion.mode.currentData(), "exclude")
        self.assertEqual(conversion.source_wavelength_nm.count(), 1)
        self.assertEqual(conversion.source_wavelength_nm.currentData(), 405)
        conversion.mode.setCurrentIndex(conversion.mode.findData("calculate"))
        with self.assertRaises(ValueError):
            editor.to_config()
        conversion.slope.setText("3.4")
        conversion.intercept.setText("-0.01")
        settings = MeasurementSettings.from_dict(editor.to_config())
        self.assertEqual(settings.saturation[450], SaturationSettings("calculate", 405, 3.4, -0.01))
        reopened = self.editor(settings)
        self.assertEqual(reopened.to_config(), editor.to_config())
        conversion.mode.setCurrentIndex(conversion.mode.findData("exclude"))
        self.assertNotIn("saturation", editor.to_config())
        conversion.mode.setCurrentIndex(conversion.mode.findData("calculate"))
        self.assertEqual(editor.to_config(), settings.to_dict())

    def test_read_only_editor_disables_conversion_and_cutoff_edits(self):
        editor = self.editor(MeasurementSettings("switch_on_concentration", 450, cutoff_concentration=100), editable=False)
        self.assertTrue(editor.cutoff_concentration.isReadOnly())
        self.assertFalse(editor.mode.isEnabled())
        self.assertFalse(editor.saturation[450].isEnabled())

    def test_reference_is_excluded_from_conversion_targets_and_sources(self):
        settings = MeasurementSettings(
            "switch_on_concentration", 450, reference_wavelength_nm=620, cutoff_concentration=100,
            saturation={450: SaturationSettings("calculate", 405, 3.4, 0)},
        )
        editor = self.editor(settings, wavelengths_nm=(405, 450, 620))
        self.assertEqual(editor.to_config(), settings.to_dict())
        self.assertTrue(editor.saturation[620].isHidden())
        self.assertEqual(editor.saturation[450].source_wavelength_nm.count(), 1)
        self.assertEqual(editor.saturation[450].source_wavelength_nm.currentData(), 405)
        editor.reference_wavelength_nm.setCurrentIndex(editor.reference_wavelength_nm.findData(405))
        self.assertTrue(editor.saturation[405].isHidden())
        self.assertFalse(editor.saturation[620].isHidden())
        self.assertEqual(editor.saturation[450].source_wavelength_nm.count(), 1)
        self.assertEqual(editor.saturation[450].source_wavelength_nm.currentData(), 620)
        self.assertEqual(editor.to_config()["saturation"]["450"]["source_wavelength_nm"], 620)
        editor.reference_wavelength_nm.setCurrentIndex(editor.reference_wavelength_nm.findData(None))
        self.assertFalse(editor.saturation[405].isHidden())
        self.assertEqual(editor.saturation[450].source_wavelength_nm.count(), 2)


if __name__ == "__main__":
    unittest.main()
