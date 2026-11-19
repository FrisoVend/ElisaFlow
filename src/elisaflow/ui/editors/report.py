"""Edit QC, calculation and report settings for a SKAX protocol."""

from PySide6.QtWidgets import (
    QCheckBox,
    QLineEdit,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from elisaflow.configuration.assay_settings import ReportSettings
from elisaflow.ui.editors.fields import _number
from elisaflow.ui.theme import FormSection, helper_text


class ReportEditor(QWidget):
    def __init__(self, settings: ReportSettings, editable: bool) -> None:
        super().__init__()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(16)
        layout.addWidget(helper_text("Aantal decimalen in de PDF: 0 geeft gehele getallen. De berekening en QC-beoordeling gebruiken de onafgeronde waarden."))
        decimal_section = FormSection("Decimalen")
        layout.addWidget(decimal_section)
        form = self._decimal_form = decimal_section.form
        self.decimals = {}
        for name, label in (
            ("result_decimals", "Resultaten:"),
            ("extinction_decimals", "Extincties:"),
            ("calibration_decimals", "Kalibratieconcentraties:"),
            ("percentage_decimals", "Percentages:"),
        ):
            field = QSpinBox()
            field.setButtonSymbols(QSpinBox.ButtonSymbols.PlusMinus)
            field.setFixedWidth(100)
            field.setRange(0, ReportSettings.MAX_DECIMALS)
            field.setValue(getattr(settings, name))
            form.addRow(label, field)
            self.decimals[name] = field
        warning_section = FormSection("PDF-waarschuwingen")
        layout.addWidget(warning_section)
        form = self._warning_form = warning_section.form
        form.addRow(helper_text("Waarden buiten deze grenzen worden rood weergegeven in de PDF. Dit verandert de QC-status niet."))
        self.cv_warning_percent = QLineEdit(str(settings.cv_warning_percent))
        self.cv_warning_percent.setToolTip("CV waarden boven dit percentage worden roodgekleurd.")
        form.addRow("CV waarschuwing (%):", self.cv_warning_percent)
        self.recovery_warning_percent = QLineEdit(str(settings.recovery_warning_percent))
        self.recovery_warning_percent.setToolTip("Toegestane afwijking van 100% terugvinding: bij 20% worden waarden onder 80% of boven 120% rood.")
        form.addRow("Terugvinding afwijking (%):", self.recovery_warning_percent)
        self.show_result_cv = QCheckBox("CV van resultaten tonen")
        self.show_result_cv.setChecked(settings.show_result_cv)
        form.addRow(self.show_result_cv)
        layout.addStretch()
        self.setEnabled(editable)

    def set_calibration(self, calibration: bool) -> None:
        self._decimal_form.setRowVisible(self.decimals["calibration_decimals"], calibration)
        self._warning_form.setRowVisible(self.recovery_warning_percent, calibration)

    def to_config(self) -> dict:
        return ReportSettings(
            **{name: field.value() for name, field in self.decimals.items()},
            show_result_cv=self.show_result_cv.isChecked(),
            cv_warning_percent=_number(self.cv_warning_percent, "CV-waarschuwing"),
            recovery_warning_percent=_number(self.recovery_warning_percent, "Afwijking van de terugvinding"),
        ).to_dict()
