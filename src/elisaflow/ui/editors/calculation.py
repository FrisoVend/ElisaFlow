"""Own patient and control formula drafts and their result settings."""

from PySide6.QtWidgets import (
    QCheckBox,
    QLineEdit,
    QVBoxLayout,
    QWidget,
)

from elisaflow.configuration.assay_settings import MeasurementSettings, ResultSettings
from elisaflow.models import ResponseScale
from elisaflow.ui.editors.fields import _response_scale_selector, _result_config
from elisaflow.ui.editors.formula import FormulaEditor
from elisaflow.ui.theme import FormSection, editor_tabs, helper_text


class ControlCalculationEditor(QWidget):
    """Own one control's inheritance, formula and exported result settings."""

    def __init__(
        self, name: str, override: dict | None, result: dict | None,
        control_names: tuple[str, ...], references: dict[str, str], sample_result: dict,
        *, sample_response_scale: str, default_response_scale: ResponseScale,
    ):
        super().__init__()
        self._name = name
        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(16)
        section = FormSection(f"Resultaatberekening voor {name}")
        layout.addWidget(section)
        form = section.form
        self.use_sample_result = QCheckBox("Zelfde formule, invoerschaal, resultaatlabel en eenheid als patiënten")
        self.use_sample_result.setChecked(override is None and result is None)
        self._sample_response_scale = sample_response_scale
        self._custom_response_scale = (override or {}).get("response_scale", default_response_scale)
        self._custom_result = result is not None
        self.use_sample_result.toggled.connect(self._toggle_sample_result)
        form.addRow(self.use_sample_result)
        self.result_label = QLineEdit((result or sample_result).get("label", "Resultaat"))
        self.result_unit = QLineEdit((result or sample_result).get("unit", ""))
        self.result_unit.setPlaceholderText("Geen eenheid")
        form.addRow("Resultaatlabel:", self.result_label)
        form.addRow("Resultaateenheid:", self.result_unit)
        self.response_scale = _response_scale_selector(self._custom_response_scale)
        self.response_scale.setToolTip(f"{self.response_scale.toolTip()} De keuze geldt voor deze formule en alle gebruikte controlegemiddelden.")
        form.addRow("Responsschaal:", self.response_scale)
        self._sample_label = "Resultaat"
        self._sample_unit = ""
        self.formula = FormulaEditor(
            {"expression": override["expression"]} if override is not None and override.get("expression") is not None else None,
            control_names, optional=True,
            sample_name=name, references=references, title=f"Formule voor {name}", empty_text="Patientformule", clear_text="Herstellen",
        )
        form.addRow("Formule:", self.formula)
        self._update_result_fields()
        layout.addStretch()

    def set_sample_result(self, label: str, unit: str, response_scale: str) -> None:
        self._sample_label = label
        self._sample_unit = unit
        self._sample_response_scale = response_scale
        self._update_result_fields()

    def _toggle_sample_result(self, use_sample: bool) -> None:
        if use_sample:
            self._custom_response_scale = self.response_scale.currentData()
        else:
            self.response_scale.setCurrentIndex(self.response_scale.findData(self._custom_response_scale))
        self._custom_result = not use_sample
        self._update_result_fields()

    def _update_result_fields(self) -> None:
        use_sample = self.use_sample_result.isChecked()
        self.result_label.setEnabled(not use_sample)
        self.result_unit.setEnabled(not use_sample)
        self.formula.setEnabled(not use_sample)
        self.response_scale.setEnabled(not use_sample)
        if use_sample:
            self.result_label.setText(self._sample_label)
            self.result_unit.setText(self._sample_unit)
            self.formula.clear()
            self.response_scale.setCurrentIndex(self.response_scale.findData(self._sample_response_scale))

    def to_config(self) -> tuple[dict | None, dict | None]:
        """Export the override with its used references and independent result settings."""
        if self.use_sample_result.isChecked():
            return None, None
        try:
            formula = self.formula.to_config()
        except ValueError as error:
            raise ValueError(f"Afwijkende formule voor controle {self._name!r}: {error}") from error
        override = {"response_scale": self.response_scale.currentData(), **(formula or {})}
        label = self.result_label.text().strip()
        unit = self.result_unit.text().strip()
        result = None
        if self._custom_result or (label, unit) != (self._sample_label, self._sample_unit):
            try:
                result = ResultSettings(label, unit or None).to_dict()
            except ValueError as error:
                raise ValueError(f"Resultaatinstellingen voor controle {self._name!r}: {error}") from error
        return override, result


class FormulaCalculationEditor(QWidget):
    """Keep formula drafts and control overrides independent of calibration."""

    def __init__(self, config: dict, control_names: tuple[str, ...], editable: bool, *, has_blank: bool = False) -> None:
        super().__init__()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(16)
        calculation = config.get("calculation", {})
        if calculation.get("method") != "formula":
            calculation = {"expression": "sample"}
        self.response_scale = _response_scale_selector(calculation.get("response_scale", ResponseScale.default_for_blank(has_blank)))
        self.response_scale.setToolTip(
            f"{self.response_scale.toolTip()} De keuze geldt voor de patiëntformule en alle controlegemiddelden waarnaar deze verwijst."
        )
        self.response_scale.setEnabled(editable)
        result = config.get("result", {}) if config.get("calculation", {}).get("method") == "formula" else {}
        self.result_label = QLineEdit(result.get("label", "Resultaat"))
        self.result_unit = QLineEdit(result.get("unit", ""))
        self.tabs = editor_tabs()
        patients = QWidget()
        patient_layout = QVBoxLayout(patients)
        patient_layout.setContentsMargins(16, 16, 16, 16)
        patient_section = FormSection("Resultaatberekening")
        patient_layout.addWidget(patient_section)
        patient_layout.addStretch()
        patient_form = patient_section.form
        patient_form.addRow("Resultaatlabel:", self.result_label)
        patient_form.addRow("Resultaateenheid:", self.result_unit)
        patient_form.addRow("Responsschaal:", self.response_scale)
        self.tabs.addTab(patients, "Patiënten")
        controls = QWidget()
        control_layout = QVBoxLayout(controls)
        control_layout.setContentsMargins(0, 0, 0, 0)
        control_layout.setSpacing(0)
        self.control_tabs = editor_tabs()
        self.control_tabs.setVisible(bool(control_names))
        control_layout.addWidget(self.control_tabs)
        if not control_names:
            control_layout.addWidget(helper_text("Geen controles gedefinieerd in het SKAX-protocol."))
        self.tabs.addTab(controls, "Controles")
        layout.addWidget(self.tabs)
        self.formula = FormulaEditor(calculation, control_names, sample_name="Patient", title="Resultaatformule")
        patient_form.addRow("Formule:", self.formula)
        self.overrides = {}
        sample_result = {"label": self.result_label.text(), "unit": self.result_unit.text()}
        for name in control_names:
            override = calculation.get("control_overrides", {}).get(name)
            editor = ControlCalculationEditor(
                name, override, calculation.get("control_results", {}).get(name), control_names, self.formula.references, sample_result,
                sample_response_scale=self.response_scale.currentData(), default_response_scale=ResponseScale.default_for_blank(has_blank),
            )
            self.overrides[name] = editor
            self.control_tabs.addTab(editor, name)
        self.result_label.textChanged.connect(self._update_control_results)
        self.result_unit.textChanged.connect(self._update_control_results)
        self.response_scale.currentIndexChanged.connect(self._update_control_results)
        self._update_control_results()
        patients.setEnabled(editable)
        for editor in self.overrides.values():
            editor.setEnabled(editable)

    def set_measurement(self, settings: MeasurementSettings, wavelengths: tuple[int, ...]) -> None:
        """Formula drafts do not depend on the measurement selection."""

    def _update_control_results(self) -> None:
        label = self.result_label.text().strip() or "Resultaat"
        unit = self.result_unit.text().strip()
        for card in self.overrides.values():
            card.set_sample_result(label, unit, self.response_scale.currentData())

    def to_config(self) -> dict:
        calculation = {"method": "formula", **self.formula.to_config(), "response_scale": self.response_scale.currentData()}
        references = dict(calculation.get("references", {}))
        overrides = {}
        control_results = {}
        for name, card in self.overrides.items():
            override, control_result = card.to_config()
            if override is not None:
                references.update(override.pop("references", {}))
                overrides[name] = override
            if control_result is not None:
                control_results[name] = control_result
        if references:
            calculation["references"] = references
        if overrides:
            calculation["control_overrides"] = overrides
        if control_results:
            calculation["control_results"] = control_results
        return {"calculation": calculation, "result": _result_config(self.result_label, self.result_unit)}
