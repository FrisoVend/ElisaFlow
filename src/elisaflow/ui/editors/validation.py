"""Edit final-result criteria and measurement-specific extinction criteria."""

from math import isfinite

from PySide6.QtCore import QSignalBlocker, Qt, Signal
from PySide6.QtWidgets import (
    QHeaderView,
    QLabel,
    QLineEdit,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from elisaflow.configuration.assay_settings import MeasurementSettings
from elisaflow.models import ResponseScale
from elisaflow.ui.editors.fields import _response_scale_selector
from elisaflow.ui.theme import configure_table, editor_tabs, helper_text


class ValidationScopeEditor(QTableWidget):
    errors_changed = Signal()

    def __init__(self, config: dict, protocol, editable: bool, *, metrics=("extinction", "result")):
        super().__init__(0, 1 + 2 * len(metrics) + ("extinction" in metrics))
        self.setObjectName("validationCriteria")
        self.metrics = metrics
        self.visible_metrics = metrics
        self.fields = {}
        self.response_scales = {}
        self._editable = editable
        self._columns = {}
        self._scale_column = None
        self._targets = {}
        self._metric_labels = {"extinction": "Extinctie", "result": "Resultaat"}
        self.error_message = ""
        configure_table(self)
        targets = [("blank", None, "Blanco")]
        targets.extend(("controls", name, name) for name in protocol.controls)
        targets.extend(("standards", str(value), f"{value:g} ({', '.join(names)})") for value, names in protocol.standards.items())
        self._set_headers()
        self.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        for column in range(1, self.columnCount()):
            self.horizontalHeader().setSectionResizeMode(column, QHeaderView.ResizeMode.ResizeToContents)
        for row, (collection, key, label) in enumerate(targets):
            rules = config.get("blank", {}) if collection == "blank" else config.get(collection, {}).get(key, {})
            self.insertRow(row)
            item = QTableWidgetItem(label)
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.setItem(row, 0, item)
            for metric in metrics:
                if collection == "blank" and metric == "result":
                    continue
                bounds = []
                for offset, bound in enumerate(("min", "max")):
                    value = rules.get(metric, {}).get(bound)
                    field = QLineEdit("" if value is None else str(value))
                    field.setPlaceholderText("Geen")
                    field.setAlignment(Qt.AlignmentFlag.AlignLeft)
                    field.setReadOnly(not editable)
                    field.setMinimumWidth(112)
                    field.setAccessibleName(f"{label}: {self._metric_labels[metric]}, {bound}")
                    self.setCellWidget(row, self._columns[metric][offset], field)
                    bounds.append(field)
                    field.textChanged.connect(self._refresh_errors)
                identity = (collection, key, metric)
                self.fields[identity] = tuple(bounds)
                self._targets[identity] = (row, label)
                if metric == "extinction":
                    default_scale = ResponseScale.RAW if collection == "blank" else ResponseScale.default_for_blank(protocol.has_blank)
                    selector = _response_scale_selector(rules.get(metric, {}).get("response_scale", default_scale), raw_only=collection == "blank")
                    selector.setEnabled(editable and collection != "blank")
                    selector.setAccessibleName(f"{label}: extinctieschaal")
                    selector.setToolTip(
                        f"{selector.toolTip()} De keuze geldt voor de extinctiegrenzen en hun CV, vóór kalibratievoorbewerking."
                    )
                    self.setCellWidget(row, self._scale_column, selector)
                    self.response_scales[collection, key] = selector
        self._refresh_errors()

    def _set_headers(self) -> None:
        headers = ["Monster"]
        for metric in self.metrics:
            self._columns[metric] = (len(headers), len(headers) + 1)
            headers.extend(f"{self._metric_labels[metric]}\n{bound}" for bound in ("Minimum", "Maximum"))
            if metric == "extinction":
                self._scale_column = len(headers)
                headers.append("Extinctie\nSchaal")
        self.setHorizontalHeaderLabels(headers)

    def set_extinction_label(self, label: str) -> None:
        if self._metric_labels["extinction"] == label:
            return
        self._metric_labels["extinction"] = label
        self._set_headers()
        self._refresh_errors()

    def set_metrics(self, metrics: tuple[str, ...]) -> None:
        if self.visible_metrics == metrics:
            return
        self.visible_metrics = metrics
        for metric, columns in self._columns.items():
            for column in columns:
                self.setColumnHidden(column, metric not in metrics)
        if self._scale_column is not None:
            self.setColumnHidden(self._scale_column, "extinction" not in metrics)
        for (collection, _), selector in self.response_scales.items():
            selector.setEnabled(self._editable and collection != "blank" and "extinction" in metrics)
        for key, fields in self.fields.items():
            for field in fields:
                field.setEnabled(key[2] in metrics)
        for row in range(self.rowCount()):
            self.setRowHidden(row, not any(self._targets[key][0] == row and key[2] in metrics for key in self.fields))
        self._refresh_errors()

    def copy_extinction_from(self, source: ValidationScopeEditor) -> None:
        for key, destination in self.fields.items():
            if key[2] == "extinction" and key in source.fields:
                for field, original in zip(destination, source.fields[key], strict=True):
                    field.setText(original.text())
                selector = self.response_scales[key[:2]]
                selector.setCurrentIndex(selector.findData(source.response_scales[key[:2]].currentData()))

    def _limits(self, fields: tuple[QLineEdit, QLineEdit]) -> dict[str, float]:
        limits = {}
        for name, field in zip(("min", "max"), fields, strict=True):
            text = field.text().strip()
            if not text:
                continue
            try:
                value = float(text.replace(",", "."))
            except ValueError as error:
                raise ValueError("Voer een getal in.") from error
            if not isfinite(value):
                raise ValueError("Voer een eindig getal in.")
            limits[name] = value
        if "min" in limits and "max" in limits and limits["min"] > limits["max"]:
            raise ValueError("Het minimum mag niet groter zijn dan het maximum.")
        return limits

    def _refresh_errors(self) -> None:
        messages = []
        for key, fields in self.fields.items():
            error = ""
            if key[2] in self.visible_metrics:
                try:
                    self._limits(fields)
                except ValueError as caught:
                    error = str(caught)
                    messages.append(f"{self._targets[key][1]} – {self._metric_labels[key[2]]}: {error}")
            for bound, field in zip(("ondergrens", "bovengrens"), fields, strict=True):
                invalid = bool(error)
                if field.property("invalid") != invalid:
                    field.setProperty("invalid", invalid)
                    field.style().unpolish(field)
                    field.style().polish(field)
                    field.update()
                field.setToolTip(error or f"Leeg = geen {bound}. De grens zelf is toegestaan.")
        self.error_message = messages[0] if messages else ""
        self.errors_changed.emit()

    def to_config(self) -> dict:
        self._refresh_errors()
        if self.error_message:
            for key, fields in self.fields.items():
                if fields[0].property("invalid"):
                    self.scrollToItem(self.item(self._targets[key][0], 0))
                    fields[0].setFocus()
                    break
            raise ValueError(self.error_message)
        config = {}
        for (collection, key, metric), fields in self.fields.items():
            if metric not in self.visible_metrics:
                continue
            limits = self._limits(fields)
            if limits:
                target = config.setdefault("blank", {}) if collection == "blank" else config.setdefault(collection, {}).setdefault(key, {})
                target[metric] = (
                    {**limits, "response_scale": self.response_scales[collection, key].currentData()} if metric == "extinction" else limits
                )
        return config


class ValidationEditor(QWidget):
    def __init__(self, config: dict, protocol, editable: bool):
        super().__init__()
        self._protocol = protocol
        self._editable = editable
        self._switching = False
        self._single_wavelength = None
        self._preferred_wavelength = None
        self._wavelengths = ()
        self._settings = MeasurementSettings()
        self.editors = {}
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(12)
        self.hint = helper_text()
        layout.addWidget(self.hint)
        self.tabs = editor_tabs()
        layout.addWidget(self.tabs, 1)
        self.error = QLabel()
        self.error.setWordWrap(True)
        self.error.setObjectName("validationError")
        layout.addWidget(self.error)
        self.add_scope(None, config.get("result", {}))
        for wavelength, rules in config.get("wavelengths", {}).items():
            self.add_scope(int(wavelength), rules)
        self.tabs.addTab(self.editors[None], "Resultaat")
        self._visible_scopes = (None,)
        self.tabs.tabBar().hide()
        self.tabs.currentChanged.connect(self._show_context)
        self._show_context()

    def add_scope(self, wavelength: int | None, config: dict) -> None:
        if wavelength not in self.editors:
            metrics = ("extinction", "result") if wavelength is None else ("extinction",)
            editor = ValidationScopeEditor(config, self._protocol, self._editable, metrics=metrics)
            editor.errors_changed.connect(self._show_context)
            self.editors[wavelength] = editor

    def _show_context(self) -> None:
        editor = self.tabs.currentWidget()
        if not self._switching:
            text = "Alle QC-criteria voor deze meting."
        elif editor is self.editors[None]:
            text = "Criteria voor het gerapporteerde resultaat."
        else:
            text = "Extinctiecriteria voor deze meting."
        self.hint.setText(text)
        message = editor.error_message if editor is not None else ""
        self.error.setText(message)
        self.error.setVisible(bool(message))

    def _update_labels(self) -> None:
        for wavelength, editor in self.editors.items():
            measurement = wavelength if wavelength is not None else self._single_wavelength
            label = "Extinctie" if measurement is None else f"Extinctie ({self._settings.wavelength_label(measurement)})"
            editor.set_extinction_label(label)

    def set_measurement(self, settings: MeasurementSettings, wavelengths: tuple[int, ...]) -> None:
        result = self.editors[None]
        for wavelength in wavelengths:
            self.add_scope(wavelength, {})
        if settings.switching and not self._switching:
            source = self._single_wavelength or settings.preferred_wavelength_nm
            has_extinction = any(field.text().strip() for key, fields in result.fields.items() if key[2] == "extinction" for field in fields)
            if source is not None and (self._single_wavelength is not None or has_extinction):
                self.add_scope(source, {})
                self.editors[source].copy_extinction_from(result)
        single = wavelengths[0] if len(wavelengths) == 1 else self._preferred_wavelength
        if not settings.switching:
            if self._switching and single in self.editors:
                result.copy_extinction_from(self.editors[single])
            elif self._single_wavelength is not None and single != self._single_wavelength:
                self.editors[self._single_wavelength].copy_extinction_from(result)
                if single in self.editors:
                    result.copy_extinction_from(self.editors[single])
            self._single_wavelength = single
        self._preferred_wavelength = settings.preferred_wavelength_nm or self._preferred_wavelength
        self._switching = settings.switching
        self._settings = settings
        self._wavelengths = wavelengths
        result.set_metrics(("result",) if settings.switching else ("extinction", "result"))
        scopes = (None, *wavelengths) if settings.switching else (None,)
        if scopes != self._visible_scopes:
            current = self.tabs.currentWidget()
            with QSignalBlocker(self.tabs):
                self.tabs.clear()
                for wavelength in scopes:
                    self.tabs.addTab(self.editors[wavelength], "Resultaat" if wavelength is None else f"Meting {settings.wavelength_label(wavelength)}")
                if self.tabs.indexOf(current) >= 0:
                    self.tabs.setCurrentWidget(current)
            self._visible_scopes = scopes
        elif settings.switching:
            for index, wavelength in enumerate(wavelengths, 1):
                self.tabs.setTabText(index, f"Meting {settings.wavelength_label(wavelength)}")
        self.tabs.tabBar().setVisible(settings.switching)
        self._update_labels()
        self._show_context()

    def to_config(self) -> dict:
        config = {}
        active = (None, *self._wavelengths) if self._switching else (None,)
        for wavelength in active:
            editor = self.editors[wavelength]
            try:
                criteria = editor.to_config()
            except ValueError:
                self.tabs.setCurrentWidget(editor)
                self._show_context()
                raise
            if criteria:
                if wavelength is None:
                    config["result"] = criteria
                else:
                    config.setdefault("wavelengths", {})[str(wavelength)] = criteria
        return config
