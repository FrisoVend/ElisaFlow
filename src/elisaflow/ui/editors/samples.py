"""Edit QC, calculation and report settings for a SKAX protocol."""

from collections.abc import Mapping

from PySide6.QtWidgets import (
    QLineEdit,
    QVBoxLayout,
    QWidget,
)

from elisaflow.ui.editors.fields import _ReplicateFilterEditor
from elisaflow.ui.theme import editor_tabs, form_layout, helper_text


class PatientEditor(QWidget):
    def __init__(self, config: dict, editable: bool, *, has_blank: bool = False) -> None:
        super().__init__()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        self.replicate_filter = _ReplicateFilterEditor(config.get("replicate_filter"), has_blank=has_blank)
        layout.addWidget(self.replicate_filter)
        layout.addStretch()
        self.setEnabled(editable)

    def to_config(self) -> dict:
        rule = self.replicate_filter.to_config()
        return {"replicate_filter": rule} if rule is not None else {}


class SampleEditor(QWidget):
    def __init__(
        self,
        config: dict,
        editable: bool,
        *,
        concentration: float | None = None,
        control: bool = False,
        has_blank: bool = False,
    ) -> None:
        super().__init__()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(16)
        form = form_layout()
        layout.addLayout(form)
        self.concentration = None
        if concentration is not None:
            self.concentration = QLineEdit(str(concentration))
            self.concentration.setReadOnly(True)
            form.addRow("Concentratie:", self.concentration)
        self.lot_number = QLineEdit(config.get("lot_number", ""))
        form.addRow("Lotnummer:", self.lot_number)
        self.replicate_filter = None
        if control or concentration is not None:
            self.replicate_filter = _ReplicateFilterEditor(config.get("replicate_filter"), has_blank=has_blank)
            layout.addWidget(self.replicate_filter)
        layout.addStretch()
        self.setEnabled(editable)

    def to_config(self) -> dict:
        config = {}
        lot_number = self.lot_number.text().strip()
        if lot_number:
            config["lot_number"] = lot_number
        if self.replicate_filter is not None:
            rule = self.replicate_filter.to_config()
            if rule is not None:
                config["replicate_filter"] = rule
        return config


class SampleCollection(QWidget):
    def __init__(
        self,
        samples: dict,
        editable: bool,
        *,
        names: Mapping[str | float, str],
        standard: bool = False,
        has_blank: bool = False,
    ) -> None:
        super().__init__()
        self._configured = set(samples)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.tabs = editor_tabs()
        self.tabs.setVisible(bool(names))
        layout.addWidget(self.tabs)
        self.editors = {}
        if not names:
            label = "standaarden" if standard else "controles"
            layout.addWidget(helper_text(f"Geen {label} gedefinieerd in het SKAX-protocol."))
        for key, name in names.items():
            editor = SampleEditor(
                samples.get(str(key), {}),
                editable,
                concentration=key if standard else None,
                control=not standard,
                has_blank=has_blank,
            )
            self.editors[key] = editor
            self.tabs.addTab(editor, name)

    def to_config(self) -> dict:
        samples = {}
        for key, editor in self.editors.items():
            config = editor.to_config()
            # Merely opening a protocol sample does not add empty QC settings.
            if config or str(key) in self._configured:
                samples[str(key)] = config
        return samples
