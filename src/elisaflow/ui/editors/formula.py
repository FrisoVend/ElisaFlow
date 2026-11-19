"""One expression editor for results and optional preprocessing."""

from copy import deepcopy

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from elisaflow.formatting import format_formula
from elisaflow.formulas.expressions import Formula
from elisaflow.formulas.names import FormulaNames
from elisaflow.ui.theme import PageHeader, apply_theme, form_layout, helper_text, page_layout, surface


def _formula_config(formula: Formula, control_names: tuple[str, ...]) -> dict:
    for control in formula.control_names:
        if control not in control_names:
            raise ValueError(
                f"Controle {control!r} ontbreekt in het SKAX-protocol. "
                "Verwijder de referentie uit de formule en voeg een beschikbare controle in."
            )
    used = formula.used_references
    return {"expression": formula.expression, **({"references": used} if used else {})}


class FormulaInput(QWidget):
    """Insert the sample value or SKAX control means into a formula."""

    def __init__(self, expression: str, control_names: tuple[str, ...], references: dict[str, str], *, sample_name: str):
        super().__init__()
        self._control_names = control_names
        self._names = FormulaNames(sample_name, references)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.expression = QLineEdit(self._names.format(expression))
        form = form_layout()
        form.addRow("Formule:", self.expression)
        layout.addLayout(form)
        operands = QHBoxLayout()
        self.operand = QComboBox()
        self.operand.addItem(f"{self._names.sample_name} (meetwaarde)")
        for name in control_names:
            self.operand.addItem(f"{name} (controlegemiddelde)", name)
        self.operand.setPlaceholderText("Selecteer een waarde")
        self.operand.setCurrentIndex(-1)
        self.operand.setAccessibleName("Waarde voor de formule")
        self.insert_button = QPushButton("Invoegen")
        self.insert_button.setEnabled(False)
        self.operand.currentIndexChanged.connect(lambda index: self.insert_button.setEnabled(index >= 0))
        self.insert_button.clicked.connect(self._insert_operand)
        operands.addWidget(self.operand, 1)
        operands.addWidget(self.insert_button)
        layout.addLayout(operands)

    def _insert_operand(self) -> None:
        if self.operand.currentIndex() < 0:
            return
        control = self.operand.currentData()
        self.expression.insert(self._names.sample_symbol if control is None else self._names.control_symbol(control))
        self.expression.setFocus()

    def to_config(self) -> dict:
        return _formula_config(self._names.parse(self.expression.text()), self._control_names)


class FormulaDialog(QDialog):
    """Validate a private equation draft before committing it to the assay editor."""

    def __init__(
        self, title: str, config: dict | None, control_names: tuple[str, ...], references: dict[str, str], parent=None, *, sample_name: str,
    ):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(680, 380)
        apply_theme(self)
        self._config = None
        layout = page_layout(self)
        layout.addWidget(PageHeader(title))
        card = surface()
        content = QVBoxLayout(card)
        content.setContentsMargins(20, 20, 20, 20)
        content.setSpacing(16)
        self.input = FormulaInput((config or {}).get("expression", ""), control_names, dict(references), sample_name=sample_name)
        content.addWidget(self.input)
        content.addWidget(helper_text(
            "Kies de meetwaarde of een controle en klik op Invoegen. Gebruik getallen, +, -, *, / en haakjes. "
            "Controles gebruiken oorspronkelijke controlegemiddelden."
        ))
        layout.addWidget(card)
        layout.addStretch()
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel)
        save = buttons.button(QDialogButtonBox.StandardButton.Save)
        save.setText("Opslaan")
        save.setProperty("primary", True)
        buttons.button(QDialogButtonBox.StandardButton.Cancel).setText("Annuleer")
        buttons.accepted.connect(self._save)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _save(self) -> None:
        try:
            config = self.input.to_config()
        except ValueError as error:
            QMessageBox.warning(self, "Formule is ongeldig", str(error))
            self.input.expression.setFocus()
            return
        self._config = config
        self.accept()

    def to_config(self) -> dict:
        if self._config is None:
            raise ValueError("Sla eerst een geldige formule op.")
        return deepcopy(self._config)


class FormulaEditor(QWidget):
    """A compact saved equation with editing delegated to its own dialog."""

    def __init__(
        self, config: dict | None, control_names: tuple[str, ...], *, sample_name: str, optional: bool = False,
        references: dict[str, str] | None = None, title: str = "Formule bewerken", empty_text: str = "Geen voorbewerking",
        clear_text: str = "Verwijderen",
    ):
        super().__init__()
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._control_names = control_names
        self._sample_name = sample_name
        self._title = title
        self._empty_text = empty_text
        if config is None and not optional:
            config = {"expression": "sample"}
        self._config = deepcopy(config)
        self.references = references if references is not None else dict((config or {}).get("references", {}))
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(10)
        self.summary = QLineEdit()
        self.summary.setReadOnly(True)
        self.summary.setAccessibleName("Formule")
        layout.addWidget(self.summary, 1)
        self.edit_button = QPushButton("Bewerken…")
        self.edit_button.setAccessibleName(title)
        self.edit_button.clicked.connect(self._edit)
        layout.addWidget(self.edit_button)
        self.clear_button = QPushButton(clear_text) if optional else None
        if self.clear_button is not None:
            self.clear_button.clicked.connect(self.clear)
            layout.addWidget(self.clear_button)
        self._refresh()

    def _refresh(self) -> None:
        self.summary.setText(
            format_formula(Formula(self._config["expression"], self.references), self._sample_name) if self._config is not None else self._empty_text
        )
        self.summary.setCursorPosition(0)
        self.summary.setToolTip(self.summary.text())
        if self.clear_button is not None:
            self.clear_button.setEnabled(self._config is not None)

    def _edit(self) -> None:
        dialog = FormulaDialog(self._title, self._config, self._control_names, self.references, self, sample_name=self._sample_name)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self._config = dialog.to_config()
            self.references.update(self._config.get("references", {}))
            self._refresh()

    def clear(self) -> None:
        if self.clear_button is not None:
            self._config = None
            self._refresh()

    def to_config(self) -> dict | None:
        if self._config is None:
            return None
        return _formula_config(Formula(self._config["expression"], self.references), self._control_names)
