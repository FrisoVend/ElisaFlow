"""Assay management page: view, change, or create an assay."""

from contextlib import nullcontext
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QWidget,
)

from elisaflow.assays.access_control import MIN_PASSWORD_LENGTH, AccessControl
from elisaflow.assays.management import AssayManagement
from elisaflow.formatting import format_protocol_name
from elisaflow.ui.assay_editor import AssayEditor
from elisaflow.ui.theme import PageHeader, helper_text, page_layout


class AssaySettingsWidget(QWidget):
    def __init__(self, assays: AssayManagement, access: AccessControl, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._assays = assays
        self._access = access
        self.setObjectName("appContent")
        layout = page_layout(self)
        layout.addWidget(PageHeader("Assay beheer", "Bekijk protocollen en beheer de instellingen per assay."))

        toolbar = QHBoxLayout()
        toolbar.setSpacing(12)
        self.search = QLineEdit()
        self.search.setPlaceholderText("Zoek een assay…")
        self.search.setAccessibleName("Assays zoeken")
        self.search.setClearButtonEnabled(True)
        toolbar.addWidget(self.search, 1)
        self.new_button = QPushButton("Nieuwe assay")
        self.new_button.clicked.connect(self._new_assay)
        toolbar.addWidget(self.new_button)
        layout.addLayout(toolbar)
        self.count = helper_text()
        layout.addWidget(self.count)

        self.assay_list = QListWidget()
        self._populate_assay_list(assays.list_assays())
        self.assay_list.itemSelectionChanged.connect(self._selection_changed)
        self.assay_list.itemDoubleClicked.connect(lambda _item: self._open_assay())
        layout.addWidget(self.assay_list)
        self.empty_state = QLabel()
        self.empty_state.setObjectName("emptyState")
        self.empty_state.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.empty_state.setWordWrap(True)
        layout.addWidget(self.empty_state, 1)
        self.search.textChanged.connect(self._filter_assays)

        buttons = QHBoxLayout()
        buttons.setSpacing(10)
        self.password_button = QPushButton("Wachtwoord beheren")
        self.password_button.clicked.connect(self._change_password)
        buttons.addWidget(self.password_button)
        buttons.addStretch()
        self.view_button = QPushButton("Bekijk assay")
        self.view_button.setProperty("primary", True)
        self.view_button.setEnabled(False)
        self.view_button.clicked.connect(self._open_assay)

        self.edit_button = QPushButton("Verander assay")
        self.edit_button.setEnabled(False)
        self.edit_button.clicked.connect(lambda: self._open_assay(read_only=False))

        buttons.addWidget(self.view_button)
        buttons.addWidget(self.edit_button)
        layout.addLayout(buttons)
        self._filter_assays()
        self._selection_changed()

    def _request_edit_password(self) -> str | None:
        if not self._access.has_password():
            QMessageBox.information(self, "Wachtwoord ontbreekt", "Stel eerst een wachtwoord in via Wachtwoord beheren.")
            return None
        password, accepted = QInputDialog.getText(self, "Assay bewerken", "Wachtwoord:", QLineEdit.EchoMode.Password)
        return password if accepted else None

    def _change_password(self) -> None:
        try:
            current_password = None
            title = "Wachtwoord instellen"
            if self._access.has_password():
                title = "Wachtwoord wijzigen"
                current_password, accepted = QInputDialog.getText(self, title, "Huidig wachtwoord:", QLineEdit.EchoMode.Password)
                if not accepted:
                    return
                if not self._access.verify_password(current_password):
                    raise PermissionError("Het huidige wachtwoord is onjuist.")
            password, accepted = QInputDialog.getText(
                self, title, f"Nieuw wachtwoord (minimaal {MIN_PASSWORD_LENGTH} tekens):", QLineEdit.EchoMode.Password,
            )
            if not accepted:
                return
            confirmation, accepted = QInputDialog.getText(self, title, "Herhaal het nieuwe wachtwoord:", QLineEdit.EchoMode.Password)
            if not accepted:
                return
            if password != confirmation:
                raise ValueError("De nieuwe wachtwoorden komen niet overeen.")
            self._access.set_new_password(password, current_password=current_password)
            QMessageBox.information(self, title, "Het wachtwoord is opgeslagen.")
        except (OSError, ValueError) as error:
            QMessageBox.warning(self, "Wachtwoord kon niet opgeslagen worden", str(error))

    def _populate_assay_list(self, names: tuple[str, ...]) -> None:
        self.assay_list.clear()
        for name in names:
            item = QListWidgetItem(format_protocol_name(name))
            item.setData(Qt.ItemDataRole.UserRole, name)
            self.assay_list.addItem(item)

    def _filter_assays(self) -> None:
        query = self.search.text().strip().casefold()
        visible = 0
        for index in range(self.assay_list.count()):
            item = self.assay_list.item(index)
            matches = query in item.text().casefold()
            item.setHidden(not matches)
            visible += matches
        current = self.assay_list.currentItem()
        if current is not None and current.isHidden():
            self.assay_list.setCurrentRow(-1)
        total = self.assay_list.count()
        self.count.setText(f"{visible} van {total} assays" if query else "1 assay" if total == 1 else f"{total} assays")
        self.assay_list.setVisible(visible > 0)
        self.empty_state.setVisible(visible == 0)
        self.empty_state.setText(
            "Geen assays gevonden.\nPas de zoekopdracht aan." if total else "Nog geen assays beschikbaar.\nVoeg een assay toe via Nieuwe assay."
        )

    def _selection_changed(self) -> None:
        selected = self.assay_list.currentItem() is not None
        self.view_button.setEnabled(selected)
        self.edit_button.setEnabled(selected)

    def _open_assay(self, *, read_only: bool = True) -> None:
        selected = self.assay_list.currentItem()
        if selected is None:
            return
        try:
            authorization = nullcontext()
            if not read_only:
                password = self._request_edit_password()
                if password is None:
                    return
                authorization = self._access.edit_assays(password)
            with authorization:
                if not read_only:
                    self._access.require_assay_editor()
                settings = self._assays.load_assay(selected.data(Qt.ItemDataRole.UserRole))
                editor = AssayEditor(settings, self._assays, self._access, parent=self, read_only=read_only)
                editor.exec()
        except (OSError, ValueError) as error:
            QMessageBox.warning(self, "Assay kon niet geopend worden", str(error))
            return

    def _new_assay(self) -> None:
        try:
            password = self._request_edit_password()
            if password is None:
                return
            with self._access.edit_assays(password):
                self._access.require_assay_editor()
                filename, _ = QFileDialog.getOpenFileName(
                    self,
                    "Selecteer een SkanIt-protocol",
                    "",
                    "SkanIt-protocollen (*.skax)",
                )
                if not filename:
                    return
                path = Path(filename)
                editor = AssayEditor(None, self._assays, self._access, parent=self, protocol_path=path)
                if editor.exec() == QDialog.DialogCode.Accepted:
                    names = self._assays.list_assays()
                    self._populate_assay_list(names)
                    self.search.clear()
                    self._filter_assays()
                    self.assay_list.setCurrentRow(names.index(path.resolve().name))
        except (OSError, ValueError) as error:
            QMessageBox.warning(self, "Assay kon niet worden aangemaakt", str(error))
