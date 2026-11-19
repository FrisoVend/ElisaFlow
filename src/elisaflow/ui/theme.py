"""The desktop design system: one palette, stylesheet and layout rhythm."""

from string import Template

from PySide6.QtCore import QLibraryInfo, QSize, Qt, QTranslator
from PySide6.QtGui import QColor, QFont, QIcon, QPalette, QPixmap
from PySide6.QtWidgets import (
    QApplication,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QScrollArea,
    QSizePolicy,
    QTableWidget,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from elisaflow.branding import ACCENT, BORDER, BRAND, HEADER, ICON_PATH, INK, LOGO_PATH, MUTED, STRIPE

_STYLESHEET = Template("""
QWidget {
    background-color: white;
    color: $ink;
    font-family: "Segoe UI";
    font-size: 14px;
}
QMainWindow, QDialog, QWidget#appContent {
    background-color: $stripe;
}
QLabel, QCheckBox, QStackedWidget, QTabWidget {
    background-color: transparent;
}
QCheckBox {
    spacing: 8px;
}
QCheckBox::indicator {
    width: 16px;
    height: 16px;
    border: 1px solid $muted;
    border-radius: 3px;
    background-color: white;
}
QCheckBox::indicator:hover, QCheckBox:focus::indicator {
    border-color: $accent;
}
QCheckBox::indicator:checked {
    background-color: $brand;
    border-color: $brand;
    image: url("$checkmark");
}
QCheckBox::indicator:disabled {
    background-color: $stripe;
    border-color: $border;
}
QCheckBox::indicator:checked:disabled {
    background-color: $muted;
    border-color: $muted;
}
QWidget#surface {
    background-color: white;
    border: 1px solid $border;
    border-radius: 8px;
}
QWidget#pageHeader {
    background-color: transparent;
}
QLabel#pageTitle {
    color: $ink;
    font-size: 24px;
    font-weight: 600;
}
QLabel#pageSubtitle, QLabel#helperText {
    color: $muted;
    font-size: 13px;
}
QLabel#sectionTitle {
    font-size: 15px;
    font-weight: 600;
    padding-top: 8px;
    padding-bottom: 4px;
}
QLabel#statusTitle {
    font-size: 17px;
    font-weight: 600;
}
QLabel#statusTitle[state="success"] { color: $brand; }
QLabel#statusTitle[state="error"] { color: #B22222; }
QLabel#emptyState {
    color: $muted;
    padding: 32px;
}
QPushButton {
    background-color: white;
    color: $brand;
    border: 1px solid $border;
    border-radius: 6px;
    padding: 8px 14px;
    min-height: 20px;
    font-weight: 600;
}
QPushButton:hover, QPushButton:focus {
    background-color: $header;
    border-color: $accent;
}
QPushButton:pressed {
    border-color: $brand;
}
QPushButton[primary="true"] {
    background-color: $brand;
    color: white;
    border-color: $brand;
}
QPushButton[primary="true"]:hover, QPushButton[primary="true"]:focus {
    background-color: $ink;
    border-color: $accent;
}
QPushButton[utility="true"] {
    background-color: transparent;
    color: $muted;
    border-color: transparent;
}
QPushButton[utility="true"]:hover, QPushButton[utility="true"]:focus {
    background-color: transparent;
    color: $brand;
    border-color: transparent;
}
QPushButton[utility="true"]:pressed {
    color: $brand;
    border-color: transparent;
}
QPushButton:disabled, QPushButton[primary="true"]:disabled, QPushButton[utility="true"]:disabled {
    background-color: $header;
    color: $muted;
    border-color: $border;
}
QLineEdit, QComboBox, QSpinBox {
    background-color: white;
    border: 1px solid $border;
    border-radius: 5px;
    padding: 7px 9px;
    min-height: 20px;
    selection-background-color: $brand;
    selection-color: white;
}
QComboBox {
    padding-right: 24px;
}
QComboBox::drop-down {
    subcontrol-origin: padding;
    subcontrol-position: top right;
    width: 24px;
    border: none;
    background-color: transparent;
}
QComboBox::down-arrow {
    image: url("$chevron");
    width: 12px;
    height: 12px;
}
QSpinBox {
    padding-right: 30px;
}
QSpinBox QLineEdit {
    border: none;
    border-radius: 0;
    padding: 0;
    min-height: 0;
    background-color: transparent;
}
QSpinBox::up-button, QSpinBox::down-button {
    subcontrol-origin: border;
    width: 26px;
    border: none;
    border-left: 1px solid $border;
    background-color: $stripe;
}
QSpinBox::up-button {
    subcontrol-position: top right;
    border-top-right-radius: 4px;
    border-bottom: 1px solid $border;
}
QSpinBox::down-button {
    subcontrol-position: bottom right;
    border-bottom-right-radius: 4px;
}
QSpinBox::up-button:hover, QSpinBox::down-button:hover {
    background-color: $header;
}
QSpinBox::up-arrow {
    image: url("$plus");
    width: 12px;
    height: 12px;
}
QSpinBox::down-arrow {
    image: url("$minus");
    width: 12px;
    height: 12px;
}
QLineEdit:focus, QComboBox:focus, QSpinBox:focus {
    border-color: $accent;
}
QLineEdit[invalid="true"] {
    border-color: #B22222;
    background-color: #FFF5F5;
}
QLabel#validationError {
    color: #B22222;
}
QLineEdit:read-only {
    background-color: $stripe;
    color: $ink;
}
QLineEdit:disabled, QComboBox:disabled, QSpinBox:disabled {
    background-color: $stripe;
    color: $muted;
}
QAbstractItemView {
    background-color: white;
    alternate-background-color: $stripe;
    border: 1px solid $border;
    selection-background-color: $header;
    selection-color: $ink;
    border-radius: 6px;
}
QListWidget::item {
    padding: 12px;
    border-bottom: 1px solid $stripe;
}
QListWidget::item:selected {
    background-color: $header;
    color: $brand;
    border-left: 3px solid $accent;
    padding-left: 9px;
}
QListWidget::item:hover:!selected {
    background-color: $stripe;
}
QListWidget#editorNavigation {
    background-color: transparent;
    border: none;
}
QListWidget#editorNavigation::item {
    border: none;
    border-radius: 5px;
    margin-bottom: 4px;
    padding: 11px 12px;
}
QListWidget#editorNavigation::item:selected {
    background-color: $header;
    color: $brand;
    border: none;
}
QTableView {
    border: none;
    border-radius: 0;
    padding-bottom: 16px;
}
QTableView::item {
    padding: 0 12px;
    border: none;
    border-bottom: 1px solid $stripe;
}
QTableView QLineEdit, QTableView QComboBox {
    margin: 6px 8px;
}
QTableView QLineEdit:read-only {
    border-color: transparent;
    background-color: transparent;
}
QAbstractItemView::item:disabled {
    color: $muted;
}
QComboBox QAbstractItemView {
    background-color: white;
    color: $ink;
    padding: 4px;
}
QHeaderView {
    border: none;
    border-radius: 0;
    background-color: $stripe;
}
QHeaderView::section {
    background-color: $stripe;
    color: $muted;
    font-weight: 600;
    border: none;
    border-bottom: 1px solid $border;
    padding: 10px 12px;
}
QTableView#validationCriteria QHeaderView::section {
    padding-left: 20px;
}
QTableView#validationCriteria QHeaderView::section:first {
    padding-left: 12px;
}
QTableCornerButton::section {
    background-color: $header;
    border: none;
}
QTabWidget::pane {
    background-color: white;
    border: none;
    border-top: 1px solid $border;
    border-radius: 0;
    top: -1px;
}
QTabWidget#mainTabs::pane {
    background-color: $stripe;
}
QTabBar {
    background-color: transparent;
    qproperty-drawBase: 0;
}
QTabBar::tab {
    background-color: transparent;
    color: $muted;
    border-bottom: 3px solid transparent;
    padding: 11px 16px;
    font-weight: 400;
}
QTabBar::tab:hover {
    background-color: transparent;
    color: $brand;
}
QTabBar::tab:selected {
    background-color: transparent;
    color: $brand;
    border-bottom-color: $accent;
    font-weight: 600;
}
QScrollArea {
    background-color: white;
    border: none;
}
QGroupBox {
    background-color: white;
    border: none;
    border-top: 1px solid $border;
    border-radius: 0;
    margin-top: 28px;
    padding: 16px 0 0;
    font-weight: 600;
}
QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 0;
    padding: 0;
}
QProgressBar {
    background-color: $header;
    border: none;
    border-radius: 3px;
    max-height: 6px;
}
QProgressBar::chunk {
    background-color: $accent;
    border-radius: 3px;
}
QScrollBar:vertical {
    background-color: $stripe;
    width: 10px;
    margin: 0;
}
QScrollBar::handle:vertical {
    background-color: $border;
    min-height: 28px;
    border-radius: 4px;
    margin: 2px;
}
QScrollBar::handle:vertical:hover {
    background-color: $muted;
}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0;
}
QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical {
    background-color: $stripe;
}
QScrollBar:horizontal {
    background-color: $stripe;
    height: 10px;
    margin: 0;
}
QScrollBar::handle:horizontal {
    background-color: $border;
    min-width: 28px;
    border-radius: 4px;
    margin: 2px;
}
QScrollBar::handle:horizontal:hover { background-color: $muted; }
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal { width: 0; }
QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal { background-color: $stripe; }
QToolTip {
    background-color: $header;
    color: $ink;
    border: 1px solid $border;
    padding: 4px;
}
""").substitute(
    ink=INK,
    muted=MUTED,
    brand=BRAND,
    accent=ACCENT,
    border=BORDER,
    header=HEADER,
    stripe=STRIPE,
    chevron=LOGO_PATH.with_name("chevron_down.svg").as_posix(),
    checkmark=LOGO_PATH.with_name("checkmark.svg").as_posix(),
    plus=LOGO_PATH.with_name("plus.svg").as_posix(),
    minus=LOGO_PATH.with_name("minus.svg").as_posix(),
)


def apply_theme(window: QWidget) -> None:
    """Initialize the shared application theme once, including dialogs and popups."""
    application = QApplication.instance()
    window.setWindowIcon(QIcon(str(ICON_PATH)))
    if application.property("elisaflowTheme"):
        application.setPalette(application.property("elisaflowPalette"))
        return
    application.setWindowIcon(window.windowIcon())
    translator = QTranslator(application)
    translator.load("qtbase_nl", QLibraryInfo.path(QLibraryInfo.LibraryPath.TranslationsPath))
    application.installTranslator(translator)
    application.styleHints().setColorScheme(Qt.ColorScheme.Light)
    if application.style().objectName().lower() != "fusion":
        application.setStyle("Fusion")
    application.setFont(QFont("Segoe UI", 10))
    palette = QPalette()
    for role, colour in (
        (QPalette.ColorRole.Window, STRIPE),
        (QPalette.ColorRole.WindowText, INK),
        (QPalette.ColorRole.Base, "#FFFFFF"),
        (QPalette.ColorRole.AlternateBase, STRIPE),
        (QPalette.ColorRole.Text, INK),
        (QPalette.ColorRole.Button, "#FFFFFF"),
        (QPalette.ColorRole.ButtonText, INK),
        (QPalette.ColorRole.Highlight, BRAND),
        (QPalette.ColorRole.HighlightedText, "#FFFFFF"),
        (QPalette.ColorRole.PlaceholderText, MUTED),
        (QPalette.ColorRole.Accent, ACCENT),
        (QPalette.ColorRole.Light, "#FFFFFF"),
        (QPalette.ColorRole.Midlight, HEADER),
        (QPalette.ColorRole.Mid, BORDER),
        (QPalette.ColorRole.Dark, MUTED),
        (QPalette.ColorRole.Shadow, INK),
        (QPalette.ColorRole.BrightText, "#FFFFFF"),
        (QPalette.ColorRole.ToolTipBase, HEADER),
        (QPalette.ColorRole.ToolTipText, INK),
        (QPalette.ColorRole.Link, BRAND),
        (QPalette.ColorRole.LinkVisited, BRAND),
    ):
        palette.setColor(role, QColor(colour))
    for role in (
        QPalette.ColorRole.WindowText,
        QPalette.ColorRole.Text,
        QPalette.ColorRole.ButtonText,
    ):
        palette.setColor(QPalette.ColorGroup.Disabled, role, QColor(MUTED))
    application.setPalette(palette)
    application.setProperty("elisaflowPalette", palette)
    application.setStyleSheet(_STYLESHEET)
    application.setProperty("elisaflowTheme", True)


def form_layout(parent: QWidget | None = None, *, inset: bool = False) -> QFormLayout:
    """Create a top-aligned form with consistent spacing and narrow-window wrapping."""
    form = QFormLayout(parent)
    margin = 16 if inset else 0
    form.setContentsMargins(margin, margin, margin, margin)
    form.setHorizontalSpacing(24)
    form.setVerticalSpacing(12)
    form.setLabelAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
    form.setFormAlignment(Qt.AlignmentFlag.AlignTop)
    form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
    form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
    return form


class FormSection(QGroupBox):
    """A named set of related fields with a divider beneath its heading."""

    def __init__(self, title: str) -> None:
        super().__init__(title)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        self.form = form_layout(self)


def page_layout(parent: QWidget) -> QVBoxLayout:
    layout = QVBoxLayout(parent)
    layout.setContentsMargins(24, 20, 24, 20)
    layout.setSpacing(20)
    return layout


def surface() -> QWidget:
    widget = QWidget()
    widget.setObjectName("surface")
    widget.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
    return widget


def scroll_page(widget: QWidget) -> QScrollArea:
    """Give a form one scroll container; tables and tab containers scroll themselves."""
    scroll = QScrollArea()
    scroll.setWidgetResizable(True)
    scroll.setWidget(widget)
    return scroll


def editor_tabs() -> QTabWidget:
    """Use the same native tab bar and full-width divider at every level."""
    tabs = QTabWidget()
    tabs.tabBar().setExpanding(False)
    return tabs


def configure_table(table: QTableWidget) -> None:
    """Use the same header alignment, row rhythm and form controls in every table."""
    table.verticalHeader().hide()
    table.verticalHeader().setDefaultSectionSize(52)
    table.verticalHeader().setMinimumSectionSize(52)
    table.horizontalHeader().setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
    table.horizontalHeader().setMinimumSectionSize(120)
    table.setAlternatingRowColors(False)
    table.setShowGrid(False)
    table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
    table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
    table.setFocusPolicy(Qt.FocusPolicy.NoFocus)


def helper_text(text: str = "") -> QLabel:
    label = QLabel(text)
    label.setObjectName("helperText")
    label.setWordWrap(True)
    label.setTextFormat(Qt.TextFormat.PlainText)
    return label


class PageHeader(QWidget):
    def __init__(self, title: str, subtitle: str = "") -> None:
        super().__init__()
        self.setObjectName("pageHeader")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(20)
        text = QVBoxLayout()
        text.setSpacing(6)
        text.setAlignment(Qt.AlignmentFlag.AlignVCenter)
        title_label = QLabel(title)
        title_label.setObjectName("pageTitle")
        title_label.setTextFormat(Qt.TextFormat.PlainText)
        title_label.setWordWrap(True)
        text.addWidget(title_label)
        if subtitle:
            subtitle_label = QLabel(subtitle)
            subtitle_label.setObjectName("pageSubtitle")
            subtitle_label.setTextFormat(Qt.TextFormat.PlainText)
            subtitle_label.setWordWrap(True)
            text.addWidget(subtitle_label)
        layout.addLayout(text, 1)

        self.logo = QLabel()
        self.logo.setAccessibleName("Reinier Diagnostiek")
        self.logo.setFixedSize(48, 64)
        self.logo.setAlignment(Qt.AlignmentFlag.AlignCenter)
        pixmap = QPixmap(str(LOGO_PATH))
        if pixmap.isNull():
            self.logo.setText("Reinier\nDiagnostiek")
        else:
            ratio = self.devicePixelRatioF()
            pixmap = pixmap.scaled(
                QSize(round(48 * ratio), round(64 * ratio)),
                Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation,
            )
            pixmap.setDevicePixelRatio(ratio)
            self.logo.setPixmap(pixmap)
        layout.addWidget(self.logo)
