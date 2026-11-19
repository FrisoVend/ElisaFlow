"""Column definitions and grouped headers shared by report tables."""

from dataclasses import dataclass
from xml.sax.saxutils import escape

from reportlab.lib.colors import Color
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.platypus import Paragraph


@dataclass(frozen=True)
class TableColumn:
    key: str
    label: str
    group: str | None = None
    weight: float = 1
    border_after: bool = False
    shared: bool = False
    header_alignment: int = TA_LEFT
    minimum_width_label: str | None = None
    numeric: bool = False


@dataclass(frozen=True)
class TableHeaders:
    """Header rows, styling and the column ranges shared by each group."""

    rows: list[list[Paragraph]]
    commands: list[tuple]
    groups: tuple[tuple[int, int], ...]


def paragraph_width(cell: Paragraph) -> float:
    return stringWidth(cell.getPlainText(), cell.style.fontName, cell.style.fontSize)


def table_headers(
    columns: list[TableColumn],
    style: ParagraphStyle,
    *,
    rows: int,
    border: Color,
) -> TableHeaders:
    styles = {
        TA_LEFT: style,
        TA_CENTER: ParagraphStyle(f"{style.name}Metric", parent=style, alignment=TA_CENTER),
    }
    headings, subheadings, commands = [], [], []
    groups = []
    for index, column in enumerate(columns):
        heading = column.label
        if column.group is not None:
            heading = ""
            if index == 0 or column.group != columns[index - 1].group:
                heading = column.group
                end = index
                while end + 1 < len(columns) and columns[end + 1].group == column.group:
                    end += 1
                groups.append((index, end + 1))
                if end > index:
                    commands.extend(
                        [
                            ("SPAN", (index, 0), (end, 0)),
                            ("LINEBELOW", (index, 0), (end, 0), 0.25, border),
                        ]
                    )
        elif rows == 2:
            commands.extend(
                [
                    ("SPAN", (index, 0), (index, 1)),
                    ("VALIGN", (index, 0), (index, 1), "MIDDLE"),
                ]
            )
        alignment = TA_CENTER if column.group else column.header_alignment
        headings.append(Paragraph(escape(heading), styles[alignment]))
        subheadings.append(Paragraph(escape(column.label) if column.group else "", styles[TA_CENTER]))
        if column.border_after:
            commands.append(("LINEAFTER", (index, 0), (index, -1), 0.35, border))
    return TableHeaders([headings, *([subheadings] if rows == 2 else [])], commands, tuple(groups))
