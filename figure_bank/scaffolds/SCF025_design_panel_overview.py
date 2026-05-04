"""SCF025 — Design Panel Overview Table.

A structured table figure for design panel comparisons. Five fixed zones:

  1. Header     — title + subtitle + tagline (h=58)
  2. Pill strip — small caller-provided label + N evenly-spaced colored pills (h=62)
  3. Col header — dark bar with 8 caller-provided column labels (h=20)
  4. Data rows  — N rows at auto-computed height (from text-line count)
  5. Footer     — bar-axis tick labels, legend swatches, strategy note, callout (h=71)

ViewBox width is 960. Height computed at render time:
  viewbox_h = 144 + sum(row_heights) + 71

Row height rule (from total text lines):
  total_lines = 1 (name) + len(description) + (1 if annotation else 0)
  4+ lines → 60 px | 3 lines → 50 px | 2 lines → 40 px

The scaffold owns all geometry. Caller provides all text labels.
Bar axis ticks are locked at 0/20/40/60/80/90/100% (3.27 px/%).

Coverage: FIG167 (design panel overview). Any 8-column structured design table.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from ._runtime import SlotOverflowError, SlotTypeError, _escape, svg_close, svg_open


BARCODE = "SCF025"
DESCRIPTION = (
    "Design panel overview table — variable-height, 8 columns, pill strip, range bars"
)


# ---------------------------------------------------------------------------
# Locked geometry constants
# ---------------------------------------------------------------------------

_VB_W = 960

_HEADER_H = 58          # header zone
_PILL_TOP = 62          # pill strip starts (4 px gap after header)
_PILL_STRIP_H = 62      # pill strip height (ends at y=124)
_COL_HDR_Y = 124        # column header bar top
_COL_HDR_H = 20         # column header bar height
_ROWS_TOP = 144         # data rows start here  (= 58 + 4 + 62 + 20)
_FOOTER_H = 71          # footer zone height

# Column mid-x positions for header bar labels (locked, matches FIG167)
_COL_MIDS = (50, 185, 283, 335, 449, 523, 567, 762)

# Outcome bar geometry
_BAR_X = 608
_BAR_W = 327
_BAR_SCALE = 3.27        # px per %
_BAR_TICK_XS = (608, 673, 739, 804, 870, 902, 935)   # 0/20/40/60/80/90/100%

# Category badge geometry (col 3)
_BADGE_X0 = 301          # first badge in a row
_BADGE_DX = 36           # x step between badges in same row
_BADGE_W = 32
_BADGE_H = 13
_BADGE_ROW_DY = 19       # y step between badge rows

# Modifier badge geometry (col 5)
_MOD_X = 511
_MOD_W = 24
_MOD_H = 14
_MOD_CX = _MOD_X + _MOD_W // 2  # = 523

# Tier / level column center-x (col 6)
_LEVEL_CX = 567


# ---------------------------------------------------------------------------
# Public dataclasses
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Pill:
    """One key-finding pill in the strip above the table.

    Args:
        number:  Max 2 chars, e.g. "①"
        title:   Max 40 chars — bold first line inside the pill
        lines:   Exactly 2 body lines, max 55 chars each
        fill:    Background hex color, e.g. "#E8F5E9"
        stroke:  Border hex color, e.g. "#43A047"
    """
    number: str
    title: str
    lines: tuple[str, str]
    fill: str
    stroke: str


@dataclass(frozen=True)
class Badge:
    """A colored pill or plain-text label for category / modifier columns.

    Args:
        text:       Displayed text, max 8 chars
        color:      Hex → colored pill background; None → plain text
        text_color: Pill text color (default white)
    """
    text: str
    color: Optional[str] = None
    text_color: str = "white"


@dataclass(frozen=True)
class LegendEntry:
    """One swatch + label in the footer legend."""
    color: str    # swatch hex
    label: str    # max 40 chars


@dataclass(frozen=True)
class Row:
    """One data row in the table.

    Row height is computed from text content:
      total_lines = 1 + len(description) + (1 if annotation else 0)
      4+ → 60 px  |  3 → 50 px  |  2 → 40 px

    Args:
        letter:          Single character shown large in the left color block
        letter_color:    Left block background hex
        row_color:       Row background hex (alternating or category-coded)
        panel_sublabel:  Small text under the letter (max 14 chars)
        panel_type:      Tiny label below sublabel, h≥60 only (max 12 chars)
        name:            Bold panel name, max 32 chars
        description:     2–3 body lines at 8.5 pt gray, max 45 chars each
        annotation:      1 accent-colored note at 8 pt ("" → omitted), max 45 chars
        annotation_color: Hex for annotation text (default dark gray)
        n:               Integer count shown in column 2
        category_badges: Column 3 badges, 1–4 entries (2 per visual row)
        variant_lines:   Column 4 text lines, 1–3, max 40 chars each
        modifier:        Column 5 badge or None → "—"
        level:           Column 6 main text, max 10 chars
        level_sublabel:  Column 6 small gray sublabel ("" → omitted), max 14 chars
        outcome_range:   (min_pct, max_pct) for the filled bar, or None → no fill
        outcome_color:   Hex for filled bar range
        outcome_uncertain: True → show "?" centered on the filled bar
    """
    letter: str
    letter_color: str
    row_color: str
    panel_sublabel: str
    panel_type: str
    name: str
    description: tuple[str, ...]
    annotation: str
    annotation_color: str = "#546E7A"
    n: int = 0
    category_badges: tuple[Badge, ...] = ()
    variant_lines: tuple[str, ...] = ()
    modifier: Optional[Badge] = None
    level: str = ""
    level_sublabel: str = ""
    outcome_range: Optional[tuple[float, float]] = None
    outcome_color: str = "#1565C0"
    outcome_uncertain: bool = False


# ---------------------------------------------------------------------------
# Input validation helpers
# ---------------------------------------------------------------------------

def _check(name: str, text: str, max_chars: int) -> None:
    if len(text) > max_chars:
        raise SlotOverflowError(name, len(text), max_chars, unit=" chars")


def _row_height(row: Row) -> int:
    total_lines = 1 + len(row.description) + (1 if row.annotation else 0)
    if total_lines >= 4:
        return 60
    elif total_lines == 3:
        return 50
    else:
        return 40


# ---------------------------------------------------------------------------
# Zone 1: Header
# ---------------------------------------------------------------------------

def _render_header(title: str, subtitle: str, tagline: str) -> str:
    _check("title", title, 70)
    _check("subtitle", subtitle, 110)
    _check("tagline", tagline, 110)
    parts = [
        f'<rect width="{_VB_W}" height="58" fill="#F8F9FA"/>',
        f'<text x="480" y="28" text-anchor="middle" font-size="15.5" '
        f'font-weight="bold" fill="#1A237E" data-fig-role="title">{_escape(title)}</text>',
        f'<text x="480" y="45" text-anchor="middle" font-size="10" '
        f'fill="#546E7A" data-fig-role="subtitle">{_escape(subtitle)}</text>',
        f'<text x="480" y="58" text-anchor="middle" font-size="10" '
        f'fill="#546E7A">{_escape(tagline)}</text>',
    ]
    return "".join(parts)


# ---------------------------------------------------------------------------
# Zone 2: Pill strip
# ---------------------------------------------------------------------------

def _render_pill_strip(pills_header: str, pills: tuple[Pill, ...]) -> str:
    _check("pills_header", pills_header, 60)
    n = len(pills)
    if n < 1:
        raise SlotTypeError("pills must contain at least 1 Pill")
    if n > 8:
        raise SlotTypeError(f"pills: max 8 pills, got {n}")

    # Evenly distribute pills across x=20–940 (total 920 px), 8 px gaps between
    total_gap = (n - 1) * 8
    pill_w = (920 - total_gap) // n

    label_y = _PILL_TOP + 11      # = 73
    pill_top = _PILL_TOP + 15     # = 77

    parts = [
        f'<text x="20" y="{label_y}" font-size="7.5" font-weight="bold" '
        f'fill="#90A4AE" letter-spacing="0.8">{_escape(pills_header)}</text>',
    ]
    for i, pill in enumerate(pills):
        _check(f"pills[{i}].number", pill.number, 2)
        _check(f"pills[{i}].title", pill.title, 40)
        _check(f"pills[{i}].lines[0]", pill.lines[0], 55)
        _check(f"pills[{i}].lines[1]", pill.lines[1], 55)
        px = 20 + i * (pill_w + 8)
        cx = px + pill_w // 2
        parts += [
            f'<rect x="{px}" y="{pill_top}" width="{pill_w}" height="38" rx="4" '
            f'fill="{pill.fill}" stroke="{pill.stroke}" stroke-width="1.2"/>',
            f'<text x="{cx}" y="{pill_top + 13}" text-anchor="middle" '
            f'font-size="9" font-weight="bold" fill="{pill.stroke}">'
            f'{_escape(pill.number)} {_escape(pill.title)}</text>',
            f'<text x="{cx}" y="{pill_top + 26}" text-anchor="middle" '
            f'font-size="8.5" fill="{pill.stroke}">{_escape(pill.lines[0])}</text>',
            f'<text x="{cx}" y="{pill_top + 36}" text-anchor="middle" '
            f'font-size="8" fill="{pill.stroke}">{_escape(pill.lines[1])}</text>',
        ]
    return "".join(parts)


# ---------------------------------------------------------------------------
# Zone 3: Column header bar
# ---------------------------------------------------------------------------

def _render_col_header(col_headers: tuple[str, ...]) -> str:
    if len(col_headers) != 8:
        raise SlotTypeError(f"col_headers must have exactly 8 entries, got {len(col_headers)}")
    for i, h in enumerate(col_headers):
        _check(f"col_headers[{i}]", h, 30)

    label_y = _COL_HDR_Y + 14   # = 138
    parts = [
        f'<rect x="10" y="{_COL_HDR_Y}" width="940" height="{_COL_HDR_H}" '
        f'fill="#37474F" rx="2"/>',
    ]
    for i, (cx, text) in enumerate(zip(_COL_MIDS, col_headers)):
        parts.append(
            f'<text x="{cx}" y="{label_y}" text-anchor="middle" '
            f'font-size="8.5" font-weight="bold" fill="#ECEFF1">{_escape(text)}</text>'
        )
    return "".join(parts)


# ---------------------------------------------------------------------------
# Zone: target-zone overlay (drawn under rows, over header bar)
# ---------------------------------------------------------------------------

def _render_target_zone(
    target_zone: Optional[tuple[float, float]],
    target_zone_label: str,
    rows_bottom: int,
    footer_y: int,
) -> str:
    if target_zone is None:
        return ""
    lo_pct, hi_pct = target_zone
    x1 = int(_BAR_X + lo_pct * _BAR_SCALE)
    x2 = int(_BAR_X + hi_pct * _BAR_SCALE)
    zone_w = x2 - x1
    zone_top = _COL_HDR_Y           # overlay starts at column header top
    zone_h = rows_bottom - zone_top  # ends at last row's bottom
    badge_y = footer_y + 10         # target zone badge in footer

    _check("target_zone_label", target_zone_label, 12)

    badge_w = max(56, len(target_zone_label) * 7 + 8)
    badge_cx = x1 + zone_w // 2

    parts = [
        # Semi-transparent fill band
        f'<rect x="{x1}" y="{zone_top}" width="{zone_w}" height="{zone_h}" '
        f'fill="#C8E6C9" opacity="0.30"/>',
        # Dashed boundary lines
        f'<line x1="{x1}" y1="{zone_top}" x2="{x1}" y2="{rows_bottom}" '
        f'stroke="#43A047" stroke-width="1" stroke-dasharray="4,3" opacity="0.55"/>',
        f'<line x1="{x2}" y1="{zone_top}" x2="{x2}" y2="{rows_bottom}" '
        f'stroke="#43A047" stroke-width="1" stroke-dasharray="4,3" opacity="0.55"/>',
        # Badge in footer
        f'<rect x="{badge_cx - badge_w//2}" y="{badge_y}" width="{badge_w}" height="14" '
        f'rx="3" fill="#C8E6C9"/>',
        f'<text x="{badge_cx}" y="{badge_y + 11}" text-anchor="middle" '
        f'font-size="8" font-weight="bold" fill="#1B5E20">{_escape(target_zone_label)}</text>',
    ]
    return "".join(parts)


# ---------------------------------------------------------------------------
# Zone 4: Data rows
# ---------------------------------------------------------------------------

def _render_row(row: Row, row_y: int, row_h: int) -> str:
    _check(f"row[{row.letter}].name", row.name, 32)
    for di, d in enumerate(row.description):
        _check(f"row[{row.letter}].description[{di}]", d, 45)
    if row.annotation:
        _check(f"row[{row.letter}].annotation", row.annotation, 45)
    _check(f"row[{row.letter}].panel_sublabel", row.panel_sublabel, 14)
    _check(f"row[{row.letter}].panel_type", row.panel_type, 12)
    _check(f"row[{row.letter}].level", row.level, 10)
    _check(f"row[{row.letter}].level_sublabel", row.level_sublabel, 14)

    # --- Vertical position tables keyed by row height ---
    if row_h >= 60:
        p = {
            "letter":    row_y + 28,
            "sublabel":  row_y + 43,
            "type_lbl":  row_y + 54,
            "name":      row_y + 17,
            "desc":      [row_y + 31, row_y + 44, row_y + 57],
            "annot":     row_y + 56,
            "n":         row_y + 34,
            "badge0":    row_y + 9,
            "badge1":    row_y + 28,
            "var":       [row_y + 19, row_y + 31, row_y + 43],
            "mod_pill":  row_y + 13,
            "level":     row_y + 28,
            "level_sub": row_y + 41,
            "bar":       row_y + 19,
            "bar_tick1": row_y + 31,
            "bar_tick2": row_y + 36,
            "bar_lbl":   row_y + 53,
        }
    elif row_h >= 50:
        p = {
            "letter":    row_y + 27,
            "sublabel":  row_y + 42,
            "type_lbl":  None,
            "name":      row_y + 18,
            "desc":      [row_y + 32, row_y + 45],
            "annot":     row_y + 46,
            "n":         row_y + 29,
            "badge0":    row_y + 9,
            "badge1":    row_y + 28,
            "var":       [row_y + 19, row_y + 31, row_y + 43],
            "mod_pill":  row_y + 18,
            "level":     row_y + 26,
            "level_sub": row_y + 39,
            "bar":       row_y + 15,
            "bar_tick1": row_y + 27,
            "bar_tick2": row_y + 32,
            "bar_lbl":   row_y + 46,
        }
    else:  # h=40
        p = {
            "letter":    row_y + 24,
            "sublabel":  None,
            "type_lbl":  None,
            "name":      row_y + 17,
            "desc":      [row_y + 31],
            "annot":     None,
            "n":         row_y + 24,
            "badge0":    row_y + 9,
            "badge1":    None,
            "var":       [row_y + 14, row_y + 26],
            "mod_pill":  row_y + 13,
            "level":     row_y + 20,
            "level_sub": row_y + 33,
            "bar":       row_y + 14,
            "bar_tick1": row_y + 26,
            "bar_tick2": row_y + 31,
            "bar_lbl":   row_y + 36,
        }

    parts = []

    # Row background + separator line
    parts.append(
        f'<rect x="10" y="{row_y}" width="940" height="{row_h}" '
        f'fill="{row.row_color}"/>'
    )
    parts.append(
        f'<line x1="10" y1="{row_y + row_h}" x2="950" y2="{row_y + row_h}" '
        f'stroke="#ECEFF1" stroke-width="1"/>'
    )

    # Letter block
    parts.append(
        f'<rect x="10" y="{row_y}" width="80" height="{row_h}" '
        f'fill="{row.letter_color}"/>'
    )
    parts.append(
        f'<text x="50" y="{p["letter"]}" text-anchor="middle" '
        f'font-size="22" font-weight="bold" fill="white">{_escape(row.letter)}</text>'
    )
    if row.panel_sublabel and p["sublabel"] is not None:
        parts.append(
            f'<text x="50" y="{p["sublabel"]}" text-anchor="middle" '
            f'font-size="7.5" fill="white" opacity="0.85">{_escape(row.panel_sublabel)}</text>'
        )
    if row.panel_type and p["type_lbl"] is not None:
        parts.append(
            f'<text x="50" y="{p["type_lbl"]}" text-anchor="middle" '
            f'font-size="7" fill="white" opacity="0.70">{_escape(row.panel_type)}</text>'
        )

    # Name zone
    parts.append(
        f'<text x="92" y="{p["name"]}" font-size="10" font-weight="bold" '
        f'fill="{row.letter_color}">{_escape(row.name)}</text>'
    )
    # Clip desc positions when annotation is present — last slot reserved for annotation
    active_desc_slots = p["desc"]
    if row.annotation and p["annot"] is not None and len(active_desc_slots) > 1:
        active_desc_slots = active_desc_slots[:-1]
    for di, (desc_text, desc_y) in enumerate(zip(row.description, active_desc_slots)):
        parts.append(
            f'<text x="92" y="{desc_y}" font-size="8.5" '
            f'fill="#455A64">{_escape(desc_text)}</text>'
        )
    if row.annotation and p["annot"] is not None:
        parts.append(
            f'<text x="92" y="{p["annot"]}" font-size="8" '
            f'fill="{row.annotation_color}">{_escape(row.annotation)}</text>'
        )

    # n column
    parts.append(
        f'<text x="283" y="{p["n"]}" text-anchor="middle" '
        f'font-size="14" font-weight="bold" fill="{row.letter_color}">{row.n}</text>'
    )

    # Category badges (col 3, up to 4 badges in 2 rows of 2)
    for bi, badge in enumerate(row.category_badges[:4]):
        brow = bi // 2
        bcol = bi % 2
        bx = _BADGE_X0 + bcol * _BADGE_DX
        by_key = "badge0" if brow == 0 else "badge1"
        by = p[by_key]
        if by is None:
            break  # row too short for second badge row
        _check(f"row[{row.letter}].category_badges[{bi}].text", badge.text, 8)
        if badge.color is not None:
            # Pill badge
            bcx = bx + _BADGE_W // 2
            parts += [
                f'<rect x="{bx}" y="{by}" width="{_BADGE_W}" height="{_BADGE_H}" '
                f'rx="6" fill="{badge.color}"/>',
                f'<text x="{bcx}" y="{by + 10}" text-anchor="middle" '
                f'font-size="8" font-weight="bold" fill="{badge.text_color}">'
                f'{_escape(badge.text)}</text>',
            ]
        else:
            # Plain text badge
            bcx = bx + _BADGE_W // 2
            parts.append(
                f'<text x="{bcx}" y="{by + 10}" text-anchor="middle" '
                f'font-size="7.5" fill="#546E7A">{_escape(badge.text)}</text>'
            )

    # Variant text lines (col 4, x=394)
    for vi, (var_text, var_y) in enumerate(zip(row.variant_lines[:3], p["var"])):
        _check(f"row[{row.letter}].variant_lines[{vi}]", var_text, 30)
        weight = "bold" if vi == 0 else "normal"
        fill = row.letter_color if vi == 0 else ("#37474F" if vi == 1 else "#78909C")
        fs = 8 if vi == 0 else 7.5
        parts.append(
            f'<text x="394" y="{var_y}" font-size="{fs}" font-weight="{weight}" '
            f'fill="{fill}">{_escape(var_text)}</text>'
        )

    # Modifier (col 5, x=511)
    mod_cx = _MOD_CX     # = 523
    if row.modifier is not None:
        mod = row.modifier
        _check(f"row[{row.letter}].modifier.text", mod.text, 8)
        if mod.color is not None:
            # Pill modifier badge
            mod_top = p["mod_pill"]
            parts += [
                f'<rect x="{_MOD_X}" y="{mod_top}" width="{_MOD_W}" height="{_MOD_H}" '
                f'rx="7" fill="{mod.color}"/>',
                f'<text x="{mod_cx}" y="{mod_top + 11}" text-anchor="middle" '
                f'font-size="9" font-weight="bold" fill="{mod.text_color}">'
                f'{_escape(mod.text)}</text>',
            ]
        else:
            # Plain text modifier
            parts.append(
                f'<text x="{mod_cx}" y="{p["level"]}" text-anchor="middle" '
                f'font-size="10" fill="#546E7A">{_escape(mod.text)}</text>'
            )
    else:
        parts.append(
            f'<text x="{mod_cx}" y="{p["level"]}" text-anchor="middle" '
            f'font-size="9.5" fill="#90A4AE">—</text>'
        )

    # Level column (col 6)
    if row.level:
        parts.append(
            f'<text x="{_LEVEL_CX}" y="{p["level"]}" text-anchor="middle" '
            f'font-size="10" font-weight="bold" fill="{row.letter_color}">'
            f'{_escape(row.level)}</text>'
        )
    if row.level_sublabel and p["level_sub"] is not None:
        parts.append(
            f'<text x="{_LEVEL_CX}" y="{p["level_sub"]}" text-anchor="middle" '
            f'font-size="7.5" fill="#90A4AE">{_escape(row.level_sublabel)}</text>'
        )

    # Outcome bar (col 7)
    bar_y = p["bar"]
    tick1 = p["bar_tick1"]
    tick2 = p["bar_tick2"]
    bar_lbl_y = p["bar_lbl"]

    # Background track
    parts.append(
        f'<rect x="{_BAR_X}" y="{bar_y}" width="{_BAR_W}" height="12" '
        f'rx="6" fill="#ECEFF1"/>'
    )

    # Filled range
    if row.outcome_range is not None:
        lo, hi = row.outcome_range
        fill_x = int(_BAR_X + lo * _BAR_SCALE)
        fill_w = max(4, int((hi - lo) * _BAR_SCALE))
        parts.append(
            f'<rect x="{fill_x}" y="{bar_y}" width="{fill_w}" height="12" '
            f'rx="6" fill="{row.outcome_color}" opacity="0.80"/>'
        )
        if row.outcome_uncertain:
            unc_cx = fill_x + fill_w // 2
            parts.append(
                f'<text x="{unc_cx}" y="{bar_y + 10}" text-anchor="middle" '
                f'font-size="9" fill="white" font-weight="bold">?</text>'
            )
        # Value labels below tick area (min and max)
        lbl_lo = f"{lo:.0f}%"
        lbl_hi = f"{hi:.0f}%"
        parts += [
            f'<text x="{fill_x}" y="{bar_lbl_y}" text-anchor="middle" '
            f'font-size="7.5" fill="#78909C">{_escape(lbl_lo)}</text>',
            f'<text x="{fill_x + fill_w}" y="{bar_lbl_y}" text-anchor="middle" '
            f'font-size="7.5" fill="#78909C">{_escape(lbl_hi)}</text>',
        ]

    # Tick marks (always shown)
    for tx in _BAR_TICK_XS:
        color = "#43A047" if tx in (870, 902) else "#B0BEC5"
        parts.append(
            f'<line x1="{tx}" y1="{tick1}" x2="{tx}" y2="{tick2}" '
            f'stroke="{color}" stroke-width="0.8"/>'
        )

    return "".join(parts)


# ---------------------------------------------------------------------------
# Zone 5: Footer
# ---------------------------------------------------------------------------

def _render_footer(
    legend_entries: tuple[LegendEntry, ...],
    strategy_annotation: str,
    total_callout: str,
    footer_y: int,
) -> str:
    _check("strategy_annotation", strategy_annotation, 160)
    _check("total_callout", total_callout, 130)

    tick_y = footer_y + 6          # axis tick label y
    swatch_y = footer_y + 22       # legend swatch top
    hr_y = footer_y + 42           # horizontal rule
    strategy_y = footer_y + 55     # strategy annotation
    callout_y = footer_y + 67      # total callout

    parts = []

    # Bar-axis tick labels (always rendered; target-zone ticks bold green)
    _TICK_LABELS = ("0%", "20%", "40%", "60%", "80%", "90%", "100%")
    for tx, lbl in zip(_BAR_TICK_XS, _TICK_LABELS):
        color = "#2E7D32" if tx in (870, 902) else "#90A4AE"
        weight = "bold" if tx in (870, 902) else "normal"
        parts.append(
            f'<text x="{tx}" y="{tick_y}" text-anchor="middle" '
            f'font-size="7.5" fill="{color}" font-weight="{weight}">{lbl}</text>'
        )

    # Legend swatches (left-aligned, up to 6 entries; 160 px spacing)
    swatch_x = 20
    for ei, entry in enumerate(legend_entries[:6]):
        _check(f"legend_entries[{ei}].label", entry.label, 40)
        parts += [
            f'<rect x="{swatch_x}" y="{swatch_y}" width="12" height="12" '
            f'rx="2" fill="{entry.color}"/>',
            f'<text x="{swatch_x + 17}" y="{swatch_y + 10}" '
            f'font-size="8.5" fill="#37474F">{_escape(entry.label)}</text>',
        ]
        swatch_x += 160

    # Horizontal rule + strategy annotation
    parts += [
        f'<line x1="20" y1="{hr_y}" x2="940" y2="{hr_y}" '
        f'stroke="#B0BEC5" stroke-width="0.5"/>',
        f'<text x="20" y="{strategy_y}" font-size="8.5" '
        f'fill="#546E7A" font-style="italic">{_escape(strategy_annotation)}</text>',
    ]

    # Total callout
    parts.append(
        f'<text x="480" y="{callout_y}" text-anchor="middle" '
        f'font-size="8" fill="#90A4AE">{_escape(total_callout)}</text>'
    )

    return "".join(parts)


# ---------------------------------------------------------------------------
# Public API: render() and preview()
# ---------------------------------------------------------------------------

def render(
    *,
    title: str,
    subtitle: str,
    tagline: str,
    pills_header: str,
    pills: tuple[Pill, ...],
    col_headers: tuple[str, ...],
    rows: tuple[Row, ...],
    target_zone: Optional[tuple[float, float]] = None,
    target_zone_label: str = "▲ TARGET",
    legend_entries: tuple[LegendEntry, ...] = (),
    strategy_annotation: str = "",
    total_callout: str = "",
) -> str:
    """Render an SCF025 figure.

    Args:
        title:               Main figure title, max 70 chars
        subtitle:            Second header line, max 110 chars
        tagline:             Third header line, max 110 chars
        pills_header:        Small label above the pill strip, max 60 chars
        pills:               3–6 colored finding pills
        col_headers:         Exactly 8 column header labels (left to right), max 22 each
        rows:                3–10 data rows; row height auto-computed from text content
        target_zone:         (min_pct, max_pct) light-green band overlay, or None
        target_zone_label:   Badge label for target zone (max 12 chars)
        legend_entries:      Up to 6 swatch+label entries for the footer legend
        strategy_annotation: Italic footer note, max 150 chars
        total_callout:       Centered footer metadata line, max 110 chars

    Returns:
        SVG string.
    """
    if not rows:
        raise SlotTypeError("rows must contain at least 1 Row")
    if len(rows) > 10:
        raise SlotTypeError(f"rows: max 10 rows, got {len(rows)}")

    row_heights = [_row_height(r) for r in rows]
    rows_total = sum(row_heights)
    vb_h = _ROWS_TOP + rows_total + _FOOTER_H
    rows_bottom = _ROWS_TOP + rows_total
    footer_y = rows_bottom  # footer starts immediately after last row

    parts = [svg_open((0, 0, _VB_W, vb_h))]
    parts.append(_render_header(title, subtitle, tagline))
    parts.append(_render_pill_strip(pills_header, pills))
    parts.append(_render_col_header(col_headers))

    # Draw data rows, accumulating y
    row_y = _ROWS_TOP
    for row, row_h in zip(rows, row_heights):
        parts.append(_render_row(row, row_y, row_h))
        row_y += row_h

    # Target-zone overlay (drawn after rows so it's above them visually)
    parts.append(
        _render_target_zone(target_zone, target_zone_label, rows_bottom, footer_y)
    )

    parts.append(
        _render_footer(legend_entries, strategy_annotation, total_callout, footer_y)
    )
    parts.append(svg_close())
    return "".join(parts)


def preview(show_outlines: bool = True) -> str:
    """Render a generic preview using placeholder content only.

    No domain-specific vocabulary appears here — all text is generic.
    """
    _pills = (
        Pill(
            number="①",
            title="Finding One",
            lines=("First body line of finding one", "Second body line here"),
            fill="#E8F5E9",
            stroke="#43A047",
        ),
        Pill(
            number="②",
            title="Finding Two",
            lines=("Body text for finding two.", "Additional context line"),
            fill="#FFEBEE",
            stroke="#EF5350",
        ),
        Pill(
            number="③",
            title="Finding Three",
            lines=("Third pill first line.", "Third pill second line"),
            fill="#E3F2FD",
            stroke="#1E88E5",
        ),
        Pill(
            number="④",
            title="Finding Four",
            lines=("Fourth finding body text.", "Supporting context here"),
            fill="#F3E5F5",
            stroke="#AB47BC",
        ),
        Pill(
            number="⑤",
            title="Finding Five",
            lines=("Fifth finding body text.", "Additional supporting note"),
            fill="#FFF8E1",
            stroke="#FFB300",
        ),
    )

    _col_headers = (
        "PANEL",
        "NAME & DESCRIPTION",
        "N",
        "CATEGORY",
        "VARIANT",
        "MOD",
        "LEVEL",
        "OUTCOME (%)",
    )

    def _row(letter: str, color: str, bg: str, n_lines: int) -> Row:
        # n_lines = total content lines (name + desc + annotation).
        # For n_lines=4: 2 desc + 1 annotation (not 3 desc + 1 annotation, which would overflow).
        n_desc = max(1, n_lines - 1) if n_lines < 4 else 2
        annot = f"Annotation note for panel {letter}." if n_lines >= 4 else ""
        descs = tuple(f"Description line {i+1} for panel {letter}." for i in range(n_desc))
        return Row(
            letter=letter,
            letter_color=color,
            row_color=bg,
            panel_sublabel=f"{n_lines} items",
            panel_type="Type A" if n_lines >= 4 else "",
            name=f"Panel {letter} — Generic Design Group",
            description=descs,
            annotation=annot,
            annotation_color="#558B2F",
            n=n_lines * 2,
            category_badges=(
                Badge("CAT1", color="#43A047"),
                Badge("CAT2", color="#1E88E5"),
                Badge("CAT3", color="#EF5350"),
            ) if n_lines >= 4 else (Badge("CAT1", color="#43A047"),),
            variant_lines=(
                f"Variant type {letter}-1",
                "Supporting context here",
                "Additional note",
            ) if n_lines >= 4 else (f"Variant type {letter}-1",),
            modifier=Badge("M", color="#1B5E20") if n_lines >= 3 else None,
            level="T1–T5" if n_lines >= 4 else "T2",
            level_sublabel="multi" if n_lines >= 4 else "",
            outcome_range=(20.0, 75.0) if n_lines >= 4 else (5.0, 30.0),
            outcome_color="#1565C0",
            outcome_uncertain=n_lines < 3,
        )

    _rows = (
        _row("A", "#1565C0", "#FFFFFF", 4),
        _row("B", "#1B5E20", "#F9FBE7", 4),
        _row("C", "#388E3C", "#FFFFFF", 4),
        _row("D", "#00695C", "#E0F2F1", 3),
        _row("E", "#E65100", "#FFFFFF", 2),
    )

    _legend = (
        LegendEntry("#43A047", "Category A — permissive"),
        LegendEntry("#EF5350", "Category B — restrictive"),
        LegendEntry("#1B5E20", "Modifier T"),
        LegendEntry("#EF9A9A", "Modifier G"),
    )

    return render(
        title="Generic Design Panel Overview (SCF025 Preview)",
        subtitle="N constructs · Generic target · Preview scaffold — no domain vocabulary",
        tagline="Goal: demonstrate all 5 zones with generic placeholder content only",
        pills_header="KEY DESIGN PRINCIPLES — APPLIED",
        pills=_pills,
        col_headers=_col_headers,
        rows=_rows,
        target_zone=(80.0, 90.0),
        target_zone_label="▲ TARGET",
        legend_entries=_legend,
        strategy_annotation=(
            "Strategy: Panels A–C establish performance ceiling  →  "
            "Panels D–E apply controlled modulation  →  "
            "Target: 80–90% outcome range"
        ),
        total_callout=(
            "Total: N constructs across 5 panels  ·  Generic scaffold preview  ·  SCF025"
        ),
    )
