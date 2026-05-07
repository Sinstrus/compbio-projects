"""SCF032 — Single-Construct Architecture Callout Diagram.

One DNA backbone drawn proportionally to scale. Colored ORF blocks sit below
the backbone at stacked vertical layers — VP1 (layer 0, closest), VP2 (layer 1),
VP3 (layer 2), AAP (layer 3), etc. Region bands (VR-V, VR-VIII) are subtle
background highlights spanning the full construct zone; their muted labels
stagger upward inside the zone (greedy L-to-R, no tick lines).

Up to 18 callout annotation boxes (3 rows × 6 columns, all ABOVE the construct)
connect to specific backbone positions via thin pointer lines ending in a small
filled-circle bulb. Box width is variable: max(40, text_px × 1.3) so each box
just fits its label with 30% padding — no fixed-width columns.

Coordinate system: the horizontal display range [0, bp_view_end] maps linearly
to [BODY_X, BODY_RIGHT]. bp_view_end defaults to max(block.bp_end) + 2% padding,
so the bounding box is always used efficiently with no blank backbone on the right.
Pass an explicit bp_view_end to override. All callout bp_pos values must fall
within [0, bp_view_end].

    px = BODY_X + (bp / bp_view_end) × BODY_W

bp_total records the true construct length for reference; all block/region bp
positions are validated against it. bp_view_end ≤ bp_total.

## Pointed blocks (Block.pointed=True)

Directional ORF and promoter elements can be drawn with a forward-pointing tip
instead of a flat right edge. The shape has five vertices (ABCDE clockwise):

    A ──────── B
    |           ╲
    |            C  ← tip; x = px_x + px_w  (the element's bp_end pixel)
    |           ╱
    E ──────── D

B and D are receded inward by _BLOCK_POINT_NUDGE (5 px) from the tip C.
The tip C sits exactly at the element's bp_end pixel position — the element
does NOT extend beyond its allocated bp-space. Adjacent same-layer blocks are
therefore never entered by the tip; any gap between bp_end and the next
bp_start is preserved as-is. For very narrow blocks the nudge is auto-capped
at px_w / 3 to keep the shape valid. The inside label is centered in the
rectangular body (A→B/D), not the full shape including the tip.

Layout (y-axis, 960 × 400 viewBox):

    y=30    TITLE (22pt bold, centered)
    y=52    subtitle (12pt, optional)

    y=68    ┌──────────┐  ┌──────────┐  …  row 3 above (farthest, 6 slots)
    y=96    └──────────┘  └──────────┘
    y=110   ┌──────────┐  ┌──────────┐      row 2 above (middle, 6 slots)
    y=138   └──────────┘  └──────────┘
    y=152   ┌──────────┐  ┌──────────┐      row 1 above (nearest, 6 slots)
    y=180   └──────────┘  └──────────┘
                 │                │           pointer lines (diagonal OK)
                 ●                ●           bulbs at y=236 (construct top)

    y=244   ─────────────────────────────── DNA backbone ([0, bp_view_end])

    y=248   ┌──────────────────────────────►  layer 0 (e.g. Rep78 / VP1, pointed)
    y=268
    y=272   ┌────────────────────►            layer 1 (e.g. Rep52 / VP2, pointed)
    y=292
    y=296   ┌─────────────────►               layer 2 (e.g. p40 / VP3, pointed)
    y=316
    y=320   ┌─────────►                       layer 3 (e.g. AAP, pointed)
    y=340
    y=344   ┌─────────►                       layer 4
    y=364
           VR labels stagger upward inside zone (row 0 y≈369 → row 1 y≈357 → row 2 y≈345)
    y=374   ─── construct zone bottom ────────────────────────────────

ViewBox: 960 × 400.

Column numbers 1–6 set left-to-right box ordering within each row. Actual
x-positions are computed by the packing algorithm: boxes are placed as tightly
as possible with a minimum gap = 25% of the widest box, centered at the mean
bp_px of the row's annotation targets. Pointer lines are always diagonal from
the packed box center to bp_px(bp_pos).
Callout box width = max(40, text_px × 1.3) — scales to content.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

from ._runtime import (
    SlotOverflowError,
    SlotTypeError,
    measure_text_width,
    svg_close,
    svg_open,
)


# ── Module identity ────────────────────────────────────────────────────────────

BARCODE = "SCF032"
DESCRIPTION = (
    "Single-construct architecture callout — one bp-proportional backbone, ORF "
    "blocks at stacked layers (pointed=True for directional arrow tip), region "
    "bands with greedy-upward-stagger labels inside zone (no tick lines), up to "
    "18 callouts above (3 rows × 6 cols) with variable-width boxes and diagonal "
    "pointer lines; bp_view_end auto-fits bounding box to data (960 × 400)"
)
VIEWBOX = (0, 0, 960, 400)


# ── Layout constants (locked) ──────────────────────────────────────────────────

# Title / subtitle
_TITLE_Y = 30
_TITLE_FONT = 22
_TITLE_FONT_W = "bold"
_TITLE_MAX_PX = 840.0
_SUBTITLE_Y = 52
_SUBTITLE_FONT = 12
_SUBTITLE_MAX_PX = 840.0

# Construct body (horizontal)
_BODY_X = 50        # left edge (px)
_BODY_W = 860       # width (px)
_BODY_RIGHT = _BODY_X + _BODY_W  # 910

# bp-to-pixel mapping — uses [bp_view_start, bp_view_end] as the display range
def _bp_px(bp: int, bp_view_end: int, bp_view_start: int = 0) -> float:
    """Map a bp position to pixel x in [BODY_X, BODY_RIGHT].

    BODY_X corresponds to bp_view_start; BODY_RIGHT to bp_view_end.
    """
    return _BODY_X + (bp - bp_view_start) / (bp_view_end - bp_view_start) * _BODY_W

# 6-column grid for callout text boxes
_N_COLS = 6
_COL_CENTERS: Tuple[int, ...] = tuple(
    round(_BODY_X + _BODY_W / _N_COLS * i + _BODY_W / _N_COLS / 2)
    for i in range(_N_COLS)
)  # (122, 265, 408, 552, 695, 838)

# DNA backbone
_BACKBONE_Y = 244
_BACKBONE_STROKE = "#37474F"
_BACKBONE_STROKE_W = 1.6

# ORF block layers — 5 supported (layer 0 = topmost/closest to backbone)
_BLOCK_H = 20           # px height of each ORF block
_LAYER_GAP = 4          # px vertical gap between layers
_LAYER_TOP_0 = 248      # top y of layer 0 (backbone + 4px)
_LAYER_TOPS: Tuple[int, ...] = tuple(
    _LAYER_TOP_0 + i * (_BLOCK_H + _LAYER_GAP) for i in range(5)
)  # (248, 272, 296, 320, 344)
_MAX_LAYER = 4          # layers 0–4 (5 total)
_BLOCK_POINT_NUDGE = 5    # px right extension for pointed block tip (effective right = px_x + px_w + nudge)
_BLOCK_LABEL_MAX_CHARS = 12
_INSIDE_FONT = 9
_INSIDE_FONT_W = "bold"
_BLOCK_ELEM_STROKE = "#37474F"

# Construct zone (region blocks span this y range)
_CONSTRUCT_TOP = 236    # = _ABOVE_BULB_Y — where bulbs sit / region band starts
_CONSTRUCT_BOT = 374    # layer 4 bottom (364) + 10 px padding

# Above callout zone — 3 rows, all above the construct
_ABOVE_BULB_Y = 236      # bulb center = construct top
_ABOVE_ROW1_TOP = 150    # nearest row top   (30 px box; bottom unchanged at 180)
_ABOVE_ROW1_BOT = 180    # nearest row bottom
_ABOVE_ROW2_TOP = 108    # middle row top    (30 px box; 12 px gap; bottom unchanged at 138)
_ABOVE_ROW2_BOT = 138    # middle row bottom
_ABOVE_ROW3_TOP = 66     # farthest row top  (30 px box; 12 px gap; bottom unchanged at 96)
_ABOVE_ROW3_BOT = 96     # farthest row bottom

_ROW_TOPS = {1: _ABOVE_ROW1_TOP, 2: _ABOVE_ROW2_TOP, 3: _ABOVE_ROW3_TOP}
_ROW_BOTS = {1: _ABOVE_ROW1_BOT, 2: _ABOVE_ROW2_BOT, 3: _ABOVE_ROW3_BOT}

# Callout box geometry — WIDTH IS VARIABLE (text_px × _CALLOUT_WIDTH_PAD)
_CALLOUT_H = 30
_CALLOUT_RX = 5
_CALLOUT_FONT = 10
_CALLOUT_MAX_TEXT_PX = 135.0   # max text pixel width (box_w ≤ 175 px at 1.3×)
_CALLOUT_WIDTH_PAD = 1.3       # box_w = max(40, text_px × this)
_CALLOUT_MIN_GAP = 10.0        # minimum px clearance between adjacent boxes in same row
_CALLOUT_STROKE = "#607D8B"
_CALLOUT_FILL = "white"
_CALLOUT_TEXT_COLOR = "#37474F"

# Pointer / bulb
_POINTER_STROKE = "#607D8B"
_POINTER_STROKE_W = 1.2
_BULB_R = 3
_BULB_FILL = "#607D8B"

# Region block labels — staggered upward inside the construct zone (no tick lines)
# Row 0 sits at the natural bottom-of-zone position; additional rows float upward.
# Greedy L-to-R assignment ensures same-row labels never overlap horizontally.
_REGION_LABEL_FONT = 11
_REGION_LABEL_MAX_PX = 110.0
_RLABEL_BASE_Y = _CONSTRUCT_BOT - 6     # row 0 baseline: inside zone near bottom
_RLABEL_ROW_STEP = 15                   # px upward per additional stagger row
_RLABEL_ROW_Y = tuple(
    _RLABEL_BASE_Y - i * _RLABEL_ROW_STEP for i in range(3)
)  # (369, 357, 345) — 3 rows maximum
_RLABEL_GAP = 4.0    # min horizontal gap (px) between same-row label boxes

# Right-padding fraction applied when auto-computing bp_view_end from block data
_AUTO_VIEW_PAD = 0.02   # 2% of max block bp_end


# ── Dataclasses ────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class Block:
    """One ORF or feature block, positioned proportionally on the construct.

    label     — displayed inside the block if it fits; max 12 chars. Label
                threshold is text-width-based: block px_w must be ≥ text_w + 6.
    bp_start  — base pair start position (0-indexed, within [0, bp_total]).
    bp_end    — base pair end position (> bp_start).
    color     — hex fill (e.g. "#1565C0"). Text color auto-selected for contrast.
    layer     — vertical layer: 0 = topmost (closest to backbone), 1, 2, 3, 4.
                Blocks on the same layer should not overlap in bp coordinates.
    pointed   — if True, the right end is drawn as a small forward-pointing tip
                (_BLOCK_POINT_NUDGE px past the rectangle right edge). The true
                effective right boundary is px_x + px_w + _BLOCK_POINT_NUDGE;
                callers must account for this when placing adjacent same-layer blocks.
    """
    label: str
    bp_start: int
    bp_end: int
    color: str
    layer: int
    pointed: bool = False


@dataclass(frozen=True)
class RegionBlock:
    """Subtle background band aligned to construct bp coordinates.

    Spans the full construct-zone height (construct top → construct bottom).
    A muted 8pt label is drawn near the BOTTOM of the band.
    label    — max pixel width 62 px at 8pt; muted label at zone bottom.
    bp_start — base pair start (same coordinate as Block).
    bp_end   — base pair end; must be > bp_start + 1.
    color    — light tint fill (e.g. "#E3F2FD"). No border.
    Region blocks must NOT overlap each other.
    """
    label: str
    bp_start: int
    bp_end: int
    color: str


@dataclass(frozen=True)
class Callout:
    """Annotation callout box with a pointer line (always above the construct).

    text   — annotation text; pixel-width validated at 10pt (max 135 px).
    column — 1–6: left-to-right ordering index within each row. Box x-positions
             are computed by the row packing algorithm (tightly packed with 25%
             of widest box as min gap), not anchored to a fixed column grid.
    bp_pos — base pair position where the pointer line touches the construct top.
             Must be within [0, bp_view_end]. Pointer is always diagonal from
             the packed box center to bp_px(bp_pos).
    row    — 1 (nearest to construct), 2 (middle), or 3 (farthest from construct).
    Each (column, row) pair must be unique within a render call.
    Box width = max(40, text_px × 1.3) — sized to content, not fixed.
    """
    text: str
    column: int
    bp_pos: int
    row: int


# ── Internal helpers ───────────────────────────────────────────────────────────

def _esc(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _auto_text_color(hex_fill: str) -> str:
    r = int(hex_fill[1:3], 16)
    g = int(hex_fill[3:5], 16)
    b = int(hex_fill[5:7], 16)
    brightness = 0.299 * r + 0.587 * g + 0.114 * b
    return "white" if brightness < 128 else "#1A1A1A"


def _baseline_centered(center_y: float, font_size: float) -> float:
    return center_y + font_size * 0.72 / 2


def _pointed_block_path(px_x: float, px_w: float, layer_top: int) -> str:
    """SVG path for a directional block with a pointed right end.

    Five vertices (clockwise, ABCDE):
        A  top-left    (px_x, layer_top)                     — rounded corner
        B  top-right   (px_x + px_w - nudge, layer_top)      — receded from tip
        C  tip         (px_x + px_w, mid_y)                  — TRUE right edge
        D  bottom-right (px_x + px_w - nudge, layer_top + h) — receded from tip
        E  bottom-left  (px_x, layer_top + h)                — rounded corner

    C sits at px_x + px_w — the element's allocated right edge in bp-space.
    B and D are receded inward by _BLOCK_POINT_NUDGE, creating the point.
    The element does NOT extend beyond its bp_end pixel position.
    """
    rx = 3
    h = _BLOCK_H
    mid_y = layer_top + h / 2
    nudge = min(_BLOCK_POINT_NUDGE, px_w / 3)  # cap for very narrow blocks
    b_x = px_x + px_w - nudge
    c_x = px_x + px_w
    return (
        f"M {px_x+rx:.1f},{layer_top} "
        f"L {b_x:.1f},{layer_top} "
        f"L {c_x:.1f},{mid_y:.1f} "
        f"L {b_x:.1f},{layer_top+h} "
        f"L {px_x+rx:.1f},{layer_top+h} "
        f"A {rx},{rx},0,0,1,{px_x:.1f},{layer_top+h-rx:.1f} "
        f"L {px_x:.1f},{layer_top+rx:.1f} "
        f"A {rx},{rx},0,0,1,{px_x+rx:.1f},{layer_top} "
        f"Z"
    )


def _compute_bp_view_end(
    blocks: tuple,
    region_blocks: tuple,
    bp_total: int,
    bp_view_end_override: Optional[int],
) -> int:
    """Return the display right edge in bp.

    If bp_view_end_override is given, use it directly (still validated in
    _validate). Otherwise auto-compute as max(block.bp_end) + 2% padding,
    clamped to bp_total. This ensures no blank backbone extends past the
    rightmost ORF block.
    """
    if bp_view_end_override is not None:
        return bp_view_end_override
    candidates = [b.bp_end for b in blocks]
    if not candidates:
        return bp_total
    max_bp = max(candidates)
    padded = int(max_bp * (1.0 + _AUTO_VIEW_PAD))
    return min(padded, bp_total)


# ── Validation ─────────────────────────────────────────────────────────────────

def _validate(
    title: str,
    bp_total: int,
    bp_view_start: int,
    bp_view_end: int,
    blocks: tuple,
    region_blocks: tuple,
    callouts: tuple,
    subtitle: Optional[str],
) -> None:
    # title
    if not isinstance(title, str) or not title.strip():
        raise SlotTypeError("title must be a non-empty str")
    w = measure_text_width(title, _TITLE_FONT, _TITLE_FONT_W)
    if w > _TITLE_MAX_PX:
        raise SlotOverflowError("title", w, _TITLE_MAX_PX, unit=" px")

    # subtitle
    if subtitle is not None:
        if not isinstance(subtitle, str):
            raise SlotTypeError("subtitle must be str or None")
        w = measure_text_width(subtitle, _SUBTITLE_FONT, "normal")
        if w > _SUBTITLE_MAX_PX:
            raise SlotOverflowError("subtitle", w, _SUBTITLE_MAX_PX, unit=" px")

    # bp_total
    if not isinstance(bp_total, int) or bp_total < 1:
        raise SlotTypeError("bp_total must be a positive int")

    # bp_view_start
    if not isinstance(bp_view_start, int) or bp_view_start < 0:
        raise SlotTypeError(f"bp_view_start must be a non-negative int, got {bp_view_start}")

    # bp_view_end
    if not isinstance(bp_view_end, int) or bp_view_end <= bp_view_start or bp_view_end > bp_total:
        raise SlotTypeError(
            f"bp_view_end must be an int in ({bp_view_start}, {bp_total}], "
            f"got {bp_view_end}"
        )

    # blocks
    if not blocks:
        raise SlotTypeError("blocks must contain at least one Block")
    for bi, b in enumerate(blocks):
        if not isinstance(b, Block):
            raise SlotTypeError(f"blocks[{bi}] must be a Block")
        if len(b.label) > _BLOCK_LABEL_MAX_CHARS:
            raise SlotOverflowError(
                f"blocks[{bi}].label", len(b.label), _BLOCK_LABEL_MAX_CHARS, unit=" chars"
            )
        if not (0 <= b.bp_start < b.bp_end <= bp_total):
            raise SlotTypeError(
                f"blocks[{bi}]: bp_start={b.bp_start} bp_end={b.bp_end} "
                f"out of range [0, {bp_total}]"
            )
        if b.layer not in range(_MAX_LAYER + 1):
            raise SlotTypeError(
                f"blocks[{bi}].layer must be 0–{_MAX_LAYER}, got {b.layer}"
            )

    # region_blocks
    sorted_rb = sorted(region_blocks, key=lambda r: r.bp_start)
    for ri, rb in enumerate(sorted_rb):
        if not isinstance(rb, RegionBlock):
            raise SlotTypeError(f"region_blocks[{ri}] must be a RegionBlock")
        w_lbl = measure_text_width(rb.label, _REGION_LABEL_FONT, "normal")
        if w_lbl > _REGION_LABEL_MAX_PX:
            raise SlotOverflowError(
                f"region_blocks[{ri}].label", w_lbl, _REGION_LABEL_MAX_PX, unit=" px"
            )
        if not (0 <= rb.bp_start < rb.bp_end <= bp_total):
            raise SlotTypeError(
                f"region_blocks[{ri}]: bp_start={rb.bp_start} bp_end={rb.bp_end} "
                f"out of range [0, {bp_total}]"
            )
        if ri > 0 and sorted_rb[ri - 1].bp_end > rb.bp_start:
            raise SlotTypeError(
                f"region_blocks[{ri - 1}] and [{ri}] overlap: "
                f"[{sorted_rb[ri-1].bp_start},{sorted_rb[ri-1].bp_end}) overlaps "
                f"[{rb.bp_start},{rb.bp_end})"
            )

    # callouts — bp_pos validated against bp_view_end (display range), not bp_total
    if len(callouts) > 18:
        raise SlotTypeError(f"max 18 callouts (3 rows × 6 cols), got {len(callouts)}")
    seen: set = set()
    for ci, c in enumerate(callouts):
        if not isinstance(c, Callout):
            raise SlotTypeError(f"callouts[{ci}] must be a Callout")
        if c.column not in range(1, _N_COLS + 1):
            raise SlotTypeError(f"callouts[{ci}].column must be 1–{_N_COLS}, got {c.column}")
        if c.row not in (1, 2, 3):
            raise SlotTypeError(f"callouts[{ci}].row must be 1, 2, or 3, got {c.row}")
        if not (bp_view_start <= c.bp_pos <= bp_view_end):
            raise SlotTypeError(
                f"callouts[{ci}].bp_pos={c.bp_pos} out of display range "
                f"[{bp_view_start}, {bp_view_end}]"
            )
        w = measure_text_width(c.text, _CALLOUT_FONT, "normal")
        if w > _CALLOUT_MAX_TEXT_PX:
            raise SlotOverflowError(
                f"callouts[{ci}].text", w, _CALLOUT_MAX_TEXT_PX, unit=" px"
            )
        key = (c.column, c.row)
        if key in seen:
            raise SlotTypeError(
                f"duplicate callout slot (column={c.column}, row={c.row})"
            )
        seen.add(key)


# ── Rendering ──────────────────────────────────────────────────────────────────

def _render_title_block(title: str, subtitle: Optional[str]) -> str:
    parts = [
        f'<text x="480" y="{_TITLE_Y}" text-anchor="middle" '
        f'font-size="{_TITLE_FONT}" font-weight="bold" fill="#1A1A1A">'
        f"{_esc(title)}</text>"
    ]
    if subtitle:
        parts.append(
            f'<text x="480" y="{_SUBTITLE_Y}" text-anchor="middle" '
            f'font-size="{_SUBTITLE_FONT}" fill="#546E7A">'
            f"{_esc(subtitle)}</text>"
        )
    return "".join(parts)


def _place_region_labels(
    region_blocks: tuple, bp_view_end: int, bp_view_start: int = 0
) -> list:
    """Greedy L-to-R row assignment so same-row labels never overlap."""
    result = []
    occupied: list = [[], []]

    for rb in sorted(region_blocks, key=lambda r: r.bp_start):
        center_x = (
            _bp_px(rb.bp_start, bp_view_end, bp_view_start)
            + _bp_px(rb.bp_end, bp_view_end, bp_view_start)
        ) / 2
        text_px = measure_text_width(rb.label, _REGION_LABEL_FONT, "normal")
        half_w = text_px / 2 + _RLABEL_GAP
        iv = (center_x - half_w, center_x + half_w)

        assigned = len(occupied) - 1
        for r, slots in enumerate(occupied):
            collides = any(iv[0] < slot[1] and iv[1] > slot[0] for slot in slots)
            if not collides:
                slots.append(iv)
                assigned = r
                break

        result.append((center_x, assigned, rb.label))

    return result


def _render_region_blocks(
    region_blocks: tuple, bp_view_end: int, bp_view_start: int = 0
) -> str:
    parts = []
    zone_h = _CONSTRUCT_BOT - _CONSTRUCT_TOP

    for rb in sorted(region_blocks, key=lambda r: r.bp_start):
        px_x = _bp_px(rb.bp_start, bp_view_end, bp_view_start)
        px_w = _bp_px(rb.bp_end, bp_view_end, bp_view_start) - px_x
        parts.append(
            f'<rect x="{px_x:.1f}" y="{_CONSTRUCT_TOP}" '
            f'width="{px_w:.1f}" height="{zone_h}" '
            f'rx="4" fill="{rb.color}" stroke="none"/>'
        )

    for center_x, row_idx, label in _place_region_labels(region_blocks, bp_view_end, bp_view_start):
        lbl_y = _RLABEL_ROW_Y[row_idx]
        parts.append(
            f'<text x="{center_x:.1f}" y="{lbl_y}" '
            f'text-anchor="middle" font-size="{_REGION_LABEL_FONT}" '
            f'fill="#78909C">{_esc(label)}</text>'
        )

    return "".join(parts)


def _render_backbone() -> str:
    return (
        f'<line x1="{_BODY_X}" y1="{_BACKBONE_Y}" '
        f'x2="{_BODY_RIGHT}" y2="{_BACKBONE_Y}" '
        f'stroke="{_BACKBONE_STROKE}" stroke-width="{_BACKBONE_STROKE_W}"/>'
    )


def _render_blocks(blocks: tuple, bp_view_end: int, bp_view_start: int = 0) -> str:
    # Render higher-numbered layers first (behind), layer 0 last (front)
    sorted_blocks = sorted(blocks, key=lambda b: -b.layer)
    parts = []
    for b in sorted_blocks:
        px_x = _bp_px(b.bp_start, bp_view_end, bp_view_start)
        px_w = _bp_px(b.bp_end, bp_view_end, bp_view_start) - px_x
        layer_top = _LAYER_TOPS[b.layer]
        block_cy = layer_top + _BLOCK_H / 2

        if b.pointed:
            d = _pointed_block_path(px_x, px_w, layer_top)
            parts.append(
                f'<path d="{d}" fill="{b.color}" stroke="{_BLOCK_ELEM_STROKE}" '
                f'stroke-width="0.7"/>'
            )
        else:
            parts.append(
                f'<rect x="{px_x:.1f}" y="{layer_top}" '
                f'width="{px_w:.1f}" height="{_BLOCK_H}" '
                f'rx="3" fill="{b.color}" stroke="{_BLOCK_ELEM_STROKE}" '
                f'stroke-width="0.7"/>'
            )

        if b.label:
            text_w = measure_text_width(b.label, _INSIDE_FONT, _INSIDE_FONT_W)
            # Threshold on px_w: label must fit in the rectangular body (before the tip)
            nudge_actual = min(_BLOCK_POINT_NUDGE, px_w / 3) if b.pointed else 0.0
            body_w = px_w - nudge_actual
            if body_w >= text_w + 6:
                tc = _auto_text_color(b.color)
                lbl_y = _baseline_centered(block_cy, _INSIDE_FONT)
                cx = px_x + body_w / 2  # centered in the rectangular body
                parts.append(
                    f'<text x="{cx:.1f}" y="{lbl_y:.1f}" '
                    f'text-anchor="middle" font-size="{_INSIDE_FONT}" '
                    f'font-weight="{_INSIDE_FONT_W}" fill="{tc}">'
                    f"{_esc(b.label)}</text>"
                )
    return "".join(parts)


def _pack_callout_row(row_callouts: list, bp_view_end: int, bp_view_start: int = 0) -> list:
    """Compute center_x for each callout box in one row.

    Each box starts at its ideal position — centered directly above its
    bp_pos pointer target. Boxes are displaced from ideal only as much as
    necessary to eliminate overlaps, using iterative forward/backward passes.
    This keeps each pointer line as short and vertical as possible while
    guaranteeing _CALLOUT_MIN_GAP clearance between all adjacent boxes.

    Algorithm:
      1. ideal[i] = bp_px(bp_pos[i]), clamped within body bounds
      2. Forward pass: push right if box overlaps the preceding box
      3. Backward pass: push left if the last box exceeds BODY_RIGHT (shifts
         the whole group) then ensures each box doesn't overlap the next
      4. Re-clamp left if the group underflows BODY_X
      Repeat 4 times (convergence is guaranteed for 1-D non-overlapping sets).

    Returns list of (Callout, center_x) sorted by column (left to right).
    """
    ordered = sorted(row_callouts, key=lambda c: c.column)
    n = len(ordered)
    widths = [
        max(40.0, measure_text_width(c.text, _CALLOUT_FONT, "normal") * _CALLOUT_WIDTH_PAD)
        for c in ordered
    ]
    half = [w / 2 for w in widths]

    # Ideal: box centered above its pointer target, clamped within body
    centers = [
        max(_BODY_X + half[i],
            min(_BODY_RIGHT - half[i],
                _bp_px(c.bp_pos, bp_view_end, bp_view_start)))
        for i, c in enumerate(ordered)
    ]

    for _ in range(4):
        # Forward pass — push each box right of its left neighbor
        for i in range(1, n):
            need = centers[i - 1] + half[i - 1] + _CALLOUT_MIN_GAP + half[i]
            if centers[i] < need:
                centers[i] = need

        # If the rightmost box overflows, shift the whole group left
        overflow = centers[-1] + half[-1] - _BODY_RIGHT
        if overflow > 0:
            for i in range(n):
                centers[i] -= overflow

        # Backward pass — push each box left of its right neighbor
        for i in range(n - 2, -1, -1):
            need = centers[i + 1] - half[i + 1] - _CALLOUT_MIN_GAP - half[i]
            if centers[i] > need:
                centers[i] = need

        # If the leftmost box underflows, shift the whole group right
        underflow = _BODY_X - (centers[0] - half[0])
        if underflow > 0:
            for i in range(n):
                centers[i] += underflow

    return [(c, cx) for c, cx in zip(ordered, centers)]


def _render_callouts(callouts: tuple, bp_view_end: int, bp_view_start: int = 0) -> tuple:
    """Returns (pointers_svg, boxes_svg, bulbs_svg)."""
    by_row: dict = {}
    for c in callouts:
        by_row.setdefault(c.row, []).append(c)

    packed_cx: dict = {}
    for row_id, row_cs in by_row.items():
        for c, cx in _pack_callout_row(row_cs, bp_view_end, bp_view_start):
            packed_cx[(c.column, c.row)] = cx

    pointer_parts: list[str] = []
    box_parts: list[str] = []
    bulb_parts: list[str] = []

    for c in callouts:
        col_x = packed_cx[(c.column, c.row)]
        term_x = _bp_px(c.bp_pos, bp_view_end, bp_view_start)
        box_top = _ROW_TOPS[c.row]
        box_bot = _ROW_BOTS[c.row]
        box_ctr = (box_top + box_bot) / 2

        text_px = measure_text_width(c.text, _CALLOUT_FONT, "normal")
        box_w = max(40.0, text_px * _CALLOUT_WIDTH_PAD)
        bx = col_x - box_w / 2

        # Pointer: diagonal from packed box center to bp_pos bulb
        pointer_parts.append(
            f'<line x1="{col_x:.1f}" y1="{box_bot}" '
            f'x2="{term_x:.1f}" y2="{_ABOVE_BULB_Y}" '
            f'stroke="{_POINTER_STROKE}" stroke-width="{_POINTER_STROKE_W}"/>'
        )

        # Box + centered text
        box_parts.append(
            f'<rect x="{bx:.1f}" y="{box_top}" '
            f'width="{box_w:.1f}" height="{_CALLOUT_H}" '
            f'rx="{_CALLOUT_RX}" fill="{_CALLOUT_FILL}" '
            f'stroke="{_CALLOUT_STROKE}" stroke-width="1.0"/>'
        )
        ty = _baseline_centered(box_ctr, _CALLOUT_FONT)
        box_parts.append(
            f'<text x="{col_x:.1f}" y="{ty:.1f}" '
            f'text-anchor="middle" font-size="{_CALLOUT_FONT}" '
            f'fill="{_CALLOUT_TEXT_COLOR}">'
            f"{_esc(c.text)}</text>"
        )

        # Bulb at pointer terminus on construct top
        bulb_parts.append(
            f'<circle cx="{term_x:.1f}" cy="{_ABOVE_BULB_Y}" r="{_BULB_R}" '
            f'fill="{_BULB_FILL}"/>'
        )

    return "".join(pointer_parts), "".join(box_parts), "".join(bulb_parts)


# ── Public API ─────────────────────────────────────────────────────────────────

def render(
    *,
    title: str,
    bp_total: int,
    blocks: tuple,
    region_blocks: tuple = (),
    callouts: tuple = (),
    subtitle: Optional[str] = None,
    bp_view_end: Optional[int] = None,
    bp_view_start: int = 0,
) -> str:
    """Render SCF032 to an SVG string.

    Parameters
    ----------
    title         : Figure title. Pixel-width validated at 22pt bold.
    bp_total      : True construct length in bp. Block/region positions are
                    validated against this.
    blocks        : Tuple of Block entries (ORF blocks with bp coordinates).
    region_blocks : Optional RegionBlock entries (e.g. VR-V, VR-VIII). Labels
                    stagger upward inside the construct zone. Must not overlap.
    callouts      : Up to 18 Callout entries (3 rows × 6 cols, all above the
                    construct). Box width adapts to text content (×1.3).
    subtitle      : Optional subtitle line.
    bp_view_end   : Right edge of the display in bp. Defaults to
                    max(block.bp_end) + 2% padding so no blank backbone extends
                    past the rightmost ORF. Pass explicitly to zoom or extend.
                    Callout bp_pos values must be ≤ bp_view_end.
    bp_view_start : Left edge of the display in bp (default 0). BODY_X maps to
                    this position; the backbone is drawn from bp_view_start to
                    bp_view_end. Use to crop blank upstream backbone — e.g. set
                    to min(block.bp_start) minus a small margin to put the first
                    element near the left edge without showing empty backbone.
                    All block.bp_start and callout.bp_pos must be ≥ bp_view_start.

    Returns
    -------
    str  SVG document (960 × 400 viewBox).
    """
    bve = _compute_bp_view_end(blocks, region_blocks, bp_total, bp_view_end)
    bvs = bp_view_start
    _validate(title, bp_total, bvs, bve, blocks, region_blocks, callouts, subtitle)

    ptrs, boxes, bulbs = _render_callouts(callouts, bve, bvs)

    return (
        svg_open(VIEWBOX)
        + _render_title_block(title, subtitle)
        + _render_region_blocks(region_blocks, bve, bvs)
        + _render_backbone()
        + _render_blocks(blocks, bve, bvs)
        + ptrs
        + boxes
        + bulbs
        + svg_close()
    )


# ── Preview ────────────────────────────────────────────────────────────────────

def preview(show_outlines: bool = True) -> str:
    """Filler-data render demonstrating all SCF032 features.

    Uses AVD002 (WT AAV9 repcap, 7104 bp) as the reference construct so that
    VR-IV/VR-V/VR-VIII appear at their true native sizes. bp_view_end is left at
    default so the scaffold auto-fits to VP3 end (bp 4589) + 2% padding ≈ 4681.

    Feature coords (0-indexed, from AVD002 GenBank):
      Rep78   496  – 2362  layer 0
      Rep52   1168 – 2362  layer 1
      p40     1875 – 2028  layer 2  (ends ~350 bp before VP1 ATG)
      VP1     2378 – 4589  layer 0
      VP2     2789 – 4589  layer 1
      AAP     2904 – 3498  layer 3
      VP3     2984 – 4589  layer 2
      VR-IV   3731 – 3758  (27 bp — correctly tiny)
      VR-V    3839 – 3893  (54 bp — not inflated by any insert)
      VR-VIII 4118 – 4157  (39 bp)

    With auto bp_view_end = int(4589 × 1.02) = 4681:
      px = 50 + bp / 4681 × 860
    VP3 end (4589) → x ≈ 893; bounding box usage ≈ 98%.
    VR-IV → x ≈ 687–693 (6 px wide); VR-V → x ≈ 706–717 (11 px); VR-VIII → x ≈ 757–764 (7 px).
    """
    BP = 7104  # AVD002 total construct length

    blocks = (
        Block("Rep78", 496,  2362, "#1565C0", layer=0, pointed=True),
        Block("Rep52", 1168, 2362, "#5C8AC0", layer=1, pointed=True),
        Block("p40",   1875, 2028, "#37474F", layer=2, pointed=True),
        Block("VP1",   2378, 4589, "#283593", layer=0, pointed=True),
        Block("VP2",   2789, 4589, "#3949AB", layer=1, pointed=True),
        Block("AAP",   2904, 3498, "#7986CB", layer=3, pointed=True),
        Block("VP3",   2984, 4589, "#5C6BC0", layer=2, pointed=True),
    )

    # True native VR sizes from AVD002 — all correctly small relative to the construct.
    # Four bands: VR-IV, VR-V, VR-VI, VR-VIII — labels stagger in 2 rows below the zone.
    region_blocks = (
        RegionBlock("VR-IV",  3731, 3758, "#E8F5E9"),
        RegionBlock("VR-V",   3839, 3893, "#FFF8E1"),
        RegionBlock("VR-VI",  3956, 3995, "#F3E5F5"),
        RegionBlock("VR-VIII", 4118, 4157, "#E3F2FD"),
    )

    # Callouts span all 3 rows. VR callouts use cols 4-6 (where the VRs map at
    # this scale) so pointers are short diagonals anchored near the bands.
    # All bp_pos ≤ auto bp_view_end ≈ 4681.
    callouts = (
        Callout("Rep78 start",  column=1, bp_pos=496,  row=1),
        Callout("p40 promoter", column=2, bp_pos=1875, row=2),
        Callout("VP1 ATG",      column=3, bp_pos=2378, row=1),
        Callout("VR-IV",        column=4, bp_pos=3744, row=3),
        Callout("VR-V",         column=5, bp_pos=3866, row=2),
        Callout("VR-VIII",      column=6, bp_pos=4137, row=1),
    )

    svg = render(
        title="SCF032 Preview — Construct Architecture Callout",
        subtitle="Proportional bp scale · stacked isoform layers · variable-width callout boxes",
        bp_total=BP,
        blocks=blocks,
        region_blocks=region_blocks,
        callouts=callouts,
        # bp_view_end left at default → auto-fits to VP3 end + 2% padding
    )

    if show_outlines:
        overlays = []
        # Column grid guide lines
        for cx in _COL_CENTERS:
            overlays.append(
                f'<line x1="{cx}" y1="62" x2="{cx}" y2="380" '
                f'stroke="#3b7dd8" stroke-width="0.4" stroke-dasharray="2 3" '
                f'data-slot-outline="col-{cx}"/>'
            )
        # Construct zone boundary
        overlays.append(
            f'<rect x="{_BODY_X}" y="{_CONSTRUCT_TOP}" '
            f'width="{_BODY_W}" height="{_CONSTRUCT_BOT - _CONSTRUCT_TOP}" '
            f'fill="none" stroke="#3b7dd8" stroke-width="0.5" '
            f'stroke-dasharray="3 2" data-slot-outline="construct-zone"/>'
        )
        svg = svg.replace("</svg>", "".join(overlays) + "</svg>")

    return svg


if __name__ == "__main__":
    import sys
    sys.stdout.write(preview())
