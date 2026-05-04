"""Matplotlib panel sizer — deterministic fontsize/rotation choice for SCF### body slots.

See SPEC_matplotlib_panel_sizer.md for the design rationale.

Public API:
    size_axis_labels(...) -> PanelSizing             # Phase 1: pure function
    clearance_px_from_subplots(...) -> float         # Phase 1: helper
    PanelRenderer(...)                               # Phase 2: bundled wrapper
    auto_place_legend(ax, ...) -> PlacementResult    # Phase 3: legend overlap resolver
    place_corner_text(ax, text, ...) -> PlacementResult  # Phase 3: corner badge placer
    add_headroom(ax, ...)                            # Phase 3: ylim escape hatch
    probe_legend_headroom(ax, ...) -> str            # Phase 3: two-pass shared-axis helper
    headroom_ylim_top(top, label, ...) -> float      # Phase 3: two-pass numeric converter
    max_headroom_label(labels) -> str                # Phase 3: two-pass reconciler

Strategy ladder (tried in order until min_padding is satisfied):
    1. default   — target_pt, rotation 0
    2. shrunk    — binary-search smaller pt at rotation 0
    3. rotated   — try rotation 30°, then 90°, repeating shrink at each
    4. flipped   — recommend_horizontal_bars=True (gated by axis_priority != 'y')
    5. failed    — return best-effort min_pt at rotation 90 with warnings

Each rotation is gated by `tick_label_clearance_px` (when supplied): a rotation
that would push the label past the available bottom (axis="x") or left (axis="y")
margin is disqualified — even if its horizontal fit at that rotation is fine.
This prevents the FIG177 bug where rot=90 fit horizontally but got clipped vertically.

Measurement is via render-then-extent (matplotlib `text.get_window_extent()`),
not the 0.5×fontsize×nchars approximation.
"""

from __future__ import annotations

import base64
import io
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Sequence

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


@dataclass(frozen=True)
class PanelSizing:
    fontsize_pt: float
    rotation_deg: int
    horizontal_alignment: str
    recommend_horizontal_bars: bool
    achieved_padding_ratio: float
    strategy: str
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class PlacementResult:
    """Outcome of auto_place_legend / place_corner_text.

    Attributes:
        loc: Chosen position. For legends, a matplotlib `loc` string
            (e.g. ``'upper right'``). For corner text, a ``(x_axes, y_axes)``
            tuple in ``transAxes`` coords.
        headroom_applied: ``''`` (none), ``'half_decade'``, or ``'full_decade'``.
        overlapped: True if every candidate position still overlapped data
            after exhausting the headroom ladder. The artist is left at the
            last-tried position; caller may want to escalate manually.
        overlap_artists: Names of data artists involved in the residual
            overlap (e.g. ``('Line2D',)``, ``('BarContainer',)``). Empty
            tuple on success.
    """
    loc: object
    headroom_applied: str
    overlapped: bool
    overlap_artists: tuple[str, ...] = ()


def clearance_px_from_subplots(
    *,
    figsize_inches: tuple[float, float],
    dpi: int,
    axis: Literal["x", "y"] = "x",
    subplots_bottom: float = 0.13,
    subplots_left: float = 0.16,
    axis_label_pt: float = 0.0,
    axis_label_padding_factor: float = 1.5,
) -> float:
    """Helper: compute available tick-label clearance from subplots_adjust margins.

    SCF### scaffolds render the axis label as SVG text *outside* the matplotlib
    panel — so by default `axis_label_pt=0` (the full margin is tick-label
    budget). Pass a non-zero `axis_label_pt` only if the matplotlib output
    itself contains an in-panel axis label (e.g. via `ax.set_xlabel(...)`).

    Returns clearance in PNG pixels at the caller's render dpi.
    """
    if axis == "x":
        total_margin_px = figsize_inches[1] * dpi * subplots_bottom
    else:
        total_margin_px = figsize_inches[0] * dpi * subplots_left
    axis_label_px = axis_label_pt * dpi / 72.0 * axis_label_padding_factor
    return max(0.0, total_margin_px - axis_label_px)


def _measure(
    labels: Sequence[str],
    fontsize_pt: float,
    rotation_deg: int,
    font_family: str,
    dpi: int,
) -> list[tuple[float, float]]:
    """Return (width_px, height_px) per label as it would render."""
    fig = plt.figure(figsize=(10, 10), dpi=dpi)
    ax = fig.add_subplot(111)
    objs = []
    for lab in labels:
        t = ax.text(0.5, 0.5, lab, fontsize=fontsize_pt,
                    rotation=rotation_deg, family=font_family,
                    transform=ax.transAxes)
        objs.append(t)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    out = []
    for o in objs:
        bb = o.get_window_extent(r)
        out.append((bb.width, bb.height))
    plt.close(fig)
    return out


def _check_label_distinctiveness(labels: list[str]) -> str | None:
    """Return a warning string if labels share a dominant common prefix."""
    if len(labels) < 2:
        return None
    common = os.path.commonprefix(labels)
    mean_len = sum(len(lbl) for lbl in labels) / len(labels)
    if mean_len > 0 and len(common) >= 3 and len(common) / mean_len > 0.5:
        return (f"low prefix distinctiveness: '{common}' is "
                f"{len(common)/mean_len:.0%} of mean label length "
                f"({len(labels)} labels)")
    return None


def _try_fit(
    labels: Sequence[str],
    fontsize_pt: float,
    rotation_deg: int,
    font_family: str,
    dpi: int,
    per_label_budget_px: float,
    min_padding: float,
    *,
    axis: Literal["x", "y"] = "x",
    tick_label_clearance_px: float | None = None,
) -> tuple[bool, float, float]:
    """Try fontsize/rotation. Return (fits, padding_ratio, max_fit_dim_px).

    Two checks must pass for `fits=True`:
      1. **Padding fit** — the dimension parallel to the axis (width for x,
         height for y) leaves ≥ min_padding gap between adjacent labels.
      2. **Clearance fit** — the dimension perpendicular to the axis (height
         for x, width for y) fits within `tick_label_clearance_px` if supplied.

    matplotlib's `get_window_extent()` returns the rotation-aware AABB, so
    `bb.width` / `bb.height` already reflect the rotated label's footprint.
    """
    dims = _measure(labels, fontsize_pt, rotation_deg, font_family, dpi)
    if not dims:
        return True, 1.0, 0.0
    max_w = max(w for w, _ in dims)
    max_h = max(h for _, h in dims)

    if axis == "x":
        fit_dim = max_w        # parallel to axis: budgeted by per_label_budget
        clearance_dim = max_h  # perpendicular: must fit in bottom margin
    else:
        fit_dim = max_h
        clearance_dim = max_w

    # Clearance check (perpendicular)
    if (tick_label_clearance_px is not None
            and clearance_dim > tick_label_clearance_px):
        return False, -1.0, fit_dim

    # Padding check (parallel)
    if fit_dim == 0:
        return True, 1.0, 0.0
    padding = (per_label_budget_px - fit_dim) / fit_dim
    fits = padding >= min_padding
    return fits, padding, fit_dim


def _binary_search_fontsize(
    labels: Sequence[str],
    rotation_deg: int,
    min_pt: float,
    max_pt: float,
    font_family: str,
    dpi: int,
    per_label_budget_px: float,
    min_padding: float,
    *,
    axis: Literal["x", "y"] = "x",
    tick_label_clearance_px: float | None = None,
) -> tuple[float, float] | None:
    """Find largest fontsize in [min_pt, max_pt] that fits with min_padding.

    Returns (fontsize_pt, padding) or None if no size in range fits.
    """
    fits_max, pad_max, _ = _try_fit(
        labels, max_pt, rotation_deg, font_family, dpi,
        per_label_budget_px, min_padding,
        axis=axis, tick_label_clearance_px=tick_label_clearance_px,
    )
    if fits_max:
        return (max_pt, pad_max)

    fits_min, pad_min, _ = _try_fit(
        labels, min_pt, rotation_deg, font_family, dpi,
        per_label_budget_px, min_padding,
        axis=axis, tick_label_clearance_px=tick_label_clearance_px,
    )
    if not fits_min:
        return None

    lo, hi = min_pt, max_pt
    best_pt, best_pad = min_pt, pad_min
    for _ in range(20):
        if hi - lo < 0.1:
            break
        mid = (lo + hi) / 2
        fits, pad, _ = _try_fit(
            labels, mid, rotation_deg, font_family, dpi,
            per_label_budget_px, min_padding,
            axis=axis, tick_label_clearance_px=tick_label_clearance_px,
        )
        if fits:
            best_pt, best_pad = mid, pad
            lo = mid
        else:
            hi = mid

    return (round(best_pt, 1), best_pad)


def size_axis_labels(
    *,
    body_w_px: int,
    body_h_px: int,
    figsize_inches: tuple[float, float],
    dpi: int,
    labels: Sequence[str],
    axis: Literal["x", "y"] = "x",
    axis_priority: Literal["x", "y", "either"] = "y",
    max_pt: float = 11.0,
    target_pt: float = 10.0,
    min_pt: float = 6.0,
    min_padding: float = 0.20,
    allow_rotation: bool = True,
    allow_horizontal_bar_flip: bool = True,
    font_family: str = "DejaVu Sans",
    reserved_axis_fraction: float = 0.85,
    tick_label_clearance_px: float | None = None,
) -> PanelSizing:
    """Pick fontsize/rotation/orientation for a tick-label set.

    Pre-condition: figsize aspect should match body box aspect (caller's job).
    Pass `tick_label_clearance_px` to gate rotations against bottom-margin
    clipping — typically computed by `clearance_px_from_subplots()` from the
    same subplots_adjust margins the caller will pass to matplotlib.

    See SPEC_matplotlib_panel_sizer.md for full algorithm.
    """
    if min_padding < 0:
        raise ValueError(f"min_padding must be ≥ 0, got {min_padding}")
    if target_pt < min_pt:
        raise ValueError(f"target_pt {target_pt} < min_pt {min_pt}")
    if target_pt > max_pt:
        raise ValueError(f"target_pt {target_pt} > max_pt {max_pt}")

    warnings_list: list[str] = []
    dist = _check_label_distinctiveness(list(labels))
    if dist:
        warnings_list.append(dist)

    if len(labels) == 0:
        return PanelSizing(
            fontsize_pt=target_pt, rotation_deg=0, horizontal_alignment="center",
            recommend_horizontal_bars=False, achieved_padding_ratio=1.0,
            strategy="default", warnings=tuple(warnings_list),
        )

    axis_dim_inches = figsize_inches[0] if axis == "x" else figsize_inches[1]
    axis_dim_png_px = axis_dim_inches * dpi * reserved_axis_fraction
    per_label_budget_px = axis_dim_png_px / len(labels)

    fit_kwargs = dict(
        font_family=font_family, dpi=dpi,
        per_label_budget_px=per_label_budget_px, min_padding=min_padding,
        axis=axis, tick_label_clearance_px=tick_label_clearance_px,
    )
    bsearch_kwargs = dict(
        font_family=font_family, dpi=dpi,
        per_label_budget_px=per_label_budget_px, min_padding=min_padding,
        axis=axis, tick_label_clearance_px=tick_label_clearance_px,
    )

    # Strategy 1: default at target_pt, rotation 0
    fits, pad, _ = _try_fit(labels, target_pt, 0, **fit_kwargs)
    if fits:
        return PanelSizing(
            fontsize_pt=target_pt, rotation_deg=0, horizontal_alignment="center",
            recommend_horizontal_bars=False, achieved_padding_ratio=pad,
            strategy="default", warnings=tuple(warnings_list),
        )

    # Strategy 2: shrink at rotation 0
    result = _binary_search_fontsize(labels, 0, min_pt, target_pt, **bsearch_kwargs)
    if result:
        pt, pad = result
        return PanelSizing(
            fontsize_pt=pt, rotation_deg=0, horizontal_alignment="center",
            recommend_horizontal_bars=False, achieved_padding_ratio=pad,
            strategy="shrunk", warnings=tuple(warnings_list),
        )

    if not allow_rotation:
        _, pad, _ = _try_fit(labels, min_pt, 0, **{**fit_kwargs, "min_padding": 0.0})
        warnings_list.append(
            f"min_padding {min_padding} unmet at min_pt {min_pt}; rotation disallowed"
        )
        return PanelSizing(
            fontsize_pt=min_pt, rotation_deg=0, horizontal_alignment="center",
            recommend_horizontal_bars=False, achieved_padding_ratio=pad,
            strategy="failed", warnings=tuple(warnings_list),
        )

    # Strategy 3: rotation 30
    fits, pad, _ = _try_fit(labels, target_pt, 30, **fit_kwargs)
    if fits:
        return PanelSizing(
            fontsize_pt=target_pt, rotation_deg=30, horizontal_alignment="right",
            recommend_horizontal_bars=False, achieved_padding_ratio=pad,
            strategy="rotated", warnings=tuple(warnings_list),
        )
    result = _binary_search_fontsize(labels, 30, min_pt, target_pt, **bsearch_kwargs)
    if result:
        pt, pad = result
        return PanelSizing(
            fontsize_pt=pt, rotation_deg=30, horizontal_alignment="right",
            recommend_horizontal_bars=False, achieved_padding_ratio=pad,
            strategy="rotated", warnings=tuple(warnings_list),
        )

    # Strategy 4: rotation 90
    fits, pad, _ = _try_fit(labels, target_pt, 90, **fit_kwargs)
    if fits:
        return PanelSizing(
            fontsize_pt=target_pt, rotation_deg=90, horizontal_alignment="center",
            recommend_horizontal_bars=False, achieved_padding_ratio=pad,
            strategy="rotated", warnings=tuple(warnings_list),
        )
    result = _binary_search_fontsize(labels, 90, min_pt, target_pt, **bsearch_kwargs)
    if result:
        pt, pad = result
        return PanelSizing(
            fontsize_pt=pt, rotation_deg=90, horizontal_alignment="center",
            recommend_horizontal_bars=False, achieved_padding_ratio=pad,
            strategy="rotated", warnings=tuple(warnings_list),
        )

    # Strategy 5: horizontal bar flip (gated)
    if (axis == "x" and allow_horizontal_bar_flip
            and axis_priority in ("x", "either")):
        flipped = size_axis_labels(
            body_w_px=body_w_px, body_h_px=body_h_px,
            figsize_inches=figsize_inches, dpi=dpi,
            labels=labels, axis="y", axis_priority=axis_priority,
            max_pt=max_pt, target_pt=target_pt, min_pt=min_pt,
            min_padding=min_padding, allow_rotation=allow_rotation,
            allow_horizontal_bar_flip=False,
            font_family=font_family,
            reserved_axis_fraction=reserved_axis_fraction,
            tick_label_clearance_px=tick_label_clearance_px,
        )
        return PanelSizing(
            fontsize_pt=flipped.fontsize_pt,
            rotation_deg=flipped.rotation_deg,
            horizontal_alignment=flipped.horizontal_alignment,
            recommend_horizontal_bars=True,
            achieved_padding_ratio=flipped.achieved_padding_ratio,
            strategy="flipped",
            warnings=tuple(warnings_list) + flipped.warnings,
        )

    # Strategy 6: failed — best-effort
    _, pad, _ = _try_fit(labels, min_pt, 90, **{**fit_kwargs, "min_padding": 0.0})
    warnings_list.append(
        f"min_padding {min_padding} unmet at all strategies; "
        f"using min_pt={min_pt} rotation=90"
    )
    return PanelSizing(
        fontsize_pt=min_pt, rotation_deg=90, horizontal_alignment="center",
        recommend_horizontal_bars=False, achieved_padding_ratio=pad,
        strategy="failed", warnings=tuple(warnings_list),
    )


# ============================================================================
# Phase 3 — overlap detection + auto-placement helpers
# ============================================================================


def _data_bboxes(ax: plt.Axes, renderer) -> list[tuple[str, object]]:
    """Collect (kind_label, Bbox) for each plotted data primitive on ``ax``.

    Walks Line2D markers + segments, patch artists (bars, fill regions),
    and PathCollection offsets (scatter). Used by overlap detectors to
    flag legends/text that cover data.
    """
    from matplotlib.transforms import Bbox

    out: list[tuple[str, object]] = []

    for line in ax.get_lines():
        try:
            xy = line.get_xydata()
        except Exception:
            continue
        if len(xy) == 0:
            continue
        try:
            disp = ax.transData.transform(xy)
        except Exception:
            continue
        # Marker bboxes only — segment AABBs are over-conservative on steep
        # slopes (a thin diagonal line generates a tall axis-aligned box that
        # claims the empty space above the line). The actual *data* is at
        # the markers; covering the connector between two markers is OK.
        ms_pt = (line.get_markersize() or 6) * 0.5
        ms_px = max(ms_pt * line.figure.dpi / 72.0, 3.0)
        for x, y in disp:
            out.append(("Line2D", Bbox.from_extents(
                x - ms_px, y - ms_px, x + ms_px, y + ms_px,
            )))

    for patch in ax.patches:
        try:
            bb = patch.get_window_extent(renderer)
        except Exception:
            continue
        if bb.width > 0 and bb.height > 0:
            out.append(("Patch", bb))

    for coll in ax.collections:
        try:
            offsets = coll.get_offsets()
        except Exception:
            continue
        if offsets is None or len(offsets) == 0:
            continue
        try:
            disp = coll.get_transform().transform(offsets)
        except Exception:
            continue
        for x, y in disp:
            out.append(("PathCollection", Bbox.from_extents(
                x - 4, y - 4, x + 4, y + 4,
            )))

    return out


def add_headroom(
    ax: plt.Axes,
    *,
    log_decades: float = 0.5,
    linear_fraction: float = 0.25,
    shared_panels: Sequence[plt.Axes] | None = None,
    force: bool = False,
) -> None:
    """Extend ``ax``'s upper y-limit to make room above data.

    The first call captures the original ``ylim`` as ``_original_ylim_top`` /
    ``_original_ylim_bot`` on the axes, so subsequent calls compute new_top
    relative to the **original** range — passing larger ``log_decades`` later
    monotonically grows the headroom, it does not compound. Step IDs gate
    redundant work; ``force=True`` recomputes anyway.

    For shared-axis scaffolds, pass every panel in ``shared_panels`` so all
    panels get the same new top.
    """
    targets = [ax]
    if shared_panels:
        targets.extend(shared_panels)

    for a in targets:
        if not hasattr(a, "_original_ylim_top"):
            bot, top = a.get_ylim()
            a._original_ylim_top = top
            a._original_ylim_bot = bot

    is_log = ax.get_yscale() == "log"
    step_id = f"log:{log_decades}" if is_log else f"lin:{linear_fraction}"
    if getattr(ax, "_headroom_applied", None) == step_id and not force:
        return

    base_top = ax._original_ylim_top
    base_bot = ax._original_ylim_bot
    if is_log:
        new_top = base_top * (10 ** log_decades)
    else:
        new_top = base_top + (base_top - base_bot) * linear_fraction

    for a in targets:
        a.set_ylim(top=new_top)
        a._headroom_applied = step_id


_DEFAULT_LEGEND_CANDIDATES: tuple[str, ...] = (
    "upper left", "upper center", "upper right",
    "center right", "lower right", "lower center",
    "lower left", "center left",
)

# (name, (x_axes, y_axes), va, ha) — clockwise from upper right.
_DEFAULT_CORNER_ANCHORS: tuple[tuple[str, tuple[float, float], str, str], ...] = (
    ("upper right",  (0.96, 0.96), "top",    "right"),
    ("lower right",  (0.96, 0.04), "bottom", "right"),
    ("lower left",   (0.04, 0.04), "bottom", "left"),
    ("upper left",   (0.04, 0.96), "top",    "left"),
    ("upper center", (0.50, 0.96), "top",    "center"),
    ("center right", (0.96, 0.50), "center", "right"),
    ("lower center", (0.50, 0.04), "bottom", "center"),
    ("center left",  (0.04, 0.50), "center", "left"),
)


_HEADROOM_LADDER: tuple[tuple[dict | None, str], ...] = (
    (None, ""),
    ({"log_decades": 0.5, "linear_fraction": 0.25}, "half_decade"),
    ({"log_decades": 1.0, "linear_fraction": 0.50}, "full_decade"),
)


def _overlap_area(bb_a, bb_b) -> float:
    """Return overlap area between two Bbox objects (0 if disjoint)."""
    if not bb_a.overlaps(bb_b):
        return 0.0
    x_overlap = min(bb_a.x1, bb_b.x1) - max(bb_a.x0, bb_b.x0)
    y_overlap = min(bb_a.y1, bb_b.y1) - max(bb_a.y0, bb_b.y0)
    return max(0.0, x_overlap) * max(0.0, y_overlap)


# Ordered from least to most demanding (exhausted = full_decade applied, still overlapping).
_HEADROOM_ORDER: tuple[str, ...] = ("", "half_decade", "full_decade", "exhausted")


def max_headroom_label(labels) -> str:
    """Return the most demanding headroom label from a collection.

    Use this in the two-pass shared-axis workflow to reconcile the per-panel
    probe results into a single unified headroom level for the whole grid::

        labels = [probe_legend_headroom(ax_a, **kw), probe_legend_headroom(ax_b, **kw)]
        needed = max_headroom_label(labels)
        shared_top = headroom_ylim_top(base_top, needed, yscale='log')
    """
    best = ""
    for label in labels:
        if label not in _HEADROOM_ORDER:
            raise ValueError(f"Unknown headroom label: {label!r}")
        if _HEADROOM_ORDER.index(label) > _HEADROOM_ORDER.index(best):
            best = label
    return best


def headroom_ylim_top(current_top: float, label: str, *, yscale: str = "log") -> float:
    """Convert a headroom label to a new numeric ylim top value.

    Converts the string returned by :func:`probe_legend_headroom` into a
    multiplier, letting the caller build a shared ylim for a grid of panels
    *before* rendering any of them. This is the key step in the two-pass
    shared-axis workflow::

        new_top = headroom_ylim_top(base_ylim[1], needed_label, yscale='log')
        shared_ylim = (base_ylim[0], new_top)
        # then render all panels with shared_ylim and allow_headroom=False

    ``'exhausted'`` is treated as ``'full_decade'`` — it means full-decade
    headroom was applied but overlap remained; the caller should still extend
    the ylim and accept the residual overlap.
    """
    _log_factors = {"": 1.0, "half_decade": 10 ** 0.5, "full_decade": 10.0, "exhausted": 10.0}
    _lin_factors = {"": 1.0, "half_decade": 1.25, "full_decade": 1.50, "exhausted": 1.50}
    factors = _log_factors if yscale == "log" else _lin_factors
    if label not in factors:
        raise ValueError(f"Unknown headroom label: {label!r}")
    return current_top * factors[label]


def probe_legend_headroom(
    ax: plt.Axes,
    *,
    candidates: Sequence[str] | None = None,
    **legend_kwargs,
) -> str:
    """Detect which headroom level yields a clean legend placement.

    Runs the same candidate → headroom escalation ladder as
    :func:`auto_place_legend` but is *non-destructive*: it restores ``ylim``
    and removes the test legend when done.

    Returns the headroom label of the first clean placement:
    ``''`` (no headroom needed), ``'half_decade'``, ``'full_decade'``, or
    ``'exhausted'`` (full-decade applied, still overlapping).

    Typical two-pass shared-axis usage::

        # Pass 1 — probe each legend-bearing panel independently
        fig_a, ax_a = renderer.fig_ax()
        _plot_data(ax_a, ...)
        renderer.size_x_labels(ax_a, ...); renderer.size_y_labels(ax_a, ...)
        label_a = probe_legend_headroom(ax_a, fontsize=7, ncol=3)
        plt.close(fig_a)

        # Reconcile across all panels, compute shared ylim
        needed = max_headroom_label([label_a, ...])  # '' if no panel needs headroom
        shared_top = headroom_ylim_top(base_ylim[1], needed, yscale='log')

        # Pass 2 — render all panels with pre-extended ylim; disable re-escalation
        render_panel(ax, ylim=(base_ylim[0], shared_top), allow_headroom=False)
    """
    legend_kwargs.pop("loc", None)
    cands = list(candidates) if candidates is not None else list(_DEFAULT_LEGEND_CANDIDATES)

    orig_ylim = ax.get_ylim()
    orig_hd_applied = getattr(ax, "_headroom_applied", None)
    orig_hd_top = getattr(ax, "_original_ylim_top", None)
    orig_hd_bot = getattr(ax, "_original_ylim_bot", None)

    try:
        for headroom_kwargs, hd_label in _HEADROOM_LADDER:
            if headroom_kwargs is not None:
                add_headroom(ax, **headroom_kwargs)

            for loc in cands:
                old = ax.get_legend()
                if old is not None:
                    old.remove()
                legend = ax.legend(loc=loc, **legend_kwargs)
                ax.figure.canvas.draw()
                renderer = ax.figure.canvas.get_renderer()

                try:
                    leg_bb = legend.get_window_extent(renderer)
                except Exception:
                    continue

                data = _data_bboxes(ax, renderer)
                if sum(_overlap_area(leg_bb, bb) for _, bb in data) == 0.0:
                    return hd_label

        return "exhausted"
    finally:
        old = ax.get_legend()
        if old is not None:
            old.remove()
        ax.set_ylim(*orig_ylim)
        if orig_hd_applied is None:
            ax.__dict__.pop("_headroom_applied", None)
        else:
            ax._headroom_applied = orig_hd_applied
        if orig_hd_top is None:
            ax.__dict__.pop("_original_ylim_top", None)
        else:
            ax._original_ylim_top = orig_hd_top
        if orig_hd_bot is None:
            ax.__dict__.pop("_original_ylim_bot", None)
        else:
            ax._original_ylim_bot = orig_hd_bot


def auto_place_legend(
    ax: plt.Axes,
    *,
    candidates: Sequence[str] | None = None,
    shared_panels: Sequence[plt.Axes] | None = None,
    allow_headroom: bool = True,
    **legend_kwargs,
) -> PlacementResult:
    """Pick a legend ``loc`` that doesn't overlap plotted data.

    For each candidate ``loc``, measures total overlap area against data
    bboxes. Returns the first candidate with zero overlap. If none are
    clean at the current headroom level, escalates (half-decade →
    full-decade) and retries. If the full ladder is exhausted with no
    clean candidate, returns the position with the *minimum* overlap
    area (least-bad), with ``overlapped=True``.

    ``allow_headroom=False`` disables the headroom escalation steps and
    keeps ``ylim`` untouched — use this for shared-axis-locked figures
    where each panel must keep its caller-set ``ylim`` so the grid stays
    visually aligned.

    ``legend_kwargs`` are forwarded to ``ax.legend`` verbatim — only ``loc``
    is auto-resolved (passing ``loc=`` in kwargs is silently overridden).
    """
    legend_kwargs.pop("loc", None)
    cands = list(candidates) if candidates is not None else list(_DEFAULT_LEGEND_CANDIDATES)

    best_overall: tuple[float, str, str, tuple[str, ...]] | None = None  # (area, loc, hd_label, artists)

    ladder = _HEADROOM_LADDER if allow_headroom else _HEADROOM_LADDER[:1]
    for headroom_kwargs, hd_label in ladder:
        if headroom_kwargs is not None:
            add_headroom(ax, shared_panels=shared_panels, **headroom_kwargs)

        for loc in cands:
            old = ax.get_legend()
            if old is not None:
                old.remove()
            legend = ax.legend(loc=loc, **legend_kwargs)
            ax.figure.canvas.draw()
            renderer = ax.figure.canvas.get_renderer()

            try:
                leg_bb = legend.get_window_extent(renderer)
            except Exception:
                continue

            data = _data_bboxes(ax, renderer)
            total_area = sum(_overlap_area(leg_bb, bb) for _, bb in data)
            artists = tuple(sorted({
                k for k, bb in data if leg_bb.overlaps(bb)
            }))

            if total_area == 0.0:
                return PlacementResult(
                    loc=loc, headroom_applied=hd_label,
                    overlapped=False, overlap_artists=(),
                )
            if best_overall is None or total_area < best_overall[0]:
                best_overall = (total_area, loc, hd_label, artists)

    # All exhausted — restore the best-area candidate.
    assert best_overall is not None
    _, best_loc, best_hd, best_artists = best_overall
    old = ax.get_legend()
    if old is not None:
        old.remove()
    ax.legend(loc=best_loc, **legend_kwargs)
    return PlacementResult(
        loc=best_loc, headroom_applied=best_hd,
        overlapped=True, overlap_artists=best_artists,
    )


def place_corner_text(
    ax: plt.Axes,
    text: str,
    *,
    prefer: str = "upper right",
    shared_panels: Sequence[plt.Axes] | None = None,
    allow_headroom: bool = True,
    **text_kwargs,
) -> PlacementResult:
    """Place ``text`` at a corner of ``ax`` that doesn't overlap data.

    Tries ``prefer`` first, then the remaining 7 anchor positions clockwise.
    Anchor list = 4 corners + 4 edge midpoints in ``transAxes`` coords. If
    all overlap, escalates through the half/full-decade headroom ladder.

    ``allow_headroom=False`` disables the headroom escalation steps —
    use it for shared-axis figures where ``ylim`` is locked.

    ``text_kwargs`` are forwarded to ``ax.text`` (font, color, bbox, zorder…).
    Only position (``transform``, ``x``, ``y``, ``va``, ``ha``) is
    auto-resolved.
    """
    for k in ("transform", "x", "y", "va", "ha",
              "verticalalignment", "horizontalalignment"):
        text_kwargs.pop(k, None)

    by_name = {name: (anc, va, ha) for name, anc, va, ha in _DEFAULT_CORNER_ANCHORS}
    if prefer not in by_name:
        raise ValueError(
            f"Unknown corner: {prefer!r}. Choose from {list(by_name.keys())}"
        )

    ordered: list[tuple[str, tuple[float, float], str, str]] = []
    found = False
    for name, anc, va, ha in _DEFAULT_CORNER_ANCHORS:
        if name == prefer:
            found = True
        if found:
            ordered.append((name, anc, va, ha))
    for name, anc, va, ha in _DEFAULT_CORNER_ANCHORS:
        if name == prefer:
            break
        ordered.append((name, anc, va, ha))

    best_overall: tuple[float, tuple[float, float], str, str, str, tuple[str, ...]] | None = None
    # (area, anc, va, ha, hd_label, artists)

    ladder = _HEADROOM_LADDER if allow_headroom else _HEADROOM_LADDER[:1]
    for headroom_kwargs, hd_label in ladder:
        if headroom_kwargs is not None:
            add_headroom(ax, shared_panels=shared_panels, **headroom_kwargs)

        for name, (x_axes, y_axes), va, ha in ordered:
            existing = getattr(ax, "_auto_corner_text", None)
            if existing is not None and existing in ax.texts:
                existing.remove()
            ax._auto_corner_text = None

            t = ax.text(
                x_axes, y_axes, text,
                transform=ax.transAxes, va=va, ha=ha,
                **text_kwargs,
            )
            ax._auto_corner_text = t
            ax.figure.canvas.draw()
            renderer = ax.figure.canvas.get_renderer()

            try:
                t_bb = t.get_window_extent(renderer)
            except Exception:
                continue

            data = _data_bboxes(ax, renderer)
            total_area = sum(_overlap_area(t_bb, bb) for _, bb in data)
            artists = tuple(sorted({
                k for k, bb in data if t_bb.overlaps(bb)
            }))

            if total_area == 0.0:
                return PlacementResult(
                    loc=(x_axes, y_axes), headroom_applied=hd_label,
                    overlapped=False, overlap_artists=(),
                )
            if best_overall is None or total_area < best_overall[0]:
                best_overall = (total_area, (x_axes, y_axes), va, ha, hd_label, artists)

    # All exhausted — restore the best-area anchor.
    assert best_overall is not None
    _, best_anc, best_va, best_ha, best_hd, best_artists = best_overall
    existing = getattr(ax, "_auto_corner_text", None)
    if existing is not None and existing in ax.texts:
        existing.remove()
    t = ax.text(
        best_anc[0], best_anc[1], text,
        transform=ax.transAxes, va=best_va, ha=best_ha,
        **text_kwargs,
    )
    ax._auto_corner_text = t
    return PlacementResult(
        loc=best_anc, headroom_applied=best_hd,
        overlapped=True, overlap_artists=best_artists,
    )


# ============================================================================
# Phase 2 — PanelRenderer class wrapper
# ============================================================================


_DEFAULT_MARGINS = {"left": 0.16, "right": 0.97, "top": 0.95, "bottom": 0.13}


class PanelRenderer:
    """Bundles figure setup, sizing, collision detection, and PNG export.

    Replaces the per-script `_render_clean_panel(plot_fn, ..., margins=...)`
    pattern in build_retrofit_*.py. Caller plots data on the returned axes,
    sizes labels via `size_x_labels()` / `size_y_labels()`, optionally checks
    for collisions, then saves with `save_as_data_uri()`.

    Usage:
        r = PanelRenderer(
            body_w_px=429, body_h_px=298,
            figsize_inches=(2.86, 1.99),
            margins={"left": 0.16, "right": 0.97, "top": 0.95, "bottom": 0.13},
            axis_priority="y",
        )
        fig, ax = r.fig_ax()
        ax.bar(x, y, color="...")
        ax.set_xticks(x)
        r.size_x_labels(ax, labels=["1M", "310K", "98K", "31K"])
        warnings = r.check_collisions(ax)
        uri = r.save_as_data_uri(fig)
    """

    def __init__(
        self,
        *,
        body_w_px: int,
        body_h_px: int,
        figsize_inches: tuple[float, float],
        dpi: int = 200,
        margins: dict | None = None,
        axis_priority: Literal["x", "y", "either"] = "y",
        font_family: str = "DejaVu Sans",
    ):
        self.body_w_px = body_w_px
        self.body_h_px = body_h_px
        self.figsize_inches = figsize_inches
        self.dpi = dpi
        self.margins = dict(_DEFAULT_MARGINS) if margins is None else dict(margins)
        self.axis_priority = axis_priority
        self.font_family = font_family

    def fig_ax(self) -> tuple[plt.Figure, plt.Axes]:
        """Create a fresh figure + axes at the configured figsize/dpi."""
        fig, ax = plt.subplots(figsize=self.figsize_inches, dpi=self.dpi)
        return fig, ax

    def _clearance(self, axis: Literal["x", "y"]) -> float:
        """Compute tick-label clearance from this renderer's margins."""
        if axis == "x":
            return clearance_px_from_subplots(
                figsize_inches=self.figsize_inches, dpi=self.dpi,
                axis="x", subplots_bottom=self.margins.get("bottom", 0.13),
            )
        return clearance_px_from_subplots(
            figsize_inches=self.figsize_inches, dpi=self.dpi,
            axis="y", subplots_left=self.margins.get("left", 0.16),
        )

    def size_x_labels(
        self,
        ax: plt.Axes,
        *,
        labels: Sequence[str],
        **sizing_kwargs,
    ) -> PanelSizing:
        """Compute sizing for x-axis tick labels and apply to `ax`.

        Caller is responsible for setting tick positions via `ax.set_xticks()`
        before calling this. `sizing_kwargs` are forwarded to size_axis_labels
        (e.g. `min_padding=0.15` to relax breathing room).
        """
        sizing = size_axis_labels(
            body_w_px=self.body_w_px, body_h_px=self.body_h_px,
            figsize_inches=self.figsize_inches, dpi=self.dpi,
            labels=labels, axis="x", axis_priority=self.axis_priority,
            tick_label_clearance_px=self._clearance("x"),
            font_family=self.font_family,
            **sizing_kwargs,
        )
        kwargs: dict = dict(
            fontsize=sizing.fontsize_pt,
            rotation=sizing.rotation_deg,
            ha=sizing.horizontal_alignment,
        )
        if sizing.rotation_deg != 0:
            kwargs["rotation_mode"] = "anchor"
        ax.set_xticklabels(labels, **kwargs)
        return sizing

    def size_y_labels(
        self,
        ax: plt.Axes,
        *,
        labels: Sequence[str],
        **sizing_kwargs,
    ) -> PanelSizing:
        """Compute sizing for y-axis tick labels and apply to `ax`."""
        sizing = size_axis_labels(
            body_w_px=self.body_w_px, body_h_px=self.body_h_px,
            figsize_inches=self.figsize_inches, dpi=self.dpi,
            labels=labels, axis="y", axis_priority=self.axis_priority,
            tick_label_clearance_px=self._clearance("y"),
            font_family=self.font_family,
            **sizing_kwargs,
        )
        ax.set_yticklabels(labels, fontsize=sizing.fontsize_pt)
        return sizing

    def check_collisions(
        self,
        ax: plt.Axes,
        *,
        include_data: bool = True,
    ) -> list[str]:
        """Detect bbox overlaps on ``ax``.

        Always compares tick labels, legend, and ``ax.texts`` against each
        other. With ``include_data=True`` (default), also compares legend
        and ``ax.texts`` against plotted data primitives (Line2D markers/
        segments, patches, scatter offsets) — catches the legend-covers-
        line and corner-text-covers-datapoint cases.

        Returns human-readable collision descriptions. Empty list means
        no overlaps detected.
        """
        fig = ax.figure
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()

        items: list[tuple[str, str, object]] = []

        for lbl in ax.get_xticklabels():
            if lbl.get_text():
                items.append(("xtick", lbl.get_text(), lbl.get_window_extent(renderer)))
        for lbl in ax.get_yticklabels():
            if lbl.get_text():
                items.append(("ytick", lbl.get_text(), lbl.get_window_extent(renderer)))

        legend = ax.get_legend()
        if legend is not None:
            try:
                items.append(("legend", "", legend.get_window_extent(renderer)))
            except Exception:
                pass

        for txt in ax.texts:
            if txt.get_text():
                items.append(("text", txt.get_text(), txt.get_window_extent(renderer)))

        collisions = []
        for i, (kind_a, text_a, bb_a) in enumerate(items):
            for kind_b, text_b, bb_b in items[i + 1:]:
                if bb_a.overlaps(bb_b):
                    a_repr = f"{kind_a}({text_a!r})" if text_a else kind_a
                    b_repr = f"{kind_b}({text_b!r})" if text_b else kind_b
                    collisions.append(f"{a_repr} overlaps {b_repr}")

        if include_data:
            data = _data_bboxes(ax, renderer)
            seen_pairs: set[tuple[str, str]] = set()
            for kind, text, bb in items:
                if kind not in ("legend", "text"):
                    continue
                for data_kind, data_bb in data:
                    if not bb.overlaps(data_bb):
                        continue
                    a_repr = f"{kind}({text!r})" if text else kind
                    b_repr = f"data({data_kind})"
                    pair = (a_repr, b_repr)
                    if pair in seen_pairs:
                        continue
                    seen_pairs.add(pair)
                    collisions.append(f"{a_repr} overlaps {b_repr}")

        return collisions

    def auto_place_legend(
        self,
        ax: plt.Axes,
        *,
        candidates: Sequence[str] | None = None,
        shared_panels: Sequence[plt.Axes] | None = None,
        allow_headroom: bool = True,
        **legend_kwargs,
    ) -> PlacementResult:
        """Pick a legend ``loc`` that doesn't overlap data. See
        :func:`auto_place_legend` for the full algorithm."""
        return auto_place_legend(
            ax, candidates=candidates, shared_panels=shared_panels,
            allow_headroom=allow_headroom, **legend_kwargs,
        )

    def place_corner_text(
        self,
        ax: plt.Axes,
        text: str,
        *,
        prefer: str = "upper right",
        shared_panels: Sequence[plt.Axes] | None = None,
        allow_headroom: bool = True,
        **text_kwargs,
    ) -> PlacementResult:
        """Place ``text`` at a non-overlapping corner. See
        :func:`place_corner_text` for the full algorithm."""
        return place_corner_text(
            ax, text, prefer=prefer, shared_panels=shared_panels,
            allow_headroom=allow_headroom, **text_kwargs,
        )

    def add_headroom(
        self,
        ax: plt.Axes,
        *,
        log_decades: float = 0.5,
        linear_fraction: float = 0.25,
        shared_panels: Sequence[plt.Axes] | None = None,
        force: bool = False,
    ) -> None:
        """Extend ``ax``'s upper y-limit. See :func:`add_headroom`."""
        add_headroom(
            ax,
            log_decades=log_decades, linear_fraction=linear_fraction,
            shared_panels=shared_panels, force=force,
        )

    def probe_legend_headroom(
        self,
        ax: plt.Axes,
        *,
        candidates: Sequence[str] | None = None,
        **legend_kwargs,
    ) -> str:
        """Non-destructive headroom probe for two-pass shared-axis grids.
        See :func:`probe_legend_headroom` for full docs."""
        return probe_legend_headroom(ax, candidates=candidates, **legend_kwargs)

    def save_as_data_uri(
        self,
        fig: plt.Figure,
        path: Path | None = None,
    ) -> str:
        """Apply margins, save figure, return base64 data URI.

        If `path` is given, also writes the PNG to disk (useful for caching
        and inspection). Always closes the figure after save.
        """
        fig.subplots_adjust(**self.margins)
        if path is not None:
            fig.savefig(path, dpi=self.dpi)
            data = path.read_bytes()
        else:
            buf = io.BytesIO()
            fig.savefig(buf, format='png', dpi=self.dpi)
            data = buf.getvalue()
        plt.close(fig)
        return f"data:image/png;base64,{base64.b64encode(data).decode('ascii')}"
