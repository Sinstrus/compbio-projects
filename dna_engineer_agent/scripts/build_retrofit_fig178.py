"""FIG178 retrofit — wrap EXP26000387 absolute-RLU 2×2 onto SCF024.

This is a *shared-y-axis* figure: the visual story (cross-panel magnitude
comparison from 1K to 1M RLU) depends on every panel using the SAME log y
range. Independent matplotlib panels rendered into 4 separate body slots
will only LOOK shared if we coordinate them precisely:

  1. **Identical ylim**: every panel forces (1e3, 2e6) — same min/max.
  2. **Identical yscale**: every panel forces 'log'.
  3. **Identical y-ticks + labels**: 1K / 10K / 100K / 1M, formatted the
     same way, so the gridline pixel locations are identical inside each
     PNG.
  4. **Identical subplots_adjust margins**: same left/top/right/bottom
     across all 4 PNGs, so the plot area sits at the same pixel offsets
     within each body slot. This is what makes the eye perceive a single
     shared y-axis across the 2×2 grid.

Trade-off: each panel renders its OWN y-tick labels (1K/10K/100K/1M show
on all 4, including right-column panels). Acceptable cost for the visual
clarity of full-height panels.

This figure uses **SCF024** (tall variant) rather than SCF016 because
SCF016's 142-px-tall body compressed the 3-decade log Y axis to ~47
px/decade, halving the original ~100 px/decade and visually merging
lines that should be distinguishable. SCF024's 960×800 viewBox gives
each panel a ~298-px-tall body — restoring the resolution.

**Sizing**: x/y tick label fontsize and rotation are now picked
programmatically by `figure_bank.matplotlib_panel.PanelRenderer`. One
renderer instance owns the shared config (figsize/margins/dpi) for all
4 panels — this is the shared-axis lock encoded as data, not as four
copies of the same hand-tune.

Layout:
  Panel A (top-left)     ALPL+ Lysate, 5×8 fraction lines
  Panel B (top-right)    ALPL+ VOY-1631 dose-response bars
  Panel C (bottom-left)  HEK   Lysate, 5×8 fraction lines
  Panel D (bottom-right) HEK   VOY-1631 dose-response bars

Source: EXP26000387 NOTES.md "Run 4 — Full Fraction Walk" plate tables.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Sequence

import numpy as np
import matplotlib
matplotlib.use('Agg')

sys.path.insert(0, str(Path.home() / "projects"))

import matplotlib.pyplot as plt

from figure_bank.matplotlib_panel import (
    PanelRenderer,
    probe_legend_headroom,
    headroom_ylim_top,
)
from figure_bank.scaffolds import SCF024_plot_grid_4panel_tall as scf024


FIGURE_BANK = Path.home() / "projects" / "figure_bank"
PANELS_DIR = FIGURE_BANK / "_retrofit_panels"
PANELS_DIR.mkdir(exist_ok=True)


# === Data: EXP26000387 Run 4 (2026-04-22, AutoGain 204) ======================
# Plate 1 = HEK293T (ALPL−); Plate 2 = HEK293T-ALPL Rivera (ALPL+).
# Cols 1, 3, 5, 7, 9 = neat conditions; Rows A-H = F1-F8 fractions.

_FRACTIONS = ("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8")

# Plate 2 — ALPL+ (HEK293T-ALPL Rivera), neat readings only.
_ALPL_LYS = {
    "WT":          [37307, 46679, 46300, 1045031, 695920, 273526, 56288, 20996],
    "VHH3-Mos":    [ 1373,  3167,  2187,   17811,  51085,  41249, 20397,  5610],
    "VHH3-Tc":     [ 3453,  3771,  5041,   87078, 117757,  51515, 25352,  4357],
    "VHH3-VP2n":   [13195, 18417, 21903,  421217, 450159, 186839, 28490, 15094],
    "VHH3-IR":     [ 2955,  6412,  5025,   87717, 116886,  80095, 25077,  8599],
}

# Plate 1 — ALPL− (HEK293T), neat readings only.
_HEK_LYS = {
    "WT":          [60426, 76301, 71656, 1274152, 991885, 421604, 137664, 56833],
    "VHH3-Mos":    [ 2082,  3433,  2925,   34128,  69805,  43407,  39561, 10322],
    "VHH3-Tc":     [ 4353,  8299,  6974,  115145, 142671,  81535,  35325,  9332],
    "VHH3-VP2n":   [26051, 32521, 30400,  585180, 489019, 206278,  58034, 40470],
    "VHH3-IR":     [ 6805, 14308,  8765,  155500, 156871, 171450,  73776, 30329],
}

_COLORS = {
    "WT":          "#888888",
    "VHH3-Mos":    "#1f77b4",
    "VHH3-Tc":     "#2ca02c",
    "VHH3-VP2n":   "#d62728",
    "VHH3-IR":     "#9467bd",
}

# VOY-1631 (col 11). Doses A=1.0e6, B=3.1e5, C=9.8e4, D=3.1e4 VG/cell.
# Short labels so the bar tick row fits cleanly without rotation.
_VOY_DOSE_LABELS = ("1M", "310K", "98K", "31K")
_VOY_ALPL = (493314, 469416, 339036, 217838)
_VOY_HEK  = (144112,  77472,  26598,  10158)


# === Shared-axis lock =======================================================

_YLIM = (1e3, 2e6)
_YTICKS = [1e3, 1e4, 1e5, 1e6]
_YLABELS = ["1K", "10K", "100K", "1M"]
# Identical margins across all 4 panels so plot areas pixel-align.
_PANEL_MARGINS = {"left": 0.16, "right": 0.97, "top": 0.95, "bottom": 0.13}


def _apply_shared_y(ax) -> None:
    """The lock — call last in each panel plot fn.

    Sets y-scale, ylim, and yticks (positions). Y-tick *labels* are applied
    later by `renderer.size_y_labels()` so the label fontsize is picked by
    the sizer, not hand-coded here.
    """
    ax.set_yscale('log')
    ax.set_ylim(*_YLIM)
    ax.set_yticks(_YTICKS)
    ax.grid(True, which='major', axis='y', alpha=0.3, linewidth=0.5)
    ax.grid(True, which='minor', axis='y', alpha=0.12, linewidth=0.3)
    ax.set_axisbelow(True)


# === Per-panel plot fns =====================================================

def _plot_lysate_lines(ax, data: dict[str, list[float]]) -> None:
    """Line chart: 5 conditions × 8 fractions.

    Plots data + sets x positions only; tick labels (and their fontsize) are
    applied by the renderer in build_fig178(). Legends, when shown, are
    placed by `renderer.auto_place_legend(...)` after this returns.
    """
    x = np.arange(len(_FRACTIONS))
    for cond, ys in data.items():
        ys_safe = [max(y, 1.0) for y in ys]
        ax.plot(x, ys_safe, marker='o', linewidth=1.2, markersize=3.5,
                color=_COLORS[cond], label=cond)
    ax.set_xticks(x)
    ax.set_xlim(-0.4, 7.4)
    _apply_shared_y(ax)


def _plot_voy_bars(ax, ratios: tuple[float, ...]) -> None:
    """Bar chart: 4 dose VOY-1631 readings.

    Plots data + sets x positions only; tick labels applied by renderer.
    Value annotations above bars stay hand-tuned (free-form text, not ticks).
    """
    x = np.arange(len(_VOY_DOSE_LABELS))
    ax.bar(x, ratios, color="#2e7d32", edgecolor="#1b5e20",
           linewidth=0.8, width=0.65)
    for xi, yi in zip(x, ratios):
        if yi >= 1e3:
            kval = yi / 1e3
            label = f"{kval:.0f}K" if kval >= 100 else f"{kval:.1f}K"
            ax.text(xi, yi * 1.4, label,
                    ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax.set_xticks(x)
    _apply_shared_y(ax)


# === Render helper ==========================================================

def _render_panel(
    renderer: PanelRenderer,
    plot_fn,
    x_labels: Sequence[str],
    out_path: Path,
    show_legend: bool = False,
    ylim_override: tuple[float, float] | None = None,
) -> str:
    """Render one panel via the shared renderer; size labels + place legend.

    NOT using bbox_inches='tight' — that would crop differently per panel
    (panel A has a legend, panel B's bars don't) and shift the pixel
    location of the plot area, breaking visual alignment across the grid.
    The renderer's locked margins (set on construction) keep all 4 panels
    pixel-aligned.

    ``ylim_override`` is applied after plot_fn() so the two-pass shared-axis
    workflow can supply a pre-computed unified ylim (e.g. extended by half a
    decade to create legend headroom) without desyncing the grid.  When set,
    ``allow_headroom=False`` is still passed to auto_place_legend because the
    caller has already baked in the needed headroom.
    """
    fig, ax = renderer.fig_ax()
    plot_fn(ax)
    if ylim_override is not None:
        ax.set_ylim(*ylim_override)
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylabel('')

    x_sizing = renderer.size_x_labels(ax, labels=list(x_labels))
    y_sizing = renderer.size_y_labels(ax, labels=_YLABELS)

    legend_result = None
    if show_legend:
        legend_result = renderer.auto_place_legend(
            ax,
            allow_headroom=False,   # ylim pre-extended by two-pass probe; no further growth
            fontsize=7, ncol=3, frameon=True,
            framealpha=0.9, edgecolor='#bbb',
            handlelength=1.0, columnspacing=0.6,
            handletextpad=0.3, borderpad=0.3,
        )

    collisions = renderer.check_collisions(ax)

    if (x_sizing.warnings or y_sizing.warnings or collisions
            or (legend_result and legend_result.overlapped)):
        print(f"  ! {out_path.name}")
        for w in x_sizing.warnings:
            print(f"      x-sizing: {w}")
        for w in y_sizing.warnings:
            print(f"      y-sizing: {w}")
        for c in collisions:
            print(f"      collision: {c}")
        if legend_result and legend_result.overlapped:
            print(f"      legend STILL overlaps {legend_result.overlap_artists} "
                  f"after compass search at {legend_result.loc!r}")
    elif legend_result is not None:
        print(f"  OK    {out_path.name}  legend@{legend_result.loc}")

    return renderer.save_as_data_uri(fig, out_path)


# === Main ===================================================================

_LEGEND_KWARGS = dict(
    fontsize=7, ncol=3, frameon=True,
    framealpha=0.9, edgecolor='#bbb',
    handlelength=1.0, columnspacing=0.6,
    handletextpad=0.3, borderpad=0.3,
)


def _probe_shared_ylim(renderer: PanelRenderer) -> tuple[float, float]:
    """Pass 1 — detect legend headroom needed for panel A, return unified ylim.

    Panel A (ALPL+ lysate) is the only panel with a legend. We plot it in a
    temporary figure, probe which headroom level yields a clean placement,
    then compute a new shared top that is applied to ALL 4 panels so the
    grid stays visually aligned.
    """
    fig, ax = renderer.fig_ax()
    _plot_lysate_lines(ax, _ALPL_LYS)
    renderer.size_x_labels(ax, labels=list(_FRACTIONS))
    renderer.size_y_labels(ax, labels=_YLABELS)
    needed = probe_legend_headroom(ax, **_LEGEND_KWARGS)
    plt.close(fig)
    new_top = headroom_ylim_top(_YLIM[1], needed, yscale="log")
    shared_ylim = (_YLIM[0], new_top)
    print(f"  Probe → panel A legend needs {needed!r} headroom "
          f"→ shared ylim top {new_top:.3g} (was {_YLIM[1]:.3g})")
    return shared_ylim


def build_fig178() -> None:
    base = "FIG178-exp26000387-attempt4-akaluc-absolute-rlu"

    # SINGLE renderer config = the shared-axis lock. All 4 panels share figsize,
    # margins, dpi, axis_priority — so the plot area sits at the same pixel
    # offsets within each body slot.
    renderer = PanelRenderer(
        body_w_px=429, body_h_px=298,           # SCF024 body
        figsize_inches=(2.86, 1.99),            # matches body aspect 1.44
        dpi=200,
        margins=_PANEL_MARGINS,
        axis_priority="y",                      # log Y is the story
    )

    # Two-pass: probe panel A headroom first, extend shared ylim uniformly.
    shared_ylim = _probe_shared_ylim(renderer)

    a_uri = _render_panel(
        renderer,
        lambda ax: _plot_lysate_lines(ax, _ALPL_LYS),
        x_labels=_FRACTIONS,
        out_path=PANELS_DIR / f"{base}-pA.png",
        show_legend=True,
        ylim_override=shared_ylim,
    )
    b_uri = _render_panel(
        renderer,
        lambda ax: _plot_voy_bars(ax, _VOY_ALPL),
        x_labels=_VOY_DOSE_LABELS,
        out_path=PANELS_DIR / f"{base}-pB.png",
        ylim_override=shared_ylim,
    )
    c_uri = _render_panel(
        renderer,
        lambda ax: _plot_lysate_lines(ax, _HEK_LYS),
        x_labels=_FRACTIONS,
        out_path=PANELS_DIR / f"{base}-pC.png",
        ylim_override=shared_ylim,
    )
    d_uri = _render_panel(
        renderer,
        lambda ax: _plot_voy_bars(ax, _VOY_HEK),
        x_labels=_VOY_DOSE_LABELS,
        out_path=PANELS_DIR / f"{base}-pD.png",
        ylim_override=shared_ylim,
    )

    panels = [
        scf024.Panel(
            panel_id="A",
            panel_title="ALPL+ — Lysate fractions (5×8)",
            x_axis_label="Iodixanol fraction",
            y_axis_label="AkaLuc RLU (log)",
            image_href=a_uri,
        ),
        scf024.Panel(
            panel_id="B",
            panel_title="ALPL+ — VOY-1631 dose-response",
            x_axis_label="VG/cell (dose)",
            y_axis_label="AkaLuc RLU (log)",
            image_href=b_uri,
        ),
        scf024.Panel(
            panel_id="C",
            panel_title="ALPL− — Lysate fractions (5×8)",
            x_axis_label="Iodixanol fraction",
            y_axis_label="AkaLuc RLU (log)",
            image_href=c_uri,
        ),
        scf024.Panel(
            panel_id="D",
            panel_title="ALPL− — VOY-1631 dose-response",
            x_axis_label="VG/cell (dose)",
            y_axis_label="AkaLuc RLU (log)",
            image_href=d_uri,
        ),
    ]
    svg = scf024.render(
        title="EXP26000387 — Run 4 AkaLuc Absolute RLU (shared log Y, 1K-1M)",
        subtitle="Top: ALPL+ HEK293T-Rivera. Bottom: ALPL− HEK293T. Single replicate per condition.",
        panels=panels,
    )
    out = FIGURE_BANK / f"{base}.svg"
    out.write_text(svg)
    print(f"  OK    {out.name}")


if __name__ == "__main__":
    build_fig178()
