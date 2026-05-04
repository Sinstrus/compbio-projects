"""FIG168-173 — VR4 Splice Donor Screen: Scatter Plots.

Six scatter figures (one per combination of panel × feature):
  FIG168: FP VP3% vs H-bond count
  FIG169: HEK VP3% vs H-bond count
  FIG170: FP VP3% vs ΔG (x-axis inverted; more stable →)
  FIG171: HEK VP3% vs ΔG (x-axis inverted)
  FIG172: FP VP3% vs MaxEntScan score
  FIG173: HEK VP3% vs MaxEntScan score

All use SCF013 (1-panel, 890×316 body). Color = GTAAG presence (red/blue).
Marker shape = exonic_ctx (o/s/^/D). Sigmoid baseline for non-GTAAG CAG/AAG.

Data: experiments/EXP26000390-vr4-splicing-screen-rtpcr/donor_feature_table.csv
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import NamedTuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd

sys.path.insert(0, str(Path.home() / "projects"))

from figure_bank.scaffolds import SCF013_plot_grid_1panel as scf013
from figure_bank.matplotlib_panel import PanelRenderer, auto_place_legend

FIGURE_BANK = Path.home() / "projects" / "figure_bank"
EXP_DIR = Path.home() / "projects" / "experiments" / "EXP26000390-vr4-splicing-screen-rtpcr"
PANELS_DIR = FIGURE_BANK / "_retrofit_panels"
PANELS_DIR.mkdir(exist_ok=True)

# --- Data -----------------------------------------------------------------

df_all = pd.read_csv(EXP_DIR / "donor_feature_table.csv", encoding="utf-8-sig")
df_all["vp3_pct_wt"] = df_all["vp3_pct_wt"].clip(lower=0)

# --- Shared PanelRenderer -------------------------------------------------

renderer = PanelRenderer(
    body_w_px=890,
    body_h_px=316,
    figsize_inches=(5.93, 2.11),
    dpi=150,
    margins={"left": 0.10, "right": 0.99, "top": 0.97, "bottom": 0.13},
    axis_priority="y",
)

# --- Encoding constants ---------------------------------------------------

GTAAG_COLOR = "#C0392B"    # red — GTAAG present
BASE_COLOR   = "#2B7BB9"   # blue — GTAAG absent

CTX_MARKERS = {
    "CAG": "o",
    "AAG": "s",
    "ACG": "^",
    "AGG": "D",
    "GGG": "P",
}

def _ctx_marker(ctx: str) -> str:
    return CTX_MARKERS.get(ctx, "o")


# --- Plot body ------------------------------------------------------------

class _FigSpec(NamedTuple):
    fig_id: str
    panel_letter: str
    x_col: str
    x_label_svg: str    # SVG slot text (may use HTML entities)
    x_label_ax: str     # matplotlib axis label (plain text, not used — scaffold owns)
    invert_x: bool
    ref_lines: list     # list of (x_value, label_str, color)
    title: str
    subtitle: str
    out_stem: str


SPECS = [
    _FigSpec(
        fig_id="FIG168",
        panel_letter="A",
        x_col="hb_count",
        x_label_svg="H-bond count (U1:donor)",
        x_label_ax="H-bond count",
        invert_x=False,
        ref_lines=[],
        title="FP Panel: VP3% vs U1 snRNA H-bond Count (n=32)",
        subtitle="GTAAG=red, no-GTAAG=blue. Markers: o=CAG, s=AAG, ^=ACG, D=AGG.",
        out_stem="FIG168-fp-vp3-vs-hbcount-scatter",
    ),
    _FigSpec(
        fig_id="FIG169",
        panel_letter="A",
        x_col="hb_count",
        x_label_svg="H-bond count (U1:donor)",
        x_label_ax="H-bond count",
        invert_x=False,
        ref_lines=[],
        title="HEK Panel: VP3% vs U1 snRNA H-bond Count (n=41)",
        subtitle="GTAAG=red, no-GTAAG=blue. Markers: o=CAG, s=AAG, ^=ACG, D=AGG.",
        out_stem="FIG169-hek-vp3-vs-hbcount-scatter",
    ),
    _FigSpec(
        fig_id="FIG170",
        panel_letter="A",
        x_col="delta_g_kcal",
        x_label_svg="ΔG U1:donor (kcal/mol, more stable →)",
        x_label_ax="dG U1:donor (kcal/mol)",
        invert_x=True,
        ref_lines=[],
        title="FP Panel: VP3% vs U1:Donor Duplex Stability (n=32)",
        subtitle="X-axis inverted: more stable →. GTAAG=red, non-GTAAG=blue.",
        out_stem="FIG170-fp-vp3-vs-delta-g-scatter",
    ),
    _FigSpec(
        fig_id="FIG171",
        panel_letter="A",
        x_col="delta_g_kcal",
        x_label_svg="ΔG U1:donor (kcal/mol, more stable →)",
        x_label_ax="dG U1:donor (kcal/mol)",
        invert_x=True,
        ref_lines=[],
        title="HEK Panel: VP3% vs U1:Donor Duplex Stability (n=41)",
        subtitle="X-axis inverted: more stable →. GTAAG=red, non-GTAAG=blue.",
        out_stem="FIG171-hek-vp3-vs-delta-g-scatter",
    ),
    _FigSpec(
        fig_id="FIG172",
        panel_letter="A",
        x_col="maxent_score",
        x_label_svg="MaxEntScan score (bits)",
        x_label_ax="MaxEntScan score (bits)",
        invert_x=False,
        ref_lines=[],
        title="FP Panel: VP3% vs MaxEntScan Score (n=32)",
        subtitle="GTAAG=red, non-GTAAG=blue. Markers: o=CAG, s=AAG, ^=ACG, D=AGG.",
        out_stem="FIG172-fp-vp3-vs-maxent-scatter",
    ),
    _FigSpec(
        fig_id="FIG173",
        panel_letter="A",
        x_col="maxent_score",
        x_label_svg="MaxEntScan score (bits)",
        x_label_ax="MaxEntScan score (bits)",
        invert_x=False,
        ref_lines=[],
        title="HEK Panel: VP3% vs MaxEntScan Score (n=41)",
        subtitle="GTAAG=red, non-GTAAG=blue. Markers: o=CAG, s=AAG, ^=ACG, D=AGG.",
        out_stem="FIG173-hek-vp3-vs-maxent-scatter",
    ),
]

PANEL_FILTER = {
    "FIG168": "FP", "FIG169": "HEK",
    "FIG170": "FP", "FIG171": "HEK",
    "FIG172": "FP", "FIG173": "HEK",
}


def _plot_scatter(ax: plt.Axes, df: pd.DataFrame, spec: _FigSpec) -> None:
    # Plot each construct as a colored, shaped scatter point
    for _, row in df.iterrows():
        color = GTAAG_COLOR if row["gtaag_present"] else BASE_COLOR
        marker = _ctx_marker(row["exonic_ctx"])
        ax.scatter(row[spec.x_col], row["vp3_pct_wt"],
                   color=color, marker=marker, s=28, zorder=4,
                   edgecolors="none", alpha=0.88)

    # Axis limits + grid
    ax.set_ylim(0, 105)
    ax.yaxis.set_major_locator(plt.MultipleLocator(25))
    ax.grid(axis="y", color="#EEEEEE", lw=0.5, zorder=0)
    ax.set_axisbelow(True)

    if spec.invert_x:
        ax.invert_xaxis()

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # Text-clean: scaffold owns title/label text
    ax.set_title("")
    ax.set_xlabel("")
    ax.set_ylabel("")

    # Y tick labels via PanelRenderer
    renderer.size_y_labels(ax, labels=["0%", "25%", "50%", "75%", "100%"])

    # Legend: GTAAG color handles + context marker handles, combined
    h_gtaag = mlines.Line2D([], [], color=GTAAG_COLOR, marker="o", ls="none",
                             markersize=5, label="GTAAG")
    h_base  = mlines.Line2D([], [], color=BASE_COLOR,  marker="o", ls="none",
                             markersize=5, label="No GTAAG")
    ctx_handles = []
    for ctx, mk in CTX_MARKERS.items():
        present = df["exonic_ctx"].isin([ctx]).any()
        if present:
            ctx_handles.append(
                mlines.Line2D([], [], color="#555", marker=mk, ls="none",
                              markersize=5, label=ctx)
            )

    all_handles = [h_gtaag, h_base] + ctx_handles
    result = auto_place_legend(
        ax,
        handles=all_handles,
        fontsize=6,
        framealpha=0.85,
        edgecolor="#ccc",
        handlelength=1.0,
        handletextpad=0.4,
        borderpad=0.5,
        labelspacing=0.3,
    )
    if result.overlapped:
        print(f"  ! legend overlap: {result.overlap_artists}")


def build_scatter(spec: _FigSpec) -> None:
    panel_filter = PANEL_FILTER[spec.fig_id]
    df = df_all[df_all["panel"] == panel_filter].copy()

    out_panel = PANELS_DIR / f"{spec.out_stem}.png"
    fig, ax = renderer.fig_ax()
    _plot_scatter(ax, df, spec)
    for c in renderer.check_collisions(ax):
        print(f"  ! collision {spec.fig_id}: {c}")
    uri = renderer.save_as_data_uri(fig, out_panel)

    svg = scf013.render(
        title=spec.title,
        subtitle=spec.subtitle,
        panels=[scf013.Panel(
            panel_id=spec.panel_letter,
            panel_title=f"VP3 % of WT vs {spec.x_label_ax}",
            x_axis_label=spec.x_label_svg,
            y_axis_label="VP3 % of WT",
            image_href=uri,
        )],
    )
    out = FIGURE_BANK / f"{spec.out_stem}.svg"
    out.write_text(svg, encoding="utf-8")
    print(f"  OK    {out.name}")


def build_all() -> None:
    for spec in SPECS:
        build_scatter(spec)


if __name__ == "__main__":
    build_all()
