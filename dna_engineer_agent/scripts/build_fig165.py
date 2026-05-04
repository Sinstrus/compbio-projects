"""FIG165 — +6=G Mismatch Penalty: Consistent 10–15 pp Cost Across Contexts.

SCF014 (2-panel).
  Panel A: CAG +5=A series (AVD153/156/155/154) sorted ascending VP3%.
           Shows +6=G consistently lowest VP3 within tier.
  Panel B: Cross-series +6=G vs +6=T pairs (T4 ACG + T5-T6 CAG).

Data: experiments/EXP26000390-vr4-splicing-screen-rtpcr/donor_feature_table.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path.home() / "projects"))

from figure_bank.scaffolds import SCF014_plot_grid_2panel as scf014
from figure_bank.matplotlib_panel import PanelRenderer

FIGURE_BANK = Path.home() / "projects" / "figure_bank"
EXP_DIR = Path.home() / "projects" / "experiments" / "EXP26000390-vr4-splicing-screen-rtpcr"
PANELS_DIR = FIGURE_BANK / "_retrofit_panels"
PANELS_DIR.mkdir(exist_ok=True)

# --- Data -----------------------------------------------------------------

df = pd.read_csv(EXP_DIR / "donor_feature_table.csv", encoding="utf-8-sig")
df["vp3_pct_wt"] = df["vp3_pct_wt"].clip(lower=0)
df = df.set_index("AVD_construct")

# Panel A: CAG +5=A, sorted ascending by VP3% (worst to best)
PA_AVDS  = ["AVD153", "AVD156", "AVD155", "AVD154"]
PA_VP3   = [df.loc[a, "vp3_pct_wt"] for a in PA_AVDS]
PA_PLUS6 = [df.loc[a, "plus6"] for a in PA_AVDS]
PA_LABELS = [f"+6={p}" for p in PA_PLUS6]

# Colors: G=red (kinetic trap), T=dark green (best despite model-worst), others=steelblue
def _pa_color(plus6):
    if plus6 == "G":
        return "#C0392B"
    if plus6 == "T":
        return "#1B6B2B"
    return "#2B7BB9"

PA_COLORS = [_pa_color(p) for p in PA_PLUS6]

# Panel B: cross-series pairs
# Pair 1 (T4 ACG): AVD143 (+6=G, 1.59%) vs AVD144 (+6=T, 19.41%)
# Pair 2 (T5-T6 CAG): AVD165 (+6=G, 44.84%) vs AVD163 (+6=T, 56.13%)
PB_AVDS   = ["AVD143", "AVD144", "AVD165", "AVD163"]
PB_VP3    = [df.loc[a, "vp3_pct_wt"] for a in PB_AVDS]
PB_PLUS6  = [df.loc[a, "plus6"] for a in PB_AVDS]
PB_COLORS = ["#C0392B" if p == "G" else "#1B6B2B" for p in PB_PLUS6]
PB_X      = [0, 0.65, 1.55, 2.2]   # pair gap between index 1 and 2
PB_LABELS = ["+6=G\nACG", "+6=T\nACG", "+6=G\nCAG", "+6=T\nCAG"]

# --- PanelRenderer --------------------------------------------------------

renderer = PanelRenderer(
    body_w_px=426,
    body_h_px=320,
    figsize_inches=(2.84, 2.13),
    dpi=150,
    margins={"left": 0.18, "right": 0.98, "top": 0.97, "bottom": 0.22},
    axis_priority="y",
)


def _make_panel_a(ax: plt.Axes) -> None:
    x = np.arange(len(PA_AVDS))
    ax.bar(x, PA_VP3, width=0.55, color=PA_COLORS, edgecolor="none", zorder=3)

    # Value labels above bars
    for xi, yv in zip(x, PA_VP3):
        ax.text(xi, yv + 1.5, f"{yv:.1f}%", ha="center", va="bottom",
                fontsize=7, color="#333", fontweight="bold")

    ax.set_xlim(-0.6, len(PA_AVDS) - 0.4)
    ax.set_ylim(0, 70)
    ax.yaxis.set_major_locator(plt.MultipleLocator(25))
    ax.grid(axis="y", color="#EEEEEE", lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_xticks(x)
    ax.set_title("")
    ax.set_xlabel("")
    ax.set_ylabel("")

    sizing = renderer.size_x_labels(ax, labels=PA_LABELS, max_pt=9, target_pt=8, min_pt=6)
    if sizing.warnings:
        for w in sizing.warnings:
            print(f"  ! panelA x-sizing: {w}")
    ax.set_yticks([0, 25, 50])
    renderer.size_y_labels(ax, labels=["0%", "25%", "50%"])
    ax.set_ylim(0, 70)


def _make_panel_b(ax: plt.Axes) -> None:
    for xi, yv, col in zip(PB_X, PB_VP3, PB_COLORS):
        ax.bar(xi, yv, width=0.55, color=col, edgecolor="none", zorder=3)

    # Value labels
    for xi, yv in zip(PB_X, PB_VP3):
        ax.text(xi, yv + 1.5, f"{yv:.1f}%", ha="center", va="bottom",
                fontsize=7, color="#333", fontweight="bold")

    ax.set_xlim(-0.5, max(PB_X) + 0.5)
    ax.set_ylim(0, 70)
    ax.yaxis.set_major_locator(plt.MultipleLocator(25))
    ax.grid(axis="y", color="#EEEEEE", lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_xticks(PB_X)
    ax.set_title("")
    ax.set_xlabel("")
    ax.set_ylabel("")

    sizing = renderer.size_x_labels(ax, labels=PB_LABELS, max_pt=9, target_pt=8, min_pt=6)
    if sizing.warnings:
        for w in sizing.warnings:
            print(f"  ! panelB x-sizing: {w}")
    ax.set_yticks([0, 25, 50])
    renderer.size_y_labels(ax, labels=["0%", "25%", "50%"])
    ax.set_ylim(0, 70)


# --- Build ----------------------------------------------------------------

def build_fig165() -> None:
    out_a = PANELS_DIR / "FIG165-kinetic-trap-pA.png"
    fig_a, ax_a = renderer.fig_ax()
    _make_panel_a(ax_a)
    for c in renderer.check_collisions(ax_a):
        print(f"  ! collision A: {c}")
    uri_a = renderer.save_as_data_uri(fig_a, out_a)

    out_b = PANELS_DIR / "FIG165-kinetic-trap-pB.png"
    fig_b, ax_b = renderer.fig_ax()
    _make_panel_b(ax_b)
    for c in renderer.check_collisions(ax_b):
        print(f"  ! collision B: {c}")
    uri_b = renderer.save_as_data_uri(fig_b, out_b)

    svg = scf014.render(
        title="+6=G Mismatch Penalty: Consistent 10–15 pp Cost Across Contexts",
        subtitle="+6=G mismatches U1 snRNA at pos 6 (corrected map); costs 10–15 pp VP3 per tier.",
        panels=[
            scf014.Panel(
                panel_id="A",
                panel_title="CAG +5=A Series (T5-T7, sorted)",
                x_axis_label="+6 position",
                y_axis_label="VP3 % of WT",
                image_href=uri_a,
            ),
            scf014.Panel(
                panel_id="B",
                panel_title="+6=G vs +6=T: ACG and CAG",
                x_axis_label="Exonic context / +6",
                y_axis_label="VP3 % of WT",
                image_href=uri_b,
            ),
        ],
    )
    out = FIGURE_BANK / "FIG165-plus6g-mismatch-penalty.svg"
    out.write_text(svg, encoding="utf-8")
    print(f"  OK    {out.name}")


if __name__ == "__main__":
    build_fig165()
