"""FIG164 — Controlled Series: Position-Specific Effects Beyond H-bond Tier.

SCF014 (2-panel).
  Panel A: CAG +6 series (AVD153/155/156/154) — same +5=A context, varying +6.
           Shows +6=G gives worst VP3 despite maximum H-bond contribution.
  Panel B: ACG +5 series (AVD145/143/144/146) — same ACG context, varying +5.
           Shows +5=G (completing GTAAG motif) gives 93.56% vs <20% otherwise.

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

# Panel A: CAG +6 series, sorted ascending by VP3%
PANEL_A_AVDS = ["AVD153", "AVD156", "AVD155", "AVD154"]   # +6=G,C,C,T sorted by VP3
PANEL_A_LABELS = ["+6=G", "+6=C", "+6=C", "+6=T"]
PANEL_A_VP3 = [df.loc[a, "vp3_pct_wt"] for a in PANEL_A_AVDS]
PANEL_A_COLORS = ["#C0392B" if df.loc[a, "plus6"] == "G" else "#2B7BB9" for a in PANEL_A_AVDS]

# Panel B: ACG +5 series, sorted ascending by VP3%
PANEL_B_AVDS = ["AVD145", "AVD143", "AVD144", "AVD146"]   # 0.19, 1.59, 19.41, 93.56
PANEL_B_LABELS = ["+5=C", "+5=A\n+6=G", "+5=A\n+6=T", "+5=G\nGTAAG"]
PANEL_B_VP3 = [df.loc[a, "vp3_pct_wt"] for a in PANEL_B_AVDS]
PANEL_B_COLORS = ["#FF8C00" if df.loc[a, "gtaag_present"] else "#2B7BB9" for a in PANEL_B_AVDS]

# --- PanelRenderer --------------------------------------------------------

# SCF014 per-panel body: 426 × 320 px
renderer = PanelRenderer(
    body_w_px=426,
    body_h_px=320,
    figsize_inches=(2.84, 2.13),
    dpi=150,
    margins={"left": 0.18, "right": 0.98, "top": 0.97, "bottom": 0.22},
    axis_priority="y",
)


def _plot_bar_panel(ax: plt.Axes, labels: list, vp3: list, colors: list,
                    ylim: float = 105) -> None:
    x = np.arange(len(labels))
    ax.bar(x, vp3, width=0.55, color=colors, edgecolor="none", zorder=3)

    # Value labels above bars
    for xi, yv in zip(x, vp3):
        ax.text(xi, yv + 1.5, f"{yv:.1f}%", ha="center", va="bottom",
                fontsize=7, color="#333", fontweight="bold")

    ax.set_xlim(-0.6, len(labels) - 0.4)
    ax.set_ylim(0, ylim)
    ax.yaxis.set_major_locator(plt.MultipleLocator(25))
    ax.grid(axis="y", color="#EEEEEE", lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    ax.set_xticks(x)
    ax.set_title("")
    ax.set_xlabel("")
    ax.set_ylabel("")


def _make_panel_a(ax: plt.Axes) -> None:
    _plot_bar_panel(ax, PANEL_A_LABELS, PANEL_A_VP3, PANEL_A_COLORS, ylim=65)
    sizing = renderer.size_x_labels(ax, labels=PANEL_A_LABELS,
                                    max_pt=9, target_pt=8, min_pt=6)
    if sizing.warnings:
        for w in sizing.warnings:
            print(f"  ! panelA x-sizing: {w}")
    ax.set_yticks([0, 25, 50])
    renderer.size_y_labels(ax, labels=["0%", "25%", "50%"])
    ax.set_ylim(0, 65)


def _make_panel_b(ax: plt.Axes) -> None:
    _plot_bar_panel(ax, PANEL_B_LABELS, PANEL_B_VP3, PANEL_B_COLORS, ylim=105)
    # Star callout on AVD146 — placed well above the value label (which sits at VP3+1.5)
    ax.text(3, PANEL_B_VP3[3] + 8, "★ GTAAG", ha="center", va="bottom",
            fontsize=7, color="#FF8C00", fontweight="bold")
    sizing = renderer.size_x_labels(ax, labels=PANEL_B_LABELS,
                                    max_pt=9, target_pt=8, min_pt=6)
    if sizing.warnings:
        for w in sizing.warnings:
            print(f"  ! panelB x-sizing: {w}")
    ax.set_yticks([0, 25, 50, 75, 100])
    renderer.size_y_labels(ax, labels=["0%", "25%", "50%", "75%", "100%"])


# --- Build ----------------------------------------------------------------

def build_fig164() -> None:
    # Panel A
    out_a = PANELS_DIR / "FIG164-controlled-series-pA.png"
    fig_a, ax_a = renderer.fig_ax()
    _make_panel_a(ax_a)
    for c in renderer.check_collisions(ax_a):
        print(f"  ! collision A: {c}")
    uri_a = renderer.save_as_data_uri(fig_a, out_a)

    # Panel B
    out_b = PANELS_DIR / "FIG164-controlled-series-pB.png"
    fig_b, ax_b = renderer.fig_ax()
    _make_panel_b(ax_b)
    for c in renderer.check_collisions(ax_b):
        print(f"  ! collision B: {c}")
    uri_b = renderer.save_as_data_uri(fig_b, out_b)

    svg = scf014.render(
        title="Controlled Series: Position-Specific Effects Beyond H-bond Tier",
        subtitle="Same exonic context and thermodynamic tier; single-position changes dominate VP3 outcome.",
        panels=[
            scf014.Panel(
                panel_id="A",
                panel_title="CAG +6 Series (+5=A, T5-T7)",
                x_axis_label="+6 position",
                y_axis_label="VP3 % of WT",
                image_href=uri_a,
            ),
            scf014.Panel(
                panel_id="B",
                panel_title="ACG +5 Series (T3-T4, GTAAG)",
                x_axis_label="+5 position",
                y_axis_label="VP3 % of WT",
                image_href=uri_b,
            ),
        ],
    )
    out = FIGURE_BANK / "FIG164-controlled-series-comparisons.svg"
    out.write_text(svg, encoding="utf-8")
    print(f"  OK    {out.name}")


if __name__ == "__main__":
    build_fig164()
