#!/usr/bin/env python3
"""
Cryptic Splice Site Visualization (v2) — MaxEntScan Only

Generates a publication-quality figure comparing cryptic splice sites in the
447 bp intron body (GGGGS5 + VHH3 + GGGGS1) before and after codon optimization.

Top panel:  BEFORE optimization (original VHH3)
Bottom panel: AFTER optimization (context-optimized VHH3)

All splice sites plotted upward (taller = stronger). Donors and acceptors
distinguished by color and marker shape (donors: circles, acceptors: triangles).

Flanking clouds show the 41 intended donor/acceptor scores from the AVD059-099
intron-retention series, emphasizing that intended sites are far stronger than
any residual cryptic site.

Usage:
    python scripts/tools/visualize_cryptic_splice_sites.py
"""

import csv
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# Add tools directory to path for imports
TOOLS_DIR = Path(__file__).parent
sys.path.insert(0, str(TOOLS_DIR))

from build_splice_cassettes import (
    score_donor_maxentscan,
    score_acceptor_maxentscan,
)

# ============================================================================
# Constants
# ============================================================================

# BEFORE: Original GGGGS from AVD006 canonical (varied-codon)
GGGGS5_BEFORE = "GGTGGAGGCGGATCTGGAGGCGGTGGTTCAGGCGGTGGAGGAAGTGGTGGCGGAGGTTCTGGTGGAGGCGGTTCT"
GGGGS1_BEFORE = "GGTGGAGGCGGTTCT"

# AFTER: GT/AG-free GGGGS (from optimize_intron_body_codons.py, final after all stages)
GGGGS5_AFTER = "GGCGGGGGCGGATCTGGGGGCGGGGGCTCCGGCGGCGGGGGATCGGGGGGGGGCGGCTCTGGCGGGGGGGGATCG"
GGGGS1_AFTER = "GGGGGCGGCGGATCG"

# BEFORE: AVD006 VHH3 (canonical, pre-optimization)
VHH3_BEFORE = (
    "GAGGTGCAACTGGTTGAAAGCGGCGGAGGACTTGTTCAACCCGGCGGCAGCCTTAGGCTTTCTTGCGCTGCCAGC"
    "GGCTTCACCTTTAGCACCGCCGACATGGGCTGGTTTAGGCAAGCTCCCGGAAAAGGCAGGGAACTTGTTGCCGCT"
    "GTGAGCGGCAGCGGCTTCAGCACCTACTCTGATAGCGTTGAGGGCAGGTTCACCATCAGCAGGGACAACGCCAAG"
    "AGGATGGTGTACCTGCAGATGAACAGCTTGAGGGCCGAGGACACCGCCGTGTACTACTGCGCCAAGGCCACAATT"
    "AGCCTGTACTACGCCATGGATGTGTGGGGACAGGGCACCACCGTGACCGTGAGCAGC"
)

# AFTER: DNAchisel + MaxEntScan optimized VHH3 (from optimize_intron_body_codons.py Stages 1-2, seed=42)
# Updated 2026-02-12: relaxed do-no-harm safe_ceiling=0.0, pos396 GTA→GTC
VHH3_AFTER = (
    "GAGGTCCAACTCGTCGAATCCGGCGGGGGACTCGTGCAACCCGGCGGCTCCCTGCGCCTA"
    "TCATGCGCGGCATCCGGATTCACATTTTCAACCGCCGATATGGGATGGTTCCGGCAGGCA"
    "CCCGGGAAGGGACGGGAACTCGTCGCCGCCGTTTCCGGATCCGGATTTTCAACATATTCG"
    "GATTCCGTGGAGGGACGATTCACCATATCCCGGGACAATGCAAAACGGATGGTCTACTTG"
    "CAAATGAATTCATTGCGCGCGGAGGACACCGCCGTCTACTATTGCGCGAAGGCAACAATT"
    "TCATTATATTATGCCATGGACGTCTGGGGGCAGGGAACCACCGTCACCGTTTCCTCC"
)

# Region boundaries in the 447 bp body
GGGGS5_END = len(GGGGS5_BEFORE)      # 75
VHH3_END = GGGGS5_END + 357          # 432
BODY_LEN = GGGGS5_END + 357 + 15     # 447

# Summary CSV with 41 intended donors/acceptors
SUMMARY_CSV = (
    Path(__file__).parent.parent.parent
    / "projects" / "AVD_VHH_Display_ALPL" / "plasmids" / "intron_retention"
    / "AVD059-099_intron_retention_summary.csv"
)

OUTPUT_DIR = TOOLS_DIR / "output"


# ============================================================================
# Scoring helpers
# ============================================================================

def find_dinuc_positions(seq: str, dinuc: str) -> List[int]:
    """Find all 0-indexed positions of a dinucleotide."""
    return [i for i in range(len(seq) - 1) if seq[i:i + 2] == dinuc]


def score_gt_donors(seq: str) -> List[Dict]:
    """Score all GT dinucleotides as potential splice donors using MaxEntScan."""
    gt_positions = find_dinuc_positions(seq, "GT")
    results = []
    for pos in gt_positions:
        s9 = pos - 3
        e9 = pos + 6
        if 0 <= s9 and e9 <= len(seq):
            niner = seq[s9:e9]
            score = score_donor_maxentscan(niner)
        else:
            score = None
        results.append({"pos": pos, "score": score})
    return results


def score_ag_acceptors(seq: str) -> List[Dict]:
    """Score AG dinucleotides using MaxEntScan acceptor scoring.

    For a cryptic AG at 0-indexed position P in the intron body, the 23-mer
    is body[P-18 : P+5] (18 nt intronic context + AG + 3 nt exonic context).
    Valid when P >= 18 and P + 5 <= len(seq).
    """
    ag_positions = find_dinuc_positions(seq, "AG")
    results = []
    for pos in ag_positions:
        start = pos - 18
        end = pos + 5  # AG occupies pos,pos+1; +3 exonic = pos+2..pos+4 → end = pos+5
        if start >= 0 and end <= len(seq):
            mer23 = seq[start:end]
            score = score_acceptor_maxentscan(mer23)
        else:
            score = None
        results.append({"pos": pos, "score": score})
    return results


def load_intended_scores(csv_path: str):
    """Load intended donor and acceptor scores from the AVD summary CSV.

    Returns (donor_scores, acceptor_scores) as lists of floats.
    - donor 15bp column → 9-mer = donor[6:15]
    - acceptor 30bp column → 23-mer = acceptor[1:24]
    """
    donor_scores = []
    acceptor_scores = []

    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Donor: extract 9-mer from 15bp donor column
            donor_15 = row["donor"]
            if len(donor_15) >= 15:
                niner = donor_15[6:15]
                ds = score_donor_maxentscan(niner)
                if ds > -100:  # filter invalid (no GT at expected position)
                    donor_scores.append(ds)

            # Acceptor: extract 23-mer from 30bp acceptor column
            acceptor_30 = row["acceptor"]
            if len(acceptor_30) >= 24:
                mer23 = acceptor_30[1:24]
                acs = score_acceptor_maxentscan(mer23)
                if acs > -100:  # filter invalid (no AG at expected position)
                    acceptor_scores.append(acs)

    return donor_scores, acceptor_scores


# ============================================================================
# Color helpers
# ============================================================================

def get_donor_color(score):
    """Color GT donors by MaxEntScan score strength."""
    if score is None:
        return "#999999"
    if score >= 3.0:
        return "#d62728"   # red — Weak or above
    if score >= 0.0:
        return "#ff7f0e"   # orange — Very Weak
    return "#999999"       # gray — Not Detected


def get_acceptor_color(score):
    """Color AG acceptors by MaxEntScan acceptor score strength."""
    if score is None:
        return "#999999"
    if score >= 3.0:
        return "#1f77b4"   # blue — Weak or above
    if score >= 0.0:
        return "#aec7e8"   # light blue — Very Weak
    return "#999999"       # gray — Not Detected


# ============================================================================
# Plotting
# ============================================================================

def plot_panel(
    ax,
    label: str,
    gt_data: List[Dict],
    ag_data: List[Dict],
    intended_donor_scores: List[float],
    intended_acceptor_scores: List[float],
    y_top: float,
    y_bot: float,
    rng: np.random.Generator,
):
    """Plot one panel (Before or After) with cryptic lollipops and intended clouds."""

    # --- Region shading for intron body (x 0..447) ---
    ax.axvspan(0, GGGGS5_END, alpha=0.12, color="#5b9bd5", zorder=0)
    ax.axvspan(GGGGS5_END, VHH3_END, alpha=0.08, color="#ffd966", zorder=0)
    ax.axvspan(VHH3_END, BODY_LEN, alpha=0.12, color="#70ad47", zorder=0)

    # Region labels
    label_y = y_top - 0.5
    ax.text(GGGGS5_END / 2, label_y, "GGGGS5", ha="center", va="top",
            fontsize=8, fontstyle="italic", color="#2e75b6", fontweight="bold")
    ax.text((GGGGS5_END + VHH3_END) / 2, label_y, "VHH3", ha="center", va="top",
            fontsize=8, fontstyle="italic", color="#bf8f00", fontweight="bold")
    ax.text((VHH3_END + BODY_LEN) / 2, label_y, "G1", ha="center", va="top",
            fontsize=7, fontstyle="italic", color="#548235", fontweight="bold")

    # --- Intended Donor cloud (left zone) ---
    cloud_donor_x_center = -35
    ax.axvspan(-55, -15, alpha=0.10, color="#70ad47", zorder=0)
    jx = rng.uniform(-8, 8, size=len(intended_donor_scores))
    ax.scatter(
        cloud_donor_x_center + jx,
        intended_donor_scores,
        c="#2ca02c", alpha=0.45, s=30, marker="o", zorder=4,
        edgecolors="white", linewidths=0.4,
    )
    median_d = float(np.median(intended_donor_scores))
    ax.hlines(median_d, -55, -15, colors="#2ca02c", linewidths=1, linestyles="--", zorder=5)
    n_d = len(intended_donor_scores)
    ax.text(-35, y_top - 0.3, f"Intended\nDonors (n={n_d})", ha="center", va="top",
            fontsize=6.5, color="#2ca02c", fontweight="bold")

    # --- Intended Acceptor cloud (right zone) ---
    cloud_acc_x_center = 482
    ax.axvspan(462, 502, alpha=0.10, color="#9467bd", zorder=0)
    jx_a = rng.uniform(-8, 8, size=len(intended_acceptor_scores))
    ax.scatter(
        cloud_acc_x_center + jx_a,
        intended_acceptor_scores,
        c="#9467bd", alpha=0.45, s=30, marker="^", zorder=4,
        edgecolors="white", linewidths=0.4,
    )
    median_a = float(np.median(intended_acceptor_scores))
    ax.hlines(median_a, 462, 502, colors="#9467bd", linewidths=1, linestyles="--", zorder=5)
    n_a = len(intended_acceptor_scores)
    ax.text(482, y_top - 0.3, f"Intended\nAcceptors (n={n_a})", ha="center", va="top",
            fontsize=6.5, color="#9467bd", fontweight="bold")

    # --- Zero line and threshold ---
    ax.axhline(y=0, color="black", linewidth=0.5, zorder=1)
    ax.axhline(y=3.0, color="gray", linewidth=0.8, linestyle="--", alpha=0.5, zorder=1)

    # --- GT donors: upward lollipops ---
    for d in gt_data:
        score = d["score"]
        if score is None:
            continue
        color = get_donor_color(score)
        ax.plot([d["pos"], d["pos"]], [0, score], color=color, linewidth=1.5,
                zorder=2, solid_capstyle="round")
        ax.plot(d["pos"], score, "o", color=color, markersize=5.5, zorder=3,
                markeredgecolor="white", markeredgewidth=0.3)

    # --- AG acceptors: upward lollipops (triangles) ---
    for d in ag_data:
        score = d["score"]
        if score is None:
            continue
        color = get_acceptor_color(score)
        ax.plot([d["pos"], d["pos"]], [0, score], color=color, linewidth=1.5,
                zorder=2, solid_capstyle="round")
        ax.plot(d["pos"], score, "^", color=color, markersize=5.5, zorder=3,
                markeredgecolor="white", markeredgewidth=0.3)

    # --- Panel title ---
    gt_count = len(gt_data)
    ag_count = len(ag_data)
    ax.set_title(f"{label}   (GT: {gt_count}, AG: {ag_count})",
                 fontsize=11, fontweight="bold", loc="left")


# ============================================================================
# Main
# ============================================================================

def main():
    print("Generating Cryptic Splice Site Visualization (v2, MaxEntScan only)")
    print("=" * 65)

    # Build full 447 bp bodies
    body_before = GGGGS5_BEFORE + VHH3_BEFORE + GGGGS1_BEFORE
    body_after = GGGGS5_AFTER + VHH3_AFTER + GGGGS1_AFTER

    assert len(body_before) == BODY_LEN, f"Before body: {len(body_before)} bp (expected {BODY_LEN})"
    assert len(body_after) == BODY_LEN, f"After body: {len(body_after)} bp (expected {BODY_LEN})"

    # Verify GT/AG counts
    print(f"  BEFORE: GT={body_before.count('GT')}, AG={body_before.count('AG')}")
    print(f"  AFTER:  GT={body_after.count('GT')}, AG={body_after.count('AG')}")

    # Score cryptic donors (MaxEntScan)
    print("\nScoring GT donors (MaxEntScan)...")
    gt_donors_before = score_gt_donors(body_before)
    gt_donors_after = score_gt_donors(body_after)

    for tag, data in [("BEFORE", gt_donors_before), ("AFTER", gt_donors_after)]:
        scored = [d for d in data if d["score"] is not None]
        if scored:
            scores = [d["score"] for d in scored]
            print(f"  {tag}: {len(data)} GT sites, scores {min(scores):.2f} to {max(scores):.2f}")

    # Score cryptic acceptors (MaxEntScan)
    print("\nScoring AG acceptors (MaxEntScan)...")
    ag_acc_before = score_ag_acceptors(body_before)
    ag_acc_after = score_ag_acceptors(body_after)

    for tag, data in [("BEFORE", ag_acc_before), ("AFTER", ag_acc_after)]:
        scored = [d for d in data if d["score"] is not None]
        if scored:
            scores = [d["score"] for d in scored]
            print(f"  {tag}: {len(data)} AG sites, scores {min(scores):.2f} to {max(scores):.2f}")

    # Load intended donor/acceptor scores from AVD059-099
    print(f"\nLoading intended splice site scores from {SUMMARY_CSV.name}...")
    if not SUMMARY_CSV.exists():
        print(f"  ERROR: CSV not found at {SUMMARY_CSV}")
        return 1

    intended_donor_scores, intended_acceptor_scores = load_intended_scores(str(SUMMARY_CSV))
    print(f"  Intended donors:    n={len(intended_donor_scores)}, "
          f"median={np.median(intended_donor_scores):.2f}, "
          f"range=[{min(intended_donor_scores):.2f}, {max(intended_donor_scores):.2f}]")
    print(f"  Intended acceptors: n={len(intended_acceptor_scores)}, "
          f"median={np.median(intended_acceptor_scores):.2f}, "
          f"range=[{min(intended_acceptor_scores):.2f}, {max(intended_acceptor_scores):.2f}]")

    # ---- Compute shared y-axis limits (everything upward) ----
    all_gt = ([d["score"] for d in gt_donors_before if d["score"] is not None]
              + [d["score"] for d in gt_donors_after if d["score"] is not None])
    all_ag = ([d["score"] for d in ag_acc_before if d["score"] is not None]
              + [d["score"] for d in ag_acc_after if d["score"] is not None])

    all_scores = all_gt + all_ag + intended_donor_scores + intended_acceptor_scores
    y_max = max(all_scores) + 2.0
    y_min = min(all_scores) - 1.0

    # ---- Create figure ----
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(15, 8), sharex=True,
        gridspec_kw={"hspace": 0.32},
    )

    rng = np.random.default_rng(42)

    for ax in (ax1, ax2):
        ax.set_ylim(y_min, y_max)
        ax.set_xlim(-65, 515)

    # Plot panels
    plot_panel(ax1, "Before Optimization",
               gt_donors_before, ag_acc_before,
               intended_donor_scores, intended_acceptor_scores,
               y_max, y_min, rng)

    plot_panel(ax2, "After Optimization",
               gt_donors_after, ag_acc_after,
               intended_donor_scores, intended_acceptor_scores,
               y_max, y_min, rng)

    # Axis labels
    ax2.set_xlabel("Nucleotide position in intron body (bp)", fontsize=10)
    for ax in (ax1, ax2):
        ax.set_ylabel("MaxEntScan Score", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Threshold annotation
    for ax in (ax1, ax2):
        ax.text(0.60, 3.3, '"Weak" splice site threshold', fontsize=6.5,
                color="gray", alpha=0.6, ha="left",
                transform=ax.get_yaxis_transform())

    # Max-score annotations
    max_before = max(d["score"] for d in gt_donors_before if d["score"] is not None)
    max_after = max(d["score"] for d in gt_donors_after if d["score"] is not None)
    ax1.annotate(
        f"Max cryptic GT: {max_before:.1f}",
        xy=(0.35, 0.92), xycoords="axes fraction",
        fontsize=8, color="#ff7f0e", ha="center",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#ff7f0e", alpha=0.8),
    )
    after_note = " (all < 3)" if max_after < 3.0 else ""
    after_color = "#999999" if max_after < 3.0 else "#ff7f0e"
    ax2.annotate(
        f"Max cryptic GT: {max_after:.1f}{after_note}",
        xy=(0.35, 0.92), xycoords="axes fraction",
        fontsize=8, color=after_color, ha="center",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=after_color, alpha=0.8),
    )

    # ---- Simplified Legend ----
    legend_elements = [
        # Donors (upward)
        Line2D([0], [0], marker="o", color="#d62728",
               label="GT donor: Weak+ (\u2265 3)", markersize=6, linestyle="-", linewidth=1),
        Line2D([0], [0], marker="o", color="#ff7f0e",
               label="GT donor: Very Weak (0\u20133)", markersize=6, linestyle="-", linewidth=1),
        Line2D([0], [0], marker="o", color="#999999",
               label="GT donor: Not Detected (< 0)", markersize=6, linestyle="-", linewidth=1),
        # Acceptors (triangles)
        Line2D([0], [0], marker="^", color="#1f77b4",
               label="AG acceptor: Weak+ (\u2265 3)", markersize=6, linestyle="-", linewidth=1),
        Line2D([0], [0], marker="^", color="#aec7e8",
               label="AG acceptor: Very Weak (0\u20133)", markersize=6, linestyle="-", linewidth=1),
        Line2D([0], [0], marker="^", color="#999999",
               label="AG acceptor: Not Detected (< 0)", markersize=6, linestyle="-", linewidth=1),
        # Clouds
        Line2D([0], [0], marker="o", color="#2ca02c", linestyle="None",
               label=f"Intended donor (n={len(intended_donor_scores)})", markersize=6, alpha=0.6),
        Line2D([0], [0], marker="^", color="#9467bd", linestyle="None",
               label=f"Intended acceptor (n={len(intended_acceptor_scores)})", markersize=6, alpha=0.6),
    ]
    fig.legend(
        handles=legend_elements, loc="lower center", ncol=4, fontsize=7.5,
        frameon=True, bbox_to_anchor=(0.5, -0.01),
    )

    fig.suptitle(
        "Cryptic Splice Site Reduction by VHH3 Codon Optimization",
        fontsize=13, fontweight="bold", y=0.98,
    )

    fig.subplots_adjust(bottom=0.13, top=0.93)

    # ---- Save ----
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "cryptic_splice_before_after.png"
    fig.savefig(str(output_path), dpi=300, bbox_inches="tight", facecolor="white")
    print(f"\nSaved figure to: {output_path}")

    try:
        plt.show(block=False)
    except Exception:
        pass

    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
