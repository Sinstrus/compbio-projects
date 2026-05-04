"""Build FIG167 via SCF025 — Wave 2 VR4 Intron Retention Design Panel Overview.

Replaces the former 435-line inline SVG (FIG167-wave2-vr4-ir-panel-overview.svg)
with a scaffold-rendered output. All domain-specific text (column headers,
pill content, row descriptions) lives here; SCF025 owns geometry only.

Source: EXP26000390.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path.home() / "projects"))

from figure_bank.scaffolds.SCF025_design_panel_overview import (
    Badge,
    LegendEntry,
    Pill,
    Row,
    render,
)

FIGURE_BANK = Path.home() / "projects" / "figure_bank"


# ---------------------------------------------------------------------------
# Key-findings pills (domain-specific content lives here, not in the scaffold)
# ---------------------------------------------------------------------------

PILLS = (
    Pill(
        number="①",
        title="Exonic −2 = A: Permissive",
        lines=(
            "CAG / AAG → functional; all others <5%",
            "34× difference from 2-nt exon change",
        ),
        fill="#E8F5E9",
        stroke="#43A047",
    ),
    Pill(
        number="②",
        title="+6 = G → Splicing Penalty",
        lines=(
            "+6=G is a mismatch (0 HB, corrected map)",
            "Costs 10–15 pp within matched tier",
        ),
        fill="#FFEBEE",
        stroke="#EF5350",
    ),
    Pill(
        number="③",
        title="GTAAG-like > canonical GTAAGT",
        lines=(
            "GT|AGGT: PABPC4 64.7%, H3-3B 46.2%",
            "GT|AAGT canonical: ATP6V0B 40.8%",
        ),
        fill="#E3F2FD",
        stroke="#1E88E5",
    ),
    Pill(
        number="④",
        title="GTAAG rescues exon context",
        lines=(
            "AVD146: ACG + GTAAG + ESE → 93.6%",
            "Despite non-permissive −2 position",
        ),
        fill="#F3E5F5",
        stroke="#AB47BC",
    ),
    Pill(
        number="⑤",
        title="No tier bottleneck at T7",
        lines=(
            "T7 + GTAAG achieves 40–51% splicing",
            "Higher HB count → higher baseline",
        ),
        fill="#FFF8E1",
        stroke="#FFB300",
    ),
)


# ---------------------------------------------------------------------------
# Column headers (caller provides all 8 labels — scaffold places them)
# ---------------------------------------------------------------------------

COL_HEADERS = (
    "PANEL",
    "NAME & DESIGN PRINCIPLE",
    "n",
    "EXONIC",
    "INTRONIC MOTIF",
    "+6",
    "TIER",
    "PREDICTED VP3% (% of WT)",
)


# ---------------------------------------------------------------------------
# Data rows (7 panels)
# ---------------------------------------------------------------------------

ROWS = (
    Row(
        letter="A",
        letter_color="#1565C0",
        row_color="#FFFFFF",
        panel_sublabel="6 constructs",
        panel_type="HEK Ports",
        name="HEK Port Validations",
        description=(
            "HEK ports (n=6): v4 acceptor + VHH27.",
            "Empirical bridge: do natural donors hold up?",
        ),
        annotation="PABPC4 · H3-3B · ATP6V0B · ELOA · MT2A",
        annotation_color="#78909C",
        n=6,
        category_badges=(
            Badge("AAG", color="#43A047"),
            Badge("CAG", color="#43A047"),
            Badge("CCG ×1", color="#EF9A9A", text_color="#B71C1C"),
        ),
        variant_lines=(
            "GT|AGGT ×2  (GTAAG-like)",
            "GT|AAGT ×1  (canonical)",
            "GT|AACA, GT|AAGA  (×2)",
        ),
        modifier=Badge("T", color="#1B5E20"),
        level="Natural",
        level_sublabel="gene-derived",
        outcome_range=(14.0, 65.0),
        outcome_color="#1565C0",
    ),
    Row(
        letter="B",
        letter_color="#1B5E20",
        row_color="#F9FBE7",
        panel_sublabel="8 constructs",
        panel_type="Ceiling",
        name="Ceiling Donors: CAG/AAG + GTAAG",
        description=(
            "Permissive exon + canonical GTAAG together.",
            "2 contexts × T2–T5 × +6=T = 8 constructs.",
        ),
        annotation="Flagship zone; basis for E & F dial-down.",
        annotation_color="#558B2F",
        n=8,
        category_badges=(
            Badge("CAG", color="#43A047"),
            Badge("AAG", color="#43A047"),
        ),
        variant_lines=(
            "GT|AAGXXX (GTAAG)",
            "+3=A, +4=A, +5=G",
        ),
        modifier=Badge("T", color="#1B5E20"),
        level="T2–T5",
        outcome_range=(75.0, 95.0),
        outcome_color="#1B5E20",
    ),
    Row(
        letter="C",
        letter_color="#388E3C",
        row_color="#FFFFFF",
        panel_sublabel="8 constructs",
        panel_type="GTAAG-like",
        name="GTAAG-like (AGGT Pattern) Sweep",
        description=(
            "Synthetic GT|AGGTXXX (PABPC4/H3-3B).",
            "2 contexts × T2–T5 × +6=T = 8 constructs.",
        ),
        annotation="Does AGGT generalize across tiers?",
        annotation_color="#558B2F",
        n=8,
        category_badges=(
            Badge("CAG", color="#43A047"),
            Badge("AAG", color="#43A047"),
        ),
        variant_lines=(
            "GT|AGGTXXX  (GTAAG-like)",
            "PABPC4/H3-3B pattern",
        ),
        modifier=Badge("T", color="#388E3C"),
        level="T2–T5",
        outcome_range=(50.0, 80.0),
        outcome_color="#388E3C",
    ),
    Row(
        letter="D",
        letter_color="#00695C",
        row_color="#E0F2F1",
        panel_sublabel="5 constructs",
        panel_type="AVD146",
        name="AVD146 Ports + CAG/AAG Upgrade",
        description=(
            "D1: port → v4+VHH27; D2: CAG; D3: AAG",
            "D4/D5: ESE ablation (GAGG→GACC/GAAT).",
        ),
        annotation="D2: top flagship — permissive + GTAAG + ESE.",
        annotation_color="#00897B",
        n=5,
        category_badges=(
            Badge("CAG", color="#43A047"),
            Badge("AAG", color="#43A047"),
            Badge("ACG ×1", color="#EF9A9A", text_color="#B71C1C"),
        ),
        variant_lines=(
            "GT|AAGAGG (AVD146)",
            "GTAAG + GAGG ESE",
            "D4/D5: ESE ablation",
        ),
        modifier=Badge("A", color=None),
        level="T4",
        level_sublabel="D1 fixed",
        outcome_range=(75.0, 100.0),
        outcome_color="#00695C",
    ),
    Row(
        letter="E",
        letter_color="#E65100",
        row_color="#FFFFFF",
        panel_sublabel="10 constructs",
        panel_type="Dial-Down",
        name="+6 = G Dial-Down",
        description=(
            "Panel B donors, +6 flipped T→G (~10–15 pp).",
            "CAG/AAG × T2–T5 (+6=G) = 8; ESE = 10",
        ),
        annotation="Primary dial-down to 80–90% zone.",
        annotation_color="#BF360C",
        n=10,
        category_badges=(
            Badge("CAG", color="#43A047"),
            Badge("AAG", color="#43A047"),
        ),
        variant_lines=(
            "GT|AAGXXG  (+6=G)",
            "GTAAG motif, +6 flipped T→G",
            "+10–15 pp penalty vs +6=T",
        ),
        modifier=Badge("G", color="#EF5350"),
        level="T1–T5",
        outcome_range=(65.0, 82.0),
        outcome_color="#E65100",
    ),
    Row(
        letter="F",
        letter_color="#F57F17",
        row_color="#FFF8E1",
        panel_sublabel="4 constructs",
        panel_type="",
        name="T1 Dial-Down: Permissive Context",
        description=(
            "CAG/AAG + GTAAG variants × Tier 1 only.",
        ),
        annotation="W1 T1 always non-permissive; unexplored gap.",
        annotation_color="#EF6C00",
        n=4,
        category_badges=(
            Badge("CAG", color="#43A047"),
            Badge("AAG", color="#43A047"),
        ),
        variant_lines=(
            "GTAAG + GTAAG-like",
        ),
        modifier=Badge("T, G", color=None),
        level="T1",
        outcome_range=(40.0, 70.0),
        outcome_color="#F57F17",
        outcome_uncertain=True,
    ),
    Row(
        letter="G",
        letter_color="#546E7A",
        row_color="#FAFAFA",
        panel_sublabel="",
        panel_type="",
        name="GT/AG-Null Controls",
        description=(
            "GT→GC + AG→AA on top B/D candidates.",
        ),
        annotation="",
        n=3,
        category_badges=(
            Badge("CAG", color="#43A047"),
        ),
        variant_lines=(
            "GT→GC / AG→AA (SS-KO)",
        ),
        modifier=Badge("—", color=None),
        level="T3–T5",
        outcome_range=(0.0, 2.0),
        outcome_color="#546E7A",
    ),
)


# ---------------------------------------------------------------------------
# Legend and footer
# ---------------------------------------------------------------------------

LEGEND = (
    LegendEntry("#43A047", "Permissive exonic (CAG / AAG)"),
    LegendEntry("#EF9A9A", "Non-permissive"),
    LegendEntry("#1B5E20", "+6 = T (no trap)"),
    LegendEntry("#EF5350", "+6 = G (mismatch)"),
)

STRATEGY = (
    "Strategy: Panels B & D establish splicing ceiling  →  "
    "Panels E & F apply controlled dial-down  →  "
    "Target: 80–90% splicing, 10–20% retention"
)

CALLOUT = (
    "Total: 44 constructs across 7 panels  ·  All VHH27  ·  "
    "v4 acceptor  ·  AVD002 backbone  ·  FLASH synthesis (AgeI/BsrGI)"
)


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

def build_fig167() -> Path:
    svg = render(
        title="Wave 2 VR4 Intron Retention — Design Panel Overview",
        subtitle=(
            "44 constructs · VHH27 (cross-species anti-ALPL) · "
            "v4 acceptor (CTAAC + GAATCC → PNSSPSPSSSSSQES scar)"
        ),
        tagline=(
            "Goal: validate >90% splicer · identify dial-down rules "
            "to reach 80–90% (10–20% retention)"
        ),
        pills_header="WAVE 1 DESIGN PRINCIPLES — APPLIED TO WAVE 2",
        pills=PILLS,
        col_headers=COL_HEADERS,
        rows=ROWS,
        target_zone=(80.0, 90.0),
        target_zone_label="▲ TARGET",
        legend_entries=LEGEND,
        strategy_annotation=STRATEGY,
        total_callout=CALLOUT,
    )
    out = FIGURE_BANK / "FIG167-wave2-vr4-ir-panel-overview.svg"
    out.write_text(svg, encoding="utf-8")
    print(f"  OK    {out.name}")
    return out


if __name__ == "__main__":
    build_fig167()
