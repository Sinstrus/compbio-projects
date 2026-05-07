"""FIG195 — FP Intron Architecture: VHH27 Splice Donor/Acceptor Design (detailed).

Conceptual bp coordinates (proportional concept figure; actual element sizes preserved
where they dominate — GGGGS5 = 75 bp, VHH27 = 363 bp; small elements widened to be legible).

    VP1 exon (left)    0 – 100  upstream VP1 context (exonic)
    CAG exonic pad   100 – 125  last exonic codon before SD; VP1 aa456
    Splice donor     125 – 160  GT dinucleotide + donor-intronic spacer
    GGGGS5 linker    160 – 235  5× GGGGS repeats (actual 75 bp, GT/AG-depleted)
    VHH27 CDS        235 – 598  anti-ALPL nanobody (actual 363 bp)
    GGGGS1 linker    598 – 623  1× GGGGS repeat (actual 15 bp, GT/AG-depleted)
    Splice acceptor  623 – 668  β-globin PPT + branch point + AG (actual 38 bp)
    HBB exon3 pad    668 – 693  first exonic codon after SA; VP1 aa457 (actual 6 bp)
    VP1 exon (right) 693 – 793  downstream VP1 context (exonic)

Intron region band: 125–668 (SD start → SA end), light amber.

Callouts (4 total, row=1):
    column=1  "aa 456"       → bp 112  CAG pad center (last VP1 aa before SD)
    column=2  "GT donor"     → bp 142  SD block center
    column=5  "AG acceptor"  → bp 645  SA block center
    column=6  "aa 457"       → bp 680  HBB pad center (first VP1 aa after SA)

AVD313 coordinate reference (0-based positions in 7614-bp construct):
    VP1 ATG:          2378
    donor_exonic_CAG: 3743–3746  offset = 1365 nt = codon 456 of VP1
    intron_FP_Tier7:  3746–4247  (GT through bglobin_acceptor AG)
    HBB_exon3_pad:    4247–4253  encodes VP1 aa457–458 in spliced product
    GGGGS5:           3755–3830  (75 bp)
    VHH27_anti-ALPL:  3830–4193  (363 bp)
    GGGGS1:           4193–4208  (15 bp)
    bglobin_acceptor: 4209–4247  (38 bp)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path.home() / "projects"))

from figure_bank.scaffolds.SCF032_construct_callout import (
    Block,
    Callout,
    RegionBlock,
    render,
)

FIGURE_BANK = Path.home() / "projects" / "figure_bank"
OUT = FIGURE_BANK / "FIG195-fp-intron-architecture-vhh27.svg"

BP_TOTAL = 793

BLOCKS = (
    Block("VP1",    0,   100, "#283593", layer=0),   # left exon context
    Block("CAG",  100,   125, "#1565C0", layer=0),   # exonic pad (VP1 aa456); part of SD signal
    Block("SD",   125,   160, "#BF360C", layer=0),   # GT dinucleotide + donor-intronic spacer
    Block("GGGGS5", 160, 235, "#2E7D32", layer=0),   # 5× GGGGS linker (75 bp actual)
    Block("VHH27", 235,  598, "#00695C", layer=0),   # anti-ALPL nanobody CDS (363 bp actual)
    Block("GS1",  598,   623, "#388E3C", layer=0),   # 1× GGGGS linker (15 bp actual)
    Block("SA",   623,   668, "#BF360C", layer=0),   # β-globin splice acceptor (PPT + BP + AG)
    Block("HBB",  668,   693, "#4A148C", layer=0),   # HBB exon3 pad (VP1 aa457)
    Block("VP1",  693,   793, "#283593", layer=0),   # right exon context
)

REGION_BLOCKS = (
    RegionBlock("Intron", 125, 668, "#FFF8E1"),
)

# Ideal-first packing: each box starts at bp_px(bp_pos) and spreads only when
# necessary to resolve overlap. Column order (1 < 2 < 5 < 6) enforces L-to-R sort.
CALLOUTS = (
    Callout("aa 456",      column=1, bp_pos=112, row=1),
    Callout("GT donor",    column=2, bp_pos=142, row=1),
    Callout("AG acceptor", column=5, bp_pos=645, row=1),
    Callout("aa 457",      column=6, bp_pos=680, row=1),
)


def main() -> None:
    svg = render(
        title="FP Intron Architecture — VHH27 Splice Design (AVD313)",
        subtitle=(
            "GGGGS5 × VHH27 × GGGGS1 intronic cassette "
            "flanked by VP1 aa456 (CAG pad) and aa457 (HBB pad)"
        ),
        bp_total=BP_TOTAL,
        blocks=BLOCKS,
        region_blocks=REGION_BLOCKS,
        callouts=CALLOUTS,
    )
    OUT.write_text(svg, encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
