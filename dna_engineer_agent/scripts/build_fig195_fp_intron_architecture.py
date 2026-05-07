"""FIG195 — FP Intron Architecture: VHH27 Splice Donor/Acceptor Design.

Schematic bp coordinates (proportional concept figure, not base-exact):
    VP1 exon (left)    0 –  75 bp   upstream VP1 context (exonic)
    Splice donor      75 – 115 bp   SD incl. GT dinucleotide
    VHH27 body       115 – 540 bp   GGGGS5 + VHH27 CDS + GGGGS1 (~425 bp)
    Splice acceptor  540 – 625 bp   PPT + branch point + AG
    VP1 exon (right) 625 – 720 bp   downstream VP1 context (exonic)

"Intron" region band spans SD → SA (bp 75–625) in light amber.

Callouts at row=1 (packing algorithm fans them to center):
    "GT donor"     → bp 77  (just inside SD, marks the splice donor GT)
    "AG acceptor"  → bp 623 (just inside SA, marks the splice acceptor AG)

No pointed blocks: these elements are intronic; directionality arrows would
imply stop-codon-like termination rather than an in-frame insertion.
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

BP_TOTAL = 720

BLOCKS = (
    Block("VP1",   0,   75, "#283593", layer=0),  # left exon context
    Block("SD",   75,  115, "#BF360C", layer=0),  # splice donor (GT)
    Block("VHH27",115,  540, "#00695C", layer=0),  # nanobody CDS within intron
    Block("SA",  540,  625, "#BF360C", layer=0),  # splice acceptor (AG)
    Block("VP1",  625,  720, "#283593", layer=0),  # right exon context
)

REGION_BLOCKS = (
    RegionBlock("Intron", 75, 625, "#FFF8E1"),
)

# Both callouts at row=1 — packing algorithm centers the group at mean bp_pos
# (~bp 350) and fans pointers outward: GT to the left (bp 77), AG to the
# right (bp 623). Column ordering (2 < 5) enforces left-to-right box sort.
CALLOUTS = (
    Callout("GT donor",    column=2, bp_pos=77,  row=1),
    Callout("AG acceptor", column=5, bp_pos=623, row=1),
)


def main() -> None:
    svg = render(
        title="FP Intron Architecture — VHH27 Splice Design",
        subtitle=(
            "Intronic VHH27 CDS flanked by splice donor (GT) and "
            "acceptor (AG) within VP1 exonic sequence"
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
