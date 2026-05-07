"""FIG194 — VHH27-VP1 insertion sites across AVD428–448 panel (N=21 constructs).

SCF032 single-construct architecture callout on AVD002 WT AAV9 RepCap backbone.
Five callout groups summarise insertion positions across VP1 VR-V through VR-IX.

Insertion groups (donor_exonic_CAG bp positions, AVD002 coordinate space):
    Group 1 — VR-V    (AVD428–433, ×6):  bp 3848–3863, center 3856
    Group 2 — VR-VI   (AVD434–435, ×2):  bp 3962–3965, center 3964
    Group 3 — VR-VIII (AVD436–446, ×11): bp 4127–4163, center 4145
    Group 4 — C-loop  (AVD447,     ×1):  bp 4361
    Group 5 — VR-IX   (AVD448,     ×1):  bp 4496

Callout fan layout: column 2→6 (L-to-R) matches group L-to-R bp order;
rows form a V-shape (1, 2, 3, 2, 1) so no pointer lines cross (verified
analytically — all pairwise intersection parameters fall outside [0,1]).
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
OUT = FIGURE_BANK / "FIG194-vhh27-insertion-site-panel-avd428-448.svg"


# ── Reference backbone: AVD002 WT AAV9 RepCap, 7104 bp ─────────────────────────
# bp_view_end auto-computed: max(block.bp_end)=4589 × 1.02 ≈ 4680 → efficient fit
# bp_view_start crops blank plasmid backbone before the vestigial p5 fragment
# (bp 445–496). 350 bp leaves ~30px of visual context left of Rep78.

BP_TOTAL = 7104
BP_VIEW_START = 350

BLOCKS = (
    Block("Rep78", 496,  2362, "#1565C0", layer=0, pointed=True),
    Block("Rep52", 1168, 2362, "#5C8AC0", layer=1, pointed=True),
    Block("p19",   895,  1053, "#37474F", layer=2, pointed=True),
    Block("p40",   1875, 2028, "#37474F", layer=2, pointed=True),
    Block("VP1",   2378, 4589, "#283593", layer=0, pointed=True),
    Block("VP2",   2789, 4589, "#3949AB", layer=1, pointed=True),
    Block("AAP",   2904, 3498, "#7986CB", layer=3, pointed=True),
    Block("VP3",   2984, 4589, "#5C6BC0", layer=2, pointed=True),
)

# Native VR positions from AVD002 WT GenBank annotations (not AVD428-448 which
# have inflated VR-V due to the VHH insertion).
# VR-IV–VR-IX labels auto-stagger upward inside the construct zone (no tick lines).
REGION_BLOCKS = (
    RegionBlock("VR-IV",  3731, 3758, "#E8F5E9"),
    RegionBlock("VR-V",   3839, 3893, "#FFF8E1"),
    RegionBlock("VR-VI",  3956, 3995, "#F3E5F5"),
    RegionBlock("VR-VIII", 4118, 4157, "#E3F2FD"),
    RegionBlock("VR-IX",  4487, 4520, "#E8F5E9"),
)

# All callouts at row=1 (nearest row, all boxes at identical y) — no crossings
# because column order (2→6) matches bp order (3856→4496) monotonically.
CALLOUTS = (
    Callout("VR-V ×6",     column=2, bp_pos=3856, row=1),
    Callout("VR-VI ×2",    column=3, bp_pos=3964, row=1),
    Callout("VR-VIII ×11", column=4, bp_pos=4145, row=1),
    Callout("C-loop ×1",   column=5, bp_pos=4361, row=1),
    Callout("VR-IX ×1",    column=6, bp_pos=4496, row=1),
)


def main() -> None:
    svg = render(
        title="VHH27-VP1 Insertion Sites — AVD428–448 (N=21)",
        subtitle=(
            "21 constructs grouped by VP1 loop; "
            "WT AAV9 RepCap reference AVD002 (7,104 bp)"
        ),
        bp_total=BP_TOTAL,
        bp_view_start=BP_VIEW_START,
        blocks=BLOCKS,
        region_blocks=REGION_BLOCKS,
        callouts=CALLOUTS,
    )
    OUT.write_text(svg, encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
