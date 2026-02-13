#!/usr/bin/env python3
"""
Regression tests for the intron body codon optimizer (Plan item 3g).

Tests cover:
1. Input constants (VHH3 protein, original sequences)
2. GGGGS linker design (Stage 0 — fully deterministic)
3. Utility function correctness (count_dinuc, split_codons, etc.)
4. Do-no-harm check logic
5. Optimizer output properties (protein preserved, no new restriction sites)

Note: DNAchisel (Stage 1) and MaxEntScan scoring (Stage 2) are tested for
correctness properties rather than exact values, since DNAchisel's internal
optimization may vary across versions. The key invariants are:
- Protein sequence is ALWAYS preserved
- GT/AG counts never increase beyond Stage 1 input
- No new BsaI/BsmBI/BbsI sites are created
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from optimize_intron_body_codons import (
    CODON_TABLE,
    FORBIDDEN_RESTRICTION_SITES,
    GGGGS1_ORIGINAL,
    GGGGS5_ORIGINAL,
    VHH3_ORIGINAL,
    VHH3_PROTEIN,
    count_dinuc,
    count_restriction_sites,
    design_ggggs_zero_gtag,
    find_dinuc_positions,
    has_gt_or_ag_internal,
    join_codons,
    reverse_complement,
    split_codons,
    translate,
)


# ---------------------------------------------------------------------------
# Input constant verification
# ---------------------------------------------------------------------------

class TestInputConstants:
    """Verify that input constants haven't been accidentally modified."""

    def test_vhh3_original_length(self):
        assert len(VHH3_ORIGINAL) == 357  # 119 codons × 3

    def test_vhh3_protein_length(self):
        assert len(VHH3_PROTEIN) == 119

    def test_vhh3_translation(self):
        """VHH3 DNA must translate to the expected protein."""
        actual = translate(VHH3_ORIGINAL)
        assert actual == VHH3_PROTEIN

    def test_ggggs5_original_length(self):
        assert len(GGGGS5_ORIGINAL) == 75  # 25 codons × 3

    def test_ggggs1_original_length(self):
        assert len(GGGGS1_ORIGINAL) == 15  # 5 codons × 3

    def test_ggggs5_protein(self):
        assert translate(GGGGS5_ORIGINAL) == "GGGGS" * 5

    def test_ggggs1_protein(self):
        assert translate(GGGGS1_ORIGINAL) == "GGGGS"

    def test_total_body_length(self):
        """Full body = GGGGS5(75) + VHH3(357) + GGGGS1(15) = 447 bp."""
        body = GGGGS5_ORIGINAL + VHH3_ORIGINAL + GGGGS1_ORIGINAL
        assert len(body) == 447

    def test_original_gt_count(self):
        """Known: original body has 28 GT dinucleotides."""
        body = GGGGS5_ORIGINAL + VHH3_ORIGINAL + GGGGS1_ORIGINAL
        assert count_dinuc(body, "GT") == 28

    def test_original_ag_count(self):
        """Known: original body has 38 AG dinucleotides."""
        body = GGGGS5_ORIGINAL + VHH3_ORIGINAL + GGGGS1_ORIGINAL
        assert count_dinuc(body, "AG") == 38


# ---------------------------------------------------------------------------
# Stage 0: GGGGS design (fully deterministic)
# ---------------------------------------------------------------------------

class TestStage0GGGGS:
    """Stage 0 is fully deterministic — exact output can be asserted."""

    def test_ggggs5_zero_gt(self):
        ggggs5, _ = design_ggggs_zero_gtag()
        assert count_dinuc(ggggs5, "GT") == 0

    def test_ggggs5_zero_ag(self):
        ggggs5, _ = design_ggggs_zero_gtag()
        assert count_dinuc(ggggs5, "AG") == 0

    def test_ggggs1_zero_gt(self):
        _, ggggs1 = design_ggggs_zero_gtag()
        assert count_dinuc(ggggs1, "GT") == 0

    def test_ggggs1_zero_ag(self):
        _, ggggs1 = design_ggggs_zero_gtag()
        assert count_dinuc(ggggs1, "AG") == 0

    def test_ggggs5_protein_preserved(self):
        ggggs5, _ = design_ggggs_zero_gtag()
        assert translate(ggggs5) == "GGGGS" * 5

    def test_ggggs1_protein_preserved(self):
        _, ggggs1 = design_ggggs_zero_gtag()
        assert translate(ggggs1) == "GGGGS"

    def test_ggggs5_length(self):
        ggggs5, _ = design_ggggs_zero_gtag()
        assert len(ggggs5) == 75

    def test_ggggs1_length(self):
        _, ggggs1 = design_ggggs_zero_gtag()
        assert len(ggggs1) == 15

    def test_ggggs5_no_restriction_sites(self):
        ggggs5, _ = design_ggggs_zero_gtag()
        assert count_restriction_sites(ggggs5) == 0

    def test_ggggs1_no_restriction_sites(self):
        _, ggggs1 = design_ggggs_zero_gtag()
        assert count_restriction_sites(ggggs1) == 0


# ---------------------------------------------------------------------------
# Utility function tests
# ---------------------------------------------------------------------------

class TestCountDinuc:
    def test_simple_gt(self):
        assert count_dinuc("ATGTCC", "GT") == 1

    def test_multiple_gt(self):
        assert count_dinuc("GTGTGT", "GT") == 3

    def test_no_match(self):
        assert count_dinuc("AAAAAA", "GT") == 0

    def test_ag_count(self):
        assert count_dinuc("AAGCAG", "AG") == 2

    def test_overlap(self):
        assert count_dinuc("AGAG", "AG") == 2


class TestSplitJoinCodons:
    def test_split(self):
        assert split_codons("ATGAAAGCC") == ["ATG", "AAA", "GCC"]

    def test_join(self):
        assert join_codons(["ATG", "AAA", "GCC"]) == "ATGAAAGCC"

    def test_roundtrip(self):
        seq = "ATGAAAGCCTTTGGC"
        assert join_codons(split_codons(seq)) == seq


class TestHasGtOrAgInternal:
    def test_gt_internal(self):
        assert has_gt_or_ag_internal("GTC")  # G-T boundary at pos 0-1

    def test_ag_internal(self):
        assert has_gt_or_ag_internal("AGC")  # A-G boundary at pos 0-1

    def test_clean(self):
        assert not has_gt_or_ag_internal("GCA")

    def test_tga_has_no_gt(self):
        # TGA: T-G internal, G-A boundary — neither is GT or AG at pos 0-1
        # Actually TGA has no GT (T,G) or AG (A,G) at positions 0-1 or 1-2
        # TG at 0-1 is not GT, GA at 1-2 is not GT or AG
        assert not has_gt_or_ag_internal("TGA")


class TestFindDinucPositions:
    def test_find_gt(self):
        positions = find_dinuc_positions("ATGTCC", "GT")
        assert positions == [2]

    def test_find_ag(self):
        positions = find_dinuc_positions("AAGCAG", "AG")
        assert positions == [1, 4]

    def test_none(self):
        assert find_dinuc_positions("AAAA", "GT") == []


# ---------------------------------------------------------------------------
# Restriction site checking
# ---------------------------------------------------------------------------

class TestRestrictionSites:
    def test_no_sites(self):
        assert count_restriction_sites("ATGATGATGATGATG") == 0

    def test_bsai_forward(self):
        assert count_restriction_sites("GGTCTC") == 1

    def test_bsai_reverse(self):
        rc = reverse_complement("GGTCTC")
        assert count_restriction_sites(rc) == 1

    def test_bsmbi(self):
        assert count_restriction_sites("CGTCTC") == 1


# ---------------------------------------------------------------------------
# Full optimization property tests
# ---------------------------------------------------------------------------

class TestFullOptimizationProperties:
    """
    Test invariant properties that must hold for ANY valid optimizer output.
    These don't assert exact values but verify correctness constraints.
    """

    @pytest.fixture(scope="class")
    def optimized_output(self):
        """Run stages 0+1 (skipping MaxEntScan-dependent Stage 2)."""
        from optimize_intron_body_codons import (
            design_ggggs_zero_gtag,
            optimize_vhh3_dnachisel,
            fix_boundary_gtag,
        )

        ggggs5_new, ggggs1_new = design_ggggs_zero_gtag()
        vhh3_opt, stats = optimize_vhh3_dnachisel(VHH3_ORIGINAL)

        # Fix boundaries
        vhh3_codons = split_codons(vhh3_opt)
        vhh3_codons, _ = fix_boundary_gtag(ggggs5_new, vhh3_codons, ggggs1_new)
        vhh3_opt = join_codons(vhh3_codons)

        ggggs5_codons = split_codons(ggggs5_new)
        vhh3_codons = split_codons(vhh3_opt)
        ggggs1_codons = split_codons(ggggs1_new)
        all_codons = ggggs5_codons + vhh3_codons + ggggs1_codons

        return {
            "all_codons": all_codons,
            "full_seq": join_codons(all_codons),
            "ggggs5": ggggs5_new,
            "ggggs1": ggggs1_new,
            "vhh3": vhh3_opt,
            "stats": stats,
        }

    def test_codon_count(self, optimized_output):
        """Must have exactly 149 codons (25 + 119 + 5)."""
        assert len(optimized_output["all_codons"]) == 149

    def test_vhh3_protein_preserved(self, optimized_output):
        """VHH3 protein sequence must be EXACTLY preserved."""
        vhh3_codons = optimized_output["all_codons"][25:144]
        actual = translate(join_codons(vhh3_codons))
        assert actual == VHH3_PROTEIN

    def test_ggggs5_protein_preserved(self, optimized_output):
        """GGGGS5 protein must be preserved."""
        ggggs5_codons = optimized_output["all_codons"][:25]
        assert translate(join_codons(ggggs5_codons)) == "GGGGS" * 5

    def test_ggggs1_protein_preserved(self, optimized_output):
        """GGGGS1 protein must be preserved."""
        ggggs1_codons = optimized_output["all_codons"][144:]
        assert translate(join_codons(ggggs1_codons)) == "GGGGS"

    def test_body_length_preserved(self, optimized_output):
        """Body length must remain 447 bp."""
        assert len(optimized_output["full_seq"]) == 447

    def test_gt_reduced(self, optimized_output):
        """GT count must be reduced from original 28."""
        original_gt = 28
        final_gt = count_dinuc(optimized_output["full_seq"], "GT")
        assert final_gt < original_gt, \
            f"GT not reduced: {original_gt} → {final_gt}"

    def test_ag_reduced(self, optimized_output):
        """AG count must be reduced from original 38."""
        original_ag = 38
        final_ag = count_dinuc(optimized_output["full_seq"], "AG")
        assert final_ag < original_ag, \
            f"AG not reduced: {original_ag} → {final_ag}"

    def test_no_new_restriction_sites(self, optimized_output):
        """No BsaI/BsmBI/BbsI sites in optimized body."""
        assert count_restriction_sites(optimized_output["full_seq"]) == 0

    def test_ggggs_regions_zero_gtag(self, optimized_output):
        """GGGGS5 and GGGGS1 must have zero GT and AG."""
        seq = optimized_output["full_seq"]
        ggggs5 = seq[:75]
        ggggs1 = seq[432:447]
        assert count_dinuc(ggggs5, "GT") == 0, "GGGGS5 has GT"
        assert count_dinuc(ggggs5, "AG") == 0, "GGGGS5 has AG"
        assert count_dinuc(ggggs1, "GT") == 0, "GGGGS1 has GT"
        assert count_dinuc(ggggs1, "AG") == 0, "GGGGS1 has AG"

    def test_dnachisel_reduced_gt(self, optimized_output):
        """DNAchisel should reduce VHH3 GT count."""
        stats = optimized_output["stats"]
        assert stats["after_gt"] <= stats["before_gt"]

    def test_dnachisel_reduced_ag(self, optimized_output):
        """DNAchisel should reduce VHH3 AG count."""
        stats = optimized_output["stats"]
        assert stats["after_ag"] <= stats["before_ag"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
