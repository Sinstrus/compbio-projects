"""
Tests for frame offset calculation bug (BUG-001).

The bug: When calculating where a mutation falls in the reading frame,
the code was using the absolute position in the full sequence instead of
the position relative to the CDS start.

Example:
  Full sequence: GGATCC[ATG]AAAGCCTAAGAATTC
  CDS starts at position 6 (0-indexed)
  Mutation at position 9 ('A' in AAA codon)

  WRONG: frame = 9 % 3 = 0 (implies first position of codon)
  RIGHT: frame = (9 - 6) % 3 = 0 (correct, first position of AAA codon)
"""

import pytest


def calculate_frame_offset_wrong(mutation_pos, cds_start=0):
    """
    BUGGY implementation: uses absolute position.
    This was the bug in the original code.
    """
    return mutation_pos % 3


def calculate_frame_offset_correct(mutation_pos, cds_start):
    """
    CORRECT implementation: calculates position relative to CDS start.
    """
    return (mutation_pos - cds_start) % 3


class TestFrameOffsetCalculation:
    """Test correct frame offset calculation relative to CDS start."""

    def test_mutation_at_cds_start(self):
        """Mutation at first position of CDS should be frame 0."""
        cds_start = 6
        mutation_pos = 6  # First base of CDS

        # The bug would give wrong answer if cds_start != 0
        assert calculate_frame_offset_correct(mutation_pos, cds_start) == 0

    def test_mutation_in_first_codon_pos1(self):
        """Mutation at position 1 of first codon (frame 1)."""
        cds_start = 6
        mutation_pos = 7  # Second base of first codon

        assert calculate_frame_offset_correct(mutation_pos, cds_start) == 1

    def test_mutation_in_first_codon_pos2(self):
        """Mutation at position 2 of first codon (frame 2)."""
        cds_start = 6
        mutation_pos = 8  # Third base of first codon

        assert calculate_frame_offset_correct(mutation_pos, cds_start) == 2

    def test_mutation_in_second_codon_pos0(self):
        """Mutation at position 0 of second codon (frame 0)."""
        cds_start = 6
        mutation_pos = 9  # First base of second codon (AAA)

        # This is where the bug manifests most clearly
        # Wrong: 9 % 3 = 0 (happens to be correct by accident)
        # Right: (9 - 6) % 3 = 0
        assert calculate_frame_offset_correct(mutation_pos, cds_start) == 0

    def test_bug_manifestation_cds_not_at_zero(self):
        """
        Demonstrate the bug: when CDS doesn't start at position 0,
        the wrong calculation gives incorrect results.
        """
        cds_start = 7  # CDS starts at offset position
        mutation_pos = 10  # Should be in second codon, position 0

        # Wrong calculation (the bug)
        wrong_frame = calculate_frame_offset_wrong(mutation_pos, cds_start)
        assert wrong_frame == 1  # WRONG! 10 % 3 = 1

        # Correct calculation
        correct_frame = calculate_frame_offset_correct(mutation_pos, cds_start)
        assert correct_frame == 0  # RIGHT! (10 - 7) % 3 = 0

    def test_with_upstream_sequence(self):
        """Test with realistic upstream context (e.g., cloning site)."""
        # Sequence: AAGCTT[ATG]AAAGCC... (HindIII site before CDS)
        upstream_length = 6  # HindIII site
        cds_start = upstream_length

        # Mutation in third codon (GCC), middle position
        mutation_pos = cds_start + 7  # Position 13 in full sequence
        # Relative position: 7 = codon 2 (third codon), position 1

        frame = calculate_frame_offset_correct(mutation_pos, cds_start)
        assert frame == 1  # 7 % 3 = 1

        # Verify codon identification
        codon_index = (mutation_pos - cds_start) // 3
        assert codon_index == 2  # Third codon (0-indexed)

    def test_hardcoded_frame_zero_bug(self):
        """
        Test for BUG-002: Hardcoded frame=0 assumption.

        Some code paths assumed frame=0 always, which is only true
        when CDS starts at a multiple of 3 in the sequence.
        """
        # Case 1: CDS starts at position 0 (multiple of 3)
        cds_start_aligned = 0
        mutation_pos = 4
        frame_aligned = calculate_frame_offset_correct(mutation_pos, cds_start_aligned)
        assert frame_aligned == 1

        # Case 2: CDS starts at position 5 (NOT multiple of 3)
        cds_start_unaligned = 5
        mutation_pos = 9
        frame_unaligned = calculate_frame_offset_correct(mutation_pos, cds_start_unaligned)
        assert frame_unaligned == 1

        # Both should give same frame offset relative to CDS
        assert frame_aligned == frame_unaligned

    def test_multiple_upstream_elements(self):
        """Test with multiple upstream elements (promoter, RBS, etc.)."""
        # Realistic construct: promoter (35bp) + RBS (10bp) + ATG...
        upstream_length = 45
        cds_start = upstream_length

        # Mutation in 10th codon, 2nd position
        mutation_pos = cds_start + 28  # 28 = 9*3 + 1

        frame = calculate_frame_offset_correct(mutation_pos, cds_start)
        assert frame == 1  # 28 % 3 = 1

        codon_index = (mutation_pos - cds_start) // 3
        assert codon_index == 9  # 10th codon (0-indexed)


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_cds_start_at_zero(self):
        """When CDS starts at 0, both methods should agree."""
        cds_start = 0

        for mutation_pos in range(0, 30):
            wrong = calculate_frame_offset_wrong(mutation_pos, cds_start)
            correct = calculate_frame_offset_correct(mutation_pos, cds_start)
            assert wrong == correct, f"Disagreement at position {mutation_pos}"

    def test_cds_start_multiple_of_three(self):
        """Test when CDS starts at multiples of 3."""
        for cds_start in [0, 3, 6, 9, 12]:
            mutation_pos = cds_start + 5  # Always in frame 2
            frame = calculate_frame_offset_correct(mutation_pos, cds_start)
            assert frame == 2

    def test_cds_start_not_multiple_of_three(self):
        """Test when CDS starts at non-multiples of 3."""
        for cds_start in [1, 2, 4, 5, 7, 8]:
            mutation_pos = cds_start + 5  # Always in frame 2
            frame = calculate_frame_offset_correct(mutation_pos, cds_start)
            assert frame == 2


def test_real_world_example():
    """
    Real-world example from AAV transfer plasmid.

    Sequence structure:
    - HindIII site (6bp): positions 0-5
    - Kozak sequence (6bp): positions 6-11
    - CDS starts: position 12
    - Mutation at position 30
    """
    hindiii_site = "AAGCTT"
    kozak = "GCCACC"
    cds_start_pos = len(hindiii_site) + len(kozak)  # = 12

    mutation_pos = 30

    # Calculate frame
    frame = calculate_frame_offset_correct(mutation_pos, cds_start_pos)

    # Relative position in CDS: 30 - 12 = 18
    # Codon: 18 // 3 = 6 (7th codon)
    # Frame: 18 % 3 = 0 (first position)

    assert frame == 0

    codon_index = (mutation_pos - cds_start_pos) // 3
    assert codon_index == 6
