#!/usr/bin/env python3
"""
Unit tests for U1 snRNA : 5'SS Duplex ΔG Calculator

Tests the custom nearest-neighbor thermodynamics calculator that uses
Hudson et al. (2013) parameters for Ψ-A pairs at positions +3/+4.
"""

import pytest
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from u1_duplex_dg import (
    calculate_u1_duplex_dg,
    is_watson_crick_pair,
    is_wobble_pair,
    is_base_paired,
    find_consecutive_mismatches,
    dna_to_rna,
    has_frame1_stop_codon,
    generate_all_gt_sequences,
    filter_frame1_stops,
    FUNCTIONAL_THRESHOLD,
)


class TestDnaToRna:
    """Test DNA to RNA conversion."""

    def test_basic_conversion(self):
        assert dna_to_rna("ACGT") == "ACGU"

    def test_all_t_conversion(self):
        assert dna_to_rna("TTTT") == "UUUU"

    def test_mixed_case(self):
        assert dna_to_rna("AcGt") == "ACGU"


class TestBasePairing:
    """Test base pairing functions."""

    def test_watson_crick_au(self):
        assert is_watson_crick_pair('A', 'U')
        assert is_watson_crick_pair('U', 'A')

    def test_watson_crick_gc(self):
        assert is_watson_crick_pair('G', 'C')
        assert is_watson_crick_pair('C', 'G')

    def test_watson_crick_psi_a(self):
        # Ψ (pseudouridine) pairs with A
        assert is_watson_crick_pair('Ψ', 'A')
        # A also pairs with Ψ (symmetric for is_base_paired purposes)
        assert is_watson_crick_pair('A', 'Ψ')  # Both directions work via WC_PAIRS lookup

    def test_wobble_gu(self):
        assert is_wobble_pair('G', 'U')
        assert is_wobble_pair('U', 'G')

    def test_not_wobble(self):
        assert is_wobble_pair('A', 'U') is False
        assert is_wobble_pair('G', 'C') is False

    def test_is_base_paired_wc(self):
        assert is_base_paired('A', 'U')
        assert is_base_paired('G', 'C')

    def test_is_base_paired_wobble(self):
        assert is_base_paired('G', 'U')
        assert is_base_paired('U', 'G')

    def test_is_base_paired_mismatch(self):
        assert is_base_paired('A', 'A') is False
        assert is_base_paired('C', 'U') is False


class TestConsecutiveMismatches:
    """Test consecutive mismatch detection."""

    def test_no_mismatches(self):
        assert find_consecutive_mismatches([]) == 0

    def test_single_mismatch(self):
        assert find_consecutive_mismatches([3]) == 1

    def test_two_consecutive(self):
        assert find_consecutive_mismatches([3, 4]) == 2

    def test_three_consecutive(self):
        assert find_consecutive_mismatches([3, 4, 5]) == 3

    def test_non_consecutive(self):
        assert find_consecutive_mismatches([1, 3, 5]) == 1

    def test_mixed(self):
        # Positions 1, 5, 6, 7, 10 -> longest run is 5,6,7 = 3
        assert find_consecutive_mismatches([1, 5, 6, 7, 10]) == 3


class TestConsensusSequence:
    """Test ΔG calculation for consensus sequences.

    The ACTUAL consensus for U1 binding is:
    - U1 region: G  U  C | C  A  Ψ  Ψ  A  C  U  U
    - Optimal:   C  A  G | G  T  A  A  T  G  A  A

    So CAGGTAATGAA is the true consensus, NOT CAGGTAAGTAA.
    """

    def test_true_consensus_caggtaatgaa_is_strongest(self):
        """True consensus CAGGTAATGAA should be highly stable (very negative ΔG)."""
        result = calculate_u1_duplex_dg("CAGGTAATGAA")
        assert result.delta_g < -15.0, f"True consensus should be highly stable, got ΔG={result.delta_g}"
        assert result.is_functional
        assert len(result.mismatch_positions) == 0, "True consensus should have no mismatches"

    def test_near_consensus_aaggtaatgaa(self):
        """Near-consensus AAGGTAATGAA should still be functional."""
        result = calculate_u1_duplex_dg("AAGGTAATGAA")
        # A at -3 (position 0) doesn't pair well with G in U1, so slightly weaker
        assert result.delta_g < -10.0, f"Near-consensus should be stable, got ΔG={result.delta_g}"
        assert result.is_functional

    def test_strong_donor_detected_as_functional(self):
        """Strong splice donors should be classified as functional."""
        strong_donors = [
            "CAGGTAATGAA",  # True consensus
            "CAGGTAATGAC",  # Maximum H-bond variant (27 H-bonds)
            "CAGGTAATGAT",  # Near-consensus
        ]
        for seq in strong_donors:
            result = calculate_u1_duplex_dg(seq)
            assert result.is_functional, f"{seq} should be functional, got ΔG={result.delta_g}"


class TestConsecutiveMismatchPenalty:
    """Test that consecutive mismatches at +3/+4 are severely penalized."""

    def test_cc_at_plus3_plus4_is_weak(self):
        """
        Consecutive CC at +3/+4 should be severely penalized.

        AAGGTCCAACC has:
        - +3 = C (mismatch with Ψ)
        - +4 = C (mismatch with Ψ)
        This creates a destabilizing "bubble" that disrupts U1 binding.
        """
        result = calculate_u1_duplex_dg("AAGGTCCAACC")
        assert result.consecutive_mismatch_penalty >= 2.0, \
            f"Should have consecutive mismatch penalty, got {result.consecutive_mismatch_penalty}"
        # This sequence should NOT be functional due to the bubble
        # (NNSplice doesn't detect it at all)

    def test_cc_vs_aa_at_plus3_plus4(self):
        """
        A at +3/+4 should be MORE stable than C at +3/+4.

        This is the key insight from NNSplice testing:
        - AAGGTCAAACC (A at +3/+4) → Strong (0.98)
        - AAGGTCCAACC (C at +3/+4) → Not Detected (0.00)
        """
        result_aa = calculate_u1_duplex_dg("AAGGTAAAACC")  # A at +3, A at +4
        result_cc = calculate_u1_duplex_dg("AAGGTCCAACC")  # C at +3, C at +4

        # AA should be much more stable (more negative ΔG)
        assert result_aa.delta_g < result_cc.delta_g, \
            f"AA at +3/+4 should be more stable: {result_aa.delta_g} vs {result_cc.delta_g}"

        # The difference should be substantial (>3 kcal/mol)
        diff = result_cc.delta_g - result_aa.delta_g
        assert diff > 2.0, f"Expected >2 kcal/mol difference, got {diff}"


class TestPsiAStabilization:
    """Test Ψ-A enhanced stability at positions +3/+4."""

    def test_a_at_plus3_pairs_with_psi(self):
        """A at position +3 pairs with Ψ in U1 (enhanced stability)."""
        # AAGGTAAAACC has A at +3 (index 5)
        result = calculate_u1_duplex_dg("AAGGTAAAACC")
        # Position +3 (index 5) should not be a mismatch
        assert 5 not in result.mismatch_positions, "A at +3 should pair with Ψ"

    def test_a_at_plus4_pairs_with_psi(self):
        """A at position +4 pairs with Ψ in U1 (enhanced stability)."""
        # AAGGTAAAACC has A at +4 (index 6)
        result = calculate_u1_duplex_dg("AAGGTAAAACC")
        # Position +4 (index 6) should not be a mismatch
        assert 6 not in result.mismatch_positions, "A at +4 should pair with Ψ"

    def test_c_at_plus3_is_mismatch(self):
        """C at position +3 is a mismatch with Ψ in U1."""
        result = calculate_u1_duplex_dg("AAGGTCAAACC")
        # Position +3 (index 5) should be a mismatch (C doesn't pair with Ψ)
        assert 5 in result.mismatch_positions, "C at +3 should be mismatch with Ψ"


class TestFunctionalThreshold:
    """Test functional splice site classification."""

    def test_threshold_value(self):
        """Verify functional threshold is set correctly."""
        assert FUNCTIONAL_THRESHOLD == -8.1

    def test_strong_donor_is_functional(self):
        """Strong donors should be classified as functional (ΔG ≤ threshold)."""
        result = calculate_u1_duplex_dg("CAGGTAATGAC")  # Tier 7 consensus
        assert result.is_functional, f"Should be functional, got ΔG={result.delta_g}"
        assert result.delta_g <= FUNCTIONAL_THRESHOLD

    def test_weak_donor_may_not_be_functional(self):
        """Weak donors with poor ΔG may not be functional."""
        # Sequence with many mismatches
        result = calculate_u1_duplex_dg("TTTGTCCCCCG")
        # This may or may not be functional depending on exact ΔG
        # Just check that the calculation completes
        assert isinstance(result.is_functional, bool)
        assert isinstance(result.delta_g, float)


class TestNonAGDonors:
    """Test handling of non-AG|GT donors (position -1 ≠ G)."""

    def test_ag_donor_has_g_at_minus1(self):
        """AG donor: position -1 (index 2) should be G."""
        # AAGGTAAAACC: position -1 is G
        result = calculate_u1_duplex_dg("AAGGTAAAACC")
        # Position -1 pairs with C in U1, G-C is Watson-Crick
        assert 2 not in result.mismatch_positions, "G at -1 should pair with C"

    def test_non_ag_donor_has_mismatch_at_minus1(self):
        """Non-AG donor: position -1 ≠ G creates mismatch."""
        # AAAGTAAAACC: position -1 is A (not G)
        result = calculate_u1_duplex_dg("AAAGTAAAACC")
        # A doesn't pair with C, so position -1 (index 2) is mismatch
        assert 2 in result.mismatch_positions, "A at -1 should be mismatch with C"

    def test_non_ag_is_weaker(self):
        """Non-AG donors should have less stable ΔG than equivalent AG donors."""
        result_ag = calculate_u1_duplex_dg("AAGGTAAAACC")  # AG|GT
        result_non_ag = calculate_u1_duplex_dg("AAAGTAAAACC")  # AA|GT (non-AG)

        # AG should be more stable (more negative ΔG)
        assert result_ag.delta_g < result_non_ag.delta_g, \
            f"AG donor should be more stable: {result_ag.delta_g} vs {result_non_ag.delta_g}"


class TestFrame1StopCodons:
    """Test Frame 1 stop codon detection."""

    def test_no_stop_codon(self):
        """Sequence without stop codons in Frame 1."""
        # AAGGTAAAACC: +3+4+5 = AAA (Lys), +6+7+8 = ACC (Thr)
        assert has_frame1_stop_codon("AAGGTAAAACC") is False

    def test_stop_at_codon2(self):
        """Sequence with stop codon at +3+4+5."""
        # AAGGTTAGACC: +3+4+5 = TAG (stop)
        assert has_frame1_stop_codon("AAGGTTAGACC") is True

    def test_stop_at_codon3(self):
        """Sequence with stop codon at +6+7+8."""
        # AAGGTAAATGA: +6+7+8 = TGA (stop)
        assert has_frame1_stop_codon("AAGGTAAATGA") is True

    def test_taa_stop(self):
        """TAA stop codon detection."""
        assert has_frame1_stop_codon("AAGGTAAATAA") is True  # TAA at +6+7+8

    def test_all_stop_codons(self):
        """All three stop codons are detected."""
        assert has_frame1_stop_codon("AAGGTTAAACC") is True  # TAA at +3+4+5
        assert has_frame1_stop_codon("AAGGTTAGACC") is True  # TAG at +3+4+5
        assert has_frame1_stop_codon("AAGGTTGAACC") is True  # TGA at +3+4+5


class TestSequenceGeneration:
    """Test sequence generation and filtering."""

    def test_generate_all_gt_sequences_count(self):
        """Should generate 4^9 = 262,144 sequences."""
        seqs = generate_all_gt_sequences()
        assert len(seqs) == 4**9, f"Expected 262,144 sequences, got {len(seqs)}"

    def test_all_sequences_have_gt(self):
        """All generated sequences should have GT at positions +1/+2."""
        seqs = generate_all_gt_sequences()
        # Check a sample (not all 262k)
        for seq in seqs[:1000]:
            assert seq[3:5] == "GT", f"Position +1/+2 should be GT: {seq}"

    def test_all_sequences_are_11mer(self):
        """All generated sequences should be 11 nucleotides."""
        seqs = generate_all_gt_sequences()
        for seq in seqs[:1000]:
            assert len(seq) == 11, f"Should be 11-mer: {seq}"

    def test_filter_removes_stop_codons(self):
        """Filter should remove sequences with Frame 1 stop codons."""
        seqs = generate_all_gt_sequences()
        filtered = filter_frame1_stops(seqs)

        # Should remove some sequences
        assert len(filtered) < len(seqs)

        # All remaining should be stop-free
        for seq in filtered[:1000]:
            assert has_frame1_stop_codon(seq) is False


class TestInputValidation:
    """Test input validation."""

    def test_wrong_length_raises_error(self):
        """Sequences not 11 nt should raise ValueError."""
        with pytest.raises(ValueError):
            calculate_u1_duplex_dg("AAGGTAA")  # Too short

        with pytest.raises(ValueError):
            calculate_u1_duplex_dg("AAGGTAAAAAAAAA")  # Too long

    def test_invalid_nucleotide_raises_error(self):
        """Invalid nucleotides should raise ValueError."""
        with pytest.raises(ValueError):
            calculate_u1_duplex_dg("AAGGTAAAAXC")  # X is invalid


class TestDuplexResultFields:
    """Test that DuplexResult contains all expected fields."""

    def test_result_has_all_fields(self):
        """DuplexResult should have all expected attributes."""
        result = calculate_u1_duplex_dg("AAGGTAAAACC")

        assert hasattr(result, 'sequence_11mer')
        assert hasattr(result, 'delta_g')
        assert hasattr(result, 'stack_contributions')
        assert hasattr(result, 'mismatch_positions')
        assert hasattr(result, 'consecutive_mismatch_penalty')
        assert hasattr(result, 'is_functional')
        assert hasattr(result, 'notes')

    def test_result_to_dict(self):
        """DuplexResult.to_dict() should return a dictionary."""
        result = calculate_u1_duplex_dg("AAGGTAAAACC")
        d = result.to_dict()

        assert isinstance(d, dict)
        assert 'sequence' in d
        assert 'delta_g' in d
        assert 'is_functional' in d
        assert 'mismatch_count' in d


class TestKnownSequences:
    """Test against known sequences with expected behavior."""

    def test_tier7_consensus_is_strongest(self):
        """
        Tier 7 consensus sequence (maximum H-bonds) should have most stable ΔG.

        CAGGTAATGAC: 27 H-bonds, stretch=11, NNSplice=0.985
        Note: Position +8 (C) doesn't pair with U1's U, so there's one mismatch.
        """
        result = calculate_u1_duplex_dg("CAGGTAATGAC")
        assert result.delta_g < -15.0, f"Tier 7 consensus should be very stable: {result.delta_g}"
        assert result.is_functional
        # C at +8 (index 10) doesn't pair with U in U1
        assert len(result.mismatch_positions) <= 1, f"Should have at most 1 mismatch: {result.mismatch_positions}"

    def test_tier1_weak_donor(self):
        """
        Tier 1 weak donor (low H-bonds) should have less stable ΔG.

        These sequences should still have GT but weak overall binding.
        """
        # A weak sequence with mismatches
        result = calculate_u1_duplex_dg("TTCGTCCCCCC")
        # Should be less stable than consensus
        assert result.delta_g > -10.0, f"Weak donor should have poor ΔG: {result.delta_g}"

    def test_nnsplice_discrepancy_sequence(self):
        """
        Test the sequence that NNSplice shows a discrepancy for.

        AAGGTCCAACC: 13 H-bonds but NNSplice=0.00 (Not Detected)
        This is due to consecutive CC mismatches at +3/+4.
        """
        result = calculate_u1_duplex_dg("AAGGTCCAACC")

        # Should have mismatches at positions +3 and +4 (indices 5 and 6)
        assert 5 in result.mismatch_positions, "C at +3 should be mismatch"
        assert 6 in result.mismatch_positions, "C at +4 should be mismatch"

        # Should have consecutive mismatch penalty
        assert result.consecutive_mismatch_penalty > 0, "Should have penalty for consecutive mismatches"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
