#!/usr/bin/env python3
"""
Unit tests for splice_donor_hbonds module.
"""

import pytest
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from splice_donor_hbonds import (
    calculate_hbonds,
    get_pair_hbonds,
    find_longest_continuous_stretch,
    is_valid_splice_site,
    generate_all_valid_sequences,
    generate_sequences_for_hbond_count,
    validate_known_sequences,
    has_frame1_stop_codon,
    get_frame1_codons,
    translate_codon,
    translate_frame1_codons,
    generate_tiered_panel,
    select_diverse_within_tier,
    calculate_sequence_diversity,
    has_consensus_motif,
    POSITION_LABELS,
    U1_PARTNERS,
    OPTIMAL_DNA,
    MAX_HBONDS,
    STOP_CODONS,
    DEFAULT_TIERS,
)


class TestGetPairHbonds:
    """Tests for get_pair_hbonds function."""

    def test_watson_crick_gc(self):
        """G-C and C-G pairs should have 3 H-bonds."""
        assert get_pair_hbonds("G", "C") == 3
        assert get_pair_hbonds("C", "G") == 3

    def test_watson_crick_au(self):
        """A-U and U-A pairs should have 2 H-bonds."""
        assert get_pair_hbonds("A", "U") == 2
        assert get_pair_hbonds("U", "A") == 2

    def test_wobble_gu(self):
        """G-U and U-G wobble pairs should have 2 H-bonds."""
        assert get_pair_hbonds("G", "U") == 2
        assert get_pair_hbonds("U", "G") == 2

    def test_pseudouridine_pairs(self):
        """Ψ pairs like U, so G-Ψ and A-Ψ should have 2 H-bonds."""
        assert get_pair_hbonds("G", "Ψ") == 2
        assert get_pair_hbonds("A", "Ψ") == 2

    def test_mismatches(self):
        """Mismatches should have 0 H-bonds."""
        assert get_pair_hbonds("A", "A") == 0
        assert get_pair_hbonds("C", "C") == 0
        assert get_pair_hbonds("A", "G") == 0
        assert get_pair_hbonds("C", "U") == 0
        assert get_pair_hbonds("U", "Ψ") == 0  # U-Ψ is not a standard pair
        assert get_pair_hbonds("C", "Ψ") == 0


class TestFindLongestContinuousStretch:
    """Tests for find_longest_continuous_stretch function."""

    def test_all_paired(self):
        """All positions paired should return full length."""
        hbonds = [3, 2, 3, 3, 2, 2, 2, 2, 3, 2, 3]
        length, start = find_longest_continuous_stretch(hbonds)
        assert length == 11
        assert start == 0

    def test_single_gap(self):
        """Gap in the middle should split the stretch."""
        hbonds = [3, 2, 3, 3, 2, 0, 2, 2, 3, 2, 3]
        length, start = find_longest_continuous_stretch(hbonds)
        assert length == 5  # Positions 0-4
        assert start == 0

    def test_gap_at_start(self):
        """Gap at start should find stretch after."""
        hbonds = [0, 2, 3, 3, 2, 2, 2, 2, 3, 2, 3]
        length, start = find_longest_continuous_stretch(hbonds)
        assert length == 10
        assert start == 1

    def test_gap_at_end(self):
        """Gap at end should find stretch before."""
        hbonds = [3, 2, 3, 3, 2, 2, 2, 2, 3, 2, 0]
        length, start = find_longest_continuous_stretch(hbonds)
        assert length == 10
        assert start == 0

    def test_multiple_gaps(self):
        """Multiple gaps should find longest stretch."""
        hbonds = [3, 2, 0, 3, 2, 2, 2, 0, 3, 2, 3]
        length, start = find_longest_continuous_stretch(hbonds)
        assert length == 4  # Positions 3-6
        assert start == 3

    def test_no_pairs(self):
        """No base pairs should return 0."""
        hbonds = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        length, start = find_longest_continuous_stretch(hbonds)
        assert length == 0
        assert start == -1


class TestCalculateHbonds:
    """Tests for calculate_hbonds function."""

    def test_optimal_sequence(self):
        """Optimal sequence should give maximum H-bonds."""
        # Optimal: C A G G T A A T G A C
        optimal_seq = "".join(OPTIMAL_DNA)
        result = calculate_hbonds(optimal_seq)
        assert result.total_hbonds == MAX_HBONDS
        assert result.num_base_pairs == 11
        assert result.longest_stretch == 11

    def test_consensus_sequence(self):
        """Consensus CAG|GTAAGTAT should have reasonable H-bonds."""
        result = calculate_hbonds("CAGGTAAGTAT")
        # Calculate expected:
        # -3: C vs G -> C-G = 3
        # -2: A vs U -> A-U = 2
        # -1: G vs C -> G-C = 3
        # +1: G vs C -> G-C = 3
        # +2: T(U) vs A -> U-A = 2
        # +3: A vs Ψ -> A-Ψ = 2
        # +4: A vs Ψ -> A-Ψ = 2
        # +5: G vs A -> mismatch = 0
        # +6: T(U) vs C -> mismatch = 0
        # +7: A vs U -> A-U = 2
        # +8: T(U) vs G -> U-G = 2
        # Total: 3+2+3+3+2+2+2+0+0+2+2 = 21
        assert result.total_hbonds == 21
        assert result.sequence == "CAGGTAAGTAT"

    def test_minimum_gt_only(self):
        """Minimal sequence with only GT should have low H-bonds."""
        # T(U) pairs: T-G wobble = 2, T-U = 0, T-C = 0, T-A = 2, T-Ψ = 0
        # Position breakdown for TTTGTTTTTAT:
        # -3: T(U)-G = 2 (wobble), -2: T(U)-U = 0, -1: T(U)-C = 0
        # +1: G-C = 3, +2: T(U)-A = 2
        # +3: T(U)-Ψ = 0, +4: T(U)-Ψ = 0, +5: T(U)-A = 2
        # +6: T(U)-C = 0, +7: A-U = 2, +8: T(U)-G = 2 (wobble)
        # Total: 2+0+0+3+2+0+0+2+0+2+2 = 13
        result = calculate_hbonds("TTTGTTTTTAT")
        assert result.total_hbonds == 13

    def test_invalid_length_short(self):
        """Too short sequence should raise ValueError."""
        with pytest.raises(ValueError, match="must be 11 nucleotides"):
            calculate_hbonds("CAGGTAAGT")  # Only 9 nt

    def test_invalid_length_long(self):
        """Too long sequence should raise ValueError."""
        with pytest.raises(ValueError, match="must be 11 nucleotides"):
            calculate_hbonds("CAGGTAAGTATAA")  # 13 nt

    def test_invalid_character(self):
        """Invalid nucleotide should raise ValueError."""
        with pytest.raises(ValueError, match="Invalid nucleotide"):
            calculate_hbonds("CAGGTAXGTAT")  # X is invalid

    def test_lowercase_accepted(self):
        """Lowercase input should be accepted and converted."""
        result = calculate_hbonds("caggtaagtat")
        assert result.sequence == "CAGGTAAGTAT"
        assert result.total_hbonds == 21

    def test_per_position_values(self):
        """Check per-position H-bond values."""
        result = calculate_hbonds("CAGGTAAGTGC")
        assert len(result.per_position) == 11
        # -3: C-G = 3, -2: A-U = 2, -1: G-C = 3
        assert result.per_position[0] == 3
        assert result.per_position[1] == 2
        assert result.per_position[2] == 3


class TestIsValidSpliceSite:
    """Tests for is_valid_splice_site function."""

    def test_valid_gt(self):
        """Sequence with GT at +1/+2 should be valid."""
        assert is_valid_splice_site("CAGGTAAGTAT") is True
        assert is_valid_splice_site("AAAGTAAAAAT") is True

    def test_invalid_gc(self):
        """Sequence with GC instead of GT should be invalid."""
        assert is_valid_splice_site("CAGGCAAGTAT") is False

    def test_invalid_at(self):
        """Sequence with AT instead of GT should be invalid."""
        assert is_valid_splice_site("CAGATAAGTAT") is False

    def test_wrong_length(self):
        """Wrong length should return False."""
        assert is_valid_splice_site("CAGGT") is False


class TestGenerateAllValidSequences:
    """Tests for generate_all_valid_sequences function."""

    def test_total_count(self):
        """Should generate 4^9 = 262,144 sequences."""
        all_sequences = generate_all_valid_sequences()
        total = sum(len(seqs) for seqs in all_sequences.values())
        assert total == 4**9  # 262,144

    def test_all_have_gt(self):
        """All generated sequences should have GT at +1/+2."""
        all_sequences = generate_all_valid_sequences()
        for hbond_count, sequences in all_sequences.items():
            for seq in sequences[:10]:  # Check first 10 per tier
                assert seq[3] == "G" and seq[4] == "T", f"Invalid: {seq}"

    def test_hbond_counts_correct(self):
        """Spot check that H-bond counts are correct."""
        all_sequences = generate_all_valid_sequences()
        for hbond_count, sequences in all_sequences.items():
            for seq in sequences[:5]:  # Check first 5 per tier
                result = calculate_hbonds(seq)
                assert result.total_hbonds == hbond_count

    def test_includes_optimal(self):
        """Should include the optimal sequence at max H-bonds."""
        all_sequences = generate_all_valid_sequences()
        optimal_seq = "".join(OPTIMAL_DNA)
        assert MAX_HBONDS in all_sequences
        assert optimal_seq in all_sequences[MAX_HBONDS]


class TestGenerateSequencesForHbondCount:
    """Tests for generate_sequences_for_hbond_count function."""

    def test_returns_correct_count(self):
        """Should return up to num_variants sequences."""
        all_sequences = generate_all_valid_sequences()
        variants = generate_sequences_for_hbond_count(20, all_sequences, num_variants=10)
        assert len(variants) <= 10

    def test_prioritizes_stretch(self):
        """With prioritize_stretch=True, first variant should have longest stretch."""
        all_sequences = generate_all_valid_sequences()
        variants = generate_sequences_for_hbond_count(20, all_sequences, num_variants=100, prioritize_stretch=True)
        if len(variants) > 1:
            first_stretch = variants[0].hbond_result.longest_stretch
            for v in variants[1:]:
                assert v.hbond_result.longest_stretch <= first_stretch

    def test_nonexistent_hbond_count(self):
        """Should return empty list for impossible H-bond counts."""
        all_sequences = generate_all_valid_sequences()
        variants = generate_sequences_for_hbond_count(100, all_sequences)  # Impossible
        assert len(variants) == 0


class TestValidateKnownSequences:
    """Tests for validate_known_sequences function."""

    def test_validation_runs(self):
        """Validation should run without errors."""
        results = validate_known_sequences()
        assert len(results) > 0

    def test_some_pass(self):
        """At least some known sequences should validate."""
        results = validate_known_sequences()
        passed_count = sum(1 for _, _, _, passed in results if passed)
        assert passed_count > 0


class TestHBondResult:
    """Tests for HBondResult data class."""

    def test_to_dict(self):
        """to_dict should return proper dictionary format."""
        result = calculate_hbonds("CAGGTAAGTAT")
        d = result.to_dict()

        assert "sequence" in d
        assert "total_hbonds" in d
        assert "hbonds_pos_m3" in d  # -3 position
        assert "hbonds_pos_p1" in d  # +1 position
        assert "longest_stretch" in d
        assert "num_base_pairs" in d


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_all_a(self):
        """All A's except GT should give specific result."""
        result = calculate_hbonds("AAAGTAAAAAA")
        # Positions: A-G=0, A-U=2, A-C=0, G-C=3, T-A=2, A-Ψ=2, A-Ψ=2, A-A=0, A-C=0, A-U=2, A-G=0
        assert result.total_hbonds == 13

    def test_all_c(self):
        """All C's except GT should give specific result."""
        result = calculate_hbonds("CCCGTCCCCCC")
        # C-G=3, C-U=0, C-C=0, G-C=3, T-A=2, C-Ψ=0, C-Ψ=0, C-A=0, C-C=0, C-U=0, C-G=3
        assert result.total_hbonds == 11

    def test_all_g(self):
        """All G's except T at +2 should give specific result."""
        result = calculate_hbonds("GGGGTGGGGGG")
        # G-G=0, G-U=2, G-C=3, G-C=3, T-A=2, G-Ψ=2, G-Ψ=2, G-A=0, G-C=3, G-U=2, G-G=0
        assert result.total_hbonds == 19

    def test_all_t(self):
        """All T's except G at +1 and A at +7 should give specific result."""
        result = calculate_hbonds("TTTGTTTTTAT")
        # T(U) pairs with U1: -3 T(U)-G=2 wobble, -2 T-U=0, -1 T-C=0,
        # +1 G-C=3, +2 T-A=2, +3 T-Ψ=0, +4 T-Ψ=0, +5 T(U)-A=2,
        # +6 T-C=0, +7 A-U=2, +8 T(U)-G=2 wobble
        # Total: 2+0+0+3+2+0+0+2+0+2+2 = 13
        assert result.total_hbonds == 13


class TestFrame1StopCodons:
    """Tests for Frame 1 stop codon detection."""

    def test_no_stop_codons(self):
        """Sequence without stop codons should return False."""
        # CAGGTAAGTAT: codon2 = AAG (Lys), codon3 = TAT (Tyr)
        assert has_frame1_stop_codon("CAGGTAAGTAT") is False

    def test_stop_codon_at_position_345(self):
        """Stop codon at +3,+4,+5 should return True."""
        # codon2 = TAA (stop)
        assert has_frame1_stop_codon("CAGGT" + "TAA" + "GGG") is True
        # codon2 = TAG (stop)
        assert has_frame1_stop_codon("CAGGT" + "TAG" + "GGG") is True
        # codon2 = TGA (stop)
        assert has_frame1_stop_codon("CAGGT" + "TGA" + "GGG") is True

    def test_stop_codon_at_position_678(self):
        """Stop codon at +6,+7,+8 should return True."""
        # codon3 = TAA (stop)
        assert has_frame1_stop_codon("CAGGTGGG" + "TAA") is True
        # codon3 = TAG (stop)
        assert has_frame1_stop_codon("CAGGTGGG" + "TAG") is True
        # codon3 = TGA (stop)
        assert has_frame1_stop_codon("CAGGTGGG" + "TGA") is True

    def test_stop_codons_at_both_positions(self):
        """Stop codons at both positions should return True."""
        assert has_frame1_stop_codon("CAGGT" + "TAA" + "TGA") is True

    def test_invalid_length(self):
        """Invalid length should raise ValueError."""
        with pytest.raises(ValueError, match="must be 11 nucleotides"):
            has_frame1_stop_codon("CAGGT")


class TestGetFrame1Codons:
    """Tests for extracting Frame 1 codons."""

    def test_codon_extraction(self):
        """Should correctly extract codons at +3-5 and +6-8."""
        codon2, codon3 = get_frame1_codons("CAGGTAAGTAT")
        assert codon2 == "AAG"  # positions 5,6,7 (0-indexed)
        assert codon3 == "TAT"  # positions 8,9,10 (0-indexed)

    def test_another_sequence(self):
        """Test with different sequence."""
        codon2, codon3 = get_frame1_codons("AAAGTGGGCCC")
        assert codon2 == "GGG"
        assert codon3 == "CCC"


class TestTranslateCodon:
    """Tests for codon translation."""

    def test_common_codons(self):
        """Test common codon translations."""
        assert translate_codon("ATG") == "M"  # Methionine (start)
        assert translate_codon("TAA") == "*"  # Stop
        assert translate_codon("TAG") == "*"  # Stop
        assert translate_codon("TGA") == "*"  # Stop
        assert translate_codon("GGG") == "G"  # Glycine
        assert translate_codon("AAA") == "K"  # Lysine

    def test_lowercase_input(self):
        """Should handle lowercase input."""
        assert translate_codon("atg") == "M"

    def test_invalid_codon(self):
        """Invalid codon should return '?'."""
        assert translate_codon("XXX") == "?"
        assert translate_codon("AT") == "?"


class TestTranslateFrame1Codons:
    """Tests for translating Frame 1 codons from 11-mer."""

    def test_translation(self):
        """Should translate both codons."""
        aa2, aa3 = translate_frame1_codons("CAGGTAAGTAT")
        # codon2 = AAG -> K (Lysine)
        # codon3 = TAT -> Y (Tyrosine)
        assert aa2 == "K"
        assert aa3 == "Y"

    def test_with_stop_codons(self):
        """Should correctly identify stop codons as '*'."""
        aa2, aa3 = translate_frame1_codons("CAGGT" + "TAA" + "GGG")
        assert aa2 == "*"
        assert aa3 == "G"


class TestCalculateSequenceDiversity:
    """Tests for sequence diversity calculation."""

    def test_identical_sequences(self):
        """Identical sequences should have 0 diversity."""
        assert calculate_sequence_diversity("CAGGTAAGTAT", "CAGGTAAGTAT") == 0

    def test_one_difference(self):
        """One position different should return 1."""
        assert calculate_sequence_diversity("CAGGTAAGTAT", "AAGGTAAGTAT") == 1

    def test_all_different(self):
        """All positions different should return 11."""
        # These differ at all 11 positions (except GT must stay)
        # Actually let's pick sequences that truly differ
        seq1 = "AAAGTAAAAAA"
        seq2 = "CCCGTCCCGGG"  # Differs at: 0,1,2,5,6,7,8,9,10 = 9 positions
        assert calculate_sequence_diversity(seq1, seq2) == 9


class TestSelectDiverseWithinTier:
    """Tests for diverse sequence selection."""

    def test_returns_correct_count(self):
        """Should return requested number of sequences."""
        sequences = ["CAGGTAAGTAT", "AAGGTAAGTAT", "TAGGTAAGTAT", "GAGGTAAGTAT", "CAGGTAAGCAT"]
        # Use "short" mode for legacy behavior (sorted selection)
        selected = select_diverse_within_tier(sequences, 3, stretch_mode="short")
        assert len(selected) == 3

    def test_fewer_than_requested(self):
        """Should return all if fewer than requested."""
        sequences = ["CAGGTAAGTAT", "AAGGTAAGTAT"]
        selected = select_diverse_within_tier(sequences, 5, stretch_mode="short")
        assert len(selected) == 2

    def test_diversity_preference(self):
        """Should prefer diverse sequences."""
        # First three are very similar, last two are different
        sequences = [
            "CAGGTAAGTAT",  # Original
            "CAGGTAAGTAT",  # Identical
            "CAGGTAAGTAT",  # Identical
            "TTTGTCCCGGG",  # Very different
            "GGGGTAAACCC",  # Very different
        ]
        selected = select_diverse_within_tier(sequences, 3, stretch_mode="short")
        # Should include diverse sequences, not just identical ones
        assert len(set(selected)) >= 2  # At least 2 unique sequences

    def test_diverse_mode_spans_stretch_range(self):
        """Diverse mode should select sequences with varying stretch lengths."""
        all_sequences = generate_all_valid_sequences()
        # Tier 5 (18-19 H-bonds) has stretch range 3-8
        tier5_seqs = []
        for hb in range(18, 20):
            if hb in all_sequences:
                tier5_seqs.extend(all_sequences[hb])
        # Filter for consensus and no stops
        tier5_seqs = [s for s in tier5_seqs if has_consensus_motif(s) and not has_frame1_stop_codon(s)]

        selected = select_diverse_within_tier(tier5_seqs, 2, stretch_mode="diverse")
        stretches = [calculate_hbonds(s).longest_stretch for s in selected]
        # The two selected should have different stretch lengths
        assert len(set(stretches)) == 2, f"Expected 2 different stretches, got {stretches}"


class TestHasConsensusMotif:
    """Tests for AG|GT consensus motif detection."""

    def test_consensus_present(self):
        """Sequences with AG|GT at -1|+1+2 should return True."""
        # G at position -1 (index 2), GT at +1+2 (indices 3,4)
        assert has_consensus_motif("CAGGTAAGTAT") is True  # CAG|GT...
        assert has_consensus_motif("AAGGTAAGTAT") is True  # AAG|GT...
        assert has_consensus_motif("ACGGTAAGTAT") is True  # ACG|GT...
        assert has_consensus_motif("ATGGTAAGTAT") is True  # ATG|GT...

    def test_consensus_missing_g_at_m1(self):
        """Sequences without G at position -1 should return False."""
        assert has_consensus_motif("CAAGTAAGTAT") is False  # CAA|GT... (A at -1)
        assert has_consensus_motif("CACGTAAGTAT") is False  # CAC|GT... (C at -1)
        assert has_consensus_motif("CATGTAAGTAT") is False  # CAT|GT... (T at -1)

    def test_wrong_length(self):
        """Wrong length should return False."""
        assert has_consensus_motif("CAGGT") is False
        assert has_consensus_motif("CAGGTAAGTATAAA") is False

    def test_missing_gt(self):
        """Sequences without GT at +1+2 should return False."""
        assert has_consensus_motif("CAGATAAGTAT") is False  # AT instead of GT
        assert has_consensus_motif("CAGGCAAGTAT") is False  # GC instead of GT


class TestGenerateTieredPanel:
    """Tests for tiered panel generation."""

    def test_panel_generation(self):
        """Should generate panel with candidates from each tier."""
        all_sequences = generate_all_valid_sequences()
        candidates = generate_tiered_panel(all_sequences)

        # Should have candidates
        assert len(candidates) > 0

        # Should cover multiple tiers
        tiers_present = set(c.tier for c in candidates)
        assert len(tiers_present) >= 5  # Most tiers should have candidates

    def test_no_frame1_stops_by_default(self):
        """Default should filter out Frame 1 stop codons."""
        all_sequences = generate_all_valid_sequences()
        candidates = generate_tiered_panel(all_sequences, avoid_frame1_stops=True)

        for c in candidates:
            assert not c.has_frame1_stop, f"Candidate {c.sequence} has Frame 1 stop"

    def test_include_stops_option(self):
        """Should be able to include sequences with stops."""
        all_sequences = generate_all_valid_sequences()
        candidates = generate_tiered_panel(all_sequences, avoid_frame1_stops=False)

        # Some candidates might have stops (but not required to)
        assert len(candidates) > 0

    def test_candidate_attributes(self):
        """Candidates should have all required attributes."""
        all_sequences = generate_all_valid_sequences()
        candidates = generate_tiered_panel(all_sequences)

        for c in candidates:
            assert len(c.sequence) == 11
            assert c.tier >= 1
            assert len(c.tier_range) == 2
            assert c.total_hbonds >= c.tier_range[0]
            assert c.total_hbonds <= c.tier_range[1]
            assert len(c.codon2_dna) == 3
            assert len(c.codon3_dna) == 3
            assert len(c.codon2_aa) == 1
            assert len(c.codon3_aa) == 1

    def test_default_tiers_total(self):
        """Default tiers should produce approximately 16 candidates."""
        all_sequences = generate_all_valid_sequences()
        candidates = generate_tiered_panel(all_sequences)

        # Default: 2+3+3+2+2+2+2 = 16
        expected_total = sum(n for _, _, n in DEFAULT_TIERS)
        # Allow some variance if tiers are sparse
        assert len(candidates) >= expected_total - 4
        assert len(candidates) <= expected_total

    def test_consensus_motif_required_by_default(self):
        """Default should require AG|GT consensus motif."""
        all_sequences = generate_all_valid_sequences()
        candidates = generate_tiered_panel(all_sequences)

        for c in candidates:
            assert has_consensus_motif(c.sequence), f"Candidate {c.sequence} lacks AG|GT consensus"

    def test_tier_specific_selection_strategy(self):
        """Tiers 5-7 should have longer stretches than tiers 1-4."""
        all_sequences = generate_all_valid_sequences()
        candidates = generate_tiered_panel(all_sequences)

        weak_tier_stretches = [c.longest_stretch for c in candidates if c.tier <= 4]
        strong_tier_stretches = [c.longest_stretch for c in candidates if c.tier >= 5]

        if weak_tier_stretches and strong_tier_stretches:
            # Strong tiers should have longer stretches on average
            avg_weak = sum(weak_tier_stretches) / len(weak_tier_stretches)
            avg_strong = sum(strong_tier_stretches) / len(strong_tier_stretches)
            assert avg_strong > avg_weak, f"Strong tier avg ({avg_strong}) should exceed weak tier avg ({avg_weak})"

    def test_tier7_has_maximum_hbonds(self):
        """Tier 7 (positive controls) should have high H-bond sequences."""
        all_sequences = generate_all_valid_sequences()
        candidates = generate_tiered_panel(all_sequences)

        tier7_candidates = [c for c in candidates if c.tier == 7]
        if tier7_candidates:
            max_hbonds = max(c.total_hbonds for c in tier7_candidates)
            # Tier 7 should have sequences with 22+ H-bonds
            assert max_hbonds >= 22, f"Tier 7 max H-bonds ({max_hbonds}) should be >= 22"

    def test_stretch_mode_short(self):
        """stretch_mode='short' should give short stretches for all tiers."""
        all_sequences = generate_all_valid_sequences()
        candidates_short = generate_tiered_panel(all_sequences, stretch_mode="short")
        candidates_diverse = generate_tiered_panel(all_sequences, stretch_mode="diverse")

        # With short stretch mode, tier 7 should have shorter stretches than diverse
        tier7_short = [c for c in candidates_short if c.tier == 7]
        tier7_diverse = [c for c in candidates_diverse if c.tier == 7]

        if tier7_short and tier7_diverse:
            max_short = max(c.longest_stretch for c in tier7_short)
            max_diverse = max(c.longest_stretch for c in tier7_diverse)
            # Short mode should have shorter max stretch
            assert max_short <= max_diverse

    def test_stretch_mode_long(self):
        """stretch_mode='long' should give long stretches for all tiers."""
        all_sequences = generate_all_valid_sequences()
        candidates_long = generate_tiered_panel(all_sequences, stretch_mode="long")
        candidates_diverse = generate_tiered_panel(all_sequences, stretch_mode="diverse")

        # With long stretch mode, tier 1 should have longer min stretch than diverse
        tier1_long = [c for c in candidates_long if c.tier == 1]
        tier1_diverse = [c for c in candidates_diverse if c.tier == 1]

        if tier1_long and tier1_diverse:
            min_long = min(c.longest_stretch for c in tier1_long)
            min_diverse = min(c.longest_stretch for c in tier1_diverse)
            # Long mode should have longer min stretch
            assert min_long >= min_diverse

    def test_diverse_mode_gradual_progression(self):
        """stretch_mode='diverse' should create gradual stretch progression across tiers."""
        all_sequences = generate_all_valid_sequences()
        candidates = generate_tiered_panel(all_sequences, stretch_mode="diverse")

        # Collect max stretch per tier
        tier_max_stretch = {}
        for c in candidates:
            if c.tier not in tier_max_stretch:
                tier_max_stretch[c.tier] = c.longest_stretch
            else:
                tier_max_stretch[c.tier] = max(tier_max_stretch[c.tier], c.longest_stretch)

        # Higher tiers should have higher max stretches (generally)
        tier_nums = sorted(tier_max_stretch.keys())
        for i in range(len(tier_nums) - 1):
            t1, t2 = tier_nums[i], tier_nums[i+1]
            # Allow some tolerance - higher tier should have >= max stretch
            assert tier_max_stretch[t2] >= tier_max_stretch[t1] - 1, \
                f"Tier {t2} max stretch ({tier_max_stretch[t2]}) should be >= tier {t1} ({tier_max_stretch[t1]})"


class TestTieredCandidateToDict:
    """Tests for TieredCandidate.to_dict() method."""

    def test_dict_has_required_fields(self):
        """to_dict should include all required fields."""
        all_sequences = generate_all_valid_sequences()
        candidates = generate_tiered_panel(all_sequences)

        if candidates:
            d = candidates[0].to_dict()
            required_fields = [
                "tier", "tier_hbonds", "sequence", "total_hbonds",
                "longest_stretch", "frame1_stop", "codon2_pos345",
                "codon3_pos678", "codon2_aa", "codon3_aa",
                "expected_splicing", "notes"
            ]
            for field in required_fields:
                assert field in d, f"Missing field: {field}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
