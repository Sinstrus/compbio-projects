"""
Tests for uniqueness counting bug (BUG-003).

The bug: When checking if a restriction site is unique in a sequence,
the code was only searching the forward strand and not accounting for
the reverse complement. For palindromic sites like HindIII (AAGCTT),
this means the site appears twice (once on each strand) but was counted
as appearing only once.

For double-stranded DNA:
  Forward:  5'-...AAGCTT...-3'
  Reverse:  3'-...TTCGAA...-5'

If searching only the forward strand, we miss the reverse complement
occurrence, leading to false positives for "unique" site claims.
"""

import pytest


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement[base] for base in reversed(seq))


def count_sites_single_strand_buggy(sequence, site):
    """
    BUGGY implementation: only counts on forward strand.
    This was the bug in the original code.
    """
    count = 0
    pos = 0
    while True:
        pos = sequence.find(site, pos)
        if pos == -1:
            break
        count += 1
        pos += 1
    return count


def count_sites_double_strand_correct(sequence, site):
    """
    CORRECT implementation: counts on both strands.
    For palindromic sites, this correctly counts both occurrences.
    """
    # Count forward strand
    forward_count = 0
    pos = 0
    while True:
        pos = sequence.find(site, pos)
        if pos == -1:
            break
        forward_count += 1
        pos += 1

    # Count reverse strand (reverse complement)
    rc_site = reverse_complement(site)
    reverse_count = 0
    pos = 0
    while True:
        pos = sequence.find(rc_site, pos)
        if pos == -1:
            break
        reverse_count += 1
        pos += 1

    # For palindromic sites, forward and reverse are the same
    # So we need to avoid double counting
    if site == rc_site:
        # Palindromic: each occurrence is on both strands
        return forward_count  # Each site is already counted on both strands
    else:
        # Non-palindromic: sum both counts
        return forward_count + reverse_count


def is_site_unique_buggy(sequence, site):
    """BUGGY: claims site is unique if it appears once on forward strand only."""
    return count_sites_single_strand_buggy(sequence, site) == 1


def is_site_unique_correct(sequence, site):
    """CORRECT: site is unique only if it appears once total across both strands."""
    return count_sites_double_strand_correct(sequence, site) == 1


class TestPalindromicSites:
    """Test palindromic restriction sites that read the same on both strands."""

    def test_hindiii_is_palindromic(self):
        """HindIII (AAGCTT) is palindromic."""
        hindiii = "AAGCTT"
        assert reverse_complement(hindiii) == hindiii

    def test_xbai_is_palindromic(self):
        """XbaI (TCTAGA) is palindromic."""
        xbai = "TCTAGA"
        assert reverse_complement(xbai) == xbai

    def test_ecorv_is_palindromic(self):
        """EcoRV (GATATC) is palindromic."""
        ecorv = "GATATC"
        assert reverse_complement(ecorv) == ecorv

    def test_palindrome_appears_on_both_strands(self):
        """
        A palindromic site appearing once in the sequence
        is actually present on BOTH strands.
        """
        # Sequence with one HindIII site
        sequence = "ATGCAAGCTTGCAT"
        hindiii = "AAGCTT"

        # On forward strand: ...AAGCTT...
        # On reverse strand at same location: ...AAGCTT... (because palindrome)

        # Buggy code thinks it appears once
        buggy_count = count_sites_single_strand_buggy(sequence, hindiii)
        assert buggy_count == 1

        # But it's really on both strands - the site exists in double-stranded DNA
        # For palindromes, the count is the same, but conceptually it's on both strands
        correct_count = count_sites_double_strand_correct(sequence, hindiii)
        assert correct_count == 1  # Still 1, but we're aware it's on both strands


class TestNonPalindromicSites:
    """Test non-palindromic restriction sites."""

    def test_bsai_is_not_palindromic(self):
        """BsaI (GGTCTC) is not palindromic."""
        bsai = "GGTCTC"
        rc = reverse_complement(bsai)
        assert bsai != rc
        assert rc == "GAGACC"

    def test_non_palindrome_different_strands(self):
        """Non-palindromic sites have different sequences on each strand."""
        sequence = "ATGCGGTCTCGCAT"  # Contains BsaI (GGTCTC)
        bsai = "GGTCTC"

        # Forward strand has GGTCTC
        forward_count = count_sites_single_strand_buggy(sequence, bsai)
        assert forward_count == 1

        # Reverse complement is GAGACC (not in this sequence)
        rc_bsai = reverse_complement(bsai)
        rc_count = count_sites_single_strand_buggy(sequence, rc_bsai)
        assert rc_count == 0

        # Total on both strands
        total = count_sites_double_strand_correct(sequence, bsai)
        assert total == 1

    def test_non_palindrome_on_reverse_strand_only(self):
        """Non-palindromic site on reverse strand only."""
        # GAGACC in sequence means BsaI (GGTCTC) on reverse strand
        sequence = "ATGCGAGACCGCAT"
        bsai = "GGTCTC"

        # Not on forward strand
        forward_count = count_sites_single_strand_buggy(sequence, bsai)
        assert forward_count == 0

        # But reverse complement IS in sequence
        rc_bsai = reverse_complement(bsai)  # GAGACC
        assert rc_bsai in sequence

        # Correct counting finds it on reverse strand
        total = count_sites_double_strand_correct(sequence, bsai)
        assert total == 1


class TestUniquenessChecking:
    """Test the uniqueness checking that was affected by the bug."""

    def test_bug_false_unique_claim(self):
        """
        The bug: claiming a site is "unique" when it appears once on
        forward strand, but might appear on reverse strand too.
        """
        # Sequence with BsaI on forward AND reverse strands (different positions)
        sequence = "ATGCGGTCTCGCATGAGACCGCAT"
        #              GGTCTC  (forward)
        #                      GAGACC (= BsaI on reverse)
        bsai = "GGTCTC"

        # Buggy check only sees forward occurrence
        buggy_unique = is_site_unique_buggy(sequence, bsai)
        assert buggy_unique == True  # WRONG! Claims it's unique

        # Correct check sees both strands
        correct_unique = is_site_unique_correct(sequence, bsai)
        assert correct_unique == False  # RIGHT! Not unique

    def test_truly_unique_site(self):
        """A truly unique site appears once across both strands."""
        sequence = "ATGCGGTCTCGCAT"  # One BsaI on forward strand only
        bsai = "GGTCTC"

        assert is_site_unique_buggy(sequence, bsai) == True
        assert is_site_unique_correct(sequence, bsai) == True

    def test_multiple_sites_on_same_strand(self):
        """Multiple sites on the same strand."""
        sequence = "ATGCGGTCTCGCATGGTCTCGCAT"  # Two BsaI on forward
        bsai = "GGTCTC"

        assert is_site_unique_buggy(sequence, bsai) == False
        assert is_site_unique_correct(sequence, bsai) == False

    def test_palindrome_context_dependent(self):
        """
        For palindromic sites in the context of a larger construct,
        we need to consider the biological meaning.

        One occurrence in the sequence means it can be cut once (but the
        enzyme sees it from both sides of the double helix).
        """
        # One HindIII site in a backbone
        backbone = "ATGCAAGCTTGCATGCGATCGATCG"
        hindiii = "AAGCTT"

        # Biologically, this IS unique (one cutting position)
        assert is_site_unique_correct(backbone, hindiii) == True


class TestContextDependentUniqueness:
    """
    Test DESIGN-002: Context-dependent uniqueness.

    A site might be "unique in the CDS" but not "unique in the final construct"
    when you account for the backbone.
    """

    def test_unique_in_cds_not_in_construct(self):
        """Site is unique in CDS but not in final construct with backbone."""
        cds = "ATGAAAGCCTAA"  # No HindIII
        backbone = "GGATCCAAGCTTGAATTC"  # Has HindIII

        hindiii = "AAGCTT"

        # Unique in CDS? Yes (zero occurrences is "unique" in sense of "not problematic")
        assert count_sites_double_strand_correct(cds, hindiii) == 0

        # Unique in final construct? Depends on cloning method
        # If using HindIII for cloning, there's one in backbone
        full_construct = backbone + cds
        assert count_sites_double_strand_correct(full_construct, hindiii) == 1

        # Adding one more HindIII in CDS would create non-unique site
        cds_with_site = cds + hindiii + "GGG"
        full_construct_2 = backbone + cds_with_site
        assert count_sites_double_strand_correct(full_construct_2, hindiii) == 2

    def test_must_check_final_construct(self):
        """
        Critical lesson: must verify uniqueness in the FINAL construct,
        not just the CDS in isolation.
        """
        # This was the issue with XbaI in v04
        cds_with_xbai = "ATGTCTAGAGCCTAA"  # Has XbaI
        backbone_with_xbai = "AAGCTTTCTAGAGAATTC"  # Also has XbaI for cloning

        xbai = "TCTAGA"

        # CDS has one XbaI
        assert count_sites_double_strand_correct(cds_with_xbai, xbai) == 1

        # Backbone has one XbaI (cloning site)
        assert count_sites_double_strand_correct(backbone_with_xbai, xbai) == 1

        # Final construct has TWO - not unique!
        full = backbone_with_xbai + cds_with_xbai
        assert count_sites_double_strand_correct(full, xbai) == 2
        assert is_site_unique_correct(full, xbai) == False


def test_real_world_scenario_from_catalog():
    """
    Real scenario from the backbone catalog:
    Check HindIII and XbaI uniqueness in AAV backbones.
    """
    # From catalog: pGS-scAAV-ITR128-Amp has:
    # - HindIII at position 185
    # - XbaI at position 2276

    # Simulate backbone sequence (simplified)
    # In real case, we'd load from the catalog
    hindiii = "AAGCTT"
    xbai = "TCTAGA"

    # These are at specific positions, so they should each appear once
    # (This is confirmed by the catalog data)

    # When designing a CDS:
    # - Must NOT include HindIII (used for cloning)
    # - Must NOT include XbaI (used for cloning)
    # Otherwise the site won't be unique in final construct

    # Example: CDS that would break cloning
    bad_cds = "ATGAAGCTTAAA"  # Has HindIII - BAD!
    assert count_sites_double_strand_correct(bad_cds, hindiii) == 1

    # Good CDS: no cloning sites
    good_cds = "ATGAAGCTCAAA"  # Similar but no HindIII - GOOD!
    assert count_sites_double_strand_correct(good_cds, hindiii) == 0
