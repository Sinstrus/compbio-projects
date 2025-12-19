"""
Tests for silent mutation classification.

Tests the logic for determining whether a nucleotide mutation is silent
(synonymous) or causes an amino acid change (non-synonymous).

This is critical for Checkpoint 8: Silent Mutation Verification.
"""

import pytest


# Standard genetic code
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate_codon(codon):
    """Translate a codon to amino acid."""
    return CODON_TABLE.get(codon.upper(), 'X')


def is_mutation_silent(original_codon, mutated_codon):
    """
    Check if a mutation is silent (synonymous).

    Args:
        original_codon: Original codon sequence (3 bases)
        mutated_codon: Mutated codon sequence (3 bases)

    Returns:
        True if mutation is silent (same amino acid), False otherwise
    """
    original_aa = translate_codon(original_codon)
    mutated_aa = translate_codon(mutated_codon)
    return original_aa == mutated_aa


def mutate_codon(codon, position, new_base):
    """
    Mutate a codon at a specific position.

    Args:
        codon: Original codon (3 bases)
        position: Position in codon (0, 1, or 2)
        new_base: New base to substitute

    Returns:
        Mutated codon
    """
    codon_list = list(codon)
    codon_list[position] = new_base
    return ''.join(codon_list)


class TestSilentMutations:
    """Test identification of silent (synonymous) mutations."""

    def test_leucine_wobble_position(self):
        """Leucine has 6 codons - test wobble position mutations."""
        # CTT -> CTC (both Leucine)
        assert is_mutation_silent('CTT', 'CTC') == True

        # CTT -> CTA (both Leucine)
        assert is_mutation_silent('CTT', 'CTA') == True

        # CTT -> CTG (both Leucine)
        assert is_mutation_silent('CTT', 'CTG') == True

    def test_arginine_wobble_position(self):
        """Arginine has 6 codons - test wobble position mutations."""
        # CGT -> CGC (both Arginine)
        assert is_mutation_silent('CGT', 'CGC') == True

        # CGT -> CGA (both Arginine)
        assert is_mutation_silent('CGT', 'CGA') == True

        # CGT -> CGG (both Arginine)
        assert is_mutation_silent('CGT', 'CGG') == True

    def test_serine_wobble_position(self):
        """Serine has 6 codons (two families) - test mutations."""
        # TCT -> TCC (both Serine)
        assert is_mutation_silent('TCT', 'TCC') == True

        # TCT -> TCA (both Serine)
        assert is_mutation_silent('TCT', 'TCA') == True

        # AGT -> AGC (both Serine, different family)
        assert is_mutation_silent('AGT', 'AGC') == True

    def test_four_fold_degenerate_codons(self):
        """Test 4-fold degenerate codons (wobble position completely free)."""
        # Alanine: GCT, GCC, GCA, GCG
        assert is_mutation_silent('GCT', 'GCC') == True
        assert is_mutation_silent('GCT', 'GCA') == True
        assert is_mutation_silent('GCT', 'GCG') == True

        # Valine: GTT, GTC, GTA, GTG
        assert is_mutation_silent('GTT', 'GTC') == True
        assert is_mutation_silent('GTT', 'GTA') == True
        assert is_mutation_silent('GTT', 'GTG') == True

        # Proline: CCT, CCC, CCA, CCG
        assert is_mutation_silent('CCT', 'CCC') == True
        assert is_mutation_silent('CCT', 'CCA') == True
        assert is_mutation_silent('CCT', 'CCG') == True

    def test_two_fold_degenerate_codons(self):
        """Test 2-fold degenerate codons (limited wobble)."""
        # Phenylalanine: TTT, TTC
        assert is_mutation_silent('TTT', 'TTC') == True

        # Tyrosine: TAT, TAC
        assert is_mutation_silent('TAT', 'TAC') == True

        # Histidine: CAT, CAC
        assert is_mutation_silent('CAT', 'CAC') == True

        # Asparagine: AAT, AAC
        assert is_mutation_silent('AAT', 'AAC') == True

        # Aspartic acid: GAT, GAC
        assert is_mutation_silent('GAT', 'GAC') == True


class TestNonSilentMutations:
    """Test identification of non-silent (non-synonymous) mutations."""

    def test_single_codon_amino_acids(self):
        """Methionine and Tryptophan have only one codon each."""
        # Any mutation in ATG (Met) is non-silent
        assert is_mutation_silent('ATG', 'ATA') == False  # Met -> Ile
        assert is_mutation_silent('ATG', 'ACG') == False  # Met -> Thr
        assert is_mutation_silent('ATG', 'AAG') == False  # Met -> Lys

        # Any mutation in TGG (Trp) is non-silent
        assert is_mutation_silent('TGG', 'TGA') == False  # Trp -> Stop
        assert is_mutation_silent('TGG', 'TAG') == False  # Trp -> Stop
        assert is_mutation_silent('TGG', 'TCG') == False  # Trp -> Ser

    def test_first_position_usually_nonsilent(self):
        """Mutations in first position usually change amino acid."""
        # AAA (Lys) -> GAA (Glu)
        assert is_mutation_silent('AAA', 'GAA') == False

        # GCT (Ala) -> ACT (Thr)
        assert is_mutation_silent('GCT', 'ACT') == False

        # CTT (Leu) -> ATT (Ile)
        assert is_mutation_silent('CTT', 'ATT') == False

    def test_second_position_always_nonsilent(self):
        """Mutations in second position ALWAYS change amino acid."""
        # AAA (Lys) -> ATA (Ile)
        assert is_mutation_silent('AAA', 'ATA') == False

        # GCT (Ala) -> GTT (Val)
        assert is_mutation_silent('GCT', 'GTT') == False

        # CTT (Leu) -> CAT (His)
        assert is_mutation_silent('CTT', 'CAT') == False

        # TGG (Trp) -> TAG (Stop)
        assert is_mutation_silent('TGG', 'TAG') == False

    def test_stop_codon_mutations(self):
        """Mutations to/from stop codons are non-silent."""
        # TAA (Stop) -> TAT (Tyr)
        assert is_mutation_silent('TAA', 'TAT') == False

        # TAG (Stop) -> TGG (Trp)
        assert is_mutation_silent('TAG', 'TGG') == False

        # TGA (Stop) -> TTA (Leu)
        assert is_mutation_silent('TGA', 'TTA') == False

        # CAA (Gln) -> TAA (Stop)
        assert is_mutation_silent('CAA', 'TAA') == False


class TestMutationByPosition:
    """Test mutation effects by codon position."""

    def test_wobble_position_systematic(self):
        """Systematically test wobble position (position 2) for all codons."""
        test_cases = [
            ('GCT', 'GCC', True),   # Ala -> Ala
            ('GCT', 'GCA', True),   # Ala -> Ala
            ('GTT', 'GTC', True),   # Val -> Val
            ('CTT', 'CTC', True),   # Leu -> Leu
            ('CCT', 'CCC', True),   # Pro -> Pro
            ('AAT', 'AAC', True),   # Asn -> Asn
            ('GAT', 'GAC', True),   # Asp -> Asp
            ('ATT', 'ATC', True),   # Ile -> Ile
            ('ATG', 'ATA', False),  # Met -> Ile (special case)
        ]

        for original, mutated, expected in test_cases:
            result = is_mutation_silent(original, mutated)
            assert result == expected, f"{original} -> {mutated}: expected {expected}, got {result}"

    def test_mutation_position_in_context(self):
        """Test mutations in the context of a full CDS."""
        # Full CDS: ATG AAA GCC TAA
        cds = "ATGAAAGCCTAA"

        # Mutate second codon (AAA = Lys) at wobble position
        # AAA -> AAG (both Lys) - SILENT
        original = cds[3:6]  # "AAA"
        mutated = mutate_codon(original, 2, 'G')  # "AAG"
        assert is_mutation_silent(original, mutated) == True

        # Mutate second codon at middle position
        # AAA -> ATA (Lys -> Ile) - NON-SILENT
        mutated2 = mutate_codon(original, 1, 'T')  # "ATA"
        assert is_mutation_silent(original, mutated2) == False


class TestRealWorldScenarios:
    """Test real-world scenarios for silent mutation design."""

    def test_removing_restriction_site_silently(self):
        """
        Common task: remove a restriction site without changing protein.

        Example: Remove BbsI (GAAGAC) from coding sequence.
        GAAGAC codes for Glu-Asp in one frame.
        """
        # Scenario: BbsI site spans two codons
        # ...GAA GAC... (Glu-Asp)
        original_codon1 = "GAA"  # Glu
        original_codon2 = "GAC"  # Asp

        # Try to mutate to remove BbsI while keeping amino acids
        # GAA -> GAG (both Glu) - SILENT
        mutated_codon1 = "GAG"
        assert is_mutation_silent(original_codon1, mutated_codon1) == True

        # Result: ...GAG GAC... (still Glu-Asp, but no BbsI)
        assert "GAAGAC" not in "GAGGAC"

    def test_multiple_mutation_strategy(self):
        """
        Sometimes need multiple mutations to remove a site.
        Each must be verified as silent.
        """
        # Remove HindIII (AAGCTT) that spans two codons
        # ...AAG CTT... (Lys-Leu)

        # Strategy 1: Mutate wobble positions
        # AAG -> AAA (both Lys)
        assert is_mutation_silent('AAG', 'AAA') == True
        # CTT -> CTC (both Leu)
        assert is_mutation_silent('CTT', 'CTC') == True

        # Result: ...AAA CTC... (still Lys-Leu, no HindIII)
        assert "AAGCTT" not in "AAACTC"

    def test_checkpoint_8_verification(self):
        """
        Checkpoint 8: Verify all mutations are silent.

        This test simulates the verification step that should happen
        after codon optimization to remove restriction sites.
        """
        mutations = [
            # (original_codon, mutated_codon, expected_silent)
            ('GAA', 'GAG', True),   # Glu -> Glu (remove BbsI)
            ('CTT', 'CTC', True),   # Leu -> Leu (remove HindIII)
            ('AAG', 'AAA', True),   # Lys -> Lys (remove HindIII)
            ('GCT', 'GCC', True),   # Ala -> Ala (general optimization)
        ]

        all_silent = all(
            is_mutation_silent(orig, mut) == expected
            for orig, mut, expected in mutations
        )

        assert all_silent, "All mutations should be silent"

    def test_accidental_nonsense_mutation(self):
        """
        Critical check: ensure we don't accidentally create stop codons.
        """
        # TAA, TAG, TGA are stop codons
        stop_codons = ['TAA', 'TAG', 'TGA']

        # Common mistake: mutating CAA (Gln) to TAA (Stop)
        assert is_mutation_silent('CAA', 'TAA') == False

        # Verify all paths to stop codons are non-silent
        sense_codons = ['CAA', 'CAG', 'TGG']  # Gln, Gln, Trp

        for sense in sense_codons:
            for stop in stop_codons:
                if sense != stop:  # Don't test stop -> stop
                    assert is_mutation_silent(sense, stop) == False


def test_translation_correctness():
    """Verify the codon table is implemented correctly."""
    test_cases = [
        ('ATG', 'M'),   # Start codon
        ('TAA', '*'),   # Stop codon
        ('TAG', '*'),   # Stop codon
        ('TGA', '*'),   # Stop codon
        ('TGG', 'W'),   # Trp (only one codon)
        ('GCT', 'A'),   # Ala
        ('TTT', 'F'),   # Phe
    ]

    for codon, expected_aa in test_cases:
        assert translate_codon(codon) == expected_aa


def test_complete_mutation_workflow():
    """
    Test a complete workflow of identifying and fixing a restriction site.

    This simulates what the agent should do:
    1. Identify restriction site in CDS
    2. Determine which codon(s) it affects
    3. Find silent mutations to remove it
    4. Verify mutations are silent
    """
    # Step 1: CDS with problematic BbsI site
    cds = "ATGGAAGACGCCTAA"
    #      ATG GAA GAC GCC TAA
    #      Met Glu Asp Ala Stop

    # Step 2: BbsI (GAAGAC) spans codons 2-3 (positions 3-9)
    site = "GAAGAC"
    site_pos = cds.find(site)
    assert site_pos == 3

    # Step 3: Extract affected codons
    codon2 = cds[3:6]   # "GAA" (Glu)
    codon3 = cds[6:9]   # "GAC" (Asp)

    # Step 4: Find silent mutations
    # GAA -> GAG (both Glu)
    new_codon2 = "GAG"
    assert is_mutation_silent(codon2, new_codon2) == True

    # GAC -> GAT (both Asp)
    new_codon3 = "GAT"
    assert is_mutation_silent(codon3, new_codon3) == True

    # Step 5: Construct new CDS
    new_cds = cds[:3] + new_codon2 + new_codon3 + cds[9:]
    assert new_cds == "ATGGAGGATGCCTAA"

    # Step 6: Verify site is removed
    assert site not in new_cds

    # Step 7: Verify protein sequence unchanged
    # (In real implementation, would translate both and compare)
    assert translate_codon("ATG") == translate_codon("ATG")  # Met
    assert translate_codon("GAA") == translate_codon("GAG")  # Glu
    assert translate_codon("GAC") == translate_codon("GAT")  # Asp
    assert translate_codon("GCC") == translate_codon("GCC")  # Ala
    assert translate_codon("TAA") == translate_codon("TAA")  # Stop
