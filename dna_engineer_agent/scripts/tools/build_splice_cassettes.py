#!/usr/bin/env python3
"""
Build Splice Cassettes for Splicing Prediction Software

Constructs full intron cassettes combining:
- Variable splice donor (11-mer from tiered panel)
- GGGGS5 linker (5x GGGGS)
- VHH3 nanobody (codon-optimized for Frame 1)
- GGGGS1 linker (1x GGGGS)
- Constant splice acceptor (human β-globin intron 2)

Output: FASTA file with all cassettes for splicing prediction tools
"""

import argparse
import csv
import os
from pathlib import Path
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

# ============================================================================
# Constants
# ============================================================================

# Stop codons
STOP_CODONS = {"TAA", "TAG", "TGA"}

# Standard genetic code
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Reverse codon table: amino acid -> list of codons
AA_TO_CODONS = {}
for codon, aa in CODON_TABLE.items():
    if aa not in AA_TO_CODONS:
        AA_TO_CODONS[aa] = []
    AA_TO_CODONS[aa].append(codon)

# VHH3 anti-ALPL protein sequence (119 aa)
VHH3_PROTEIN = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFSTADMGWFRQAPGKGR"
    "ELVAAVSGSGFSTYSDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAKATISLYY"
    "AMDVWGQGTTVTVSS"
)

# GGGGS unit
GGGGS_UNIT = "GGGGS"

# Human β-globin intron 2 splice acceptor
# From GenBank: NG_000007.3 (human β-globin gene)
# Intron 2 is ~850 bp, we use the 3' end (last ~50 bp of intron + start of exon 3)
#
# Structure: ...branch_point...polypyrimidine_tract...AG|exon3
# The acceptor region (last 50 bp of intron 2):
#   TCTATTTTCCCACCCTTAGG = 20 bp ending at the AG splice site
#
# For splice prediction, we need:
# - Branch point region (~20-50 bp upstream of AG)
# - Polypyrimidine tract (pyrimidine-rich, immediately upstream of AG)
# - AG dinucleotide (invariant 3' splice site)
# - Downstream exon sequence

# β-globin intron 2 acceptor: last 50 bp of intron + first 20 bp of exon 3
# Intron ends: ...CTAATGAAAATAAATTTCCTTTATTAGCCAGAAGTCAGATGCTAAAGGG
# The AG is at the very end: ...AGG with the G being the first nt of exon 3
# Actually the canonical sequence is: ...tctattttcccacccttagGCTGCTGGTGGTC...
#                                                        ^^ AG splice site

BETA_GLOBIN_INTRON2_ACCEPTOR_REGION = (
    # Last ~45 bp of intron 2 (including branch point and PPT)
    "CTAATGAAAATAAATTTCCTTTATTAGCCAGAAGTCAGATGCTAAAGGGACTTTGTGTAATTTTAGAAAATA"
    "ATATTTCATCTTTTATTTCCTTTACCCTTTCATCATTTTTCTTTTCTATTTTCCCACCCTTAGG"
)
# Note: The AG is at the end (last 2 nucleotides)

# Exon 3 of β-globin (first 60 bp) - starts with G (the G after the splice site AG)
BETA_GLOBIN_EXON3_START = (
    "GCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCAC"
)

# Combined acceptor for cassettes: intronic acceptor region ending in AG + exonic sequence
# CRITICAL: Intron must END with AG, exon starts with the next nucleotide (G in this case)
SPLICE_ACCEPTOR_BGLOBIN = (
    "TTCTATTTTCCCACCCTTAG" +   # Last 20 bp of intron (PPT region ending in AG)
    BETA_GLOBIN_EXON3_START    # First 60 bp of exon 3 (starts with G)
)

# Shorter β-globin acceptor (for more compact cassettes)
SPLICE_ACCEPTOR_BGLOBIN_SHORT = (
    "TTTTTTTTTTTCTATTTTCCCACCCTTAG" +  # Extended PPT + intronic (ends in AG)
    "GCTGCTGGTGGTCTACCCTTGGACCCAG"     # Exon 3 start (29 bp, starts with G)
)

# Legacy minimal acceptor (not recommended)
SPLICE_ACCEPTOR_MINIMAL = (
    "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTAACAG"  # Synthetic PPT + AG
)


# ============================================================================
# Functions
# ============================================================================

def translate(dna: str) -> str:
    """Translate DNA to protein."""
    protein = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3].upper()
        aa = CODON_TABLE.get(codon, "X")
        protein.append(aa)
    return "".join(protein)


def find_unavoidable_frame1_stops(protein: str) -> List[Tuple[int, str, str, str]]:
    """
    Find unavoidable Frame 1 stop codons based on dipeptide constraints.

    A Frame 1 stop at the boundary of protein positions i and i+1 occurs when:
    - Codon_i[1:3] + Codon_{i+1}[0] ∈ {TAA, TAG, TGA}

    This is unavoidable when ALL codons for aa_i end with a pattern that,
    combined with ALL first characters of codons for aa_{i+1}, form a stop.

    Problematic patterns:
    - "TG" + "A" = "TGA" (stop)
    - "TA" + "A" = "TAA" (stop)
    - "TA" + "G" = "TAG" (stop)

    Args:
        protein: Amino acid sequence

    Returns:
        List of (position, aa1, aa2, reason) tuples for unavoidable stops
    """
    ALL_CODONS = {
        "A": ["GCC", "GCT", "GCA", "GCG"],
        "C": ["TGC", "TGT"],
        "D": ["GAC", "GAT"],
        "E": ["GAG", "GAA"],
        "F": ["TTC", "TTT"],
        "G": ["GGC", "GGT", "GGA", "GGG"],
        "H": ["CAC", "CAT"],
        "I": ["ATC", "ATT", "ATA"],
        "K": ["AAG", "AAA"],
        "L": ["CTG", "CTC", "CTT", "CTA", "TTG", "TTA"],
        "M": ["ATG"],
        "N": ["AAC", "AAT"],
        "P": ["CCC", "CCT", "CCA", "CCG"],
        "Q": ["CAG", "CAA"],
        "R": ["CGC", "CGG", "AGA", "AGG", "CGT", "CGA"],
        "S": ["AGC", "TCC", "TCT", "AGT", "TCA", "TCG"],
        "T": ["ACC", "ACT", "ACA", "ACG"],
        "V": ["GTG", "GTC", "GTT", "GTA"],
        "W": ["TGG"],
        "Y": ["TAC", "TAT"],
    }

    unavoidable = []

    for i in range(len(protein) - 1):
        aa1, aa2 = protein[i], protein[i+1]
        codons1 = ALL_CODONS.get(aa1, [])
        codons2 = ALL_CODONS.get(aa2, [])

        if not codons1 or not codons2:
            continue

        # Check if ALL combinations produce a stop
        all_stop = True
        for c1 in codons1:
            for c2 in codons2:
                boundary = c1[1:3] + c2[0]
                if boundary not in STOP_CODONS:
                    all_stop = False
                    break
            if not all_stop:
                break

        if all_stop:
            # Find the specific constraint
            endings = set(c[1:3] for c in codons1)
            starts = set(c[0] for c in codons2)
            reason = f"All {aa1} codons end with {endings}, all {aa2} codons start with {starts}"
            unavoidable.append((i, aa1, aa2, reason))

    return unavoidable


def find_stop_codons_in_frame(dna: str, frame: int = 0) -> List[Tuple[int, str]]:
    """
    Find all stop codons in a specific reading frame.

    Args:
        dna: DNA sequence
        frame: Reading frame (0, 1, or 2)

    Returns:
        List of (position, codon) tuples
    """
    stops = []
    for i in range(frame, len(dna) - 2, 3):
        codon = dna[i:i+3].upper()
        if codon in STOP_CODONS:
            stops.append((i, codon))
    return stops


def codon_optimize_for_frame1_stop_free(protein: str, verbose: bool = False) -> str:
    """
    Codon-optimize a protein sequence to avoid stop codons in Frame 1.

    Frame 1 reads codons starting at position 1:
    - Frame 0: [0,1,2], [3,4,5], [6,7,8], ...
    - Frame 1: [1,2,3], [4,5,6], [7,8,9], ...

    Each Frame 1 codon spans two adjacent Frame 0 codons:
    - Frame 1 codon at [1,2,3] = Frame0_codon0[1:3] + Frame0_codon1[0]
    - Frame 1 codon at [4,5,6] = Frame0_codon1[1:3] + Frame0_codon2[0]

    To avoid Frame 1 stops, for each pair of adjacent codons (prev, curr):
    - prev[1:3] + curr[0] must not be TAA, TAG, or TGA

    Args:
        protein: Amino acid sequence
        verbose: Print optimization details

    Returns:
        Codon-optimized DNA sequence
    """
    # All codons for each amino acid, ordered by human codon preference
    ALL_CODONS = {
        "A": ["GCC", "GCT", "GCA", "GCG"],
        "C": ["TGC", "TGT"],
        "D": ["GAC", "GAT"],
        "E": ["GAG", "GAA"],
        "F": ["TTC", "TTT"],
        "G": ["GGC", "GGT", "GGA", "GGG"],
        "H": ["CAC", "CAT"],
        "I": ["ATC", "ATT", "ATA"],
        "K": ["AAG", "AAA"],
        "L": ["CTG", "CTC", "CTT", "CTA", "TTG", "TTA"],
        "M": ["ATG"],
        "N": ["AAC", "AAT"],
        "P": ["CCC", "CCT", "CCA", "CCG"],
        "Q": ["CAG", "CAA"],
        "R": ["CGC", "CGG", "AGA", "AGG", "CGT", "CGA"],
        "S": ["AGC", "TCC", "TCT", "AGT", "TCA", "TCG"],
        "T": ["ACC", "ACT", "ACA", "ACG"],
        "V": ["GTG", "GTC", "GTT", "GTA"],
        "W": ["TGG"],
        "Y": ["TAC", "TAT"],
        "*": ["TGA", "TAA", "TAG"],
    }

    def is_stop(codon: str) -> bool:
        return codon.upper() in STOP_CODONS

    def frame1_boundary_codon(prev: str, curr: str) -> str:
        """Frame 1 codon at boundary = prev[1:3] + curr[0]"""
        return prev[1:3] + curr[0]

    n = len(protein)
    if n == 0:
        return ""

    # Iterative approach with local repair
    # First, assign default codons
    codons = [ALL_CODONS[aa][0] for aa in protein]

    # Now fix any Frame 1 stop boundaries
    max_iterations = n * 10  # Prevent infinite loops
    iteration = 0

    while iteration < max_iterations:
        iteration += 1
        found_stop = False

        for i in range(1, n):
            prev = codons[i-1]
            curr = codons[i]
            boundary = frame1_boundary_codon(prev, curr)

            if is_stop(boundary):
                found_stop = True
                # Try to fix by changing current codon
                fixed = False
                for alt_curr in ALL_CODONS[protein[i]]:
                    if not is_stop(frame1_boundary_codon(prev, alt_curr)):
                        codons[i] = alt_curr
                        fixed = True
                        break

                if not fixed:
                    # Try changing previous codon
                    for alt_prev in ALL_CODONS[protein[i-1]]:
                        if not is_stop(frame1_boundary_codon(alt_prev, curr)):
                            # Also check the boundary before alt_prev
                            if i >= 2:
                                prev_prev = codons[i-2]
                                if is_stop(frame1_boundary_codon(prev_prev, alt_prev)):
                                    continue
                            codons[i-1] = alt_prev
                            fixed = True
                            break

                if not fixed:
                    # Try all combinations of prev and curr
                    for alt_prev in ALL_CODONS[protein[i-1]]:
                        for alt_curr in ALL_CODONS[protein[i]]:
                            if not is_stop(frame1_boundary_codon(alt_prev, alt_curr)):
                                # Check boundary before
                                if i >= 2:
                                    prev_prev = codons[i-2]
                                    if is_stop(frame1_boundary_codon(prev_prev, alt_prev)):
                                        continue
                                codons[i-1] = alt_prev
                                codons[i] = alt_curr
                                fixed = True
                                break
                        if fixed:
                            break

                # Don't warn for every unfixable stop - these are likely unavoidable
                pass

        if not found_stop:
            break

    dna = "".join(codons)

    if verbose:
        # Report any remaining Frame 1 stops
        stops = find_stop_codons_in_frame(dna, frame=1)
        if stops:
            print(f"Frame 1 stops remaining: {len(stops)}")
            for pos, codon in stops[:10]:  # Show first 10
                print(f"  Position {pos}: {codon}")
            if len(stops) > 10:
                print(f"  ... and {len(stops) - 10} more")
        else:
            print("No Frame 1 stops - optimization successful!")

    return dna


def codon_optimize_for_dual_frame(protein: str, avoid_frames: List[int] = [0, 1]) -> str:
    """
    Codon-optimize a protein sequence to avoid stop codons in multiple frames.

    For now, this focuses on Frame 1 which is the critical frame for our
    partial intron retention design.

    Args:
        protein: Amino acid sequence
        avoid_frames: List of frames that should not contain stops

    Returns:
        Codon-optimized DNA sequence
    """
    return codon_optimize_for_frame1_stop_free(protein, verbose=False)


def fix_cryptic_splice_donors(dna: str, protein: str, weaken_all_gt: bool = True) -> str:
    """
    Fix cryptic splice donor sites in DNA sequence.

    This function:
    1. Eliminates AG|GT consensus patterns (strongest cryptic donors)
    2. Optionally weakens remaining GT sites by reducing H-bonds at position +6

    A cryptic splice donor can cause aberrant splicing. The consensus is AG|GT
    where AG is at the end of one codon and GT is at the start of the next.
    Even without AG|GT, GT sites can still be recognized if the surrounding
    context is splice-site-like.

    Args:
        dna: DNA sequence (must be in-frame, length divisible by 3)
        protein: Corresponding protein sequence
        weaken_all_gt: If True, also weaken GT sites by changing position +6

    Returns:
        Tuple of (fixed DNA sequence, list of fixes made)
    """
    if len(dna) != len(protein) * 3:
        raise ValueError(f"DNA length {len(dna)} doesn't match protein length {len(protein)} * 3")

    # Codon alternatives for each amino acid
    ALL_CODONS = {
        "A": ["GCC", "GCT", "GCA", "GCG"],
        "C": ["TGC", "TGT"],
        "D": ["GAC", "GAT"],
        "E": ["GAA", "GAG"],  # GAA first to avoid AG ending
        "F": ["TTC", "TTT"],
        "G": ["GGC", "GGT", "GGA", "GGG"],
        "H": ["CAC", "CAT"],
        "I": ["ATC", "ATT", "ATA"],
        "K": ["AAG", "AAA"],
        "L": ["CTG", "CTC", "CTT", "CTA", "TTG", "TTA"],
        "M": ["ATG"],
        "N": ["AAC", "AAT"],
        "P": ["CCC", "CCT", "CCA", "CCG"],
        "Q": ["CAA", "CAG"],  # CAA first to avoid G at end (weakens +6)
        "R": ["CGC", "CGG", "AGA", "AGG", "CGT", "CGA"],
        "S": ["AGC", "TCC", "TCT", "AGT", "TCA", "TCG"],
        "T": ["ACC", "ACT", "ACA", "ACG"],
        "V": ["GTC", "GTT", "GTA", "GTG"],  # GTG last (creates AG|GT with preceding AG)
        "W": ["TGG"],
        "Y": ["TAC", "TAT"],
    }

    # Convert to list of codons
    codons = [dna[i:i+3] for i in range(0, len(dna), 3)]
    fixes_made = []

    # Pass 1: Fix AG|GT patterns at codon boundaries
    for i in range(len(codons) - 1):
        codon1 = codons[i]
        codon2 = codons[i + 1]

        # Check for AG|GT pattern
        if codon1[1:3] == "AG" and codon2[0:2] == "GT":
            aa1 = protein[i]
            aa2 = protein[i + 1]

            # Strategy 1: Change codon1 to not end in AG
            fixed = False
            for alt_codon1 in ALL_CODONS.get(aa1, []):
                if alt_codon1[1:3] != "AG":
                    # Verify it doesn't create a new problem with previous codon
                    if i > 0 and codons[i-1][1:3] == "AG" and alt_codon1[0:2] == "GT":
                        continue
                    codons[i] = alt_codon1
                    fixes_made.append((i, codon1, alt_codon1, f"AG ending -> {alt_codon1[1:3]}"))
                    fixed = True
                    break

            # Strategy 2: If can't fix codon1, try changing codon2 to not start with GT
            if not fixed:
                for alt_codon2 in ALL_CODONS.get(aa2, []):
                    if alt_codon2[0:2] != "GT":
                        codons[i + 1] = alt_codon2
                        fixes_made.append((i + 1, codon2, alt_codon2, f"GT start -> {alt_codon2[0:2]}"))
                        fixed = True
                        break

    # Pass 2: Weaken remaining GT sites by changing position +6
    # For GT at codon[i+1][0:2]:
    #   Position +3 = codon[i+1][2]
    #   Position +4 = codon[i+2][0]
    #   Position +5 = codon[i+2][1]
    #   Position +6 = codon[i+2][2]  <-- this is what we want to weaken
    # Strong +6: G (3 H-bonds with C in U1)
    # Weak +6: A, T, C (0 H-bonds - mismatch with C in U1)
    if weaken_all_gt:
        for i in range(len(codons) - 2):
            codon2 = codons[i + 1]  # Codon containing GT at start
            codon3 = codons[i + 2]  # Codon containing position +6 at end

            # Check if codon2 starts with GT (potential splice donor)
            if codon2[0:2] == "GT":
                aa3 = protein[i + 2]

                # Position +6 is the THIRD (last) nucleotide of codon3
                if codon3[2] == "G":  # Strong position +6
                    # Try to change codon3 to NOT end with G
                    for alt_codon3 in ALL_CODONS.get(aa3, []):
                        if alt_codon3[2] != "G":
                            codons[i + 2] = alt_codon3
                            fixes_made.append((i + 2, codon3, alt_codon3,
                                f"weaken +6: {codon3} -> {alt_codon3} (G -> {alt_codon3[2]})"))
                            break

    return "".join(codons), fixes_made


def design_ggggs_linker(repeats: int, avoid_frame1_stops: bool = True) -> str:
    """
    Design GGGGS linker DNA sequence.

    Args:
        repeats: Number of GGGGS repeats
        avoid_frame1_stops: If True, choose codons that avoid Frame 1 stops

    Returns:
        DNA sequence for linker
    """
    if avoid_frame1_stops:
        # Use codon optimization
        protein = GGGGS_UNIT * repeats
        return codon_optimize_for_dual_frame(protein, avoid_frames=[1])
    else:
        # Simple version using common codons
        # GGGGS = GGT GGC GGT GGC AGC (or similar)
        unit = "GGTGGCGGTGGCAGC"  # GGT GGC GGT GGC AGC
        return unit * repeats


def design_vhh3_codon_optimized(fix_cryptic_donors: bool = True) -> str:
    """
    Design VHH3 DNA sequence optimized to avoid Frame 1 stop codons
    and cryptic splice donors.

    Args:
        fix_cryptic_donors: If True, also fix AG|GT cryptic splice donor patterns

    Returns:
        Codon-optimized VHH3 DNA sequence (357 bp)
    """
    # First optimize for Frame 1 stop codon avoidance
    dna = codon_optimize_for_dual_frame(VHH3_PROTEIN, avoid_frames=[1])

    # Then fix any cryptic splice donors (AG|GT patterns)
    if fix_cryptic_donors:
        dna, fixes = fix_cryptic_splice_donors(dna, VHH3_PROTEIN)
        # Fixes will contain any changes made (e.g., GAG->GAA at position 0)

    return dna


# Upstream exonic context - use β-globin exon 2 (realistic context)
# β-globin exon 2: ~223 bp, we use last 60 bp
BETA_GLOBIN_EXON2_END = (
    "GTGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGT"
)

# Default upstream exon: last 50 bp of a generic coding sequence ending in-frame
# This provides context for splice site recognition
# CRITICAL FIX (v5): The original exon had a cryptic AAG|GT splice donor at position 10.
# The new sequence replaces AAGGT→AAGAA and GTG→ATG to eliminate all GGT patterns.
DEFAULT_UPSTREAM_EXON = (
    "ATGGCCAAGAACACCCTGATCAACATCAAGGCCGAGATCAAGAACGACGAG"  # 51 bp, no cryptic donors
)


def build_intron_cassette(
    splice_donor_11mer: str,
    ggggs5_dna: str,
    vhh3_dna: str,
    ggggs1_dna: str,
    splice_acceptor: str = None,
    upstream_exon: str = None,
    downstream_exon: str = None
) -> str:
    """
    Build a complete intron cassette for splicing prediction.

    Structure:
    [upstream_exon][donor_exonic]gt[donor_intronic][GGGGS5][VHH3][GGGGS1][acceptor_intronic]ag[downstream_exon]

    The splice donor 11-mer positions:
    - Positions -3 to -1: Last 3 nt of upstream exon (exonic)
    - Positions +1 to +8: First 8 nt of intron (intronic, starts with GT)

    Args:
        splice_donor_11mer: 11-mer splice donor sequence (positions -3 to +8)
        ggggs5_dna: GGGGS5 linker DNA
        vhh3_dna: VHH3 DNA
        ggggs1_dna: GGGGS1 linker DNA
        splice_acceptor: Splice acceptor sequence (intronic + exonic)
        upstream_exon: Upstream exonic sequence (before the -3 position)
        downstream_exon: Additional downstream exonic sequence (optional)

    Returns:
        Complete cassette DNA sequence
    """
    if splice_acceptor is None:
        splice_acceptor = SPLICE_ACCEPTOR_BGLOBIN
    if upstream_exon is None:
        upstream_exon = DEFAULT_UPSTREAM_EXON
    if downstream_exon is None:
        downstream_exon = ""  # Acceptor already includes exonic sequence

    # The splice donor 11-mer: positions -3,-2,-1 are exonic, +1 to +8 are intronic
    # Position +1,+2 MUST be GT for a valid splice donor
    donor_exonic = splice_donor_11mer[0:3]   # Positions -3, -2, -1
    donor_intronic = splice_donor_11mer[3:]  # Positions +1 to +8 (starts with GT)

    # Verify GT at +1/+2
    if donor_intronic[0:2].upper() != "GT":
        raise ValueError(f"Invalid splice donor: must have GT at +1/+2, got {donor_intronic[0:2]}")

    # Build the cassette
    # Structure: [exon][GT...intron...AG][exon]
    cassette = (
        upstream_exon +           # Upstream exon (before -3)
        donor_exonic +            # Exonic part of donor (-3 to -1)
        # --- SPLICE JUNCTION (exon|intron) ---
        donor_intronic +          # Intronic part of donor (+1 to +8, starts with GT)
        ggggs5_dna +              # GGGGS5 linker (intron body)
        vhh3_dna +                # VHH3 nanobody (intron body)
        ggggs1_dna +              # GGGGS1 linker (intron body)
        # --- SPLICE JUNCTION (intron|exon) ---
        splice_acceptor +         # Acceptor: intronic (ends with AG) + exonic
        downstream_exon           # Additional downstream exon (if any)
    )

    return cassette


def build_minimal_cassette(
    splice_donor_11mer: str,
    intron_body: str = None,
    splice_acceptor: str = None,
    upstream_exon: str = None
) -> str:
    """
    Build a minimal cassette for basic splicing prediction.

    This is useful for quick tests that focus on the splice donor strength
    without the full VHH3 insert.

    Uses the β-globin acceptor for realistic splice site recognition.
    """
    if intron_body is None:
        # Minimal intron body: 80 bp of generic sequence
        # Introns need to be >80 bp for efficient splicing
        intron_body = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"

    if splice_acceptor is None:
        splice_acceptor = SPLICE_ACCEPTOR_BGLOBIN_SHORT

    if upstream_exon is None:
        upstream_exon = DEFAULT_UPSTREAM_EXON

    donor_exonic = splice_donor_11mer[0:3]
    donor_intronic = splice_donor_11mer[3:]

    # Verify GT at +1/+2
    if donor_intronic[0:2].upper() != "GT":
        raise ValueError(f"Invalid splice donor: must have GT at +1/+2, got {donor_intronic[0:2]}")

    return (
        upstream_exon +       # Upstream exon
        donor_exonic +        # Exonic donor (-3 to -1)
        donor_intronic +      # Intronic donor (+1 to +8, starts with GT)
        intron_body +         # Minimal intron body
        splice_acceptor       # β-globin acceptor (intronic AG + exonic)
    )


def read_panel_csv(filepath: str) -> List[Dict]:
    """Read the tiered panel CSV file."""
    candidates = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            candidates.append(row)
    return candidates


def write_fasta(sequences: List[Tuple[str, str]], filepath: str) -> None:
    """
    Write sequences to FASTA format.

    Args:
        sequences: List of (header, sequence) tuples
        filepath: Output file path
    """
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)

    with open(filepath, 'w') as f:
        for header, seq in sequences:
            f.write(f">{header}\n")
            # Wrap sequence at 60 characters
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Build splice cassettes for splicing prediction software",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Build full cassettes with VHH3
  python build_splice_cassettes.py --panel tiered_panel_candidates.csv --output cassettes.fasta

  # Build minimal cassettes for quick testing
  python build_splice_cassettes.py --panel tiered_panel_candidates.csv --output minimal.fasta --minimal

  # Show VHH3 codon optimization
  python build_splice_cassettes.py --show-vhh3

  # Score cassettes with NNSplice (BDGP neural network)
  python build_splice_cassettes.py --score-nnsplice --fasta cassettes.fasta --panel panel.csv

  # Score cassettes with BOTH MaxEntScan and NNSplice (comprehensive)
  python build_splice_cassettes.py --score-all --fasta cassettes.fasta --panel panel.csv
        """
    )

    parser.add_argument("--panel", type=str,
                        help="Path to tiered panel CSV file")
    parser.add_argument("--output", type=str, default="splice_cassettes.fasta",
                        help="Output FASTA file")
    parser.add_argument("--minimal", action="store_true",
                        help="Build minimal cassettes (no VHH3)")
    parser.add_argument("--show-vhh3", action="store_true",
                        help="Show VHH3 codon optimization details")
    parser.add_argument("--score-nnsplice", action="store_true",
                        help="Score cassettes with NNSplice (BDGP) neural network model")
    parser.add_argument("--score-all", action="store_true",
                        help="Score cassettes with BOTH MaxEntScan and NNSplice (comprehensive)")
    parser.add_argument("--fasta", type=str,
                        help="FASTA file to score (for --score-nnsplice or --score-all)")
    parser.add_argument("--output-dir", type=str, default="./output",
                        help="Output directory")

    args = parser.parse_args()

    # Show VHH3 optimization
    if args.show_vhh3:
        print("VHH3 Codon Optimization")
        print("=" * 60)
        print(f"\nProtein sequence ({len(VHH3_PROTEIN)} aa):")
        print(VHH3_PROTEIN)

        # First, identify unavoidable stops (dipeptide constraints)
        print("\n" + "-" * 60)
        print("Analyzing unavoidable Frame 1 stops (dipeptide constraints)...")
        unavoidable = find_unavoidable_frame1_stops(VHH3_PROTEIN)
        if unavoidable:
            print(f"\nFound {len(unavoidable)} unavoidable stop(s):")
            for pos, aa1, aa2, reason in unavoidable:
                print(f"  Position {pos}-{pos+1}: {aa1}-{aa2} ({reason})")
        else:
            print("No unavoidable stops found.")

        print("\n" + "-" * 60)
        print("Generating codon-optimized DNA...")
        vhh3_dna = codon_optimize_for_frame1_stop_free(VHH3_PROTEIN, verbose=False)
        print(f"\nDNA sequence ({len(vhh3_dna)} bp):")
        print(vhh3_dna)

        # Verify translation
        translated = translate(vhh3_dna)
        print(f"\nTranslated ({len(translated)} aa):")
        print(translated)

        match = translated == VHH3_PROTEIN
        print(f"Translation matches original: {match}")

        # Check for Frame 1 stops
        frame1_stops = find_stop_codons_in_frame(vhh3_dna, frame=1)
        print(f"\nFrame 1 stop codons: {len(frame1_stops)}")
        if frame1_stops:
            for pos, codon in frame1_stops:
                # Calculate which protein position this corresponds to
                protein_pos = pos // 3
                print(f"  DNA position {pos}: {codon} (between protein positions {protein_pos} and {protein_pos+1})")
        else:
            print("  None - VHH3 is Frame 1 stop-free!")

        # Check Frame 0
        frame0_stops = find_stop_codons_in_frame(vhh3_dna, frame=0)
        print(f"\nFrame 0 stop codons: {len(frame0_stops)}")
        if frame0_stops:
            for pos, codon in frame0_stops:
                print(f"  Position {pos}: {codon}")
        else:
            print("  None")

        return 0

    # Score cassettes with NNSplice
    if args.score_nnsplice:
        if not args.fasta:
            print("Error: --score-nnsplice requires --fasta")
            return 1
        if not args.panel:
            print("Error: --score-nnsplice requires --panel")
            return 1

        fasta_path = os.path.join(args.output_dir, args.fasta) if not os.path.isabs(args.fasta) else args.fasta
        panel_path = os.path.join(args.output_dir, args.panel) if not os.path.isabs(args.panel) else args.panel
        output_csv = os.path.join(args.output_dir, "nnsplice_scores.csv")

        print(f"Scoring cassettes with NNSplice (BDGP)...")
        print(f"  FASTA: {fasta_path}")
        print(f"  Panel: {panel_path}")
        print(f"  Output: {output_csv}")
        print()

        score_all_cassettes_nnsplice(fasta_path, panel_path, output_csv)
        return 0

    # Comprehensive scoring with both algorithms
    if args.score_all:
        if not args.fasta:
            print("Error: --score-all requires --fasta")
            return 1
        if not args.panel:
            print("Error: --score-all requires --panel")
            return 1

        fasta_path = os.path.join(args.output_dir, args.fasta) if not os.path.isabs(args.fasta) else args.fasta
        panel_path = os.path.join(args.output_dir, args.panel) if not os.path.isabs(args.panel) else args.panel
        output_csv = os.path.join(args.output_dir, "splice_prediction_scores.csv")

        print(f"Comprehensive scoring with MaxEntScan and NNSplice...")
        print(f"  FASTA: {fasta_path}")
        print(f"  Panel: {panel_path}")
        print(f"  Output: {output_csv}")
        print()

        score_all_cassettes_comprehensive(fasta_path, panel_path, output_csv)
        return 0

    if not args.panel:
        parser.print_help()
        return 1

    # Read panel
    print(f"Reading panel from: {args.panel}")
    candidates = read_panel_csv(args.panel)
    print(f"Found {len(candidates)} candidates")

    # Design components
    print("\nDesigning intron components...")

    if args.minimal:
        print("Building MINIMAL cassettes (no VHH3)")
        sequences = []

        for c in candidates:
            seq_name = f"Tier{c['tier']}_{c['sequence']}_{c['total_hbonds']}hb_minimal"
            cassette = build_minimal_cassette(c['sequence'])
            sequences.append((seq_name, cassette))

        print(f"Built {len(sequences)} minimal cassettes")

    else:
        print("Building FULL cassettes with VHH3")

        # Generate optimized sequences
        print("  Optimizing VHH3 for Frame 1...")
        vhh3_dna = design_vhh3_codon_optimized()
        frame1_stops = find_stop_codons_in_frame(vhh3_dna, frame=1)
        print(f"    VHH3 length: {len(vhh3_dna)} bp")
        print(f"    Frame 1 stops in VHH3: {len(frame1_stops)}")

        print("  Designing GGGGS5 linker...")
        ggggs5_dna = design_ggggs_linker(5)
        print(f"    GGGGS5 length: {len(ggggs5_dna)} bp")

        print("  Designing GGGGS1 linker...")
        ggggs1_dna = design_ggggs_linker(1)
        print(f"    GGGGS1 length: {len(ggggs1_dna)} bp")

        # Build cassettes using β-globin acceptor
        print("  Using β-globin intron 2 acceptor...")
        print(f"    Acceptor length: {len(SPLICE_ACCEPTOR_BGLOBIN)} bp")

        sequences = []

        for c in candidates:
            seq_name = f"Tier{c['tier']}_{c['sequence']}_{c['total_hbonds']}hb"

            cassette = build_intron_cassette(
                splice_donor_11mer=c['sequence'],
                ggggs5_dna=ggggs5_dna,
                vhh3_dna=vhh3_dna,
                ggggs1_dna=ggggs1_dna,
                splice_acceptor=SPLICE_ACCEPTOR_BGLOBIN
            )

            sequences.append((seq_name, cassette))

        print(f"\nBuilt {len(sequences)} full cassettes")

        # Print cassette statistics
        if sequences:
            example_cassette = sequences[0][1]
            example_len = len(example_cassette)
            print(f"Cassette length: {example_len} bp")

            # Calculate positions
            upstream_len = len(DEFAULT_UPSTREAM_EXON)
            donor_exonic_len = 3
            donor_intronic_len = 8
            intron_body_len = len(ggggs5_dna) + len(vhh3_dna) + len(ggggs1_dna)

            # Find AG in acceptor (intronic part ends with AG)
            acceptor_intronic_len = 20  # "TTCTATTTTCCCACCCTTAG" (ends with AG)
            acceptor_exonic_len = len(BETA_GLOBIN_EXON3_START)

            total_intron = donor_intronic_len + intron_body_len + acceptor_intronic_len
            print(f"Total intron length: {total_intron} bp")

            # Show structure
            print(f"\nCassette structure:")
            print(f"  Upstream exon: {upstream_len} bp")
            print(f"  Donor exonic (-3 to -1): {donor_exonic_len} bp")
            print(f"  [SPLICE DONOR GT at position {upstream_len + donor_exonic_len + 1}]")
            print(f"  Donor intronic (+1 to +8): {donor_intronic_len} bp")
            print(f"  Intron body (GGGGS5+VHH3+GGGGS1): {intron_body_len} bp")
            print(f"  Acceptor intronic (PPT+AG): {acceptor_intronic_len} bp")
            print(f"  [SPLICE ACCEPTOR AG at position {upstream_len + donor_exonic_len + total_intron - 1}]")
            print(f"  Acceptor exonic (exon 3): {acceptor_exonic_len} bp")

    # Write output
    output_path = os.path.join(args.output_dir, args.output)
    write_fasta(sequences, output_path)
    print(f"\nWrote cassettes to: {output_path}")

    # Also write a summary
    summary_path = os.path.join(args.output_dir, args.output.replace('.fasta', '_summary.txt'))
    with open(summary_path, 'w') as f:
        f.write("Splice Cassette Summary\n")
        f.write("=" * 60 + "\n\n")

        for seq_name, seq in sequences:
            f.write(f"{seq_name}\n")
            f.write(f"  Length: {len(seq)} bp\n")
            # Extract donor from name
            parts = seq_name.split('_')
            if len(parts) >= 2:
                donor = parts[1]
                f.write(f"  Donor: {donor}\n")
            f.write("\n")

    print(f"Wrote summary to: {summary_path}")

    return 0


# ============================================================================
# MaxEntScan Scoring Functions
# ============================================================================

# MaxEntScan donor scoring matrices (Yeo & Burge, 2004)
# Based on the 9-mer at positions -3 to +6 (where GT is at +1/+2)
# Scores are log2 odds ratios

# Position weight matrix for human donor sites
# Derived from splice site databases
# Format: {position: {nucleotide: score}}
MAXENT_DONOR_PWM = {
    # Position -3 (index 0 in 9-mer)
    -3: {'A': 0.24, 'C': 0.31, 'G': 0.09, 'T': 0.36},
    # Position -2 (index 1 in 9-mer)
    -2: {'A': 0.63, 'C': 0.10, 'G': 0.12, 'T': 0.15},
    # Position -1 (index 2 in 9-mer) - G strongly preferred (AG|GT consensus)
    -1: {'A': 0.10, 'C': 0.03, 'G': 0.79, 'T': 0.08},
    # Position +1 (index 3 in 9-mer) - invariant G
    1: {'A': 0.00, 'C': 0.00, 'G': 1.00, 'T': 0.00},
    # Position +2 (index 4 in 9-mer) - invariant T
    2: {'A': 0.00, 'C': 0.00, 'G': 0.00, 'T': 1.00},
    # Position +3 (index 5 in 9-mer)
    3: {'A': 0.57, 'C': 0.03, 'G': 0.36, 'T': 0.04},
    # Position +4 (index 6 in 9-mer)
    4: {'A': 0.72, 'C': 0.06, 'G': 0.12, 'T': 0.10},
    # Position +5 (index 7 in 9-mer)
    5: {'A': 0.09, 'C': 0.06, 'G': 0.76, 'T': 0.09},
    # Position +6 (index 8 in 9-mer)
    6: {'A': 0.17, 'C': 0.13, 'G': 0.13, 'T': 0.57},
}

# Background nucleotide frequencies
BACKGROUND_FREQ = {'A': 0.27, 'C': 0.23, 'G': 0.23, 'T': 0.27}

import math


def score_donor_maxentscan(sequence_9mer: str) -> float:
    """
    Score a 9-mer donor sequence using MaxEntScan-like scoring.

    This is a simplified implementation based on position weight matrices.
    The actual MaxEntScan uses a more complex maximum entropy model.

    Args:
        sequence_9mer: 9-character DNA string (positions -3 to +6)

    Returns:
        MaxEntScan-like score (higher = stronger donor)
        Typical range: -20 to +12 bits
    """
    if len(sequence_9mer) != 9:
        raise ValueError(f"Sequence must be 9 nucleotides, got {len(sequence_9mer)}")

    sequence_9mer = sequence_9mer.upper()

    # Check for valid GT at positions +1/+2 (indices 3,4)
    if sequence_9mer[3:5] != 'GT':
        return -999.0  # Invalid donor

    score = 0.0
    positions = [-3, -2, -1, 1, 2, 3, 4, 5, 6]

    for i, pos in enumerate(positions):
        nt = sequence_9mer[i]
        if nt not in 'ACGT':
            return -999.0  # Invalid nucleotide

        # Log odds score: log2(P(nt|splice site) / P(nt|background))
        p_site = MAXENT_DONOR_PWM[pos].get(nt, 0.001)
        p_bg = BACKGROUND_FREQ.get(nt, 0.25)

        if p_site > 0:
            score += math.log2(p_site / p_bg)
        else:
            score += -10  # Penalty for zero probability

    return score


def extract_donor_9mer(cassette_seq: str, donor_pos: int) -> str:
    """
    Extract the 9-mer donor sequence from a cassette.

    Args:
        cassette_seq: Full cassette sequence
        donor_pos: 0-indexed position of the G in GT

    Returns:
        9-mer sequence (positions -3 to +6)
    """
    # 9-mer spans from -3 to +6 relative to splice site
    # GT is at positions +1/+2, so GT is at donor_pos and donor_pos+1
    # Position -3 is at donor_pos - 3
    start = donor_pos - 3
    end = donor_pos + 6  # +6 inclusive, so end at +7 exclusive

    if start < 0 or end > len(cassette_seq):
        raise ValueError(f"Cannot extract 9-mer: start={start}, end={end}, seq_len={len(cassette_seq)}")

    return cassette_seq[start:end]


def classify_maxent_score(score: float) -> str:
    """Classify MaxEntScan score into strength category."""
    if score >= 9.0:
        return "Strong"
    elif score >= 6.0:
        return "Moderate"
    elif score >= 3.0:
        return "Weak"
    elif score >= 0.0:
        return "Very Weak"
    else:
        return "Not Detected"


# ============================================================================
# MaxEntScan Acceptor Scoring Functions
# ============================================================================

# MaxEntScan acceptor scoring matrices (Yeo & Burge, 2004)
# Based on the 23-mer at positions -20 to +3 (AG at positions -2/-1, indices 18-19)
# Scores are log2 odds ratios
#
# Position mapping: index 0 = position -20, index 18 = position -2 (A of AG),
#                   index 19 = position -1 (G of AG), index 22 = position +3

MAXENT_ACCEPTOR_PWM = {
    # Polypyrimidine tract region (positions -20 to -4): enriched for C/T
    -20: {'A': 0.09, 'C': 0.28, 'G': 0.07, 'T': 0.56},
    -19: {'A': 0.08, 'C': 0.30, 'G': 0.07, 'T': 0.55},
    -18: {'A': 0.09, 'C': 0.29, 'G': 0.07, 'T': 0.55},
    -17: {'A': 0.09, 'C': 0.27, 'G': 0.08, 'T': 0.56},
    -16: {'A': 0.09, 'C': 0.27, 'G': 0.08, 'T': 0.56},
    -15: {'A': 0.10, 'C': 0.26, 'G': 0.08, 'T': 0.56},
    -14: {'A': 0.10, 'C': 0.26, 'G': 0.08, 'T': 0.56},
    -13: {'A': 0.09, 'C': 0.27, 'G': 0.08, 'T': 0.56},
    -12: {'A': 0.10, 'C': 0.26, 'G': 0.08, 'T': 0.56},
    -11: {'A': 0.10, 'C': 0.26, 'G': 0.09, 'T': 0.55},
    -10: {'A': 0.10, 'C': 0.25, 'G': 0.09, 'T': 0.56},
    -9:  {'A': 0.10, 'C': 0.26, 'G': 0.08, 'T': 0.56},
    -8:  {'A': 0.10, 'C': 0.24, 'G': 0.09, 'T': 0.57},
    -7:  {'A': 0.10, 'C': 0.23, 'G': 0.08, 'T': 0.59},
    -6:  {'A': 0.10, 'C': 0.24, 'G': 0.08, 'T': 0.58},
    -5:  {'A': 0.10, 'C': 0.26, 'G': 0.07, 'T': 0.57},
    -4:  {'A': 0.12, 'C': 0.30, 'G': 0.06, 'T': 0.52},
    # YAG consensus (positions -3 to -1)
    -3:  {'A': 0.08, 'C': 0.44, 'G': 0.03, 'T': 0.45},  # Y of YAG
    # Invariant AG at positions -2, -1 (from MaxEntScan: Yeo & Burge 2004)
    -2:  {'A': 0.9903, 'C': 0.0032, 'G': 0.0034, 'T': 0.0030},
    -1:  {'A': 0.0027, 'C': 0.0037, 'G': 0.9905, 'T': 0.0030},
    # Exonic positions +1 to +3
    1:   {'A': 0.23, 'C': 0.16, 'G': 0.42, 'T': 0.19},
    2:   {'A': 0.26, 'C': 0.21, 'G': 0.24, 'T': 0.29},
    3:   {'A': 0.25, 'C': 0.22, 'G': 0.24, 'T': 0.29},
}


def score_acceptor_maxentscan(sequence_23mer: str) -> float:
    """
    Score a 23-mer acceptor sequence using MaxEntScan-like scoring.

    The 23-mer spans positions -20 to +3, with the AG dinucleotide at
    indices 18-19 (positions -2 and -1).

    Args:
        sequence_23mer: 23-character DNA string (intron[-20:] + exon[:3])

    Returns:
        MaxEntScan-like score (higher = stronger acceptor)
        Typical range: ~-15 to +12
    """
    if len(sequence_23mer) != 23:
        raise ValueError(f"Sequence must be 23 nucleotides, got {len(sequence_23mer)}")

    sequence_23mer = sequence_23mer.upper()

    # Check for valid AG at indices 18-19
    if sequence_23mer[18:20] != 'AG':
        return -999.0  # Invalid acceptor

    score = 0.0
    positions = list(range(-20, 0)) + [1, 2, 3]  # -20..-1, +1, +2, +3

    for i, pos in enumerate(positions):
        nt = sequence_23mer[i]
        if nt not in 'ACGT':
            return -999.0  # Invalid nucleotide

        p_site = MAXENT_ACCEPTOR_PWM[pos].get(nt, 0.001)
        p_bg = BACKGROUND_FREQ.get(nt, 0.25)

        if p_site > 0:
            score += math.log2(p_site / p_bg)
        else:
            score += -10  # Penalty for zero probability

    return score


def classify_maxent_acceptor_score(score: float) -> str:
    """Classify MaxEntScan acceptor score into strength category."""
    if score >= 9.0:
        return "Strong"
    elif score >= 6.0:
        return "Moderate"
    elif score >= 3.0:
        return "Weak"
    elif score >= 0.0:
        return "Very Weak"
    else:
        return "Not Detected"


# ============================================================================
# NNSplice (BDGP) Scoring Functions
# ============================================================================

def get_nnsplice_binary_path() -> Path:
    """Get path to the NNSplice donor predictor binary."""
    return Path(__file__).parent / "bin" / "fa2Donorpred"


def score_single_sequence_nnsplice(
    seq_name: str,
    seq: str,
    threshold: float = 0.1
) -> List[dict]:
    """
    Score a single sequence using NNSplice (BDGP neural network).

    NNSplice must be run on one sequence at a time as it treats multi-FASTA
    files as a single concatenated sequence.

    Args:
        seq_name: Sequence name
        seq: DNA sequence
        threshold: Minimum score threshold (default 0.1)

    Returns:
        List of dicts with position, score, pattern for each hit
    """
    import subprocess
    import re
    import tempfile

    nnsplice_bin = get_nnsplice_binary_path()

    if not nnsplice_bin.exists():
        raise FileNotFoundError(
            f"NNSplice binary not found at {nnsplice_bin}. "
            "Please build it from https://github.com/bdgp/bdgp-splice-promoter-tools"
        )

    # Write single sequence to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(f">{seq_name}\n{seq}\n")
        temp_path = f.name

    try:
        result = subprocess.run(
            [str(nnsplice_bin), "-t", str(threshold), temp_path],
            capture_output=True, text=True
        )
    finally:
        import os
        os.unlink(temp_path)

    if result.returncode != 0:
        raise RuntimeError(f"NNSplice failed: {result.stderr}")

    # Parse output
    # Format:
    # Prediction for sequence_name (N bases):
    # Hit N: position: X .. Y, exon end = Z
    # Pattern: ACGTACGTACGTACG
    # Prediction: 0.985-----------------------------------------
    hits = []
    current_pattern = None
    current_exon_end = None

    lines = result.stdout.strip().split('\n')
    for line in lines:
        line = line.strip()

        # Parse hit info: "Hit N: position: X .. Y, exon end = Z"
        if line.startswith('Hit'):
            pos_match = re.search(r'exon end = (\d+)', line)
            if pos_match:
                current_exon_end = int(pos_match.group(1))
            continue

        # Parse pattern: "Pattern: ACGTACGTACGTACG"
        if line.startswith('Pattern:'):
            current_pattern = line.split('Pattern:')[1].strip()
            continue

        # Parse prediction: "Prediction: 0.985-----------------------------------------"
        if line.startswith('Prediction:'):
            pred_part = line.split('Prediction:')[1].strip()
            # Score is before the dashes
            score_str = pred_part.split('-')[0].strip()
            try:
                score = float(score_str)
            except ValueError:
                score = 0.0

            if current_exon_end is not None:
                hits.append({
                    'sequence_name': seq_name,
                    'position': current_exon_end,  # Position of last exonic nt (1-indexed)
                    'score': score,
                    'pattern': current_pattern or "",
                })

            # Reset for next hit
            current_pattern = None
            current_exon_end = None

    return hits


def score_with_nnsplice(
    fasta_path: str,
    threshold: float = 0.1
) -> List[dict]:
    """
    Score sequences using NNSplice (BDGP neural network).

    NNSplice uses a neural network trained on human splice sites to predict
    splice donor and acceptor probabilities. Unlike SpliceAI, it works on
    short sequences without requiring 10kb genomic context.

    Note: NNSplice processes sequences one at a time as it treats multi-FASTA
    files as concatenated sequences.

    Args:
        fasta_path: Path to FASTA file
        threshold: Minimum score threshold (default 0.1)

    Returns:
        List of dicts with sequence name, position, score, pattern for each hit
    """
    # Read all sequences from FASTA
    sequences = {}
    current_name = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            sequences[current_name] = ''.join(current_seq)

    # Score each sequence individually
    all_hits = []
    for seq_name, seq in sequences.items():
        hits = score_single_sequence_nnsplice(seq_name, seq, threshold)
        all_hits.extend(hits)

    return all_hits


def score_all_cassettes_nnsplice(
    fasta_path: str,
    panel_csv_path: str,
    output_csv_path: str,
    expected_donor_position: int = 54,  # 0-indexed position of GT in default cassette
) -> None:
    """
    Score all cassettes in a FASTA file using NNSplice (BDGP).

    Args:
        fasta_path: Path to FASTA file with cassettes
        panel_csv_path: Path to tiered panel CSV (for merging results)
        output_csv_path: Path to output CSV with NNSplice scores
        expected_donor_position: Expected 0-indexed position of designed splice donor
    """
    # Read cassettes from FASTA
    cassettes = {}
    current_name = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    cassettes[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            cassettes[current_name] = ''.join(current_seq)

    print(f"Loaded {len(cassettes)} cassettes from {fasta_path}")

    # Load panel CSV
    panel_data = {}
    with open(panel_csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            seq = row['sequence']
            panel_data[seq] = row

    # Calculate expected exon_end position for the designed splice donor
    # NNSplice reports 1-indexed positions
    # Upstream exon: 51 bp (positions 1-51)
    # Donor exonic (-3,-2,-1): 3 bp (positions 52-54)
    # GT at +1/+2: positions 55-56
    # So the exon_end (1-indexed position of last exonic nt before GT) = 54
    upstream_len = len(DEFAULT_UPSTREAM_EXON)  # 51 bp
    donor_exonic_len = 3  # -3, -2, -1
    expected_exon_end_1indexed = upstream_len + donor_exonic_len  # 54 (1-indexed)

    # Score with NNSplice
    print("Running NNSplice prediction...")
    try:
        nnsplice_hits = score_with_nnsplice(fasta_path, threshold=0.1)
        print(f"NNSplice found {len(nnsplice_hits)} hits")
    except Exception as e:
        print(f"Error running NNSplice: {e}")
        nnsplice_hits = []

    # Group hits by sequence name
    hits_by_seq = {}
    for hit in nnsplice_hits:
        seq_name = hit['sequence_name']
        if seq_name not in hits_by_seq:
            hits_by_seq[seq_name] = []
        hits_by_seq[seq_name].append(hit)

    # Build results
    results = []
    for name, seq in cassettes.items():
        # Extract the 11-mer donor sequence from the name
        # Name format: Tier{N}_{11mer}_{hbonds}hb[_minimal]
        parts = name.split('_')
        donor_11mer = parts[1] if len(parts) >= 2 else "UNKNOWN"

        # Get panel data for this sequence
        panel_row = panel_data.get(donor_11mer, {})

        # Get NNSplice hits for this cassette
        seq_hits = hits_by_seq.get(name, [])

        # Find the hit at the designed donor position
        designed_hit = None
        cryptic_hits = []

        for hit in seq_hits:
            # NNSplice exon_end is 1-indexed, expected_exon_end_1indexed is also 1-indexed
            if abs(hit['position'] - expected_exon_end_1indexed) <= 2:  # Allow small tolerance
                designed_hit = hit
            else:
                cryptic_hits.append(hit)

        # Determine cryptic warning
        if cryptic_hits:
            cryptic_warning = f"Warning: {len(cryptic_hits)} cryptic donor(s) at pos {[h['position'] for h in cryptic_hits]}"
        else:
            cryptic_warning = "Clean"

        # Get NNSplice score at designed position
        if designed_hit:
            nnsplice_score = designed_hit['score']
            nnsplice_detection = "Detected"
            if nnsplice_score >= 0.9:
                nnsplice_detection = "Strong"
            elif nnsplice_score >= 0.7:
                nnsplice_detection = "Moderate"
            elif nnsplice_score >= 0.4:
                nnsplice_detection = "Weak"
            else:
                nnsplice_detection = "Very Weak"
        else:
            nnsplice_score = 0.0
            nnsplice_detection = "Not Detected"

        result = {
            'tier': panel_row.get('tier', ''),
            'sequence': donor_11mer,
            'hbonds': panel_row.get('total_hbonds', ''),
            'stretch': panel_row.get('longest_stretch', ''),
            'nnsplice_score': f"{nnsplice_score:.4f}",
            'nnsplice_detection': nnsplice_detection,
            'nnsplice_cryptic_check': cryptic_warning,
            'cassette_name': name,
            'cassette_length': len(seq),
        }
        results.append(result)

        score_str = f"{nnsplice_score:.4f}" if nnsplice_score > 0 else "N/A"
        print(f"  {name}: nnsplice_score={score_str}")

    # Write results
    Path(output_csv_path).parent.mkdir(parents=True, exist_ok=True)

    fieldnames = ['tier', 'sequence', 'hbonds', 'stretch',
                  'nnsplice_score', 'nnsplice_detection', 'nnsplice_cryptic_check',
                  'cassette_name', 'cassette_length']

    with open(output_csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"\nWrote NNSplice scores to: {output_csv_path}")


def score_all_cassettes_comprehensive(
    fasta_path: str,
    panel_csv_path: str,
    output_csv_path: str,
) -> None:
    """
    Score all cassettes using both MaxEntScan and NNSplice algorithms.

    This provides comprehensive splice site scoring using:
    1. MaxEntScan - Position weight matrix / maximum entropy scoring
    2. NNSplice (BDGP) - Neural network scoring

    Args:
        fasta_path: Path to FASTA file with cassettes
        panel_csv_path: Path to tiered panel CSV (for merging results)
        output_csv_path: Path to output CSV with comprehensive scores
    """
    # Read cassettes from FASTA
    cassettes = {}
    current_name = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    cassettes[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            cassettes[current_name] = ''.join(current_seq)

    print(f"Loaded {len(cassettes)} cassettes from {fasta_path}")

    # Load panel CSV
    panel_data = {}
    with open(panel_csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            seq = row['sequence']
            panel_data[seq] = row

    # Calculate expected donor position
    upstream_len = len(DEFAULT_UPSTREAM_EXON)  # 51 bp
    donor_exonic_len = 3  # -3, -2, -1
    gt_position_0indexed = upstream_len + donor_exonic_len  # 54 (position of G in GT)
    expected_exon_end_1indexed = gt_position_0indexed  # 54 (1-indexed)

    # Run NNSplice scoring
    print("Running NNSplice prediction...")
    try:
        nnsplice_hits = score_with_nnsplice(fasta_path, threshold=0.1)
        print(f"NNSplice found {len(nnsplice_hits)} total hits")
    except Exception as e:
        print(f"Error running NNSplice: {e}")
        nnsplice_hits = []

    # Group NNSplice hits by sequence name
    hits_by_seq = {}
    for hit in nnsplice_hits:
        seq_name = hit['sequence_name']
        if seq_name not in hits_by_seq:
            hits_by_seq[seq_name] = []
        hits_by_seq[seq_name].append(hit)

    # Score each cassette with both algorithms
    print("\nScoring cassettes with MaxEntScan and NNSplice...")
    results = []

    for name, seq in cassettes.items():
        # Extract donor 11-mer from name
        parts = name.split('_')
        donor_11mer = parts[1] if len(parts) >= 2 else "UNKNOWN"

        # Get panel data
        panel_row = panel_data.get(donor_11mer, {})

        # === MaxEntScan scoring ===
        try:
            donor_9mer = extract_donor_9mer(seq, gt_position_0indexed)
            maxent_score = score_donor_maxentscan(donor_9mer)
            maxent_detection = classify_maxent_score(maxent_score)
        except Exception as e:
            maxent_score = -999.0
            maxent_detection = f"Error: {e}"
            donor_9mer = "N/A"

        # === NNSplice scoring ===
        seq_hits = hits_by_seq.get(name, [])
        designed_hit = None
        cryptic_count = 0

        for hit in seq_hits:
            if abs(hit['position'] - expected_exon_end_1indexed) <= 2:
                designed_hit = hit
            else:
                cryptic_count += 1

        if designed_hit:
            nnsplice_score = designed_hit['score']
            if nnsplice_score >= 0.9:
                nnsplice_detection = "Strong"
            elif nnsplice_score >= 0.7:
                nnsplice_detection = "Moderate"
            elif nnsplice_score >= 0.4:
                nnsplice_detection = "Weak"
            else:
                nnsplice_detection = "Very Weak"
        else:
            nnsplice_score = 0.0
            nnsplice_detection = "Not Detected"

        result = {
            'tier': panel_row.get('tier', ''),
            'sequence': donor_11mer,
            'delta_g': panel_row.get('delta_g', ''),
            'pos6_nt': panel_row.get('pos6_nt', donor_11mer[8] if len(donor_11mer) > 8 else ''),
            'total_hbonds': panel_row.get('total_hbonds', ''),
            'longest_stretch': panel_row.get('longest_stretch', ''),
            'donor_9mer': donor_9mer,
            'maxent_score': f"{maxent_score:.2f}" if maxent_score > -100 else "N/A",
            'maxent_detection': maxent_detection,
            'nnsplice_score': f"{nnsplice_score:.4f}",
            'nnsplice_detection': nnsplice_detection,
            'cryptic_sites': cryptic_count,
            'cassette_name': name,
            'cassette_length': len(seq),
        }
        results.append(result)

        print(f"  {name}: MaxEnt={maxent_score:.2f} ({maxent_detection}), NNSplice={nnsplice_score:.4f} ({nnsplice_detection})")

    # Write results
    Path(output_csv_path).parent.mkdir(parents=True, exist_ok=True)

    fieldnames = ['tier', 'sequence', 'delta_g', 'pos6_nt', 'total_hbonds', 'longest_stretch',
                  'donor_9mer', 'maxent_score', 'maxent_detection', 'nnsplice_score',
                  'nnsplice_detection', 'cryptic_sites', 'cassette_name', 'cassette_length']

    with open(output_csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"\nWrote comprehensive scores to: {output_csv_path}")

    # Also create a nicely formatted summary table
    summary_path = output_csv_path.replace('.csv', '_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("Comprehensive Splice Donor Scoring Results\n")
        f.write("=" * 120 + "\n\n")
        f.write("Algorithms used:\n")
        f.write("  - MaxEntScan: Position weight matrix scoring (Yeo & Burge, 2004)\n")
        f.write("  - NNSplice: Neural network scoring (BDGP, Reese et al., 1997)\n\n")
        f.write("-" * 120 + "\n")
        f.write(f"{'Tier':<5} {'Sequence':<12} {'ΔG':>8} {'+6':>4} {'H-bonds':>8} {'Stretch':>8} {'MaxEnt':>8} {'MaxEnt Det.':<12} {'NNSplice':>10} {'NNSplice Det.':<14}\n")
        f.write("-" * 120 + "\n")

        for r in results:
            f.write(f"{r['tier']:<5} {r['sequence']:<12} {r['delta_g']:>8} {r['pos6_nt']:>4} {r['total_hbonds']:>8} {r['longest_stretch']:>8} "
                    f"{r['maxent_score']:>8} {r['maxent_detection']:<12} {r['nnsplice_score']:>10} {r['nnsplice_detection']:<14}\n")

        f.write("-" * 120 + "\n")
        f.write(f"\nTotal cassettes scored: {len(results)}\n")

    print(f"Wrote summary to: {summary_path}")


if __name__ == "__main__":
    exit(main())
