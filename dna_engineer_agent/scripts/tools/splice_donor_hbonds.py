#!/usr/bin/env python3
"""
Splice Donor H-Bond Prediction Pipeline

Predicts splice donor strength based on hydrogen bonding between the 5' splice site
(5'SS) and U1 snRNA. Generates synthetic sequences with target H-bond counts.

Background:
The 5'SS (positions -3 to +8, where | marks exon/intron boundary) base-pairs with
U1 snRNA in an antiparallel manner:

    5'SS:   -3 -2 -1 | +1 +2 +3 +4 +5 +6 +7 +8
    U1:     11 10  9 |  8  7  6  5  4  3  2  1
    U1nt:    G  U  C |  C  A  Ψ  Ψ  A  C  U  G

H-bond rules:
- G-C or C-G: 3 H-bonds
- A-U or U-A: 2 H-bonds
- G-U or U-G wobble: 2 H-bonds
- G-Ψ or A-Ψ: 2 H-bonds (Ψ behaves like U)
- Any mismatch: 0 H-bonds

Theoretical maximum: 27 H-bonds (all positions optimal)
"""

import argparse
import csv
import itertools
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ============================================================================
# Constants
# ============================================================================

# Position labels for the 11-mer 5'SS
POSITION_LABELS = ["-3", "-2", "-1", "+1", "+2", "+3", "+4", "+5", "+6", "+7", "+8"]

# U1 snRNA partners at each position (5'SS position index 0-10)
# Note: Ψ (pseudouridine) behaves like U for base pairing
U1_PARTNERS = ["G", "U", "C", "C", "A", "Ψ", "Ψ", "A", "C", "U", "G"]

# Optimal DNA nucleotide at each position for maximum H-bonds
# (DNA input is transcribed to RNA for pairing calculation)
OPTIMAL_DNA = ["C", "A", "G", "G", "T", "A", "A", "T", "G", "A", "C"]

# Invariant positions (GT dinucleotide at +1/+2)
INVARIANT_POSITIONS = {3: "G", 4: "T"}  # indices 3 and 4 correspond to +1 and +2

# H-bond rules for base pairs
# Key: (5'SS RNA nucleotide, U1 nucleotide), Value: number of H-bonds
HBOND_RULES: Dict[Tuple[str, str], int] = {
    # Watson-Crick pairs
    ("G", "C"): 3,
    ("C", "G"): 3,
    ("A", "U"): 2,
    ("U", "A"): 2,
    # Wobble pairs
    ("G", "U"): 2,
    ("U", "G"): 2,
    # Pseudouridine pairs (Ψ behaves like U)
    ("G", "Ψ"): 2,
    ("A", "Ψ"): 2,
    # All other combinations are mismatches (0 H-bonds)
}

# DNA to RNA conversion
DNA_TO_RNA = {"A": "A", "T": "U", "G": "G", "C": "C"}

# Maximum theoretical H-bonds
MAX_HBONDS = 27

# Minimum H-bonds for sequence generation
MIN_HBONDS_GENERATE = 10

# Stop codons for Frame 1 filtering
STOP_CODONS = {"TAA", "TAG", "TGA"}

# Standard genetic code (DNA codon -> amino acid)
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

# Default tiers for panel generation
# Format: (min_hbonds, max_hbonds, num_sequences)
DEFAULT_TIERS = [
    (10, 11, 2),   # Minimal - negative control
    (12, 13, 3),   # Primary target range
    (14, 15, 3),   # Primary target range
    (16, 17, 2),   # Moderate
    (18, 19, 2),   # Moderate-strong
    (20, 21, 2),   # Strong (consensus-like)
    (22, 27, 2),   # Maximum - positive control
]

# Expected splicing efficiency by tier (for annotation)
TIER_ANNOTATIONS = {
    (10, 11): ("negative_ctrl", "~0-5%"),
    (12, 13): ("target_range", "~5-15%"),
    (14, 15): ("target_range", "~15-25%"),
    (16, 17): ("moderate", "~25-40%"),
    (18, 19): ("moderate_strong", "~40-60%"),
    (20, 21): ("strong", "~60-80%"),
    (22, 27): ("positive_ctrl", "~80-100%"),
}


# ============================================================================
# Data Classes
# ============================================================================

@dataclass
class HBondResult:
    """Result of H-bond calculation for a 5'SS sequence."""
    sequence: str  # 11-mer DNA sequence
    total_hbonds: int  # Sum of H-bonds across all positions
    per_position: List[int]  # H-bonds at each position (11 values)
    longest_stretch: int  # Length of longest continuous base-paired stretch
    stretch_start: int  # Index where longest stretch starts
    num_base_pairs: int  # Number of positions with >0 H-bonds
    notes: str = ""  # Optional notes (e.g., "consensus", "weak +3")

    def to_dict(self) -> Dict:
        """Convert to dictionary for CSV export."""
        result = {
            "sequence": self.sequence,
            "total_hbonds": self.total_hbonds,
        }
        # Add per-position H-bonds with descriptive column names
        for i, (label, hbonds) in enumerate(zip(POSITION_LABELS, self.per_position)):
            col_name = f"hbonds_pos_{label.replace('+', 'p').replace('-', 'm')}"
            result[col_name] = hbonds
        result["longest_stretch"] = self.longest_stretch
        result["stretch_start_pos"] = POSITION_LABELS[self.stretch_start] if self.stretch_start >= 0 else ""
        result["num_base_pairs"] = self.num_base_pairs
        result["notes"] = self.notes
        return result


@dataclass
class SequenceVariant:
    """A generated sequence variant with its H-bond analysis."""
    sequence: str
    hbond_result: HBondResult
    generation_method: str = "enumeration"


@dataclass
class TieredCandidate:
    """A candidate sequence for the tiered panel."""
    sequence: str
    tier: int  # Tier number (1-7)
    tier_range: Tuple[int, int]  # (min_hbonds, max_hbonds) for this tier
    total_hbonds: int
    longest_stretch: int
    has_frame1_stop: bool
    codon2_dna: str  # Codon at +3,+4,+5
    codon3_dna: str  # Codon at +6,+7,+8
    codon2_aa: str   # Amino acid for codon 2
    codon3_aa: str   # Amino acid for codon 3
    notes: str = ""
    expected_splicing: str = ""

    def to_dict(self) -> Dict:
        """Convert to dictionary for CSV export."""
        return {
            "tier": self.tier,
            "tier_hbonds": f"{self.tier_range[0]}-{self.tier_range[1]}",
            "sequence": self.sequence,
            "total_hbonds": self.total_hbonds,
            "longest_stretch": self.longest_stretch,
            "frame1_stop": "Yes" if self.has_frame1_stop else "No",
            "codon2_pos345": self.codon2_dna,
            "codon3_pos678": self.codon3_dna,
            "codon2_aa": self.codon2_aa,
            "codon3_aa": self.codon3_aa,
            "expected_splicing": self.expected_splicing,
            "notes": self.notes,
        }


# ============================================================================
# Core Functions
# ============================================================================

def get_pair_hbonds(ss_rna_nt: str, u1_nt: str) -> int:
    """
    Get the number of H-bonds for a base pair.

    Args:
        ss_rna_nt: 5'SS RNA nucleotide (A, U, G, or C)
        u1_nt: U1 snRNA nucleotide (A, U, G, C, or Ψ)

    Returns:
        Number of H-bonds (0, 2, or 3)
    """
    return HBOND_RULES.get((ss_rna_nt, u1_nt), 0)


def find_longest_continuous_stretch(hbonds: List[int]) -> Tuple[int, int]:
    """
    Find the longest continuous stretch of base pairs (positions with hbonds > 0).

    Args:
        hbonds: List of H-bond counts per position

    Returns:
        Tuple of (length of longest stretch, start index)
    """
    max_length = 0
    max_start = -1
    current_length = 0
    current_start = -1

    for i, h in enumerate(hbonds):
        if h > 0:
            if current_length == 0:
                current_start = i
            current_length += 1
        else:
            if current_length > max_length:
                max_length = current_length
                max_start = current_start
            current_length = 0

    # Check final stretch
    if current_length > max_length:
        max_length = current_length
        max_start = current_start

    return max_length, max_start


def calculate_hbonds(sequence: str) -> HBondResult:
    """
    Calculate H-bonds for an 11-mer 5'SS sequence.

    Args:
        sequence: 11-character DNA string (positions -3 to +8)
                  Input is DNA (ACGT), internally converted to RNA for pairing

    Returns:
        HBondResult with total, per-position, and stretch analysis

    Raises:
        ValueError: If sequence is not 11 nucleotides or contains invalid characters
    """
    sequence = sequence.upper()

    # Validate length
    if len(sequence) != 11:
        raise ValueError(f"Sequence must be 11 nucleotides, got {len(sequence)}: {sequence}")

    # Validate characters and convert to RNA
    rna_seq = []
    for i, nt in enumerate(sequence):
        if nt not in DNA_TO_RNA:
            raise ValueError(f"Invalid nucleotide '{nt}' at position {i} in sequence: {sequence}")
        rna_seq.append(DNA_TO_RNA[nt])

    # Calculate H-bonds at each position
    per_position = []
    for i, (ss_nt, u1_nt) in enumerate(zip(rna_seq, U1_PARTNERS)):
        hbonds = get_pair_hbonds(ss_nt, u1_nt)
        per_position.append(hbonds)

    # Find longest continuous stretch of base pairs
    longest, start = find_longest_continuous_stretch(per_position)

    # Count positions with base pairing
    num_base_pairs = sum(1 for h in per_position if h > 0)

    return HBondResult(
        sequence=sequence,
        total_hbonds=sum(per_position),
        per_position=per_position,
        longest_stretch=longest,
        stretch_start=start,
        num_base_pairs=num_base_pairs,
        notes=""
    )


def is_valid_splice_site(sequence: str) -> bool:
    """
    Check if a sequence has the invariant GT at +1/+2.

    Args:
        sequence: 11-mer DNA sequence

    Returns:
        True if positions +1/+2 are GT
    """
    if len(sequence) != 11:
        return False
    return sequence[3].upper() == "G" and sequence[4].upper() == "T"


# ============================================================================
# Frame 1 Stop Codon and Translation Functions
# ============================================================================

def has_frame1_stop_codon(sequence: str) -> bool:
    """
    Check if 11-mer has stop codon in Frame 1 (codon spans -1/+1/+2).

    Frame 1 layout:
        Positions: -3 -2 -1 | +1 +2 +3 | +4 +5 +6 | +7 +8
        Indices:    0  1  2 |  3  4  5 |  6  7  8 |  9 10
                            └─codon1─┘ └─codon2─┘ └─codon3 (partial)

    Actually for Frame 1 where codon spans -1,+1,+2:
        Codon 1: -1, +1, +2 = indices 2, 3, 4 (contains exon -1 + GT)
        Codon 2: +3, +4, +5 = indices 5, 6, 7
        Codon 3: +6, +7, +8 = indices 8, 9, 10

    Codon 1 ends in GT, so possible codons are _GT where _ is from exon.
    Stop codons (TAA, TAG, TGA) don't end in GT, so codon 1 never has a stop.

    We only need to check codons 2 and 3 (positions +3-5 and +6-8).

    Args:
        sequence: 11-character DNA string (positions -3 to +8)

    Returns:
        True if positions +3+4+5 or +6+7+8 form a stop codon
    """
    sequence = sequence.upper()
    if len(sequence) != 11:
        raise ValueError(f"Sequence must be 11 nucleotides, got {len(sequence)}")

    # Codon 2: +3, +4, +5 = indices 5, 6, 7
    codon2 = sequence[5:8]
    # Codon 3: +6, +7, +8 = indices 8, 9, 10
    codon3 = sequence[8:11]

    return codon2 in STOP_CODONS or codon3 in STOP_CODONS


def get_frame1_codons(sequence: str) -> Tuple[str, str]:
    """
    Extract Frame 1 codons from intronic portion of 11-mer.

    Args:
        sequence: 11-character DNA string

    Returns:
        Tuple of (codon2, codon3) where:
            codon2 = positions +3, +4, +5 (indices 5, 6, 7)
            codon3 = positions +6, +7, +8 (indices 8, 9, 10)
    """
    sequence = sequence.upper()
    codon2 = sequence[5:8]
    codon3 = sequence[8:11]
    return codon2, codon3


def translate_codon(codon: str) -> str:
    """
    Translate a DNA codon to amino acid.

    Args:
        codon: 3-letter DNA codon

    Returns:
        Single-letter amino acid code, or '?' if invalid
    """
    codon = codon.upper()
    return CODON_TABLE.get(codon, "?")


def translate_frame1_codons(sequence: str) -> Tuple[str, str]:
    """
    Translate Frame 1 codons from intronic portion.

    Args:
        sequence: 11-character DNA string

    Returns:
        Tuple of (aa2, aa3) for codons 2 and 3
        '*' indicates stop codon
    """
    codon2, codon3 = get_frame1_codons(sequence)
    return translate_codon(codon2), translate_codon(codon3)


# ============================================================================
# Sequence Generation
# ============================================================================

def generate_all_valid_sequences() -> Dict[int, List[str]]:
    """
    Generate all valid 11-mer 5'SS sequences with GT at +1/+2.

    Variable positions: -3, -2, -1, +3, +4, +5, +6, +7, +8 (9 positions)
    Fixed positions: +1=G, +2=T
    Total: 4^9 = 262,144 sequences

    Returns:
        Dictionary mapping H-bond count to list of sequences
    """
    sequences_by_hbonds: Dict[int, List[str]] = {}

    # Variable position indices (0-based): 0, 1, 2, 5, 6, 7, 8, 9, 10
    variable_indices = [0, 1, 2, 5, 6, 7, 8, 9, 10]
    nucleotides = ["A", "C", "G", "T"]

    # Enumerate all combinations
    for combo in itertools.product(nucleotides, repeat=9):
        # Build sequence with fixed GT at positions 3, 4
        seq_list = [""] * 11
        seq_list[3] = "G"
        seq_list[4] = "T"

        for idx, var_idx in enumerate(variable_indices):
            seq_list[var_idx] = combo[idx]

        sequence = "".join(seq_list)

        # Calculate H-bonds
        result = calculate_hbonds(sequence)

        # Group by H-bond count
        if result.total_hbonds not in sequences_by_hbonds:
            sequences_by_hbonds[result.total_hbonds] = []
        sequences_by_hbonds[result.total_hbonds].append(sequence)

    return sequences_by_hbonds


def generate_sequences_for_hbond_count(
    target: int,
    all_sequences: Dict[int, List[str]],
    num_variants: int = 10,
    prioritize_stretch: bool = True
) -> List[SequenceVariant]:
    """
    Get sequence variants for a target H-bond count.

    Args:
        target: Target H-bond count
        all_sequences: Pre-generated sequences grouped by H-bond count
        num_variants: Maximum number of variants to return
        prioritize_stretch: If True, prefer sequences with longer continuous stretches

    Returns:
        List of SequenceVariant objects
    """
    if target not in all_sequences:
        return []

    sequences = all_sequences[target]

    # Calculate H-bonds for each and sort by stretch length
    variants = []
    for seq in sequences:
        result = calculate_hbonds(seq)
        variants.append(SequenceVariant(
            sequence=seq,
            hbond_result=result,
            generation_method="enumeration"
        ))

    if prioritize_stretch:
        # Sort by longest stretch (descending), then by stretch start position
        variants.sort(key=lambda v: (-v.hbond_result.longest_stretch, v.hbond_result.stretch_start))

    return variants[:num_variants]


def select_diverse_variants(
    all_sequences: Dict[int, List[str]],
    min_hbonds: int = MIN_HBONDS_GENERATE,
    max_hbonds: int = MAX_HBONDS,
    variants_per_tier: int = 100
) -> List[SequenceVariant]:
    """
    Select diverse variants across all H-bond tiers.

    Args:
        all_sequences: Pre-generated sequences grouped by H-bond count
        min_hbonds: Minimum H-bond count to include
        max_hbonds: Maximum H-bond count to include
        variants_per_tier: Number of variants to select per H-bond tier

    Returns:
        List of SequenceVariant objects
    """
    all_variants = []

    for hbond_count in range(min_hbonds, max_hbonds + 1):
        variants = generate_sequences_for_hbond_count(
            hbond_count, all_sequences, num_variants=variants_per_tier
        )
        all_variants.extend(variants)

    return all_variants


# ============================================================================
# Tiered Panel Generation
# ============================================================================

def calculate_sequence_diversity(seq1: str, seq2: str) -> int:
    """
    Calculate Hamming distance between two sequences.

    Args:
        seq1: First sequence
        seq2: Second sequence

    Returns:
        Number of positions that differ
    """
    return sum(1 for a, b in zip(seq1.upper(), seq2.upper()) if a != b)


def select_diverse_within_tier(
    sequences: List[str],
    num_select: int,
    stretch_mode: str = "diverse"
) -> List[str]:
    """
    Select diverse sequences within a tier.

    Strategy depends on stretch_mode:
    - "diverse": Select sequences with VARYING stretch lengths (short, medium, long)
                 to create gradual progression across the panel
    - "short": Prefer shortest stretches (weak U1 binding)
    - "long": Prefer longest stretches (consensus-like, strong donors)

    Args:
        sequences: List of candidate sequences
        num_select: Number to select
        stretch_mode: "diverse", "short", or "long"

    Returns:
        List of selected sequences
    """
    if len(sequences) <= num_select:
        return sequences

    # Calculate H-bond results for all sequences
    seq_results = [(seq, calculate_hbonds(seq)) for seq in sequences]

    if stretch_mode == "diverse":
        # Group sequences by stretch length
        by_stretch: Dict[int, List[Tuple[str, HBondResult]]] = {}
        for seq, result in seq_results:
            stretch = result.longest_stretch
            if stretch not in by_stretch:
                by_stretch[stretch] = []
            by_stretch[stretch].append((seq, result))

        # Sort each group by H-bonds (descending) for quality
        for stretch in by_stretch:
            by_stretch[stretch].sort(key=lambda x: -x[1].total_hbonds)

        stretch_values = sorted(by_stretch.keys())
        min_stretch = min(stretch_values)
        max_stretch = max(stretch_values)

        # Calculate which stretch values to target based on num_select
        # We want to spread selections EVENLY across the available stretch range
        if num_select == 1:
            target_stretches = [stretch_values[len(stretch_values) // 2]]  # Middle
        elif num_select >= len(stretch_values):
            target_stretches = stretch_values  # One from each
        else:
            # Select num_select stretches evenly distributed across the range
            target_stretches = []
            for i in range(num_select):
                # Interpolate position in stretch_values list
                pos = int(i * (len(stretch_values) - 1) / (num_select - 1)) if num_select > 1 else 0
                target_stretches.append(stretch_values[pos])

        # Select sequences from target stretches
        selected = []
        indices = {s: 0 for s in stretch_values}

        # First pass: try to get one from each target stretch
        for target in target_stretches:
            if len(selected) >= num_select:
                break
            if target in by_stretch:
                idx = indices[target]
                while idx < len(by_stretch[target]):
                    candidate_seq = by_stretch[target][idx][0]
                    # Check diversity against already selected
                    if not selected or min(calculate_sequence_diversity(candidate_seq, s) for s in selected) >= 2:
                        selected.append(candidate_seq)
                        indices[target] = idx + 1
                        break
                    idx += 1
                indices[target] = idx + 1

        # If we still need more, do round-robin through all stretches
        while len(selected) < num_select:
            added_this_round = False
            for stretch in stretch_values:
                if len(selected) >= num_select:
                    break
                idx = indices[stretch]
                if idx < len(by_stretch[stretch]):
                    candidate_seq = by_stretch[stretch][idx][0]
                    if candidate_seq not in selected:
                        # Check diversity
                        if not selected or min(calculate_sequence_diversity(candidate_seq, s) for s in selected) >= 2:
                            selected.append(candidate_seq)
                            added_this_round = True
                    indices[stretch] += 1

            # If no diverse sequences found, accept any
            if not added_this_round:
                for stretch in stretch_values:
                    if len(selected) >= num_select:
                        break
                    idx = indices[stretch]
                    while idx < len(by_stretch[stretch]):
                        candidate_seq = by_stretch[stretch][idx][0]
                        if candidate_seq not in selected:
                            selected.append(candidate_seq)
                            indices[stretch] = idx + 1
                            break
                        idx += 1
                    else:
                        indices[stretch] = idx + 1
                break  # Only do one rescue pass to avoid infinite loop

            # Safety check
            if all(indices[s] >= len(by_stretch[s]) for s in stretch_values):
                break

        return selected[:num_select]

    else:
        # Legacy behavior: sort by stretch length
        prioritize_short_stretch = (stretch_mode == "short")
        if prioritize_short_stretch:
            seq_results.sort(key=lambda x: (x[1].longest_stretch, -x[1].total_hbonds))
        else:
            seq_results.sort(key=lambda x: (-x[1].longest_stretch, -x[1].total_hbonds))

        # Greedy selection for diversity
        selected = [seq_results[0][0]]  # Start with best by stretch criteria

        for seq, _ in seq_results[1:]:
            if len(selected) >= num_select:
                break

            # Calculate minimum distance to already selected sequences
            min_dist = min(calculate_sequence_diversity(seq, s) for s in selected)

            # Only add if sufficiently different (at least 2 positions different)
            if min_dist >= 2:
                selected.append(seq)

        # If we don't have enough diverse sequences, add remaining by stretch order
        if len(selected) < num_select:
            for seq, _ in seq_results:
                if seq not in selected:
                    selected.append(seq)
                    if len(selected) >= num_select:
                        break

        return selected[:num_select]


def has_consensus_motif(sequence: str) -> bool:
    """
    Check if a sequence has the AG|GT consensus motif at the splice junction.

    The canonical splice donor consensus is AG|GT where:
    - Position -1 is G (immediately before the splice site)
    - Positions +1+2 are GT (the invariant dinucleotide)

    Together, this forms AG|GT which is required for most splice prediction
    tools (SpliceBDGP, MaxEntScan) to recognize the donor.

    Args:
        sequence: 11-mer DNA sequence (positions -3 to +8)

    Returns:
        True if position -1 is G and +1+2 is GT
    """
    if len(sequence) != 11:
        return False
    # Position -1 is index 2, positions +1+2 are indices 3,4
    return (sequence[2].upper() == "G" and
            sequence[3].upper() == "G" and
            sequence[4].upper() == "T")


def generate_tiered_panel(
    all_sequences: Dict[int, List[str]],
    tiers: Optional[List[Tuple[int, int, int]]] = None,
    avoid_frame1_stops: bool = True,
    require_consensus_motif: bool = True,
    stretch_mode: str = "diverse",
) -> List[TieredCandidate]:
    """
    Generate tiered panel of candidates for empirical testing.

    v5 CHANGE: Now uses DIVERSE stretch selection by default to create gradual
    progression across the panel. Each tier contains sequences with varying
    stretch lengths instead of all-short or all-long.

    Previous issue (v4): Tiers 1-4 all had stretch=3, tiers 5-7 jumped to stretch=8+.
    This created too dramatic a transition between tiers 4 and 5.

    v5 solution: Select DIVERSE stretches within each tier for gradual progression.

    Args:
        all_sequences: Pre-generated sequences grouped by H-bond count
        tiers: List of (min_hbonds, max_hbonds, num_sequences) tuples
               Default: DEFAULT_TIERS
        avoid_frame1_stops: If True, filter out sequences with Frame 1 stop codons
        require_consensus_motif: If True, require AG|GT consensus (G at position -1)
                                 This is CRITICAL for splice prediction tools
        stretch_mode: Selection strategy for stretch lengths:
                      "diverse" = select across stretch range (default, recommended for v5)
                      "short" = force short stretch for all tiers (weak donors)
                      "long" = force long stretch for all tiers (strong donors)

    Returns:
        List of TieredCandidate objects
    """
    if tiers is None:
        tiers = DEFAULT_TIERS

    candidates = []

    for tier_num, (min_hb, max_hb, num_seqs) in enumerate(tiers, start=1):
        # Collect all sequences in this H-bond range
        tier_sequences = []
        for hb in range(min_hb, max_hb + 1):
            if hb in all_sequences:
                tier_sequences.extend(all_sequences[hb])

        if not tier_sequences:
            print(f"Warning: No sequences found for tier {tier_num} ({min_hb}-{max_hb} H-bonds)")
            continue

        # Filter 1: Require AG|GT consensus motif if requested
        # CRITICAL: Without this, splice prediction tools won't recognize the donor!
        if require_consensus_motif:
            before_count = len(tier_sequences)
            tier_sequences = [seq for seq in tier_sequences if has_consensus_motif(seq)]
            after_count = len(tier_sequences)
            if after_count == 0 and before_count > 0:
                print(f"Warning: Tier {tier_num} has no sequences with AG|GT consensus "
                      f"(filtered {before_count} -> 0)")
                continue

        # Filter 2: Remove Frame 1 stop codons if requested
        if avoid_frame1_stops:
            tier_sequences = [seq for seq in tier_sequences if not has_frame1_stop_codon(seq)]

        if not tier_sequences:
            print(f"Warning: No sequences without Frame 1 stops for tier {tier_num}")
            continue

        # Filter 3 & 4: Select diverse sequences with specified stretch mode
        selected = select_diverse_within_tier(
            tier_sequences, num_seqs, stretch_mode
        )

        # Get annotation for this tier
        tier_key = (min_hb, max_hb)
        annotation, expected_splicing = TIER_ANNOTATIONS.get(tier_key, ("", ""))

        # Create TieredCandidate objects
        for seq in selected:
            result = calculate_hbonds(seq)
            codon2, codon3 = get_frame1_codons(seq)
            aa2, aa3 = translate_frame1_codons(seq)

            candidate = TieredCandidate(
                sequence=seq,
                tier=tier_num,
                tier_range=(min_hb, max_hb),
                total_hbonds=result.total_hbonds,
                longest_stretch=result.longest_stretch,
                has_frame1_stop=has_frame1_stop_codon(seq),
                codon2_dna=codon2,
                codon3_dna=codon3,
                codon2_aa=aa2,
                codon3_aa=aa3,
                notes=annotation,
                expected_splicing=expected_splicing,
            )
            candidates.append(candidate)

    return candidates


def export_tiered_panel(candidates: List[TieredCandidate], filepath: str) -> None:
    """
    Export tiered panel candidates to CSV.

    Args:
        candidates: List of TieredCandidate objects
        filepath: Output CSV file path
    """
    if not candidates:
        print("No candidates to export")
        return

    Path(filepath).parent.mkdir(parents=True, exist_ok=True)

    fieldnames = list(candidates[0].to_dict().keys())

    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for candidate in candidates:
            writer.writerow(candidate.to_dict())

    print(f"Exported {len(candidates)} candidates to: {filepath}")


def generate_panel_summary(candidates: List[TieredCandidate]) -> str:
    """
    Generate summary of the tiered panel.

    Args:
        candidates: List of TieredCandidate objects

    Returns:
        Summary string
    """
    lines = [
        "Tiered Splice Donor Panel Summary",
        "=" * 50,
        ""
    ]

    # Group by tier
    by_tier: Dict[int, List[TieredCandidate]] = {}
    for c in candidates:
        if c.tier not in by_tier:
            by_tier[c.tier] = []
        by_tier[c.tier].append(c)

    for tier_num in sorted(by_tier.keys()):
        tier_candidates = by_tier[tier_num]
        if tier_candidates:
            c = tier_candidates[0]
            hb_range = f"{c.tier_range[0]}-{c.tier_range[1]}"
            lines.append(
                f"Tier {tier_num} ({hb_range} H-bonds): "
                f"{len(tier_candidates)} sequences - expect {c.expected_splicing} splicing"
            )

            # Show sequences
            for candidate in tier_candidates:
                lines.append(
                    f"  {candidate.sequence} | {candidate.total_hbonds} H-bonds | "
                    f"stretch={candidate.longest_stretch} | "
                    f"AA: {candidate.codon2_aa}{candidate.codon3_aa}"
                )

    lines.append("")
    lines.append(f"Total: {len(candidates)} sequences for synthesis")

    # Amino acid statistics
    lines.append("")
    lines.append("Amino acid distribution (codons 2+3):")
    lines.append("-" * 30)

    aa_counts: Dict[str, int] = {}
    for c in candidates:
        for aa in [c.codon2_aa, c.codon3_aa]:
            aa_counts[aa] = aa_counts.get(aa, 0) + 1

    for aa in sorted(aa_counts.keys()):
        lines.append(f"  {aa}: {aa_counts[aa]}")

    return "\n".join(lines)


# ============================================================================
# Validation
# ============================================================================

# Known sequences for validation
# Format: (sequence or partial, expected_hbonds_range, source, notes)
KNOWN_SEQUENCES = [
    # Consensus sequences (calculated: 21 and 18 H-bonds respectively)
    ("CAGGTAAGTAT", (20, 22), "Consensus extended", "Strong splice site"),
    ("AAGGTAAGTAT", (17, 19), "Near-consensus", "A at -3 instead of C (-3 H-bonds)"),

    # Weak sites (both calculate to 18 H-bonds)
    ("AAGGTGAGTAT", (17, 19), "Weak +3", "G at +3 (G-Ψ wobble still 2 H-bonds)"),
    ("AAGGTAAGCAT", (17, 19), "Weak +6", "C at +6 (mismatch with C)"),

    # Minimal-ish sequence (calculated: 13 H-bonds - T wobbles with G at -3 and U at +7)
    ("TTTGTTTTTAT", (12, 14), "Minimal", "T-G wobble at -3, T-U at +7/+8"),

    # Near-maximum (calculated: 22 H-bonds)
    ("CAGGTAAGTGC", (21, 23), "Near-maximum", "Suboptimal at +5/+6"),
]


def validate_known_sequences() -> List[Tuple[str, int, str, bool]]:
    """
    Validate H-bond calculations against known sequences.

    Returns:
        List of (sequence, calculated_hbonds, source, passed) tuples
    """
    results = []

    for seq, (min_exp, max_exp), source, notes in KNOWN_SEQUENCES:
        if len(seq) != 11:
            results.append((seq, -1, source, False))
            continue

        try:
            result = calculate_hbonds(seq)
            passed = min_exp <= result.total_hbonds <= max_exp
            results.append((seq, result.total_hbonds, source, passed))
        except ValueError as e:
            results.append((seq, -1, f"{source}: {e}", False))

    return results


def run_validation(output_file: Optional[str] = None) -> bool:
    """
    Run validation and optionally save results.

    Args:
        output_file: Path to save validation results (optional)

    Returns:
        True if all validations pass
    """
    results = validate_known_sequences()

    lines = [
        "Splice Donor H-Bond Validation Results",
        "=" * 50,
        ""
    ]

    all_passed = True
    for seq, hbonds, source, passed in results:
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_passed = False
        lines.append(f"{status}: {seq} -> {hbonds} H-bonds ({source})")

    lines.append("")
    lines.append(f"Overall: {'PASSED' if all_passed else 'FAILED'}")

    output_text = "\n".join(lines)
    print(output_text)

    if output_file:
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, "w") as f:
            f.write(output_text)
        print(f"\nValidation results saved to: {output_file}")

    return all_passed


# ============================================================================
# Output Functions
# ============================================================================

def export_to_csv(variants: List[SequenceVariant], filepath: str) -> None:
    """
    Export sequence variants to CSV.

    Args:
        variants: List of SequenceVariant objects
        filepath: Output CSV file path
    """
    if not variants:
        print("No variants to export")
        return

    Path(filepath).parent.mkdir(parents=True, exist_ok=True)

    # Get column names from first result
    fieldnames = list(variants[0].hbond_result.to_dict().keys())
    fieldnames.append("generation_method")

    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for variant in variants:
            row = variant.hbond_result.to_dict()
            row["generation_method"] = variant.generation_method
            writer.writerow(row)

    print(f"Exported {len(variants)} sequences to: {filepath}")


def export_full_enumeration(all_sequences: Dict[int, List[str]], filepath: str) -> None:
    """
    Export all enumerated sequences to CSV.

    Args:
        all_sequences: Dictionary mapping H-bond count to sequences
        filepath: Output CSV file path
    """
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)

    # Collect all variants
    all_variants = []
    for hbond_count in sorted(all_sequences.keys()):
        for seq in all_sequences[hbond_count]:
            result = calculate_hbonds(seq)
            all_variants.append(SequenceVariant(
                sequence=seq,
                hbond_result=result,
                generation_method="enumeration"
            ))

    export_to_csv(all_variants, filepath)


def plot_distribution(all_sequences: Dict[int, List[str]], filepath: str) -> None:
    """
    Plot histogram of sequences by H-bond count.

    Args:
        all_sequences: Dictionary mapping H-bond count to sequences
        filepath: Output image file path
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping distribution plot")
        return

    Path(filepath).parent.mkdir(parents=True, exist_ok=True)

    # Prepare data
    hbond_counts = sorted(all_sequences.keys())
    seq_counts = [len(all_sequences[h]) for h in hbond_counts]

    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))

    bars = ax.bar(hbond_counts, seq_counts, color='steelblue', edgecolor='black', alpha=0.7)

    ax.set_xlabel("Total H-bonds", fontsize=12)
    ax.set_ylabel("Number of Sequences", fontsize=12)
    ax.set_title("Distribution of 5'SS Sequences by H-bond Count\n(with invariant GT at +1/+2)", fontsize=14)

    # Add count labels on bars
    for bar, count in zip(bars, seq_counts):
        if count > 0:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 100,
                   f'{count:,}', ha='center', va='bottom', fontsize=8)

    ax.set_xticks(hbond_counts)
    ax.set_xlim(min(hbond_counts) - 0.5, max(hbond_counts) + 0.5)
    ax.grid(axis='y', alpha=0.3)

    # Add summary statistics
    total_seqs = sum(seq_counts)
    summary_text = f"Total sequences: {total_seqs:,}\nH-bond range: {min(hbond_counts)}-{max(hbond_counts)}"
    ax.text(0.02, 0.98, summary_text, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(filepath, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"Distribution plot saved to: {filepath}")


def generate_sequence_logos(
    all_sequences: Dict[int, List[str]],
    output_dir: str,
    tiers: Optional[List[int]] = None
) -> None:
    """
    Generate sequence logos for selected H-bond tiers.

    Args:
        all_sequences: Dictionary mapping H-bond count to sequences
        output_dir: Directory to save logos
        tiers: Specific H-bond counts to generate logos for (default: all)
    """
    try:
        import logomaker
        import pandas as pd
        import matplotlib.pyplot as plt
    except ImportError:
        print("logomaker not available, skipping sequence logos")
        return

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    if tiers is None:
        tiers = sorted(all_sequences.keys())

    for hbond_count in tiers:
        if hbond_count not in all_sequences or not all_sequences[hbond_count]:
            continue

        sequences = all_sequences[hbond_count]

        # Build counts matrix
        counts = {pos: {"A": 0, "C": 0, "G": 0, "T": 0} for pos in range(11)}
        for seq in sequences:
            for i, nt in enumerate(seq):
                counts[i][nt] += 1

        # Convert to DataFrame
        df = pd.DataFrame(counts).T
        df.columns = ["A", "C", "G", "T"]
        df.index = POSITION_LABELS

        # Normalize to frequencies
        df = df.div(df.sum(axis=1), axis=0)

        # Create logo
        fig, ax = plt.subplots(figsize=(14, 3))
        logo = logomaker.Logo(df, ax=ax, color_scheme='classic')
        ax.set_ylabel("Frequency")
        ax.set_title(f"5'SS Sequence Logo - {hbond_count} H-bonds (n={len(sequences):,})")

        filepath = os.path.join(output_dir, f"logo_{hbond_count}hbonds.png")
        plt.tight_layout()
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"Logo saved: {filepath}")


def generate_summary_stats(all_sequences: Dict[int, List[str]]) -> str:
    """
    Generate summary statistics for the generated sequences.

    Args:
        all_sequences: Dictionary mapping H-bond count to sequences

    Returns:
        Summary string
    """
    lines = [
        "Splice Donor H-Bond Sequence Generation Summary",
        "=" * 50,
        ""
    ]

    total_seqs = sum(len(seqs) for seqs in all_sequences.values())
    lines.append(f"Total sequences generated: {total_seqs:,}")
    lines.append(f"H-bond range: {min(all_sequences.keys())}-{max(all_sequences.keys())}")
    lines.append("")
    lines.append("Distribution by H-bond count:")
    lines.append("-" * 30)

    for hbond_count in sorted(all_sequences.keys()):
        count = len(all_sequences[hbond_count])
        pct = 100 * count / total_seqs
        lines.append(f"  {hbond_count:2d} H-bonds: {count:6,} sequences ({pct:5.2f}%)")

    # Calculate stretch statistics for each tier
    lines.append("")
    lines.append("Average longest stretch by H-bond count:")
    lines.append("-" * 30)

    for hbond_count in sorted(all_sequences.keys()):
        sequences = all_sequences[hbond_count]
        if sequences:
            stretches = [calculate_hbonds(seq).longest_stretch for seq in sequences]
            avg_stretch = sum(stretches) / len(stretches)
            max_stretch = max(stretches)
            lines.append(f"  {hbond_count:2d} H-bonds: avg={avg_stretch:.1f}, max={max_stretch}")

    return "\n".join(lines)


# ============================================================================
# CLI Interface
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Splice Donor H-Bond Prediction Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Score a single sequence
  python splice_donor_hbonds.py --score CAGGTAAGTAT

  # Validate against known sequences
  python splice_donor_hbonds.py --validate

  # Generate all sequences and export
  python splice_donor_hbonds.py --generate --output-dir ./output

  # Generate with visualization
  python splice_donor_hbonds.py --generate --output-dir ./output --plot --logos

  # Generate tiered panel for splicing efficiency screening (Frame 1)
  # v5: Uses DIVERSE stretch selection for gradual progression across panel
  python splice_donor_hbonds.py --generate-panel --frame 1 --output-dir ./output

  # Override stretch mode (not recommended - diverse is best for gradual progression)
  python splice_donor_hbonds.py --generate-panel --stretch-mode short --output-dir ./output
  python splice_donor_hbonds.py --generate-panel --stretch-mode long --output-dir ./output

  # Generate panel without AG|GT consensus filter (not recommended)
  python splice_donor_hbonds.py --generate-panel --no-consensus-filter --output-dir ./output

  # Check Frame 1 stop codons in a sequence
  python splice_donor_hbonds.py --check-stops CAGGTAAGTAT
        """
    )

    parser.add_argument("--score", type=str, metavar="SEQ",
                       help="Score a single 11-mer sequence")
    parser.add_argument("--validate", action="store_true",
                       help="Validate against known sequences")
    parser.add_argument("--generate", action="store_true",
                       help="Generate all valid sequences")
    parser.add_argument("--generate-panel", action="store_true",
                       help="Generate tiered panel for splicing efficiency screening")
    parser.add_argument("--frame", type=int, choices=[1], default=1,
                       help="Reading frame for stop codon filtering (default: 1)")
    parser.add_argument("--check-stops", type=str, metavar="SEQ",
                       help="Check Frame 1 stop codons and translate codons for a sequence")
    parser.add_argument("--include-stops", action="store_true",
                       help="Include sequences with stop codons in panel (for analysis)")
    parser.add_argument("--no-consensus-filter", action="store_true",
                       help="Don't require AG|GT consensus motif (not recommended - tools may not detect)")
    parser.add_argument("--stretch-mode", type=str, choices=["diverse", "short", "long"],
                       default="diverse",
                       help="Stretch selection mode: diverse (default, recommended), short, or long")
    parser.add_argument("--output-dir", type=str, default="./output",
                       help="Output directory (default: ./output)")
    parser.add_argument("--plot", action="store_true",
                       help="Generate distribution plot")
    parser.add_argument("--logos", action="store_true",
                       help="Generate sequence logos (requires logomaker)")
    parser.add_argument("--logo-tiers", type=str, metavar="LIST",
                       help="Comma-separated H-bond counts for logos (default: all)")
    parser.add_argument("--summary", action="store_true",
                       help="Print summary statistics")
    parser.add_argument("--export-full", action="store_true",
                       help="Export all sequences (not just selected variants)")
    parser.add_argument("--variants-per-tier", type=int, default=100,
                       help="Number of variants per H-bond tier for selected export")

    args = parser.parse_args()

    # Score single sequence
    if args.score:
        try:
            result = calculate_hbonds(args.score)
            print(f"Sequence: {result.sequence}")
            print(f"Total H-bonds: {result.total_hbonds}")
            print(f"Per position:")
            for label, hbonds in zip(POSITION_LABELS, result.per_position):
                print(f"  {label}: {hbonds}")
            print(f"Longest continuous stretch: {result.longest_stretch} (starting at {POSITION_LABELS[result.stretch_start]})")
            print(f"Number of base-paired positions: {result.num_base_pairs}/11")

            # Check if valid splice site
            if is_valid_splice_site(result.sequence):
                print("Valid splice site (GT at +1/+2): Yes")
            else:
                print("Valid splice site (GT at +1/+2): No")
        except ValueError as e:
            print(f"Error: {e}")
            return 1
        return 0

    # Check Frame 1 stop codons
    if args.check_stops:
        try:
            seq = args.check_stops.upper()
            if len(seq) != 11:
                print(f"Error: Sequence must be 11 nucleotides, got {len(seq)}")
                return 1

            result = calculate_hbonds(seq)
            has_stop = has_frame1_stop_codon(seq)
            codon2, codon3 = get_frame1_codons(seq)
            aa2, aa3 = translate_frame1_codons(seq)

            print(f"Sequence: {seq}")
            print(f"Total H-bonds: {result.total_hbonds}")
            print(f"Longest stretch: {result.longest_stretch}")
            print()
            print("Frame 1 Analysis:")
            print(f"  Codon 2 (+3,+4,+5): {codon2} -> {aa2}")
            print(f"  Codon 3 (+6,+7,+8): {codon3} -> {aa3}")
            print(f"  Has Frame 1 stop codon: {'YES' if has_stop else 'No'}")

            if has_stop:
                if codon2 in STOP_CODONS:
                    print(f"    -> Stop at codon 2 ({codon2})")
                if codon3 in STOP_CODONS:
                    print(f"    -> Stop at codon 3 ({codon3})")
        except ValueError as e:
            print(f"Error: {e}")
            return 1
        return 0

    # Validation
    if args.validate:
        output_file = os.path.join(args.output_dir, "validation", "known_sequences_check.txt")
        success = run_validation(output_file)
        return 0 if success else 1

    # Generate tiered panel
    if args.generate_panel:
        print("Generating tiered splice donor panel (v5 - diverse stretch selection)...")
        print(f"Frame: {args.frame}")
        print(f"Filter Frame 1 stop codons: {not args.include_stops}")
        print(f"Require AG|GT consensus: {not args.no_consensus_filter}")
        print(f"Stretch selection mode: {args.stretch_mode}")
        print()

        # Generate all sequences first
        print("Generating all valid 5'SS sequences (4^9 = 262,144)...")
        all_sequences = generate_all_valid_sequences()
        print(f"Generated sequences across {len(all_sequences)} H-bond tiers")
        print()

        # Generate tiered panel with diverse stretch selection (v5 default)
        candidates = generate_tiered_panel(
            all_sequences,
            tiers=DEFAULT_TIERS,
            avoid_frame1_stops=not args.include_stops,
            require_consensus_motif=not args.no_consensus_filter,
            stretch_mode=args.stretch_mode,
        )

        # Print summary
        print(generate_panel_summary(candidates))
        print()

        # Export to CSV
        os.makedirs(args.output_dir, exist_ok=True)
        csv_path = os.path.join(args.output_dir, "tiered_panel_candidates.csv")
        export_tiered_panel(candidates, csv_path)

        # Also save summary to text file
        summary_path = os.path.join(args.output_dir, "tiered_panel_summary.txt")
        with open(summary_path, "w") as f:
            f.write(generate_panel_summary(candidates))
        print(f"Summary saved to: {summary_path}")

        return 0

    # Generate sequences
    if args.generate:
        print("Generating all valid 5'SS sequences (4^9 = 262,144)...")
        all_sequences = generate_all_valid_sequences()
        print(f"Generated sequences across {len(all_sequences)} H-bond tiers")

        # Summary statistics
        if args.summary:
            print()
            print(generate_summary_stats(all_sequences))

        # Export
        os.makedirs(args.output_dir, exist_ok=True)

        if args.export_full:
            csv_path = os.path.join(args.output_dir, "all_5ss_sequences.csv")
            export_full_enumeration(all_sequences, csv_path)
        else:
            # Export selected variants
            variants = select_diverse_variants(
                all_sequences,
                variants_per_tier=args.variants_per_tier
            )
            csv_path = os.path.join(args.output_dir, "synthetic_5ss_sequences.csv")
            export_to_csv(variants, csv_path)

        # Distribution plot
        if args.plot:
            plot_path = os.path.join(args.output_dir, "hbond_distribution.png")
            plot_distribution(all_sequences, plot_path)

        # Sequence logos
        if args.logos:
            logo_dir = os.path.join(args.output_dir, "sequence_logos")
            tiers = None
            if args.logo_tiers:
                tiers = [int(x.strip()) for x in args.logo_tiers.split(",")]
            generate_sequence_logos(all_sequences, logo_dir, tiers)

        return 0

    # If no action specified, show help
    parser.print_help()
    return 0


if __name__ == "__main__":
    exit(main())
