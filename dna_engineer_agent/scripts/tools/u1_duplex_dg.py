#!/usr/bin/env python3
"""
U1 snRNA : 5' Splice Site Duplex Stability Calculator

Calculate ΔG for U1 snRNA binding to 5' splice sites using nearest-neighbor
thermodynamics with:
- Standard Turner (1998/2004) parameters for most base pairs
- Hudson et al. (2013) parameters for pseudouridine (Ψ)-A pairs at positions +3/+4

Background:
The 5'SS (positions -3 to +8) base-pairs with U1 snRNA in an antiparallel manner:

    5'SS:   -3 -2 -1 | +1 +2 +3 +4 +5 +6 +7 +8
    U1:     11 10  9 |  8  7  6  5  4  3  2  1   (U1 positions)
    U1nt:    G  U  C |  C  A  Ψ  Ψ  A  C  U  U   (U1 nucleotides)

Key insight: Positions +3 and +4 pair with pseudouridine (Ψ) in U1 snRNA.
Ψ-A pairs are ~1.7 kcal/mol MORE stable than U-A pairs (Hudson et al. 2013).

This explains why some 13 H-bond sequences are "Not Detected" by NNSplice
while others score 0.98: consecutive mismatches at +3/+4 create a destabilizing
bubble that H-bond counting alone doesn't capture.

References:
- Hudson GA et al. (2013) RNA 19(1):36-45. DOI: 10.1261/rna.035519.112
- Turner DH & Mathews DH (2010) Nucleic Acids Res 38(D):D230-D235
- Freund M et al. (2005) Nucleic Acids Res 33(17):5498-5510
"""

import argparse
import csv
import itertools
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ============================================================================
# Constants
# ============================================================================

# U1 snRNA sequence (positions 4-14, human consensus)
# 5'→3': m2,2,7G-ppp-A-U-A-C-Ψ-Ψ-A-C-C-U...
# The region that pairs with 5'SS (positions 4-14 of U1, numbered 1-11 here):
U1_SNRNA_BINDING = ["G", "U", "C", "C", "A", "Ψ", "Ψ", "A", "C", "U", "U"]

# Position labels for the 11-mer 5'SS
POSITION_LABELS = ["-3", "-2", "-1", "+1", "+2", "+3", "+4", "+5", "+6", "+7", "+8"]

# Watson-Crick base pair partners
WC_PAIRS = {
    'A': 'U', 'U': 'A',
    'G': 'C', 'C': 'G',
    'Ψ': 'A',  # Pseudouridine pairs with A
}

# Wobble pairs
WOBBLE_PAIRS = {
    ('G', 'U'), ('U', 'G'),
}

# ============================================================================
# Nearest-Neighbor Parameters
# ============================================================================

# Standard Turner (1998/2004) nearest-neighbor parameters for RNA duplexes
# Format: (5'XY/3'AB) means X pairs with B, Y pairs with A (antiparallel)
# These are ΔG°37 values in kcal/mol
# Note: For 5'SS/U1 duplex, we read the 5'SS 5'→3' and U1 3'←5'

TURNER_NN_DNA_RNA = {
    # Watson-Crick stacks (from Turner lab RNA/RNA parameters)
    # Format: tuple((5'_ss_1, 5'_ss_2), (3'_u1_1, 3'_u1_2)) for 5'SS paired to U1
    # The 5'SS is DNA (converted to RNA: T→U), U1 is RNA with Ψ

    # AA/UU stack
    ('A', 'A', 'U', 'U'): -0.93,  # AA/UU
    ('U', 'U', 'A', 'A'): -0.93,  # UU/AA

    # AU/UA stacks
    ('A', 'U', 'U', 'A'): -1.10,  # AU/UA
    ('U', 'A', 'A', 'U'): -1.33,  # UA/AU

    # GC/CG stacks (strongest)
    ('G', 'C', 'C', 'G'): -2.36,  # GC/CG
    ('C', 'G', 'G', 'C'): -3.26,  # CG/GC

    # GA/CU and AG/UC stacks
    ('G', 'A', 'C', 'U'): -2.35,  # GA/CU
    ('A', 'G', 'U', 'C'): -2.08,  # AG/UC
    ('C', 'U', 'G', 'A'): -2.08,  # CU/GA
    ('U', 'C', 'A', 'G'): -2.35,  # UC/AG

    # GU/CA and UG/AC (involves wobble)
    ('G', 'U', 'C', 'A'): -2.24,  # GU/CA
    ('U', 'G', 'A', 'C'): -2.11,  # UG/AC
    ('C', 'A', 'G', 'U'): -2.11,  # CA/GU
    ('A', 'C', 'U', 'G'): -2.24,  # AC/UG

    # GG/CC stack
    ('G', 'G', 'C', 'C'): -3.26,  # GG/CC
    ('C', 'C', 'G', 'G'): -3.26,  # CC/GG

    # Other combinations
    ('A', 'A', 'U', 'U'): -0.93,
}

# Hudson et al. (2013) nearest-neighbor parameters for Ψ-A pairs
# These are internal Ψ-A parameters measured experimentally
# Format: same as above but involving Ψ
# ΔG°37 values in kcal/mol

HUDSON_PSI_A_NN = {
    # Stacks involving Ψ-A pairs
    # From Hudson et al. Table 2 (internal Ψ-A)

    # GΨ/CÃ (where Ã means A opposite Ψ): G-C stack followed by Ψ-A
    ('G', 'A', 'C', 'Ψ'): -2.49,  # 5'SS G-A / U1 C-Ψ → G pairs with C, A pairs with Ψ

    # ΨG/ÃC: Ψ-A stack followed by G-C
    ('A', 'G', 'Ψ', 'C'): -2.74,  # 5'SS A-G / U1 Ψ-C → A pairs with Ψ, G pairs with C

    # AΨ/UÃ: A-U stack followed by Ψ-A (but we don't have this context in U1)
    ('A', 'A', 'U', 'Ψ'): -1.62,

    # ΨA/ÃU: Ψ-A stack followed by A-U
    ('A', 'A', 'Ψ', 'U'): -1.90,

    # CΨ/GÃ: C-G stack followed by Ψ-A
    ('C', 'A', 'G', 'Ψ'): -2.19,

    # ΨC/ÃG: Ψ-A stack followed by C-G
    ('A', 'C', 'Ψ', 'G'): -3.29,

    # Stacks between two Ψ-A pairs (positions +3 and +4 both have Ψ in U1)
    # This is the critical stack for positions +3/+4
    ('A', 'A', 'Ψ', 'Ψ'): -1.75,  # Estimated from Hudson data (average of flanking contexts)
}

# Mismatch penalties
# These are approximate ΔΔG penalties for mismatches
MISMATCH_SINGLE_PENALTY = 2.0  # kcal/mol penalty per single mismatch

# Consecutive mismatch penalties (bubble formation is highly destabilizing)
CONSECUTIVE_MISMATCH_PENALTY = {
    2: 4.0,   # Two consecutive mismatches
    3: 6.5,   # Three consecutive mismatches
    4: 9.0,   # Four consecutive mismatches
}

# Initiation and terminal penalties
INITIATION_DG = 4.09  # kcal/mol (duplex initiation cost)
AU_TERMINAL_PENALTY = 0.45  # kcal/mol per A-U or Ψ-A terminal pair
GC_TERMINAL_PENALTY = 0.0   # No penalty for G-C terminal

# Functional threshold
# ΔG must be more negative than this for the splice site to be functional
FUNCTIONAL_THRESHOLD = -8.1  # kcal/mol (from Freund et al. 2005)

# ============================================================================
# Data Classes
# ============================================================================

@dataclass
class DuplexResult:
    """Result of U1:5'SS duplex ΔG calculation."""
    sequence_11mer: str
    delta_g: float              # Total ΔG (kcal/mol)
    stack_contributions: List[Tuple[str, float]]  # Each stack's contribution
    mismatch_positions: List[int]  # 0-indexed positions with mismatches
    consecutive_mismatch_penalty: float
    is_functional: bool         # ΔG ≤ threshold
    notes: List[str]

    def to_dict(self) -> Dict:
        """Convert to dictionary for CSV export."""
        return {
            "sequence": self.sequence_11mer,
            "delta_g": round(self.delta_g, 2),
            "is_functional": self.is_functional,
            "mismatch_count": len(self.mismatch_positions),
            "mismatch_positions": ",".join(str(p) for p in self.mismatch_positions),
            "consecutive_mm_penalty": round(self.consecutive_mismatch_penalty, 2),
            "notes": "; ".join(self.notes) if self.notes else "",
        }


# ============================================================================
# Core Functions
# ============================================================================

def dna_to_rna(sequence: str) -> str:
    """Convert DNA sequence to RNA (T → U)."""
    return sequence.upper().replace('T', 'U')


def is_watson_crick_pair(nt1: str, nt2: str) -> bool:
    """Check if two nucleotides form a Watson-Crick pair."""
    return WC_PAIRS.get(nt1) == nt2 or WC_PAIRS.get(nt2) == nt1


def is_wobble_pair(nt1: str, nt2: str) -> bool:
    """Check if two nucleotides form a G-U wobble pair."""
    return (nt1, nt2) in WOBBLE_PAIRS or (nt2, nt1) in WOBBLE_PAIRS


def is_base_paired(nt1: str, nt2: str) -> bool:
    """Check if two nucleotides can form a base pair (WC or wobble)."""
    return is_watson_crick_pair(nt1, nt2) or is_wobble_pair(nt1, nt2)


def get_stack_dg(nt_5ss_1: str, nt_5ss_2: str, nt_u1_1: str, nt_u1_2: str) -> float:
    """
    Get ΔG for a base pair stack.

    Args:
        nt_5ss_1, nt_5ss_2: Consecutive nucleotides on 5'SS (5'→3')
        nt_u1_1, nt_u1_2: Corresponding U1 nucleotides (3'→5', so antiparallel)

    Returns:
        ΔG for this stack in kcal/mol
    """
    # Check if this involves Ψ (positions +3/+4)
    if nt_u1_1 == 'Ψ' or nt_u1_2 == 'Ψ':
        # Use Hudson parameters if available
        key = (nt_5ss_1, nt_5ss_2, nt_u1_1, nt_u1_2)
        if key in HUDSON_PSI_A_NN:
            return HUDSON_PSI_A_NN[key]
        # Try reverse
        rev_key = (nt_5ss_2, nt_5ss_1, nt_u1_2, nt_u1_1)
        if rev_key in HUDSON_PSI_A_NN:
            return HUDSON_PSI_A_NN[rev_key]

    # Use standard Turner parameters
    key = (nt_5ss_1, nt_5ss_2, nt_u1_1, nt_u1_2)
    if key in TURNER_NN_DNA_RNA:
        return TURNER_NN_DNA_RNA[key]

    # Try to approximate from individual pairs
    # This is a fallback for combinations not explicitly defined
    dg = 0.0

    # Check pair 1 (5'SS nt1 with U1 nt2)
    if is_watson_crick_pair(nt_5ss_1, nt_u1_2):
        dg += -1.5  # Approximate WC contribution
    elif is_wobble_pair(nt_5ss_1, nt_u1_2):
        dg += -1.0  # Approximate wobble contribution

    # Check pair 2 (5'SS nt2 with U1 nt1)
    if is_watson_crick_pair(nt_5ss_2, nt_u1_1):
        dg += -1.5
    elif is_wobble_pair(nt_5ss_2, nt_u1_1):
        dg += -1.0

    return dg if dg < 0 else 0.5  # Small positive if no pairing


def find_consecutive_mismatches(mismatch_positions: List[int]) -> int:
    """
    Find the longest run of consecutive mismatch positions.

    Args:
        mismatch_positions: List of 0-indexed mismatch positions

    Returns:
        Length of longest consecutive run
    """
    if not mismatch_positions:
        return 0

    mismatch_positions = sorted(mismatch_positions)
    max_run = 1
    current_run = 1

    for i in range(1, len(mismatch_positions)):
        if mismatch_positions[i] == mismatch_positions[i-1] + 1:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1

    return max_run


def calculate_u1_duplex_dg(sequence_11mer: str) -> DuplexResult:
    """
    Calculate ΔG for U1 snRNA binding to 5' splice site.

    This implements a nearest-neighbor model with:
    - Standard Turner parameters for most stacks
    - Hudson et al. (2013) parameters for Ψ-A pairs at +3/+4
    - Penalties for mismatches, especially consecutive ones

    Args:
        sequence_11mer: 11-mer 5'SS sequence (-3 to +8), DNA format (ACGT)

    Returns:
        DuplexResult with ΔG and detailed breakdown
    """
    sequence_11mer = sequence_11mer.upper()

    if len(sequence_11mer) != 11:
        raise ValueError(f"Sequence must be 11 nucleotides, got {len(sequence_11mer)}: {sequence_11mer}")

    # Validate and convert DNA to RNA
    rna_seq = []
    for i, nt in enumerate(sequence_11mer):
        if nt == 'T':
            rna_seq.append('U')
        elif nt in 'ACGU':
            rna_seq.append(nt)
        else:
            raise ValueError(f"Invalid nucleotide '{nt}' at position {i}")

    # U1 binding region partners (3'→5' orientation to match antiparallel)
    # Position 0 (-3) of 5'SS pairs with position 0 (G) of U1
    u1_partners = U1_SNRNA_BINDING  # ["G", "U", "C", "C", "A", "Ψ", "Ψ", "A", "C", "U", "U"]

    # Track pairing at each position
    pairs = []  # List of (5'SS nt, U1 nt, is_paired)
    mismatch_positions = []
    notes = []

    for i in range(11):
        ss_nt = rna_seq[i]
        u1_nt = u1_partners[i]

        is_paired = is_base_paired(ss_nt, u1_nt)
        pairs.append((ss_nt, u1_nt, is_paired))

        if not is_paired:
            mismatch_positions.append(i)

    # Calculate consecutive mismatch penalty
    max_consecutive = find_consecutive_mismatches(mismatch_positions)
    consec_penalty = CONSECUTIVE_MISMATCH_PENALTY.get(max_consecutive, 0.0)
    if max_consecutive > 4:
        consec_penalty = 9.0 + (max_consecutive - 4) * 2.5  # Linear extrapolation

    if max_consecutive >= 2:
        notes.append(f"{max_consecutive} consecutive mismatches (+{consec_penalty:.1f} kcal/mol)")

    # Start with initiation penalty
    total_dg = INITIATION_DG
    stack_contributions = []

    # Calculate stack contributions
    for i in range(10):  # 10 stacks for 11 positions
        ss_nt_1 = rna_seq[i]
        ss_nt_2 = rna_seq[i + 1]
        u1_nt_1 = u1_partners[i]
        u1_nt_2 = u1_partners[i + 1]

        pair_1_ok = pairs[i][2]
        pair_2_ok = pairs[i + 1][2]

        if pair_1_ok and pair_2_ok:
            # Both positions are paired - use nearest-neighbor
            stack_dg = get_stack_dg(ss_nt_1, ss_nt_2, u1_nt_2, u1_nt_1)
            stack_contributions.append((f"{ss_nt_1}{ss_nt_2}/{u1_nt_1}{u1_nt_2}", stack_dg))
            total_dg += stack_dg
        elif pair_1_ok or pair_2_ok:
            # One position paired, one mismatched - partial destabilization
            # Single mismatch disrupts one stack
            penalty = MISMATCH_SINGLE_PENALTY / 2  # Split penalty between adjacent stacks
            stack_contributions.append((f"{ss_nt_1}{ss_nt_2}:half_mm", penalty))
            total_dg += penalty
        else:
            # Both mismatched - this stack contributes nothing (already penalized)
            stack_contributions.append((f"{ss_nt_1}{ss_nt_2}:double_mm", 0.0))

    # Add consecutive mismatch penalty
    total_dg += consec_penalty

    # Terminal penalties
    # Position 0 (-3) is 5' terminal of duplex
    if pairs[0][2]:  # First pair exists
        if pairs[0][1] in ('U', 'A', 'Ψ'):  # A-U or Ψ-A terminal
            total_dg += AU_TERMINAL_PENALTY
            notes.append(f"5' terminal A-U/Ψ-A penalty (+{AU_TERMINAL_PENALTY})")

    # Position 10 (+8) is 3' terminal of duplex
    if pairs[10][2]:  # Last pair exists
        if pairs[10][1] in ('U', 'A', 'Ψ'):
            total_dg += AU_TERMINAL_PENALTY
            notes.append(f"3' terminal A-U/Ψ-A penalty (+{AU_TERMINAL_PENALTY})")

    # Check special positions
    # Position -1 (index 2): Critical for AG|GT consensus
    if not pairs[2][2]:
        notes.append(f"Position -1 mismatch ({rna_seq[2]} vs C): non-AG donor")

    # Positions +3/+4 (indices 5,6): Ψ-A enhanced stability
    if pairs[5][2] and rna_seq[5] == 'A':
        notes.append("+3 has Ψ-A enhanced pairing")
    if pairs[6][2] and rna_seq[6] == 'A':
        notes.append("+4 has Ψ-A enhanced pairing")

    # Determine if functional
    is_functional = total_dg <= FUNCTIONAL_THRESHOLD

    return DuplexResult(
        sequence_11mer=sequence_11mer,
        delta_g=round(total_dg, 2),
        stack_contributions=stack_contributions,
        mismatch_positions=mismatch_positions,
        consecutive_mismatch_penalty=consec_penalty,
        is_functional=is_functional,
        notes=notes
    )


# ============================================================================
# Panel Generation Functions
# ============================================================================

STOP_CODONS = {"TAA", "TAG", "TGA"}


def has_frame1_stop_codon(sequence: str) -> bool:
    """
    Check if 11-mer has stop codon in Frame 1.

    Frame 1: codon spans -1/+1/+2, so:
        Codon 2: +3, +4, +5 = indices 5, 6, 7
        Codon 3: +6, +7, +8 = indices 8, 9, 10

    Args:
        sequence: 11-character DNA string

    Returns:
        True if positions +3+4+5 or +6+7+8 form a stop codon
    """
    sequence = sequence.upper()
    codon2 = sequence[5:8]
    codon3 = sequence[8:11]
    return codon2 in STOP_CODONS or codon3 in STOP_CODONS


def generate_all_gt_sequences() -> List[str]:
    """
    Generate all 11-mer sequences with GT at +1/+2.

    Format: NNN N GT NNNNN (positions -3 to +8)
    - Positions +1, +2 must be GT (invariant)
    - All other 9 positions variable (A, C, G, T)

    Returns:
        List of 262,144 (4^9) 11-mer sequences
    """
    bases = ['A', 'C', 'G', 'T']
    sequences = []

    # Variable positions: -3, -2, -1, +3, +4, +5, +6, +7, +8 (9 positions)
    for combo in itertools.product(bases, repeat=9):
        # Format: [pos-3][pos-2][pos-1][G][T][pos+3][pos+4][pos+5][pos+6][pos+7][pos+8]
        seq = f"{combo[0]}{combo[1]}{combo[2]}GT{combo[3]}{combo[4]}{combo[5]}{combo[6]}{combo[7]}{combo[8]}"
        sequences.append(seq)

    return sequences


def filter_frame1_stops(sequences: List[str]) -> List[str]:
    """Remove sequences with stop codons in Frame 1."""
    return [seq for seq in sequences if not has_frame1_stop_codon(seq)]


def calculate_sequence_diversity(seq1: str, seq2: str) -> int:
    """Calculate Hamming distance between two sequences."""
    return sum(1 for a, b in zip(seq1.upper(), seq2.upper()) if a != b)


def generate_dg_based_panel(
    output_dir: str,
    sequences_per_tier: Optional[List[int]] = None,
    verbose: bool = True
) -> List[Dict]:
    """
    Generate tiered panel based on U1 duplex ΔG.

    This is the v6 approach: use ΔG as the primary metric for tier selection,
    not H-bond count. ΔG directly measures binding thermodynamics and accounts
    for Ψ-A enhanced stability and consecutive mismatch penalties.

    Tier definitions (ΔG-based):
        Tier 1: ΔG > -6.0 (weakest) - negative control
        Tier 2: -6.0 to -7.5 - target range
        Tier 3: -7.5 to -9.0 - target range
        Tier 4: -9.0 to -10.5 - moderate
        Tier 5: -10.5 to -12.0 - moderate-strong
        Tier 6: -12.0 to -13.5 - strong
        Tier 7: ΔG < -13.5 (strongest) - positive control

    Args:
        output_dir: Directory for output files
        sequences_per_tier: Number of sequences to select per tier (default: [2,3,3,2,2,2,2])
        verbose: Print progress

    Returns:
        List of dictionaries with panel data
    """
    if sequences_per_tier is None:
        sequences_per_tier = [2, 3, 3, 2, 2, 2, 2]  # 16 total

    # ΔG tier boundaries (ΔG_min, ΔG_max)
    # Note: ΔG is negative, so ">" means "less negative" = weaker
    tier_boundaries = [
        (1, float('inf'), -6.0),     # Tier 1: ΔG > -6.0 (weakest)
        (2, -6.0, -7.5),             # Tier 2: -6.0 to -7.5
        (3, -7.5, -9.0),             # Tier 3: -7.5 to -9.0
        (4, -9.0, -10.5),            # Tier 4: -9.0 to -10.5
        (5, -10.5, -12.0),           # Tier 5: -10.5 to -12.0
        (6, -12.0, -13.5),           # Tier 6: -12.0 to -13.5
        (7, -13.5, float('-inf')),   # Tier 7: ΔG < -13.5 (strongest)
    ]

    tier_annotations = {
        1: ("negative_ctrl", "~0-5%"),
        2: ("target_range", "~5-15%"),
        3: ("target_range", "~15-25%"),
        4: ("moderate", "~25-40%"),
        5: ("moderate_strong", "~40-60%"),
        6: ("strong", "~60-80%"),
        7: ("positive_ctrl", "~80-100%"),
    }

    # Step 1: Generate all GT sequences
    if verbose:
        print("Generating all 262,144 GT sequences...")
    all_seqs = generate_all_gt_sequences()
    if verbose:
        print(f"Generated {len(all_seqs)} sequences")

    # Step 2: Filter stop codons
    if verbose:
        print("Filtering Frame 1 stop codons...")
    filtered_seqs = filter_frame1_stops(all_seqs)
    if verbose:
        print(f"After stop codon filter: {len(filtered_seqs)} sequences")

    # Step 3: Calculate ΔG for all
    if verbose:
        print("Calculating ΔG for all sequences...")

    seq_data = []
    for i, seq in enumerate(filtered_seqs):
        if verbose and i > 0 and i % 50000 == 0:
            print(f"  Processed {i}/{len(filtered_seqs)}...")

        dg_result = calculate_u1_duplex_dg(seq)

        # Also calculate H-bonds for comparison
        # Import from splice_donor_hbonds if available
        try:
            from splice_donor_hbonds import calculate_hbonds
            hbond_result = calculate_hbonds(seq)
            total_hbonds = hbond_result.total_hbonds
            longest_stretch = hbond_result.longest_stretch
        except ImportError:
            # Fallback: count mismatches to estimate H-bonds
            total_hbonds = 11 - len(dg_result.mismatch_positions)
            longest_stretch = 0

        seq_data.append({
            'sequence': seq,
            'delta_g': dg_result.delta_g,
            'total_hbonds': total_hbonds,
            'longest_stretch': longest_stretch,
            'mismatch_count': len(dg_result.mismatch_positions),
            'consec_mm_penalty': dg_result.consecutive_mismatch_penalty,
            'is_functional': dg_result.is_functional,
        })

    if verbose:
        dg_values = [d['delta_g'] for d in seq_data]
        print(f"ΔG range: {min(dg_values):.1f} to {max(dg_values):.1f} kcal/mol")

    # Step 4: Assign to ΔG tiers and select diverse sequences
    panel = []

    for tier_num, dg_min, dg_max in tier_boundaries:
        # Filter sequences in this ΔG range
        # Note: dg_min is the "upper" bound (less negative), dg_max is "lower" (more negative)
        tier_seqs = [
            d for d in seq_data
            if (dg_max < d['delta_g'] <= dg_min) or
               (tier_num == 1 and d['delta_g'] > dg_min) or
               (tier_num == 7 and d['delta_g'] <= dg_max)
        ]

        # Correct the logic for tier boundaries
        tier_seqs = []
        for d in seq_data:
            dg = d['delta_g']
            if tier_num == 1 and dg > -6.0:
                tier_seqs.append(d)
            elif tier_num == 2 and -7.5 < dg <= -6.0:
                tier_seqs.append(d)
            elif tier_num == 3 and -9.0 < dg <= -7.5:
                tier_seqs.append(d)
            elif tier_num == 4 and -10.5 < dg <= -9.0:
                tier_seqs.append(d)
            elif tier_num == 5 and -12.0 < dg <= -10.5:
                tier_seqs.append(d)
            elif tier_num == 6 and -13.5 < dg <= -12.0:
                tier_seqs.append(d)
            elif tier_num == 7 and dg <= -13.5:
                tier_seqs.append(d)

        if not tier_seqs:
            if verbose:
                print(f"Warning: No sequences in Tier {tier_num}")
            continue

        if verbose:
            print(f"Tier {tier_num}: {len(tier_seqs)} sequences available")

        # Sort by ΔG within tier
        tier_seqs.sort(key=lambda x: x['delta_g'])

        # Select diverse sequences
        n_select = sequences_per_tier[tier_num - 1]

        if len(tier_seqs) <= n_select:
            selected = tier_seqs
        else:
            # Select evenly spaced by ΔG, then maximize diversity
            indices = []
            step = len(tier_seqs) / n_select
            for i in range(n_select):
                idx = int(i * step)
                indices.append(idx)

            # Greedy diversity selection
            selected = [tier_seqs[indices[0]]]
            candidates = [tier_seqs[i] for i in indices[1:]]

            while len(selected) < n_select and candidates:
                # Find candidate most different from selected
                best_candidate = None
                best_min_dist = -1

                for cand in candidates:
                    min_dist = min(
                        calculate_sequence_diversity(cand['sequence'], s['sequence'])
                        for s in selected
                    )
                    if min_dist > best_min_dist:
                        best_min_dist = min_dist
                        best_candidate = cand

                if best_candidate:
                    selected.append(best_candidate)
                    candidates.remove(best_candidate)
                else:
                    break

        # Get annotation
        annotation, expected_splicing = tier_annotations.get(tier_num, ("", ""))

        # Add to panel
        for seq_info in selected:
            panel.append({
                'tier': tier_num,
                'sequence': seq_info['sequence'],
                'delta_g': seq_info['delta_g'],
                'total_hbonds': seq_info['total_hbonds'],
                'longest_stretch': seq_info['longest_stretch'],
                'mismatch_count': seq_info['mismatch_count'],
                'is_functional': seq_info['is_functional'],
                'notes': annotation,
                'expected_splicing': expected_splicing,
            })

    # Save panel
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    panel_path = Path(output_dir) / "tiered_panel_dg_v6.csv"

    fieldnames = ['tier', 'sequence', 'delta_g', 'total_hbonds', 'longest_stretch',
                  'mismatch_count', 'is_functional', 'notes', 'expected_splicing']

    with open(panel_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in panel:
            writer.writerow(row)

    if verbose:
        print(f"\nSaved {len(panel)} sequences to {panel_path}")

    # Generate summary
    summary_path = Path(output_dir) / "tiered_panel_dg_v6_summary.txt"
    with open(summary_path, 'w') as f:
        f.write("Tiered Splice Donor Panel (ΔG-Based v6)\n")
        f.write("=" * 60 + "\n\n")
        f.write("Panel selection based on U1 snRNA:5'SS duplex ΔG\n")
        f.write("Uses Hudson et al. (2013) Ψ-A parameters at positions +3/+4\n\n")

        for tier_num in range(1, 8):
            tier_seqs = [p for p in panel if p['tier'] == tier_num]
            if tier_seqs:
                f.write(f"Tier {tier_num} ({tier_annotations[tier_num][0]}, expect {tier_annotations[tier_num][1]} splicing):\n")
                for s in tier_seqs:
                    f.write(f"  {s['sequence']} | ΔG={s['delta_g']:.1f} | H-bonds={s['total_hbonds']} | stretch={s['longest_stretch']}\n")
                f.write("\n")

        f.write(f"Total: {len(panel)} sequences for synthesis\n")

    if verbose:
        print(f"Saved summary to {summary_path}")

    return panel


# ============================================================================
# Panel Generation v7: ΔG + Position +6 Stratification
# ============================================================================

def generate_dg_panel_v7(
    output_dir: str,
    verbose: bool = True
) -> List[Dict]:
    """
    Generate 32-sequence panel with ΔG tiers and +6 nucleotide stratification.

    This v7 panel design allows empirical testing of whether ΔG alone or
    position +6 better predicts splicing efficiency.

    Key design:
    - 7 ΔG tiers (percentile-based boundaries)
    - Per tier: exactly 2 sequences with +6=G, plus 2-3 sequences with +6≠G
    - All sequences have AG|GT consensus (G at -1) and AA at +3/+4
    - No Frame 1 stop codons

    Args:
        output_dir: Directory for output files
        verbose: Print progress

    Returns:
        List of dictionaries with panel data (32 sequences total)
    """
    # Import here to avoid circular dependency
    try:
        from splice_donor_hbonds import calculate_hbonds
    except ImportError:
        calculate_hbonds = None

    # Tier definitions
    # Using percentile-based boundaries on the filtered AG|GT+AA sequence space
    # Total: 7 tiers × (2 G + 2-3 non-G) = ~32 sequences
    # Configuration: 4+4+5+5+5+5+4 = 32
    tier_configs = [
        # (tier_num, percentile_low, percentile_high, n_g6, n_non_g6, expected_splicing, annotation)
        (1, 85, 100, 2, 2, "~0-5%", "negative_ctrl"),      # Weakest ΔG (highest percentile = least negative)
        (2, 70, 85, 2, 2, "~5-15%", "very_weak"),
        (3, 55, 70, 2, 3, "~15-25%", "weak"),
        (4, 40, 55, 2, 3, "~25-40%", "moderate"),
        (5, 25, 40, 2, 3, "~40-60%", "moderate_strong"),   # Changed from 2 to 3 non-G
        (6, 10, 25, 2, 3, "~60-80%", "strong"),
        (7, 0, 10, 2, 2, "~80-100%", "positive_ctrl"),     # Strongest ΔG (lowest percentile = most negative)
    ]

    # Step 1: Generate all AG|GT sequences with AA at +3/+4
    # Format: NNN|GT|AA|NNNN (positions -3 to +8)
    # Fixed: -1=G (AG consensus), +1=G, +2=T, +3=A, +4=A
    # Variable: -3, -2, +5, +6, +7, +8 (6 positions)
    if verbose:
        print("Generating all AG|GT+AA sequences...")

    bases = ['A', 'C', 'G', 'T']
    all_sequences = []

    # Variable positions: -3 (idx 0), -2 (idx 1), +5 (idx 7), +6 (idx 8), +7 (idx 9), +8 (idx 10)
    for combo in itertools.product(bases, repeat=6):
        # Build 11-mer: [pos-3][pos-2][G][G][T][A][A][pos+5][pos+6][pos+7][pos+8]
        seq = f"{combo[0]}{combo[1]}GGT AA{combo[2]}{combo[3]}{combo[4]}{combo[5]}"
        seq = seq.replace(' ', '')  # Remove space used for readability
        all_sequences.append(seq)

    if verbose:
        print(f"  Generated {len(all_sequences)} AG|GT+AA sequences (4^6 = 4096)")

    # Step 2: Filter Frame 1 stop codons
    if verbose:
        print("Filtering Frame 1 stop codons...")

    filtered_seqs = filter_frame1_stops(all_sequences)

    if verbose:
        print(f"  After filter: {len(filtered_seqs)} sequences")

    # Step 3: Calculate ΔG and H-bonds for all sequences
    if verbose:
        print("Calculating ΔG for all sequences...")

    seq_data = []
    for seq in filtered_seqs:
        dg_result = calculate_u1_duplex_dg(seq)

        # Calculate H-bonds
        if calculate_hbonds:
            hbond_result = calculate_hbonds(seq)
            total_hbonds = hbond_result.total_hbonds
            longest_stretch = hbond_result.longest_stretch
        else:
            # Fallback: estimate from mismatch count
            total_hbonds = 11 - len(dg_result.mismatch_positions)
            longest_stretch = 0

        # Position +6 is index 8 (0-indexed)
        pos6_nt = seq[8]

        seq_data.append({
            'sequence': seq,
            'delta_g': dg_result.delta_g,
            'total_hbonds': total_hbonds,
            'longest_stretch': longest_stretch,
            'mismatch_count': len(dg_result.mismatch_positions),
            'pos6_nt': pos6_nt,
            'pos6_is_g': pos6_nt == 'G',
        })

    # Sort by ΔG (ascending = most negative first)
    seq_data.sort(key=lambda x: x['delta_g'])

    if verbose:
        dg_values = [d['delta_g'] for d in seq_data]
        print(f"  ΔG range: {min(dg_values):.2f} to {max(dg_values):.2f} kcal/mol")

        # Count +6=G vs non-G
        n_g6 = sum(1 for d in seq_data if d['pos6_is_g'])
        print(f"  +6=G: {n_g6}, +6≠G: {len(seq_data) - n_g6}")

    # Step 4: Calculate percentile-based ΔG boundaries
    n_seqs = len(seq_data)

    def get_percentile_dg(percentile: float) -> float:
        """Get ΔG value at given percentile (0=strongest/most negative, 100=weakest/least negative)."""
        # Percentile 0 = most negative ΔG (strongest binder)
        # Percentile 100 = least negative ΔG (weakest binder)
        idx = int(percentile / 100 * (n_seqs - 1))
        return seq_data[idx]['delta_g']

    if verbose:
        print("\nPercentile-based ΔG tier boundaries:")
        for tier_num, p_low, p_high, _, _, _, _ in tier_configs:
            dg_low = get_percentile_dg(p_low)
            dg_high = get_percentile_dg(p_high)
            print(f"  Tier {tier_num}: percentile {p_low}-{p_high}% → ΔG {dg_low:.2f} to {dg_high:.2f}")

    # Step 5: Select sequences for each tier
    panel = []
    used_sequences = set()  # Track already selected sequences to avoid duplicates

    for tier_num, p_low, p_high, n_g6_target, n_non_g6_target, expected_splicing, annotation in tier_configs:
        # Calculate ΔG boundaries for this tier
        dg_low = get_percentile_dg(p_low)
        dg_high = get_percentile_dg(p_high)

        # Filter sequences in this ΔG range
        # Note: dg_low is the boundary at the lower percentile (more negative)
        #       dg_high is at the higher percentile (less negative)
        tier_seqs = [
            d for d in seq_data
            if dg_low <= d['delta_g'] <= dg_high
        ]

        if verbose:
            print(f"\nTier {tier_num}: {len(tier_seqs)} sequences in ΔG range [{dg_low:.2f}, {dg_high:.2f}]")

        if not tier_seqs:
            if verbose:
                print(f"  Warning: No sequences found!")
            continue

        # Exclude already-used sequences
        tier_seqs = [d for d in tier_seqs if d['sequence'] not in used_sequences]

        if not tier_seqs:
            if verbose:
                print(f"  No unused sequences available!")
            continue

        # Split by +6 nucleotide
        g6_seqs = [d for d in tier_seqs if d['pos6_is_g']]
        non_g6_seqs = [d for d in tier_seqs if not d['pos6_is_g']]

        if verbose:
            print(f"  +6=G: {len(g6_seqs)}, +6≠G: {len(non_g6_seqs)}")

        # Select +6=G sequences (exactly n_g6_target, spread across ΔG range)
        selected_g6 = _select_diverse_from_tier(g6_seqs, n_g6_target)

        # Select +6≠G sequences (n_non_g6_target, prefer diverse nucleotides A, C, T)
        selected_non_g6 = _select_diverse_non_g6(non_g6_seqs, n_non_g6_target)

        if verbose:
            print(f"  Selected: {len(selected_g6)} with +6=G, {len(selected_non_g6)} with +6≠G")

        # Add to panel and mark as used
        for seq_info in selected_g6 + selected_non_g6:
            used_sequences.add(seq_info['sequence'])
            panel.append({
                'tier': tier_num,
                'sequence': seq_info['sequence'],
                'delta_g': seq_info['delta_g'],
                'total_hbonds': seq_info['total_hbonds'],
                'longest_stretch': seq_info['longest_stretch'],
                'pos6_nt': seq_info['pos6_nt'],
                'mismatch_count': seq_info['mismatch_count'],
                'expected_splicing': expected_splicing,
                'annotation': annotation,
            })

    if verbose:
        print(f"\nTotal panel size: {len(panel)} sequences")
        g6_count = sum(1 for p in panel if p['pos6_nt'] == 'G')
        print(f"  +6=G: {g6_count}, +6≠G: {len(panel) - g6_count}")

    # Step 6: Save panel to CSV
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    panel_path = Path(output_dir) / "tiered_panel_32_v7.csv"

    fieldnames = ['tier', 'sequence', 'delta_g', 'total_hbonds', 'longest_stretch',
                  'pos6_nt', 'mismatch_count', 'expected_splicing', 'annotation']

    with open(panel_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in panel:
            writer.writerow(row)

    if verbose:
        print(f"\nSaved panel to: {panel_path}")

    # Step 7: Generate summary text file
    summary_path = Path(output_dir) / "tiered_panel_32_v7_summary.txt"
    with open(summary_path, 'w') as f:
        f.write("Tiered Splice Donor Panel v7 (32 Sequences)\n")
        f.write("=" * 70 + "\n\n")
        f.write("Design: ΔG-based tiers with position +6 stratification\n")
        f.write("Constraints:\n")
        f.write("  - AG|GT consensus (G at position -1)\n")
        f.write("  - AA at positions +3/+4 (for Ψ-A enhanced stability)\n")
        f.write("  - No Frame 1 stop codons\n")
        f.write(f"\nTotal sequences: {len(panel)}\n")
        f.write(f"  +6=G: {sum(1 for p in panel if p['pos6_nt'] == 'G')}\n")
        f.write(f"  +6≠G: {sum(1 for p in panel if p['pos6_nt'] != 'G')}\n\n")

        for tier_num in range(1, 8):
            tier_seqs = [p for p in panel if p['tier'] == tier_num]
            if tier_seqs:
                f.write(f"\nTier {tier_num} ({tier_seqs[0]['annotation']}, expect {tier_seqs[0]['expected_splicing']} splicing):\n")
                f.write("-" * 70 + "\n")
                for s in tier_seqs:
                    g6_mark = "(+6=G)" if s['pos6_nt'] == 'G' else f"(+6={s['pos6_nt']})"
                    f.write(f"  {s['sequence']} | ΔG={s['delta_g']:6.2f} | H-bonds={s['total_hbonds']:2d} | stretch={s['longest_stretch']} | {g6_mark}\n")

    if verbose:
        print(f"Saved summary to: {summary_path}")

    # Step 8: Generate distribution plot
    _generate_dg_distribution_plot_v7(seq_data, panel, output_dir, verbose)

    return panel


def _select_diverse_from_tier(seqs: List[Dict], n_select: int) -> List[Dict]:
    """
    Select diverse sequences from a tier, spread across the ΔG range.

    Args:
        seqs: List of sequence data dictionaries, sorted by ΔG
        n_select: Number of sequences to select

    Returns:
        List of selected sequence dictionaries
    """
    if len(seqs) <= n_select:
        return seqs

    # Sort by ΔG
    seqs = sorted(seqs, key=lambda x: x['delta_g'])

    # Select evenly spaced across the range
    indices = []
    if n_select == 1:
        indices = [len(seqs) // 2]
    else:
        step = (len(seqs) - 1) / (n_select - 1)
        for i in range(n_select):
            idx = int(round(i * step))
            indices.append(min(idx, len(seqs) - 1))

    # Remove duplicates while preserving order
    seen = set()
    unique_indices = []
    for idx in indices:
        if idx not in seen:
            seen.add(idx)
            unique_indices.append(idx)

    selected = [seqs[i] for i in unique_indices]

    # If we need more due to duplicates, add from the middle
    while len(selected) < n_select and len(selected) < len(seqs):
        for i in range(len(seqs)):
            if seqs[i] not in selected:
                selected.append(seqs[i])
                if len(selected) >= n_select:
                    break

    return selected[:n_select]


def _select_diverse_non_g6(seqs: List[Dict], n_select: int) -> List[Dict]:
    """
    Select diverse non-G sequences at position +6, preferring diverse nucleotides (A, C, T).

    Args:
        seqs: List of sequence data dictionaries with +6≠G
        n_select: Number of sequences to select

    Returns:
        List of selected sequence dictionaries
    """
    if len(seqs) <= n_select:
        return seqs

    # Group by +6 nucleotide
    by_nt = {'A': [], 'C': [], 'T': []}
    for s in seqs:
        nt = s['pos6_nt']
        if nt in by_nt:
            by_nt[nt].append(s)

    # Sort each group by ΔG
    for nt in by_nt:
        by_nt[nt].sort(key=lambda x: x['delta_g'])

    # Select one from each nucleotide type if possible, then fill remaining
    selected = []

    # Round-robin selection from each nucleotide type
    nt_order = ['T', 'C', 'A']  # Prefer T first (wobble pair with U in U1 at +6)
    nt_indices = {nt: 0 for nt in nt_order}

    while len(selected) < n_select:
        added_this_round = False

        for nt in nt_order:
            if len(selected) >= n_select:
                break
            if nt_indices[nt] < len(by_nt[nt]):
                candidate = by_nt[nt][nt_indices[nt]]
                # Check for sequence diversity (at least 2 positions different)
                if not selected or min(
                    calculate_sequence_diversity(candidate['sequence'], s['sequence'])
                    for s in selected
                ) >= 2:
                    selected.append(candidate)
                    added_this_round = True
                nt_indices[nt] += 1

        # If couldn't add any diverse sequences, just add the next available
        if not added_this_round:
            for nt in nt_order:
                while nt_indices[nt] < len(by_nt[nt]):
                    candidate = by_nt[nt][nt_indices[nt]]
                    nt_indices[nt] += 1
                    if candidate not in selected:
                        selected.append(candidate)
                        break
            if len(selected) >= n_select:
                break

        # Safety check for infinite loop
        if all(nt_indices[nt] >= len(by_nt[nt]) for nt in nt_order):
            break

    return selected[:n_select]


def _generate_dg_distribution_plot_v7(
    all_seq_data: List[Dict],
    panel: List[Dict],
    output_dir: str,
    verbose: bool = True
) -> None:
    """
    Generate distribution plot showing all AG|GT+AA sequences and selected panel.

    Args:
        all_seq_data: All sequence data with ΔG values
        panel: Selected panel sequences
        output_dir: Output directory
        verbose: Print progress
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import numpy as np
    except ImportError:
        if verbose:
            print("matplotlib not available, skipping distribution plot")
        return

    fig, ax = plt.subplots(figsize=(14, 8))

    # Extract ΔG values
    all_dg = [d['delta_g'] for d in all_seq_data]

    # Create histogram of all sequences
    n_bins = 50
    counts, bins, patches = ax.hist(all_dg, bins=n_bins, color='lightgray',
                                     edgecolor='gray', alpha=0.7, label='All AG|GT+AA sequences')

    # Overlay selected panel sequences
    # Color by tier, shape by +6 nucleotide
    tier_colors = {
        1: '#d73027',  # Red (weakest)
        2: '#fc8d59',  # Orange
        3: '#fee08b',  # Yellow
        4: '#d9ef8b',  # Yellow-green
        5: '#91cf60',  # Light green
        6: '#1a9850',  # Green
        7: '#006837',  # Dark green (strongest)
    }

    # Plot panel sequences
    for p in panel:
        tier = p['tier']
        dg = p['delta_g']
        is_g6 = p['pos6_nt'] == 'G'

        # Find y position (histogram height at this ΔG)
        bin_idx = np.searchsorted(bins[:-1], dg, side='right') - 1
        bin_idx = max(0, min(bin_idx, len(counts) - 1))
        y_pos = counts[bin_idx] + 50  # Offset above histogram

        color = tier_colors.get(tier, 'black')
        marker = 'o' if is_g6 else '^'  # Circle for G, triangle for non-G
        size = 100 if is_g6 else 80

        ax.scatter(dg, y_pos, c=color, marker=marker, s=size,
                   edgecolors='black', linewidths=1, zorder=10)

    # Labels and title
    ax.set_xlabel('ΔG (kcal/mol)', fontsize=12)
    ax.set_ylabel('Number of Sequences', fontsize=12)
    ax.set_title('v7 Panel: ΔG Distribution with +6 Nucleotide Stratification\n'
                 '(AG|GT consensus, AA at +3/+4, no Frame 1 stops)', fontsize=14)

    # Create legend
    legend_elements = [
        mpatches.Patch(facecolor='lightgray', edgecolor='gray', alpha=0.7,
                       label=f'All sequences (n={len(all_seq_data)})'),
    ]
    for tier in range(1, 8):
        tier_seqs = [p for p in panel if p['tier'] == tier]
        if tier_seqs:
            legend_elements.append(
                mpatches.Patch(facecolor=tier_colors[tier], edgecolor='black',
                               label=f'Tier {tier} (n={len(tier_seqs)})')
            )

    # Add marker legend
    legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                                      markersize=10, label='+6=G'))
    legend_elements.append(plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='gray',
                                      markersize=10, label='+6≠G'))

    ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

    # Add annotation
    g6_count = sum(1 for p in panel if p['pos6_nt'] == 'G')
    non_g6_count = len(panel) - g6_count
    annotation_text = (
        f"Panel: {len(panel)} sequences\n"
        f"+6=G: {g6_count} (●)\n"
        f"+6≠G: {non_g6_count} (▲)"
    )
    ax.text(0.02, 0.98, annotation_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()

    plot_path = Path(output_dir) / "dg_distribution_v7.png"
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()

    if verbose:
        print(f"Saved distribution plot to: {plot_path}")


# ============================================================================
# CLI Interface
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="U1 snRNA : 5'SS Duplex ΔG Calculator with Ψ-A Parameters",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Calculate ΔG for a single sequence
  python u1_duplex_dg.py --sequence CAGGTAAGTAT

  # Score multiple sequences from a panel CSV
  python u1_duplex_dg.py --panel tiered_panel_candidates.csv --output dg_scores.csv

  # Generate new ΔG-based tiered panel (v6)
  python u1_duplex_dg.py --generate-panel --output-dir ./output

  # Generate v7 panel (32 sequences with +6 stratification)
  python u1_duplex_dg.py --generate-panel-v7 --output-dir ./output

  # Analyze full sequence space (262k sequences)
  python u1_duplex_dg.py --analyze-space --output-dir ./output

  # Compare ΔG with H-bond predictions
  python u1_duplex_dg.py --compare-metrics --panel panel.csv --output comparison.csv
        """
    )

    parser.add_argument("--sequence", "-s", type=str,
                        help="Single 11-mer sequence to analyze")
    parser.add_argument("--panel", "-p", type=str,
                        help="CSV file with 'sequence' column to analyze")
    parser.add_argument("--output", "-o", type=str,
                        help="Output CSV file path")
    parser.add_argument("--output-dir", type=str, default="./output",
                        help="Output directory for generated files")
    parser.add_argument("--generate-panel", action="store_true",
                        help="Generate new ΔG-based tiered panel (v6)")
    parser.add_argument("--generate-panel-v7", action="store_true",
                        help="Generate v7 panel (32 sequences with ΔG tiers + position +6 stratification)")
    parser.add_argument("--analyze-space", action="store_true",
                        help="Analyze full 262k sequence space")
    parser.add_argument("--compare-metrics", action="store_true",
                        help="Compare ΔG with H-bond predictions")
    parser.add_argument("--threshold", "-t", type=float, default=FUNCTIONAL_THRESHOLD,
                        help=f"ΔG threshold for functional splice site (default: {FUNCTIONAL_THRESHOLD})")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

    # Single sequence analysis
    if args.sequence:
        try:
            result = calculate_u1_duplex_dg(args.sequence)
            print(f"Sequence: {result.sequence_11mer}")
            print(f"ΔG: {result.delta_g:.2f} kcal/mol")
            print(f"Functional (ΔG ≤ {args.threshold}): {result.is_functional}")
            print(f"Mismatches at positions: {[POSITION_LABELS[i] for i in result.mismatch_positions]}")
            print(f"Consecutive mismatch penalty: {result.consecutive_mismatch_penalty:.1f} kcal/mol")
            print()
            print("Notes:")
            for note in result.notes:
                print(f"  - {note}")
            print()
            print("Stack contributions:")
            for stack, dg in result.stack_contributions:
                print(f"  {stack}: {dg:+.2f}")
        except ValueError as e:
            print(f"Error: {e}")
            return 1
        return 0

    # Panel analysis
    if args.panel:
        results = []
        with open(args.panel) as f:
            reader = csv.DictReader(f)
            for row in reader:
                seq = row.get('sequence', row.get('donor_11mer', ''))
                if seq and len(seq) == 11:
                    try:
                        result = calculate_u1_duplex_dg(seq)
                        results.append({
                            'sequence': seq,
                            'delta_g': result.delta_g,
                            'is_functional': result.is_functional,
                            'mismatch_count': len(result.mismatch_positions),
                            'consecutive_penalty': result.consecutive_mismatch_penalty,
                            **{k: v for k, v in row.items() if k != 'sequence'}
                        })
                    except ValueError as e:
                        print(f"Warning: Could not process {seq}: {e}")

        if args.output:
            Path(args.output).parent.mkdir(parents=True, exist_ok=True)
            with open(args.output, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=results[0].keys())
                writer.writeheader()
                writer.writerows(results)
            print(f"Wrote {len(results)} results to {args.output}")
        else:
            for r in results:
                status = "✓" if r['is_functional'] else "✗"
                print(f"{status} {r['sequence']}: ΔG={r['delta_g']:.1f}, functional={r['is_functional']}")

        return 0

    # Generate ΔG-based panel (v6)
    if args.generate_panel:
        print("Generating ΔG-based tiered panel (v6)...")
        panel = generate_dg_based_panel(args.output_dir, verbose=True)
        print(f"\nGenerated panel with {len(panel)} sequences")
        return 0

    # Generate v7 panel (32 sequences with +6 stratification)
    if args.generate_panel_v7:
        print("Generating v7 panel (32 sequences with ΔG + position +6 stratification)...")
        print()
        panel = generate_dg_panel_v7(args.output_dir, verbose=True)
        print(f"\n{'='*60}")
        print(f"Generated v7 panel with {len(panel)} sequences")
        print(f"Output directory: {args.output_dir}")
        print(f"Files created:")
        print(f"  - tiered_panel_32_v7.csv")
        print(f"  - tiered_panel_32_v7_summary.txt")
        print(f"  - dg_distribution_v7.png")
        return 0

    # Analyze full sequence space
    if args.analyze_space:
        print("Analyzing full sequence space (262,144 sequences)...")
        print("This may take a few minutes...\n")

        all_seqs = generate_all_gt_sequences()
        filtered_seqs = filter_frame1_stops(all_seqs)

        print(f"Total GT sequences: {len(all_seqs)}")
        print(f"After Frame 1 stop filter: {len(filtered_seqs)}")

        # Calculate ΔG distribution
        dg_values = []
        for i, seq in enumerate(filtered_seqs[:10000]):  # Sample first 10k for quick analysis
            if i > 0 and i % 2000 == 0:
                print(f"  Processed {i}/10000...")
            result = calculate_u1_duplex_dg(seq)
            dg_values.append(result.delta_g)

        print(f"\nΔG distribution (sample of 10,000):")
        print(f"  Min: {min(dg_values):.1f} kcal/mol")
        print(f"  Max: {max(dg_values):.1f} kcal/mol")
        print(f"  Mean: {sum(dg_values)/len(dg_values):.1f} kcal/mol")

        # Count functional
        n_functional = sum(1 for dg in dg_values if dg <= FUNCTIONAL_THRESHOLD)
        print(f"  Functional (ΔG ≤ {FUNCTIONAL_THRESHOLD}): {n_functional} ({100*n_functional/len(dg_values):.1f}%)")

        return 0

    # Default: show help
    parser.print_help()
    return 0


if __name__ == "__main__":
    exit(main())
