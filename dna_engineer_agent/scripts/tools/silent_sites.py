#!/usr/bin/env python3
"""
Silent Restriction Site Finder

A tool for identifying restriction enzyme recognition sites that can be
introduced into DNA sequences with minimal mutations while maintaining
the encoded protein sequence (silent mutations).

Author: Synthetic Biology Design Tool
License: MIT
"""

import argparse
import sys
import csv
from datetime import datetime
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass


# ============================================================================
# IUPAC Nucleotide Ambiguity Codes
# ============================================================================

IUPAC_CODES = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'R': {'A', 'G'},      # puRine
    'Y': {'C', 'T'},      # pYrimidine
    'S': {'G', 'C'},      # Strong (3 H-bonds)
    'W': {'A', 'T'},      # Weak (2 H-bonds)
    'K': {'G', 'T'},      # Keto
    'M': {'A', 'C'},      # aMino
    'B': {'C', 'G', 'T'}, # not A
    'D': {'A', 'G', 'T'}, # not C
    'H': {'A', 'C', 'T'}, # not G
    'V': {'A', 'C', 'G'}, # not T
    'N': {'A', 'C', 'G', 'T'}  # aNy
}


def matches_iupac(base: str, code: str) -> bool:
    """
    Check if a DNA base matches an IUPAC code

    Args:
        base: Single DNA nucleotide (A, C, G, T)
        code: IUPAC ambiguity code

    Returns:
        True if base is in the set represented by code
    """
    base = base.upper()
    code = code.upper()

    if code not in IUPAC_CODES:
        return False

    return base in IUPAC_CODES[code]


def reverse_complement_iupac(seq: str) -> str:
    """
    Get reverse complement of a DNA sequence with IUPAC codes

    Args:
        seq: DNA sequence with IUPAC codes

    Returns:
        Reverse complement sequence with IUPAC codes
    """
    # IUPAC complement table
    iupac_complement = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
        'R': 'Y',  # R (A/G) -> Y (T/C)
        'Y': 'R',  # Y (C/T) -> R (G/A)
        'S': 'S',  # S (G/C) -> S (C/G)
        'W': 'W',  # W (A/T) -> W (T/A)
        'K': 'M',  # K (G/T) -> M (C/A)
        'M': 'K',  # M (A/C) -> K (T/G)
        'B': 'V',  # B (C/G/T) -> V (G/C/A)
        'D': 'H',  # D (A/G/T) -> H (A/C/T)
        'H': 'D',  # H (A/C/T) -> D (A/G/T)
        'V': 'B',  # V (A/C/G) -> B (T/G/C)
        'N': 'N'   # N (any) -> N (any)
    }

    return ''.join(iupac_complement.get(base.upper(), base) for base in reversed(seq))


def is_palindromic(iupac_pattern: str) -> bool:
    """
    Check if an IUPAC pattern is palindromic (same as its reverse complement)

    Args:
        iupac_pattern: Pattern with IUPAC codes

    Returns:
        True if palindromic, False otherwise
    """
    return iupac_pattern.upper() == reverse_complement_iupac(iupac_pattern.upper())


def matches_iupac_pattern(dna_window: str, iupac_pattern: str) -> bool:
    """
    Check if a DNA window matches an IUPAC pattern

    Args:
        dna_window: DNA sequence
        iupac_pattern: Pattern with IUPAC codes

    Returns:
        True if the window matches the pattern
    """
    if len(dna_window) != len(iupac_pattern):
        return False

    for base, code in zip(dna_window, iupac_pattern):
        if not matches_iupac(base, code):
            return False

    return True


def count_pattern_occurrences(dna_seq: str, iupac_pattern: str) -> int:
    """
    Count occurrences of an IUPAC pattern in a DNA sequence (both strands)

    For restriction enzymes, we need to check both the forward strand and reverse complement:
    - Palindromic sites (e.g., GAATTC): appear on both strands at same position → count once
    - Non-palindromic sites (e.g., GAAGAC): can appear on either strand → count both

    Args:
        dna_seq: DNA sequence to search
        iupac_pattern: Pattern with IUPAC codes

    Returns:
        Total number of restriction sites (cutting positions)
    """
    dna_seq = dna_seq.upper()
    iupac_pattern = iupac_pattern.upper()

    # Count on forward strand
    forward_count = 0
    pattern_len = len(iupac_pattern)

    for i in range(len(dna_seq) - pattern_len + 1):
        window = dna_seq[i:i + pattern_len]
        if matches_iupac_pattern(window, iupac_pattern):
            forward_count += 1

    # Check if pattern is palindromic
    if is_palindromic(iupac_pattern):
        # For palindromes, forward strand count equals total sites
        # (same sites appear on reverse complement at same positions)
        return forward_count

    # For non-palindromes, also count on reverse complement strand
    reverse_comp = reverse_complement_iupac(dna_seq)
    reverse_count = 0

    for i in range(len(reverse_comp) - pattern_len + 1):
        window = reverse_comp[i:i + pattern_len]
        if matches_iupac_pattern(window, iupac_pattern):
            reverse_count += 1

    return forward_count + reverse_count


def count_pattern_occurrences_in_range(dna_seq: str, iupac_pattern: str, start_pos: int, end_pos: int) -> int:
    """
    Count occurrences of an IUPAC pattern that START within a specific range (both strands).
    Sites can extend beyond the range boundary, but must start within [start_pos, end_pos).

    For restriction enzymes, we need to check both the forward strand and reverse complement:
    - Palindromic sites: appear on both strands at same position → count once
    - Non-palindromic sites: can appear on either strand → count both

    Args:
        dna_seq: DNA sequence to search
        iupac_pattern: Pattern with IUPAC codes
        start_pos: Start position (0-indexed, inclusive)
        end_pos: End position (0-indexed, exclusive)

    Returns:
        Total number of restriction sites starting within the specified range
    """
    dna_seq = dna_seq.upper()
    iupac_pattern = iupac_pattern.upper()

    # Count on forward strand
    forward_count = 0
    pattern_len = len(iupac_pattern)

    # Only search positions that start within the range
    for i in range(start_pos, min(end_pos, len(dna_seq) - pattern_len + 1)):
        window = dna_seq[i:i + pattern_len]
        if matches_iupac_pattern(window, iupac_pattern):
            forward_count += 1

    # Check if pattern is palindromic
    if is_palindromic(iupac_pattern):
        # For palindromes, forward strand count equals total sites
        return forward_count

    # For non-palindromes, also count on reverse complement strand
    reverse_comp = reverse_complement_iupac(dna_seq)
    reverse_count = 0

    # For reverse complement, a site at position i in the reverse complement
    # corresponds to position (len-i-pattern_len) in the forward strand
    # We want sites that START in the range [start_pos, end_pos) in forward coordinates
    for i in range(len(reverse_comp) - pattern_len + 1):
        # Convert reverse complement position to forward position
        forward_pos = len(dna_seq) - i - pattern_len

        # Check if this position is within our range
        if start_pos <= forward_pos < end_pos:
            window = reverse_comp[i:i + pattern_len]
            if matches_iupac_pattern(window, iupac_pattern):
                reverse_count += 1

    return forward_count + reverse_count


# ============================================================================
# Genetic Code (Standard)
# ============================================================================

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
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def translate(dna_seq: str, frame: int = 0, stop_at_terminator: bool = True) -> str:
    """
    Translate DNA sequence to protein

    Args:
        dna_seq: DNA sequence string
        frame: Reading frame offset (0, 1, or 2)
        stop_at_terminator: If True, stop at first stop codon. If False, continue through stop codons (represented as '*')

    Returns:
        Protein sequence (single letter codes)
    """
    dna_seq = dna_seq.upper()[frame:]
    protein = []

    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        if len(codon) == 3:
            aa = CODON_TABLE.get(codon, 'X')  # X for unknown
            if aa == '*' and stop_at_terminator:  # Stop codon
                break
            protein.append(aa)

    return ''.join(protein)


def find_protein_in_dna(dna_seq: str, protein_seq: str) -> Optional[Tuple[int, int, int]]:
    """
    Find protein sequence within DNA by trying all reading frames

    Args:
        dna_seq: DNA sequence
        protein_seq: Target protein sequence

    Returns:
        Tuple of (start_pos, end_pos, frame) or None if not found
    """
    protein_seq = protein_seq.upper()

    for frame in [0, 1, 2]:
        # Translate through stop codons to find proteins anywhere in the sequence
        translated = translate(dna_seq, frame, stop_at_terminator=False)

        if protein_seq in translated:
            # Find the position in protein
            protein_start = translated.index(protein_seq)

            # Convert to DNA coordinates
            dna_start = frame + (protein_start * 3)
            dna_end = dna_start + (len(protein_seq) * 3)

            return (dna_start, dna_end, frame)

    return None


# ============================================================================
# Edit Distance with IUPAC Support
# ============================================================================

def iupac_edit_distance(seq1: str, seq2_iupac: str) -> Tuple[int, List[str]]:
    """
    Calculate minimum edit distance between a DNA sequence and an IUPAC pattern
    Uses dynamic programming (Levenshtein distance) with IUPAC matching

    Args:
        seq1: DNA sequence (concrete bases only)
        seq2_iupac: Target sequence with IUPAC codes

    Returns:
        Tuple of (edit_distance, list of mutation descriptions)
    """
    seq1 = seq1.upper()
    seq2_iupac = seq2_iupac.upper()

    n = len(seq1)
    m = len(seq2_iupac)

    # DP matrix: dp[i][j] = min edits to transform seq1[:i] to seq2_iupac[:j]
    dp = [[float('inf')] * (m + 1) for _ in range(n + 1)]

    # Backtracking matrix to reconstruct mutations
    backtrack = [[None] * (m + 1) for _ in range(n + 1)]

    # Base cases
    for i in range(n + 1):
        dp[i][0] = i  # Delete all from seq1
        if i > 0:
            backtrack[i][0] = ('del', i-1, seq1[i-1])

    for j in range(m + 1):
        dp[0][j] = j  # Insert all from seq2
        if j > 0:
            backtrack[0][j] = ('ins', 0, seq2_iupac[j-1])

    dp[0][0] = 0

    # Fill DP table
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Match or mismatch cost
            if matches_iupac(seq1[i-1], seq2_iupac[j-1]):
                match_cost = 0
            else:
                match_cost = 1  # Substitution

            match = dp[i-1][j-1] + match_cost
            delete = dp[i-1][j] + 1
            insert = dp[i][j-1] + 1

            min_cost = min(match, delete, insert)
            dp[i][j] = min_cost

            # Record which operation was used
            if min_cost == match:
                if match_cost == 0:
                    backtrack[i][j] = ('match', i-1, seq1[i-1])
                else:
                    backtrack[i][j] = ('sub', i-1, seq1[i-1], seq2_iupac[j-1])
            elif min_cost == delete:
                backtrack[i][j] = ('del', i-1, seq1[i-1])
            else:
                backtrack[i][j] = ('ins', i-1, seq2_iupac[j-1])

    # Reconstruct mutations
    mutations = []
    i, j = n, m

    while i > 0 or j > 0:
        if backtrack[i][j] is None:
            break

        op = backtrack[i][j]

        if op[0] == 'match':
            i -= 1
            j -= 1
        elif op[0] == 'sub':
            mutations.append(f"{op[2]}{op[1]+1}{op[3]}")
            i -= 1
            j -= 1
        elif op[0] == 'del':
            mutations.append(f"del{op[1]+1}")
            i -= 1
        elif op[0] == 'ins':
            mutations.append(f"ins{i+1}{op[2]}")
            j -= 1

    mutations.reverse()

    return (dp[n][m], mutations)


def apply_mutations(dna_seq: str, mutations: List[str]) -> str:
    """
    Apply a list of mutations to a DNA sequence

    Args:
        dna_seq: Original DNA sequence
        mutations: List of mutation strings (e.g., ['A10T', 'del5', 'ins3G'])

    Returns:
        Modified DNA sequence
    """
    seq_list = list(dna_seq.upper())
    offset = 0  # Track position shifts from indels

    for mut in sorted(mutations, key=lambda x: int(''.join(filter(str.isdigit, x)))):
        if mut.startswith('del'):
            pos = int(mut[3:]) - 1 + offset
            if 0 <= pos < len(seq_list):
                seq_list.pop(pos)
                offset -= 1
        elif mut.startswith('ins'):
            # Format: ins{pos}{base}
            pos_str = ''.join(filter(str.isdigit, mut[3:]))
            base = mut[3 + len(pos_str):]
            pos = int(pos_str) - 1 + offset
            seq_list.insert(pos, base)
            offset += 1
        else:
            # Substitution: {old}{pos}{new}
            pos = int(''.join(filter(str.isdigit, mut))) - 1 + offset
            new_base = mut[-1]
            if 0 <= pos < len(seq_list):
                seq_list[pos] = new_base

    return ''.join(seq_list)


# ============================================================================
# Restriction Enzyme Database - CURATED LIST
# ============================================================================
# Only includes reliable, non-fussy enzymes
# Excludes: nicking enzymes, mega-enzymes, highly ambiguous sites,
#           temperature-sensitive enzymes, double-cutters

RESTRICTION_ENZYMES = [
    # Standard 6-bp cutters (most common, reliable enzymes)
    ("AfeI", "AGCGCT"),
    ("AgeI", "ACCGGT"),
    ("AccI", "GTMKAC"),
    ("AciI", "CCGC"),
    ("AclI", "AACGTT"),
    ("AcuI", "CTGAAG"),
    ("AfeI", "AGCGCT"),
    ("AflII", "CTTAAG"),
    ("AflIII", "ACRYGT"),
    ("AgeI", "ACCGGT"),
    ("AhdI", "GACNNNNNGTC"),
    ("AleI", "CACNNNNGTG"),
    ("AluI", "AGCT"),
    ("AlwI", "GGATC"),
    ("AlwNI", "CAGNNNCTG"),
    ("ApaI", "GGGCCC"),
    ("ApaLI", "GTGCAC"),
    ("ApeKI", "GCWGC"),
    ("ApoI", "RAATTY"),
    ("AscI", "GGCGCGCC"),
    ("AseI", "ATTAAT"),
    ("AsiSI", "GCGATCGC"),
    ("AvaI", "CYCGRG"),
    ("AvaII", "GGWCC"),
    ("AvrII", "CCTAGG"),
    ("BaeGI", "GKGCMC"),
    ("BaeI", "ACNNNNGTAYC"),
    ("BamHI", "GGATCC"),
    ("BanI", "GGYRCC"),
    ("BanII", "GRGCYC"),
    ("BbsI", "GAAGAC"),
    ("BbvCI", "CCTCAGC"),
    ("BbvI", "GCAGC"),
    ("BccI", "CCATC"),
    ("BceAI", "ACGGC"),
    ("BcgI", "CGANNNNNNTGC"),
    ("BciVI", "GTATCC"),
    ("BclI", "TGATCA"),
    ("BcoDI", "GTCTC"),
    ("BfaI", "CTAG"),
    ("BfuAI", "ACCTGC"),
    ("BglI", "GCCNNNNNGGC"),
    ("BglII", "AGATCT"),
    ("BlpI", "GCTNAGC"),
    ("BmgBI", "CACGTC"),
    ("BmrI", "ACTGGG"),
    ("BmtI", "GCTAGC"),
    ("BpmI", "CTGGAG"),
    ("BpuEI", "CTTGAG"),
    ("Bpu10I", "CCTNAGC"),
    ("BsaAI", "YACGTR"),
    ("BsaBI", "GATNNNNATC"),
    ("BsaHI", "GRCGYC"),
    ("BsaI", "GGTCTC"),
    ("BsaJI", "CCNNGG"),
    ("BsaWI", "WCCGGW"),
    ("BsaXI", "ACNNNNNCTCC"),
    ("BseRI", "GAGGAG"),
    ("BseYI", "CCCAGC"),
    ("BsgI", "GTGCAG"),
    ("BsiEI", "CGRYCG"),
    ("BsiHKAI", "GWGCWC"),
    ("BsiWI", "CGTACG"),
    ("BslI", "CCNNNNNNNGG"),
    ("BsmAI", "GTCTC"),
    ("BsmBI", "CGTCTC"),
    ("BsmFI", "GGGAC"),
    ("BsmI", "GAATGC"),
    ("BsoBI", "CYCGRG"),
    ("BspCNI", "CTCAG"),
    ("BspDI", "ATCGAT"),
    ("BspEI", "TCCGGA"),
    ("BspHI", "TCATGA"),
    ("Bsp1286I", "GDGCHC"),
    ("BspMI", "ACCTGC"),
    ("BspQI", "GCTCTTC"),
    ("BsrBI", "CCGCTC"),
    ("BsrDI", "GCAATG"),
    ("BsrFI", "RCCGGY"),
    ("BsrGI", "TGTACA"),
    ("BsrI", "ACTGG"),
    ("BssHII", "GCGCGC"),
    ("BssSI", "CACGAG"),
    ("BstAPI", "GCANNNNNTGC"),
    ("BstBI", "TTCGAA"),
    ("BstEII", "GGTNACC"),
    ("BstNI", "CCWGG"),
    ("BstUI", "CGCG"),
    ("BstXI", "CCANNNNNNTGG"),
    ("BstYI", "RGATCY"),
    ("BstZ17I", "GTATAC"),
    ("Bsu36I", "CCTNAGG"),
    ("BtgI", "CCRYGG"),
    ("BtgZI", "GCGATG"),
    ("BtsCI", "GGATG"),
    ("BtsIMutI", "CAGTG"),
    ("BtsI", "GCAGTG"),
    ("Cac8I", "GCNNGC"),
    ("ClaI", "ATCGAT"),
    ("CspCI", "CAANNNNNGTGG"),
    ("CviKI-1", "RGCY"),
    ("CviQI", "GTAC"),
    ("DdeI", "CTNAG"),
    ("DpnI", "GATC"),
    ("DpnII", "GATC"),
    ("DraI", "TTTAAA"),
    ("DraIII", "CACNNNGTG"),
    ("DrdI", "GACNNNNNNGTC"),
    ("EaeI", "YGGCCR"),
    ("EagI", "CGGCCG"),
    ("EarI", "CTCTTC"),
    ("EciI", "GGCGGA"),
    ("Eco53kI", "GAGCTC"),
    ("EcoNI", "CCTNNNNNAGG"),
    ("EcoO109I", "RGGNCCY"),
    ("EcoP15I", "CAGCAG"),
    ("EcoRI", "GAATTC"),
    ("EcoRV", "GATATC"),
    ("Esp3I", "CGTCTC"),
    ("FatI", "CATG"),
    ("FauI", "CCCGC"),
    ("Fnu4HI", "GCNGC"),
    ("FokI", "GGATG"),
    ("FseI", "GGCCGGCC"),
    ("FspI", "TGCGCA"),
    ("HaeII", "RGCGCY"),
    ("HaeIII", "GGCC"),
    ("HgaI", "GACGC"),
    ("HhaI", "GCGC"),
    ("HincII", "GTYRAC"),
    ("HindIII", "AAGCTT"),
    ("HinfI", "GANTC"),
    ("HinP1I", "GCGC"),
    ("HpaI", "GTTAAC"),
    ("HpaII", "CCGG"),
    ("HphI", "GGTGA"),
    ("HpyAV", "CCTTC"),
    ("HpyCH4III", "ACNGT"),
    ("HpyCH4IV", "ACGT"),
    ("HpyCH4V", "TGCA"),
    ("Hpy99I", "CGWCG"),
    ("Hpy188I", "TCNGA"),
    ("Hpy166II", "GTNNAC"),
    ("Hpy188III", "TCNNGA"),
    ("I-CeuI", "TAACTATAACGGTCCTAAGGTAGCGAA"),
    ("I-SceI", "TAGGGATAACAGGGTAAT"),
    ("KasI", "GGCGCC"),
    ("KpnI", "GGTACC"),
    ("MboI", "GATC"),
    ("MboII", "GAAGA"),
    ("MfeI", "CAATTG"),
    ("MluCI", "AATT"),
    ("MluI", "ACGCGT"),
    ("MlyI", "GAGTC"),
    ("MmeI", "TCCRAC"),
    ("MnlI", "CCTC"),
    ("MscI", "TGGCCA"),
    ("MseI", "TTAA"),
    ("MslI", "CAYNNNNRTG"),
    ("MspA1I", "CMGCKG"),
    ("MspI", "CCGG"),
    ("MspJI", "CNNR"),
    ("MwoI", "GCNNNNNNNGC"),
    ("NaeI", "GCCGGC"),
    ("NarI", "GGCGCC"),
    ("Nb.BbvCI", "CCTCAGC"),
    ("Nb.BsmI", "GAATGC"),
    ("Nb.BsrDI", "GCAATG"),
    ("Nb.BssSI", "CACGAG"),
    ("Nb.BtsI", "GCAGTG"),
    ("NciI", "CCSGG"),
    ("NcoI", "CCATGG"),
    ("NdeI", "CATATG"),
    ("NgoMIV", "GCCGGC"),
    ("NheI", "GCTAGC"),
    ("NlaIII", "CATG"),
    ("NlaIV", "GGNNCC"),
    ("NmeAIII", "GCCGAG"),
    ("NotI", "GCGGCCGC"),
    ("NruI", "TCGCGA"),
    ("NsiI", "ATGCAT"),
    ("NspI", "RCATGY"),
    ("Nt.AlwI", "GGATC"),
    ("Nt.BbvCI", "CCTCAGC"),
    ("Nt.BsmAI", "GTCTC"),
    ("Nt.BspQI", "GCTCTTC"),
    ("Nt.BstNBI", "GAGTC"),
    ("Nt.CviPII", "CCD"),
    ("PacI", "TTAATTAA"),
    ("PaeR7I", "CTCGAG"),
    ("PaqCI", "CACCTGC"),
    ("PciI", "ACATGT"),
    ("PflFI", "GACNNNGTC"),
    ("PflMI", "CCANNNNNTGG"),
    ("PI-PspI", "TGGCAAACAGCTATTATGGGTATTATGGGT"),
    ("PI-SceI", "ATCTATGTCGGGTGCGGAGAAAGAGGTAAT"),
    ("PleI", "GAGTC"),
    ("PluTI", "GGCGCC"),
    ("PmeI", "GTTTAAAC"),
    ("PmlI", "CACGTG"),
    ("PpuMI", "RGGWCCY"),
    ("PshAI", "GACNNNNGTC"),
    ("PsiI", "TTATAA"),
    ("PspGI", "CCWGG"),
    ("PspOMI", "GGGCCC"),
    ("PspXI", "VCTCGAGB"),
    ("PstI", "CTGCAG"),
    ("PvuI", "CGATCG"),
    ("PvuII", "CAGCTG"),
    ("RsaI", "GTAC"),
    ("RsrII", "CGGWCCG"),
    ("SacI", "GAGCTC"),
    ("SacII", "CCGCGG"),
    ("SalI", "GTCGAC"),
    ("SapI", "GCTCTTC"),
    ("Sau3AI", "GATC"),
    ("Sau96I", "GGNCC"),
    ("SbfI", "CCTGCAGG"),
    ("ScaI", "AGTACT"),
    ("ScrFI", "CCNGG"),
    ("SexAI", "ACCWGGT"),
    ("SfaNI", "GCATC"),
    ("SfcI", "CTRYAG"),
    ("SfiI", "GGCCNNNNNGGCC"),
    ("SfoI", "GGCGCC"),
    ("SgrAI", "CRCCGGYG"),
    ("SmaI", "CCCGGG"),
    ("SmlI", "CTYRAG"),
    ("SnaBI", "TACGTA"),
    ("SpeI", "ACTAGT"),
    ("SphI", "GCATGC"),
    ("SrfI", "GCCCGGGC"),
    ("SspI", "AATATT"),
    ("StuI", "AGGCCT"),
    ("StyD4I", "CCNGG"),
    ("StyI", "CCWWGG"),
    ("SwaI", "ATTTAAAT"),
    ("TaqI", "TCGA"),
    ("TfiI", "GAWTC"),
    ("TseI", "GCWGC"),
    ("Tsp45I", "GTSAC"),
    ("TspMI", "CCCGGG"),
    ("TspRI", "NNCASTGNN"),
    ("Tth111I", "GACNNNGTC"),
    ("XbaI", "TCTAGA"),
    ("XcmI", "CCANNNNNNNNNTGG"),
    ("XhoI", "CTCGAG"),
    ("XmaI", "CCCGGG"),
    ("XmnI", "GAANNNNTTC"),
    ("ZraI", "GACGTC"),
]


# ============================================================================
# Data Structures
# ============================================================================

@dataclass
class Candidate:
    """Represents a candidate restriction site"""
    position: int           # 1-indexed position in DNA
    enzyme: str             # Enzyme name
    site_sequence: str      # Recognition sequence
    edits_required: int     # Number of mutations needed
    mutations: List[str]    # List of mutation strings
    mutation_type: str      # Silent, Non-Silent (Missense), etc.
    original_window: str    # Original DNA sequence at this position
    uniqueness_dna: str     # "Unique" or "Non-unique (X sites)" in full DNA
    uniqueness_roi: str     # "Unique" or "Non-unique (X sites)" in ROI (if ROI specified)


# ============================================================================
# Main Algorithm
# ============================================================================

def find_roi_in_dna(dna_seq: str, roi_seq: str) -> Optional[Tuple[int, int]]:
    """
    Find region of interest within DNA sequence

    Args:
        dna_seq: Full DNA sequence
        roi_seq: Region of interest sequence

    Returns:
        Tuple of (start_pos, end_pos) or None if not found
    """
    dna_seq = dna_seq.upper()
    roi_seq = roi_seq.upper()

    pos = dna_seq.find(roi_seq)
    if pos == -1:
        return None

    return (pos, pos + len(roi_seq))


def find_candidates(
    dna_seq: str,
    protein_seq: str,
    max_mutations: int,
    min_length: int,
    roi_seq: Optional[str] = None
) -> List[Candidate]:
    """
    Find all candidate restriction sites within max_mutations distance

    Args:
        dna_seq: DNA sequence
        protein_seq: Protein sequence to maintain
        max_mutations: Maximum number of edits allowed
        min_length: Minimum recognition site length
        roi_seq: Optional region of interest - only search for mutations within this region

    Returns:
        List of Candidate objects
    """
    dna_seq = dna_seq.upper()
    protein_seq = protein_seq.upper()

    # Find protein coding region in DNA
    protein_location = find_protein_in_dna(dna_seq, protein_seq)

    if protein_location is None:
        print(f"ERROR: Protein sequence not found in DNA sequence")
        print(f"Tried all 3 forward reading frames")
        sys.exit(1)

    coding_start, coding_end, frame = protein_location
    print(f"Found protein at positions {coding_start+1}-{coding_end} (frame {frame})")

    # Find ROI if specified
    roi_start = 0
    roi_end = len(dna_seq)

    if roi_seq:
        roi_seq = roi_seq.upper()
        roi_location = find_roi_in_dna(dna_seq, roi_seq)

        if roi_location is None:
            print(f"ERROR: Region of interest not found in DNA sequence")
            print(f"ROI: {roi_seq[:50]}..." if len(roi_seq) > 50 else f"ROI: {roi_seq}")
            sys.exit(1)

        roi_start, roi_end = roi_location
        print(f"Region of interest: positions {roi_start+1}-{roi_end} ({roi_end - roi_start} bp)")
        print(f"Restricting mutation search to ROI only")

    # Filter enzymes by minimum length
    filtered_enzymes = [(name, seq) for name, seq in RESTRICTION_ENZYMES
                        if len(seq) >= min_length]

    print(f"Searching {len(filtered_enzymes)} enzymes (min length: {min_length} bp)")

    candidates = []

    # Sliding window search (restricted to ROI if specified)
    for pos in range(roi_start, roi_end):
        for enzyme_name, enzyme_site in filtered_enzymes:
            site_len = len(enzyme_site)

            # Skip if window would exceed sequence
            if pos + site_len > len(dna_seq):
                continue

            # Extract window from DNA
            window = dna_seq[pos:pos + site_len]

            # Calculate edit distance
            edit_dist, mutations = iupac_edit_distance(window, enzyme_site)

            # Check if within threshold
            if edit_dist <= max_mutations:
                # Apply mutations to create modified sequence
                modified_dna = dna_seq[:pos] + apply_mutations(window, mutations) + dna_seq[pos + site_len:]

                # Check impact on protein
                mutation_type = classify_mutation_impact(
                    original_dna=dna_seq,
                    modified_dna=modified_dna,
                    coding_start=coding_start,
                    coding_end=coding_end,
                    frame=frame,
                    original_protein=protein_seq
                )

                # Count occurrences of this site in modified DNA
                site_count_dna = count_pattern_occurrences(modified_dna, enzyme_site)

                if site_count_dna == 1:
                    uniqueness_dna = "Unique"
                else:
                    uniqueness_dna = f"Non-unique ({site_count_dna} sites)"

                # Count occurrences in ROI if specified
                # Count sites that START within the ROI (they can extend beyond the boundary)
                if roi_seq:
                    site_count_roi = count_pattern_occurrences_in_range(
                        modified_dna, enzyme_site, roi_start, roi_end
                    )

                    if site_count_roi == 1:
                        uniqueness_roi = "Unique"
                    else:
                        uniqueness_roi = f"Non-unique ({site_count_roi} sites)"
                else:
                    uniqueness_roi = "N/A"

                # Create candidate
                candidate = Candidate(
                    position=pos + 1,  # 1-indexed
                    enzyme=enzyme_name,
                    site_sequence=enzyme_site,
                    edits_required=edit_dist,
                    mutations=mutations,
                    mutation_type=mutation_type,
                    original_window=window,
                    uniqueness_dna=uniqueness_dna,
                    uniqueness_roi=uniqueness_roi
                )

                candidates.append(candidate)

    return candidates


def classify_mutation_impact(
    original_dna: str,
    modified_dna: str,
    coding_start: int,
    coding_end: int,
    frame: int,
    original_protein: str
) -> str:
    """
    Classify the impact of mutations on protein sequence

    Args:
        original_dna: Original DNA sequence
        modified_dna: DNA with mutations applied
        coding_start: Start of coding region
        coding_end: End of coding region
        frame: Reading frame
        original_protein: Original protein sequence

    Returns:
        Classification string
    """
    # Extract and translate coding regions
    original_coding = original_dna[coding_start:coding_end]
    modified_coding = modified_dna[coding_start:coding_end]

    # Check for length changes (frameshifts)
    if len(modified_coding) != len(original_coding):
        return "Non-Silent (Frameshift)"

    # Translate both (using the correct frame offset)
    original_translation = translate(original_coding, frame=frame)
    modified_translation = translate(modified_coding, frame=frame)

    # Check for stop codons
    if '*' in modified_translation and '*' not in original_translation:
        return "Non-Silent (Premature Stop)"

    # Compare sequences
    if original_translation == modified_translation:
        return "Silent"
    else:
        # Count differences
        differences = sum(1 for a, b in zip(original_translation, modified_translation) if a != b)
        return f"Non-Silent (Missense, {differences} aa changed)"


# ============================================================================
# Output Formatting
# ============================================================================

def print_results(candidates: List[Candidate], show_all: bool = False, has_roi: bool = False):
    """
    Print candidates in a formatted table

    Args:
        candidates: List of Candidate objects
        show_all: If False, only show silent mutations
        has_roi: If True, show both DNA and ROI uniqueness columns
    """
    if not candidates:
        print("\nNo candidates found within mutation threshold.")
        return

    # Filter for silent if requested
    if not show_all:
        candidates = [c for c in candidates if c.mutation_type == "Silent"]

    if not candidates:
        print("\nNo SILENT candidates found. Use --show-all to see non-silent options.")
        return

    # Sort by edits required, then by position
    candidates.sort(key=lambda c: (c.edits_required, c.position))

    # Column widths
    pos_w = 8
    enzyme_w = 15
    site_w = 20
    edits_w = 6
    mut_w = 22
    type_w = 25

    if has_roi:
        uniq_dna_w = 18
        uniq_roi_w = 18
        total_width = 132
    else:
        uniq_w = 20
        total_width = 120

    # Print header
    print(f"\n{'='*total_width}")
    print(f"Found {len(candidates)} candidate{'s' if len(candidates) != 1 else ''}")
    print(f"{'='*total_width}")

    # Header
    if has_roi:
        header = (
            f"{'Position':<{pos_w}} "
            f"{'Enzyme':<{enzyme_w}} "
            f"{'Site Sequence':<{site_w}} "
            f"{'Edits':<{edits_w}} "
            f"{'Mutations':<{mut_w}} "
            f"{'Type':<{type_w}} "
            f"{'Unique(DNA)':<{uniq_dna_w}} "
            f"{'Unique(ROI)':<{uniq_roi_w}}"
        )
    else:
        header = (
            f"{'Position':<{pos_w}} "
            f"{'Enzyme':<{enzyme_w}} "
            f"{'Site Sequence':<{site_w}} "
            f"{'Edits':<{edits_w}} "
            f"{'Mutations':<{mut_w}} "
            f"{'Type':<{type_w}} "
            f"{'Uniqueness':<{uniq_w}}"
        )

    print(header)
    print("-" * total_width)

    # Rows
    for c in candidates:
        mut_str = ', '.join(c.mutations) if c.mutations else "None (exact match)"
        if len(mut_str) > mut_w - 3:
            mut_str = mut_str[:mut_w-3] + "..."

        if has_roi:
            row = (
                f"{c.position:<{pos_w}} "
                f"{c.enzyme:<{enzyme_w}} "
                f"{c.site_sequence:<{site_w}} "
                f"{c.edits_required:<{edits_w}} "
                f"{mut_str:<{mut_w}} "
                f"{c.mutation_type:<{type_w}} "
                f"{c.uniqueness_dna:<{uniq_dna_w}} "
                f"{c.uniqueness_roi:<{uniq_roi_w}}"
            )
        else:
            row = (
                f"{c.position:<{pos_w}} "
                f"{c.enzyme:<{enzyme_w}} "
                f"{c.site_sequence:<{site_w}} "
                f"{c.edits_required:<{edits_w}} "
                f"{mut_str:<{mut_w}} "
                f"{c.mutation_type:<{type_w}} "
                f"{c.uniqueness_dna:<{uniq_w}}"
            )
        print(row)

    print("=" * total_width)

    # Summary
    silent_count = sum(1 for c in candidates if c.mutation_type == "Silent")
    unique_dna_count = sum(1 for c in candidates if c.uniqueness_dna == "Unique")

    if has_roi:
        unique_roi_count = sum(1 for c in candidates if c.uniqueness_roi == "Unique")
        print(f"\nSummary: {silent_count} silent, {len(candidates) - silent_count} non-silent | DNA: {unique_dna_count} unique | ROI: {unique_roi_count} unique")
    else:
        print(f"\nSummary: {silent_count} silent, {len(candidates) - silent_count} non-silent | {unique_dna_count} unique, {len(candidates) - unique_dna_count} non-unique")


def export_to_csv(candidates: List[Candidate], output_file: str, show_all: bool = False):
    """
    Export candidates to CSV file

    Args:
        candidates: List of Candidate objects
        output_file: Path to output CSV file
        show_all: If False, only export silent mutations
    """
    # Filter for silent if requested
    if not show_all:
        candidates = [c for c in candidates if c.mutation_type == "Silent"]

    if not candidates:
        print(f"No candidates to export to CSV")
        return

    # Sort by edits required, then by position
    candidates.sort(key=lambda c: (c.edits_required, c.position))

    # Write CSV
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = [
            'Position',
            'Enzyme',
            'Site_Sequence',
            'Edits_Required',
            'Mutations',
            'Type',
            'Uniqueness_DNA',
            'Uniqueness_ROI',
            'Original_Window'
        ]

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for c in candidates:
            mut_str = ', '.join(c.mutations) if c.mutations else "None (exact match)"

            writer.writerow({
                'Position': c.position,
                'Enzyme': c.enzyme,
                'Site_Sequence': c.site_sequence,
                'Edits_Required': c.edits_required,
                'Mutations': mut_str,
                'Type': c.mutation_type,
                'Uniqueness_DNA': c.uniqueness_dna,
                'Uniqueness_ROI': c.uniqueness_roi,
                'Original_Window': c.original_window
            })

    print(f"Results exported to: {output_file}")


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Find restriction enzyme sites that can be introduced with minimal silent mutations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python silent_sites.py --dna ATGGCTAGC... --protein MAST... --mutations 2

This tool finds restriction sites that can be added to your DNA sequence
while keeping the protein sequence unchanged (silent mutations).
        """
    )

    parser.add_argument(
        '--dna',
        required=True,
        help='DNA sequence (case insensitive)'
    )

    parser.add_argument(
        '--protein',
        required=True,
        help='Protein sequence (single letter codes, case insensitive)'
    )

    parser.add_argument(
        '--mutations',
        type=int,
        default=2,
        help='Maximum number of mutations allowed (default: 2)'
    )

    parser.add_argument(
        '--min-length',
        type=int,
        default=6,
        help='Minimum restriction site length (default: 6)'
    )

    parser.add_argument(
        '--max-mutations-cap',
        type=int,
        default=3,
        help='Hard cap on mutations to prevent slow execution (default: 3)'
    )

    parser.add_argument(
        '--show-all',
        action='store_true',
        help='Show non-silent mutations as well (default: only silent)'
    )

    parser.add_argument(
        '--roi',
        type=str,
        default=None,
        help='Region of interest - only search for mutations within this DNA sequence (must be found within --dna)'
    )

    args = parser.parse_args()

    # Validate mutations cap
    if args.mutations > args.max_mutations_cap:
        print(f"ERROR: --mutations ({args.mutations}) exceeds safety cap ({args.max_mutations_cap})")
        print(f"       Use --max-mutations-cap to override (warning: may be slow)")
        sys.exit(1)

    # Run analysis
    print("\n" + "="*100)
    print("Silent Restriction Site Finder")
    print("="*100)
    print(f"DNA length: {len(args.dna)} bp")
    print(f"Protein length: {len(args.protein)} aa")
    print(f"Max mutations: {args.mutations}")
    print(f"Min site length: {args.min_length} bp")

    candidates = find_candidates(
        dna_seq=args.dna,
        protein_seq=args.protein,
        max_mutations=args.mutations,
        min_length=args.min_length,
        roi_seq=args.roi
    )

    print_results(candidates, show_all=args.show_all, has_roi=(args.roi is not None))

    # Export to CSV
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_filename = f"silent_sites_results_{timestamp}.csv"
    export_to_csv(candidates, csv_filename, show_all=args.show_all)


if __name__ == "__main__":
    main()
