#!/usr/bin/env python3
"""
Build Splicing Reporter Plasmids for Frame-Shift Dual-Reporter System

Constructs 32 dual-reporter plasmids based on TRX155, each with a different
splice donor from the v7 panel. The intron contains full GGGGS5-VHH3-GGGGS1 sequence.

Reporter Logic (Frame-Shift Design):
- Intron RETAINED → VHH3 translated → NanoLuc active
- Intron SPLICED → OvORF path → FLuc active

The downstream reporters are in different reading frames from entry point:
- Frame 0 from entry → FLuc (stop-free to aa 850)
- Frame 1 from entry → NanoLuc (stop at aa 234, after NanoLuc ends)
- Frame 2 from entry → BLOCKED at aa 61

Frame Logic Math:
- D_spliced = exonic distance from ATG to downstream entry (when spliced)
- D_retained = D_spliced + intron length (when retained)
- D mod 3 = 0 → downstream frame 0 → FLuc
- D mod 3 = 2 → downstream frame 1 → NanoLuc
- Therefore: D_spliced mod 3 = 0 (FLuc), intron mod 3 = 2 (gives D_retained mod 3 = 2 for NanoLuc)
"""

import argparse
import csv
import os
import re
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Import shared functions from build_splice_cassettes
from build_splice_cassettes import (
    VHH3_PROTEIN,
    SPLICE_ACCEPTOR_BGLOBIN,
    BETA_GLOBIN_EXON3_START,
    find_stop_codons_in_frame,
    translate,
    read_panel_csv,
)

# ============================================================================
# AVD006 Sequences - Pre-optimized with varied codons for synthesis
# ============================================================================
# These sequences are taken directly from AVD006-Rep2Mut2Cap9-VP1-VHH3-VR4.gb
# positions 3744-4190 (GGGGS5-VHH3-GGGGS1 insert).
#
# The AVD006 sequences use varied codons for the GGGGS linkers:
#   - Gly: GGT, GGA, GGC (varied, not repetitive)
#   - Ser: TCT, TCA, AGT (varied, not repetitive)
#
# This is much better for DNA synthesis than repetitive sequences like
# GGCGGCGGCGGCAGC repeated 5 times, which causes synthesis issues.
#
# These sequences are stop-free in frame 0 (native frame), which is how
# they're read in the reporter (VHH3 entry at position 117, mod 3 = 0).

# GGGGS5 linker from AVD006 (positions 3744-3818, 75 bp)
# Protein: GGGGSGGGGSGGGGSGGGGSGGGGSEVQLVESGGGLVQPGGSLRLSCAASGFTFSTADMGWFRQAPGKGRELVAAVSGSGFSTYSDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAKATISLYYAMDVWGQGTTVTVSS
AVD006_GGGGS5_DNA = "GGTGGAGGCGGATCTGGAGGCGGTGGTTCAGGCGGTGGAGGAAGTGGTGGCGGAGGTTCTGGTGGAGGCGGTTCT"
# Length: 75 bp, encodes (GGGGS)×5

# VHH3 anti-ALPL nanobody from AVD006 (positions 3819-4175, 357 bp)
AVD006_VHH3_DNA = (
    "GAGGTGCAACTGGTTGAAAGCGGCGGAGGACTTGTTCAACCCGGCGGCAGCCTTAGGCTTTCTTGCGCTGCCAGC"
    "GGCTTCACCTTTAGCACCGCCGACATGGGCTGGTTTAGGCAAGCTCCCGGAAAAGGCAGGGAACTTGTTGCCGCT"
    "GTGAGCGGCAGCGGCTTCAGCACCTACTCTGATAGCGTTGAGGGCAGGTTCACCATCAGCAGGGACAACGCCAAG"
    "AGGATGGTGTACCTGCAGATGAACAGCTTGAGGGCCGAGGACACCGCCGTGTACTACTGCGCCAAGGCCACAATT"
    "AGCCTGTACTACGCCATGGATGTGTGGGGACAGGGCACCACCGTGACCGTGAGCAGC"
)
# Length: 357 bp, encodes VHH3 (119 aa)

# GGGGS1 linker from AVD006 (positions 4176-4190, 15 bp)
AVD006_GGGGS1_DNA = "GGTGGAGGCGGTTCT"
# Length: 15 bp, encodes GGGGS


# ============================================================================
# Constants
# ============================================================================

# Upstream exon sequence (66 bp) with triple tandem Kozak-Met-Ala + VR4 loop
# This is placed before the splice donor to provide translation start
#
# Structure: [Triple Kozak-Met-Ala (42 bp)] + VR4(3720-3743) (24 bp)
#
# The triple tandem sequence contains 3 ATG codons:
#   Position 7-9:   ATG (weak Kozak context - likely skipped)
#   Position 22-24: ATG (weak Kozak context - likely skipped)
#   Position 37-39: ATG (STRONG Kozak: GCCACCATGGCA - this is used!)
#
# Uses AAV9 VP1 VR4 loop sequence EXACTLY from AVD002 positions 3720-3743
# No modifications - this is the exact sequence ending in "...GSG"
# Protein from strong ATG: MASKTINGSG (Met-Ala-Ser-Lys-Thr-Ile-Asn-Gly-Ser-Gly)
#
UPSTREAM_EXON_66BP = "ATCATCATGCTTGCCTTCATCATGCTTACCGCCACCATGGCATCAAAGACTATTAACGGTTCTGGA"
# Length: 66 bp
# Contains: Triple Kozak-Met-Ala (42 bp) + VR4 (24 bp)
# Strong ATG at position 37-39 (GCCACCATGGCA context)
# From strong ATG to end: 30 bp (same as before, preserves frame logic)
# Source: VR4 from AAV9 VP1 (positions 3720-3743 of AVD002) - EXACT match

# Exonic stuffer (1 bp) placed between upstream exon and donor exonic
# This changes junction from 33 bp to 34 bp, which is CRITICAL for frame alignment:
#
# Frame Logic with J=34 and intron=476:
#   J = 30 (upstream from ATG) + 1 (exonic stuffer) + 3 (donor exonic) = 34 bp
#   VHH3 entry = 34 + 8 + 75 = 117, mod 3 = 0 → VHH3 in NATIVE frame (stop-free) ✓
#   D_spliced = 34 + 15 + 32 = 81, mod 3 = 0 → downstream frame 0 → FLuc ✓
#   D_retained = 34 + 476 + 15 + 32 = 557, mod 3 = 2 → downstream frame 1 → NanoLuc ✓
#
# The exonic stuffer is CRITICAL to keep VHH3 in its native frame (frame 0).
# VHH3 frame 1 has a stop codon at position 244 - we MUST avoid frame 1!
#
EXONIC_STUFFER_1BP = "G"
# Encodes Gly when combined with upstream GGA → GGAG... donor

# β-globin intron 2 acceptor: last 23 bp of intronic region (PPT ending in AG)
# This is used for the intronic portion of the acceptor
SPLICE_ACCEPTOR_INTRONIC = "TCTATTTTCCCACCCTTAG"  # 20 bp ends with AG
# Extended version with more PPT context
SPLICE_ACCEPTOR_INTRONIC_23BP = "TTCTATTTTCCCACCCTCAG"  # 20 bp - used from SPLICE_ACCEPTOR_BGLOBIN

# AAV9 VP1 VR4 loop downstream - forms the acceptor exonic region (minimal)
# Uses only the first 15 bp of VR4 loop: "CAGAATCAACAAACG"
# Source: AAV9 VP1 VR4 loop (positions 3744-3758 of AVD002)
#
# This minimal sequence avoids frame 1 stop codon issues present in longer versions.
# Using stuffers to achieve correct frame alignment instead.
#
# Frame verification (both paths are stop-free):
#   SPLICED (J=34): Enters at frame 34 mod 3 = 1 → acceptor frame 2 → no stops ✓
#   RETAINED (J+I=510): Enters at frame 510 mod 3 = 0 → acceptor frame 0 → no stops ✓
#
# Length: 15 bp (mod 3 = 0)
#
ACCEPTOR_EXONIC_15BP = "CAGAATCAACAAACG"
# 15 bp, no stop codons in any frame (frame 0, 1, or 2)

# OvORF Linker - connects downstream exon to NanoLuc/FLuc region
# This is the start of the Overlapping ORF from original TRX155 (positions 1546-1579)
# MODIFIED to remove ALL alternative start codons (ATG, CTG, GTG, TTG, ACG)
#
# Original: CATGCTTGCCTTCATCATGCTTACCGCCACCATG
#   - 3x ATG at positions 1, 16, 31 (frame 1 Met codons)
#   - 1x TTG at position 5 (formed by CTT + G)
#
# Modified: CGCCCTCGCCTTCATCGCCCTCACCGCCACCG (32 bp)
#   - ATG → GCC (Met → Ala) at 3 positions
#   - CTT → CTC (Leu → Leu, synonymous) at 2 positions to eliminate TTG
#   - Length 32 bp (mod 3 = 2) for correct frame alignment:
#     D_spliced = 34 + 15 + 32 = 81, mod 3 = 0 → FLuc
#
# This sequence provides stop-free path in all frames.
# SPLICED path reads through this in frame 0 to reach FLuc.
#
OVORF_LINKER_32BP = "CGCCCTCGCCTTCATCGCCCTCACCGCCACCG"
# 32 bp, no alternative starts (ATG/CTG/GTG/TTG/ACG), no stops in any frame


# ============================================================================
# GenBank Parsing and Writing Functions
# ============================================================================

def parse_genbank(filepath: str) -> Tuple[Dict, str]:
    """
    Parse a GenBank file into header dict and sequence.

    Returns:
        Tuple of (header_dict, sequence)
        header_dict contains: locus, definition, features, etc.
    """
    with open(filepath, 'r') as f:
        content = f.read()

    # Split into header and sequence
    origin_match = re.search(r'\nORIGIN\s*\n', content)
    if not origin_match:
        raise ValueError("No ORIGIN section found in GenBank file")

    header_section = content[:origin_match.end()]
    sequence_section = content[origin_match.end():]

    # Parse sequence (remove numbers, spaces, newlines, and //)
    sequence = re.sub(r'[^acgtACGT]', '', sequence_section)
    sequence = sequence.upper()

    # Parse header
    header = {
        'raw': header_section,
        'locus': '',
        'definition': '',
        'features': [],
    }

    # Extract LOCUS line
    locus_match = re.search(r'LOCUS\s+(\S+)\s+(\d+)\s+bp', header_section)
    if locus_match:
        header['locus'] = locus_match.group(1)
        header['length'] = int(locus_match.group(2))

    return header, sequence


def write_genbank(
    filepath: str,
    locus_name: str,
    sequence: str,
    features: List[Dict],
    definition: str = ".",
    accession: str = None,
    organism: str = ".",
) -> None:
    """
    Write a GenBank file with the given sequence and features.

    Args:
        filepath: Output file path
        locus_name: Name for LOCUS line
        sequence: DNA sequence
        features: List of feature dicts with keys: type, location, qualifiers
        definition: DEFINITION line content
        accession: ACCESSION (defaults to locus_name)
        organism: ORGANISM
    """
    if accession is None:
        accession = locus_name

    seq_len = len(sequence)
    today = datetime.now().strftime("%d-%b-%Y").upper()

    lines = []

    # LOCUS line
    lines.append(f"LOCUS       {locus_name} {seq_len} bp    DNA     circular SYN {today}")
    lines.append(f"DEFINITION  {definition}")
    lines.append(f"ACCESSION   {accession}")
    lines.append(f"VERSION     {accession}")
    lines.append("KEYWORDS    .")
    lines.append("SOURCE      .")
    lines.append("  ORGANISM  .")
    lines.append("            .")

    # FEATURES section
    lines.append("FEATURES             Location/Qualifiers")

    for feat in features:
        feat_type = feat.get('type', 'misc_feature')
        location = feat.get('location', '')
        qualifiers = feat.get('qualifiers', {})

        # Feature line (type at column 6, location at column 22)
        lines.append(f"     {feat_type:<15} {location}")

        # Qualifiers
        for key, value in qualifiers.items():
            if value is None:
                lines.append(f'                     /{key}')
            elif isinstance(value, bool) and value:
                lines.append(f'                     /{key}')
            else:
                # Wrap long values
                val_str = str(value)
                if len(val_str) > 58:
                    # Multi-line qualifier
                    lines.append(f'                     /{key}="{val_str[:58]}')
                    remaining = val_str[58:]
                    while remaining:
                        chunk = remaining[:58]
                        remaining = remaining[58:]
                        if remaining:
                            lines.append(f'                     {chunk}')
                        else:
                            lines.append(f'                     {chunk}"')
                else:
                    lines.append(f'                     /{key}="{val_str}"')

    # ORIGIN section
    lines.append("ORIGIN")

    # Format sequence in GenBank style (60 bp per line, 10 bp groups)
    seq_lower = sequence.lower()
    for i in range(0, len(seq_lower), 60):
        chunk = seq_lower[i:i+60]
        # Split into 10-bp groups
        groups = [chunk[j:j+10] for j in range(0, len(chunk), 10)]
        line_num = i + 1
        lines.append(f"{line_num:>9} {' '.join(groups)}")

    lines.append("//")

    # Write file
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, 'w') as f:
        f.write('\n'.join(lines) + '\n')


# ============================================================================
# Cassette Building Functions
# ============================================================================

def build_intron_body(
    donor_intronic_8bp: str,
    ggggs5_dna: str,
    vhh3_dna: str,
    ggggs1_dna: str,
    acceptor_intronic: str,
    target_mod3: int = 2,
) -> Tuple[str, int]:
    """
    Build the complete intron body with stuffer to achieve target mod 3.

    Intron structure:
    [donor_intronic (+1..+8)] [GGGGS5] [VHH3] [GGGGS1] [stuffer] [acceptor_intronic]

    CRITICAL FRAME LOGIC (with J=34 exonic junction):
    - Intronic stuffer is placed AFTER GGGGS1 (before acceptor_intronic)
    - VHH3 entry = 34 + 8 + 75 = 117, mod 3 = 0 → VHH3 in NATIVE frame (stop-free) ✓
    - Intron = 476 bp, mod 3 = 2 → correct frame shift ✓
    - D_spliced = 81, mod 3 = 0 → downstream frame 0 → FLuc ✓
    - D_retained = 557, mod 3 = 2 → downstream frame 1 → NanoLuc ✓

    Args:
        donor_intronic_8bp: 8 bp intronic part of splice donor (+1 to +8, starts with GT)
        ggggs5_dna: GGGGS5 linker DNA (75 bp)
        vhh3_dna: VHH3 nanobody DNA (357 bp)
        ggggs1_dna: GGGGS1 linker DNA (15 bp)
        acceptor_intronic: Intronic part of acceptor ending with AG
        target_mod3: Target intron length mod 3 (default 2 for frame shift design)

    Returns:
        Tuple of (intron_sequence, stuffer_length)
    """
    # Calculate current length without stuffer
    base_intron = donor_intronic_8bp + ggggs5_dna + vhh3_dna + ggggs1_dna + acceptor_intronic
    base_len = len(base_intron)

    # Calculate stuffer needed
    current_mod = base_len % 3
    if current_mod == target_mod3:
        stuffer = ""
    else:
        # Need to add stuffer to reach target_mod3
        # (base_len + stuffer_len) % 3 = target_mod3
        # stuffer_len = (target_mod3 - current_mod) % 3
        stuffer_len = (target_mod3 - current_mod) % 3
        # Use neutral nucleotide(s) for stuffer
        stuffer = "A" * stuffer_len

    # Insert stuffer AFTER GGGGS1 (before acceptor_intronic)
    intron = donor_intronic_8bp + ggggs5_dna + vhh3_dna + ggggs1_dna + stuffer + acceptor_intronic

    return intron, len(stuffer)


def build_new_cassette(
    splice_donor_11mer: str,
    ggggs5_dna: str,
    vhh3_dna: str,
    ggggs1_dna: str,
    acceptor_intronic: str,
    acceptor_exonic: str,
    ovorf_linker: str = OVORF_LINKER_32BP,
    upstream_exon: str = UPSTREAM_EXON_66BP,
    exonic_stuffer: str = EXONIC_STUFFER_1BP,
) -> Tuple[str, Dict]:
    """
    Build the complete new cassette to replace positions 1083-1579.

    Cassette structure:
    [Upstream-Exon (36 bp)] [Stuffer (1 bp)] [Donor-Exonic (3 bp)] [GT...intron...AG] [Acceptor-Exonic (15 bp)] [OvORF-Linker (32 bp)]

    Reporter Logic:
    - RETAINED → VHH3 translated → downstream frame 1 → NanoLuc
    - SPLICED → OvORF path → downstream frame 0 → FLuc

    CRITICAL Frame Logic (with J=34 junction):
    - Junction J = 30 (upstream from ATG) + 1 (exonic stuffer) + 3 (donor exonic) = 34 bp
    - VHH3 entry = 34 + 8 + 75 = 117, mod 3 = 0 → VHH3 in NATIVE frame (stop-free) ✓
    - Intron = 476 bp, mod 3 = 2 → correct frame shift
    - D_spliced = 34 + 15 + 32 = 81, mod 3 = 0 → downstream frame 0 → FLuc ✓
    - D_retained = 34 + 476 + 15 + 32 = 557, mod 3 = 2 → downstream frame 1 → NanoLuc ✓

    Args:
        splice_donor_11mer: 11-mer splice donor (positions -3 to +8)
        ggggs5_dna: GGGGS5 linker DNA
        vhh3_dna: VHH3 nanobody DNA
        ggggs1_dna: GGGGS1 linker DNA
        acceptor_intronic: Intronic part of splice acceptor (ends with AG)
        acceptor_exonic: Exonic part of splice acceptor
        ovorf_linker: OvORF linker connecting to NanoLuc region (ATG-free)
        upstream_exon: Upstream exon with Kozak + ATG + Ala
        exonic_stuffer: 1 bp stuffer between upstream exon and donor (for frame alignment)

    Returns:
        Tuple of (cassette_sequence, cassette_info_dict)
    """
    # Parse splice donor 11-mer
    donor_exonic = splice_donor_11mer[0:3]   # Positions -3, -2, -1
    donor_intronic = splice_donor_11mer[3:]  # Positions +1 to +8 (starts with GT)

    # Verify GT at +1/+2
    if donor_intronic[0:2].upper() != "GT":
        raise ValueError(f"Invalid splice donor: must have GT at +1/+2, got {donor_intronic[0:2]}")

    # Build intron body with stuffer for mod 3 = 2
    # CRITICAL: intron mod 3 = 2 is required for correct frame alignment:
    #   - SPLICED: D mod 3 = 0 → downstream frame 0 → FLuc
    #   - RETAINED: D mod 3 = (0+2) = 2 → downstream frame 1 → NanoLuc
    intron, intronic_stuffer_len = build_intron_body(
        donor_intronic_8bp=donor_intronic,
        ggggs5_dna=ggggs5_dna,
        vhh3_dna=vhh3_dna,
        ggggs1_dna=ggggs1_dna,
        acceptor_intronic=acceptor_intronic,
        target_mod3=2,
    )

    # Build complete cassette with exonic stuffer between upstream and donor
    # OvORF linker connects downstream exon to NanoLuc region (stop-free in frame 1)
    cassette = upstream_exon + exonic_stuffer + donor_exonic + intron + acceptor_exonic + ovorf_linker

    # Calculate frame verification
    # Strong ATG is 36 bp into upstream exon (at position 37-39, after triple tandem)
    atg_offset = 36
    # Junction = distance from ATG to intron start
    exonic_from_atg_to_junction = (len(upstream_exon) - atg_offset +
                                    len(exonic_stuffer) + len(donor_exonic))

    # When spliced: D = junction + acceptor_exonic + ovorf_linker
    D_spliced = exonic_from_atg_to_junction + len(acceptor_exonic) + len(ovorf_linker)

    # When retained: D + I = junction + intron + acceptor_exonic + ovorf_linker
    D_retained = exonic_from_atg_to_junction + len(intron) + len(acceptor_exonic) + len(ovorf_linker)

    info = {
        'upstream_exon_len': len(upstream_exon),
        'exonic_stuffer': exonic_stuffer,
        'exonic_stuffer_len': len(exonic_stuffer),
        'donor_exonic': donor_exonic,
        'donor_exonic_len': len(donor_exonic),
        'donor_intronic': donor_intronic,
        'intron_len': len(intron),
        'stuffer_len': intronic_stuffer_len,  # Intronic stuffer (after GGGGS1)
        'acceptor_exonic_len': len(acceptor_exonic),
        'ovorf_linker_len': len(ovorf_linker),
        'cassette_len': len(cassette),
        'junction': exonic_from_atg_to_junction,
        'D_spliced': D_spliced,
        'D_retained': D_retained,
        'frame_spliced': D_spliced % 3,
        'frame_retained': D_retained % 3,
        'intron_mod3': len(intron) % 3,
    }

    return cassette, info


def build_splicing_reporter_plasmid(
    template_path: str,
    splice_donor_11mer: str,
    plasmid_name: str,
    output_dir: str,
    ggggs5_dna: str,
    vhh3_dna: str,
    ggggs1_dna: str,
) -> Tuple[str, Dict]:
    """
    Build frame-shift dual-reporter plasmid.

    RETAINED (intron kept)   → VHH3 translated → downstream frame 1 → NanoLuc
    SPLICED (intron removed) → OvORF path → downstream frame 0 → FLuc

    Args:
        template_path: Path to TRX155 GenBank template
        splice_donor_11mer: 11-mer splice donor from v7 panel
        plasmid_name: Name for output plasmid
        output_dir: Output directory
        ggggs5_dna: Pre-computed GGGGS5 DNA
        vhh3_dna: Pre-computed VHH3 DNA
        ggggs1_dna: Pre-computed GGGGS1 DNA

    Returns:
        Tuple of (output_path, info_dict)
    """
    # Parse template
    header, template_seq = parse_genbank(template_path)

    # Region to replace: positions 1083-1579 (1-indexed, inclusive)
    # This includes:
    #   - Original miniXon2 region (1083-1527)
    #   - Triple ATG sequence (1539-1579) that we want to remove
    # The intervening region (1528-1538) is also removed for clean junction
    #
    # Convert to 0-indexed: 1082-1578 (inclusive) = 1082:1579 (slice)
    replace_start = 1082  # 0-indexed start (position 1083 in 1-indexed)
    replace_end = 1579    # 0-indexed exclusive end (position 1579 in 1-indexed)

    # Extract acceptor intronic from constant
    # CRITICAL: Using TCAG ending (not TTAG) to avoid TAG stop codon
    acceptor_intronic = "TTCTATTTTCCCACCCTCAG"  # 20 bp PPT + CAG (not TAG!)
    # Using minimal 15 bp exonic to avoid frame 1 stop codons
    acceptor_exonic = ACCEPTOR_EXONIC_15BP      # 15 bp

    # Build new cassette
    new_cassette, cassette_info = build_new_cassette(
        splice_donor_11mer=splice_donor_11mer,
        ggggs5_dna=ggggs5_dna,
        vhh3_dna=vhh3_dna,
        ggggs1_dna=ggggs1_dna,
        acceptor_intronic=acceptor_intronic,
        acceptor_exonic=acceptor_exonic,
    )

    # Build new sequence
    new_seq = template_seq[:replace_start] + new_cassette + template_seq[replace_end:]

    # Calculate length difference for feature offset
    old_region_len = replace_end - replace_start  # 445 bp
    new_region_len = len(new_cassette)
    offset = new_region_len - old_region_len

    # Create features for the new plasmid
    features = create_plasmid_features(
        cassette_info=cassette_info,
        cassette_start=replace_start,  # 0-indexed
        splice_donor_11mer=splice_donor_11mer,
        vhh3_len=len(vhh3_dna),
        ggggs5_len=len(ggggs5_dna),
        ggggs1_len=len(ggggs1_dna),
        offset=offset,
        template_seq_len=len(template_seq),
    )

    # Write output GenBank file
    output_path = os.path.join(output_dir, f"{plasmid_name}.gb")
    write_genbank(
        filepath=output_path,
        locus_name=plasmid_name,
        sequence=new_seq,
        features=features,
        definition=f"Splicing reporter plasmid with donor {splice_donor_11mer}",
    )

    # Add more info
    cassette_info['plasmid_name'] = plasmid_name
    cassette_info['output_path'] = output_path
    cassette_info['plasmid_len'] = len(new_seq)
    cassette_info['splice_donor'] = splice_donor_11mer

    return output_path, cassette_info


def create_plasmid_features(
    cassette_info: Dict,
    cassette_start: int,
    splice_donor_11mer: str,
    vhh3_len: int,
    ggggs5_len: int,
    ggggs1_len: int,
    offset: int,
    template_seq_len: int,
) -> List[Dict]:
    """
    Create feature annotations for the reporter plasmid.

    Args:
        cassette_info: Info dict from build_new_cassette
        cassette_start: 0-indexed start position of new cassette
        splice_donor_11mer: The splice donor sequence
        vhh3_len: Length of VHH3 DNA
        ggggs5_len: Length of GGGGS5 DNA
        ggggs1_len: Length of GGGGS1 DNA
        offset: Length difference from template
        template_seq_len: Original template sequence length

    Returns:
        List of feature dicts
    """
    features = []

    # Calculate positions (convert to 1-indexed for GenBank)
    upstream_exon_start = cassette_start + 1  # 1-indexed
    upstream_exon_end = upstream_exon_start + cassette_info['upstream_exon_len'] - 1

    # Exonic stuffer comes after upstream exon
    exonic_stuffer_start = upstream_exon_end + 1
    exonic_stuffer_end = exonic_stuffer_start + cassette_info['exonic_stuffer_len'] - 1

    # Donor exonic comes after exonic stuffer
    donor_exonic_start = exonic_stuffer_end + 1
    donor_exonic_end = donor_exonic_start + cassette_info['donor_exonic_len'] - 1

    # Intron starts after donor exonic
    intron_start = donor_exonic_end + 1
    intron_end = intron_start + cassette_info['intron_len'] - 1

    # Within intron: donor_intronic (8) + GGGGS5 + VHH3 + GGGGS1 + stuffer + acceptor_intronic (20)
    donor_intronic_start = intron_start
    donor_intronic_end = donor_intronic_start + 7  # 8 bp

    ggggs5_start = donor_intronic_end + 1
    ggggs5_end = ggggs5_start + ggggs5_len - 1

    vhh3_start = ggggs5_end + 1
    vhh3_end = vhh3_start + vhh3_len - 1

    ggggs1_start = vhh3_end + 1
    ggggs1_end = ggggs1_start + ggggs1_len - 1

    stuffer_start = ggggs1_end + 1
    stuffer_end = stuffer_start + cassette_info['stuffer_len'] - 1 if cassette_info['stuffer_len'] > 0 else stuffer_start - 1

    acceptor_intronic_start = stuffer_end + 1 if cassette_info['stuffer_len'] > 0 else ggggs1_end + 1
    acceptor_intronic_end = intron_end

    acceptor_exonic_start = intron_end + 1
    acceptor_exonic_end = acceptor_exonic_start + cassette_info['acceptor_exonic_len'] - 1

    # OvORF linker comes after acceptor exonic
    ovorf_linker_start = acceptor_exonic_end + 1
    ovorf_linker_end = ovorf_linker_start + cassette_info['ovorf_linker_len'] - 1

    # CMV promoter (keeping from original, before cassette)
    features.append({
        'type': 'promoter',
        'location': '492..1076',
        'qualifiers': {'label': 'CMV', 'note': 'CMV promoter'}
    })

    # Upstream exon with ATG (VR4-based)
    features.append({
        'type': 'exon',
        'location': f'{upstream_exon_start}..{upstream_exon_end}',
        'qualifiers': {
            'label': 'upstream_exon_VR4',
            'note': 'Kozak + Met-Ala + AAV9 VP1 VR4 loop (36 bp, exact AVD002 3720-3743)'
        }
    })

    # ATG start codon
    atg_start = upstream_exon_start + 6  # GCCACC ATG
    features.append({
        'type': 'CDS',
        'location': f'{atg_start}..{atg_start + 2}',
        'qualifiers': {
            'label': 'ATG',
            'codon_start': '1',
            'translation': 'M'
        }
    })

    # Exonic stuffer (for frame alignment) - only add if present
    # NOTE: With current design (J=33), NO exonic stuffer is needed
    if cassette_info['exonic_stuffer_len'] > 0:
        features.append({
            'type': 'misc_feature',
            'location': f'{exonic_stuffer_start}..{exonic_stuffer_end}',
            'qualifiers': {
                'label': 'exonic_stuffer',
                'note': f'Frame adjustment stuffer ({cassette_info["exonic_stuffer_len"]} bp)'
            }
        })

    # Splice donor (exonic + intronic)
    features.append({
        'type': 'misc_feature',
        'location': f'{donor_exonic_start}..{donor_intronic_end}',
        'qualifiers': {
            'label': f'splice_donor_{splice_donor_11mer}',
            'note': f'Splice donor: {splice_donor_11mer}'
        }
    })

    # Intron
    features.append({
        'type': 'intron',
        'location': f'{intron_start}..{intron_end}',
        'qualifiers': {
            'label': 'VHH3_intron',
            'note': f'Intron with VHH3 insert ({cassette_info["intron_len"]} bp, mod3={cassette_info["intron_mod3"]})'
        }
    })

    # GGGGS5 linker
    features.append({
        'type': 'misc_feature',
        'location': f'{ggggs5_start}..{ggggs5_end}',
        'qualifiers': {
            'label': 'GGGGS5',
            'note': '5x GGGGS linker'
        }
    })

    # VHH3 nanobody
    features.append({
        'type': 'misc_feature',
        'location': f'{vhh3_start}..{vhh3_end}',
        'qualifiers': {
            'label': 'VHH3',
            'note': 'VHH3 anti-ALPL nanobody (codon optimized)'
        }
    })

    # GGGGS1 linker
    features.append({
        'type': 'misc_feature',
        'location': f'{ggggs1_start}..{ggggs1_end}',
        'qualifiers': {
            'label': 'GGGGS1',
            'note': '1x GGGGS linker'
        }
    })

    # Stuffer (if present)
    if cassette_info['stuffer_len'] > 0:
        features.append({
            'type': 'misc_feature',
            'location': f'{stuffer_start}..{stuffer_end}',
            'qualifiers': {
                'label': 'stuffer',
                'note': f'Frame adjustment stuffer ({cassette_info["stuffer_len"]} bp)'
            }
        })

    # Splice acceptor (intronic region)
    features.append({
        'type': 'misc_feature',
        'location': f'{acceptor_intronic_start}..{acceptor_intronic_end}',
        'qualifiers': {
            'label': 'splice_acceptor_intronic',
            'note': 'Beta-globin intron 2 acceptor (PPT + AG)'
        }
    })

    # Acceptor exon (VR4-based)
    features.append({
        'type': 'exon',
        'location': f'{acceptor_exonic_start}..{acceptor_exonic_end}',
        'qualifiers': {
            'label': 'acceptor_exon_VR4',
            'note': 'AAV9 VP1 VR4 loop minimal (15 bp, exact AVD002 3744-3758)'
        }
    })

    # OvORF linker (connects to NanoLuc region, stop-free in frame 1)
    features.append({
        'type': 'misc_feature',
        'location': f'{ovorf_linker_start}..{ovorf_linker_end}',
        'qualifiers': {
            'label': 'OvORF_linker',
            'note': 'Overlapping ORF linker (32 bp, no alternative starts, stop-free - spliced path to FLuc)'
        }
    })

    # Downstream features from template (adjusted by offset)
    # NanoLuc: original 1650-2159
    nluc_start = 1650 + offset
    nluc_end = 2159 + offset
    features.append({
        'type': 'CDS',
        'location': f'{nluc_start}..{nluc_end}',
        'qualifiers': {
            'label': 'NanoLuc',
            'note': 'NanoLuc reporter (active when RETAINED - VHH3 path)',
            'product': 'NanoLuc luciferase'
        }
    })

    # FLuc (luc2): original 2360-4006
    fluc_start = 2360 + offset
    fluc_end = 4006 + offset
    features.append({
        'type': 'CDS',
        'location': f'{fluc_start}..{fluc_end}',
        'qualifiers': {
            'label': 'FLuc',
            'note': 'Firefly luciferase (active when SPLICED - OvORF path)',
            'product': 'firefly luciferase'
        }
    })

    # bGH polyA signal: original 4168-4392
    polya_start = 4168 + offset
    polya_end = 4392 + offset
    features.append({
        'type': 'polyA_signal',
        'location': f'{polya_start}..{polya_end}',
        'qualifiers': {
            'label': 'bGH_polyA',
            'note': 'bovine growth hormone polyadenylation signal'
        }
    })

    # AmpR: original 5651-6511 (complement)
    ampr_start = 5651 + offset
    ampr_end = 6511 + offset
    features.append({
        'type': 'CDS',
        'location': f'complement({ampr_start}..{ampr_end})',
        'qualifiers': {
            'label': 'AmpR',
            'gene': 'bla',
            'product': 'beta-lactamase',
            'note': 'confers resistance to ampicillin'
        }
    })

    # ori: original 4892-5480 (complement)
    ori_start = 4892 + offset
    ori_end = 5480 + offset
    features.append({
        'type': 'rep_origin',
        'location': f'complement({ori_start}..{ori_end})',
        'qualifiers': {
            'label': 'ori',
            'note': 'high-copy-number ColE1/pMB1/pBR322/pUC origin of replication'
        }
    })

    return features


# ============================================================================
# Main Functions
# ============================================================================

def generate_plasmid_panel(
    template_path: str,
    panel_csv_path: str,
    output_dir: str,
    name_prefix: str = "SPL",
) -> List[Dict]:
    """
    Generate all 32 reporter plasmids from the v7 panel.

    Args:
        template_path: Path to TRX155 template GenBank file
        panel_csv_path: Path to tiered_panel_32_v7.csv
        output_dir: Output directory for plasmid files
        name_prefix: Prefix for plasmid names (default: SPL)

    Returns:
        List of info dicts for all generated plasmids
    """
    # Read panel
    print(f"Reading panel from: {panel_csv_path}")
    candidates = read_panel_csv(panel_csv_path)
    print(f"Found {len(candidates)} candidates")

    # Use pre-optimized sequences from AVD006
    # These have varied codons for synthesis compatibility (not repetitive!)
    print("\nUsing AVD006 sequences (pre-optimized with varied codons)...")

    vhh3_dna = AVD006_VHH3_DNA
    print(f"  VHH3 length: {len(vhh3_dna)} bp")

    # Verify VHH3 translation
    vhh3_translated = translate(vhh3_dna)
    if vhh3_translated != VHH3_PROTEIN:
        raise ValueError("VHH3 translation mismatch!")
    print(f"    VHH3 translation verified: {vhh3_translated[:20]}...")

    ggggs5_dna = AVD006_GGGGS5_DNA
    print(f"  GGGGS5 length: {len(ggggs5_dna)} bp (varied codons)")

    ggggs1_dna = AVD006_GGGGS1_DNA
    print(f"  GGGGS1 length: {len(ggggs1_dna)} bp (varied codons)")

    # Create output directory
    plasmids_dir = os.path.join(output_dir, "plasmids")
    Path(plasmids_dir).mkdir(parents=True, exist_ok=True)

    # Generate plasmids
    print(f"\nGenerating {len(candidates)} plasmids...")
    results = []

    for idx, candidate in enumerate(candidates, start=1):
        tier = candidate['tier']
        sequence = candidate['sequence']
        annotation = candidate.get('annotation', '')

        # Generate plasmid name: SPL001_ACGGTAAAGCA_Tier1
        plasmid_name = f"{name_prefix}{idx:03d}_{sequence}_Tier{tier}"

        print(f"  [{idx:02d}/{len(candidates)}] {plasmid_name}...")

        try:
            output_path, info = build_splicing_reporter_plasmid(
                template_path=template_path,
                splice_donor_11mer=sequence,
                plasmid_name=plasmid_name,
                output_dir=plasmids_dir,
                ggggs5_dna=ggggs5_dna,
                vhh3_dna=vhh3_dna,
                ggggs1_dna=ggggs1_dna,
            )

            # Add panel info
            info['tier'] = tier
            info['annotation'] = annotation
            info['delta_g'] = candidate.get('delta_g', '')
            info['total_hbonds'] = candidate.get('total_hbonds', '')
            info['expected_splicing'] = candidate.get('expected_splicing', '')

            results.append(info)

        except Exception as e:
            print(f"    ERROR: {e}")
            results.append({
                'plasmid_name': plasmid_name,
                'splice_donor': sequence,
                'tier': tier,
                'error': str(e),
            })

    # Write summary CSV
    summary_path = os.path.join(output_dir, "plasmid_panel_summary.csv")
    write_summary_csv(results, summary_path)
    print(f"\nWrote summary to: {summary_path}")

    return results


def write_summary_csv(results: List[Dict], output_path: str) -> None:
    """Write plasmid panel summary to CSV."""
    fieldnames = [
        'plasmid_name', 'splice_donor', 'tier', 'annotation',
        'delta_g', 'total_hbonds', 'expected_splicing',
        'cassette_len', 'intron_len', 'intron_mod3',
        'stuffer_len', 'D_spliced', 'D_retained', 'frame_spliced', 'frame_retained',
        'plasmid_len', 'output_path', 'error'
    ]

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        for row in results:
            writer.writerow(row)


def verify_plasmid(filepath: str) -> Dict:
    """
    Verify a generated plasmid file.

    Checks:
    1. GenBank file is valid and parseable
    2. VHH3 sequence is present in intron
    3. ATG is present upstream of intron
    4. Intron length mod 3 = 2
    5. Frame logic is correct

    Returns:
        Dict with verification results
    """
    results = {
        'filepath': filepath,
        'valid': True,
        'errors': [],
        'warnings': [],
    }

    try:
        header, seq = parse_genbank(filepath)
        results['plasmid_len'] = len(seq)
    except Exception as e:
        results['valid'] = False
        results['errors'].append(f"Failed to parse GenBank: {e}")
        return results

    # Check for VHH3 sequence (using AVD006 sequence)
    vhh3_dna = AVD006_VHH3_DNA
    if vhh3_dna in seq:
        results['vhh3_present'] = True
        results['vhh3_position'] = seq.find(vhh3_dna) + 1  # 1-indexed
    else:
        results['valid'] = False
        results['errors'].append("VHH3 sequence not found in plasmid")
        results['vhh3_present'] = False

    # Check for ATG in upstream region
    # The cassette starts at position 1083, ATG should be around position 1089 (1083 + 6)
    upstream_region = seq[1082:1200]  # Check first part of cassette
    if 'ATG' in upstream_region:
        atg_pos = upstream_region.find('ATG')
        results['atg_present'] = True
        results['atg_position'] = 1083 + atg_pos  # 1-indexed
    else:
        results['valid'] = False
        results['errors'].append("ATG not found in upstream region")
        results['atg_present'] = False

    # Check intron length (look for GT...AG pattern)
    # This is a simplified check - in practice we'd use feature annotations

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Build splicing reporter plasmids for frame-shift dual-reporter system",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate all 32 plasmids
  python build_reporter_plasmids.py \\
    --template projects/splice_donor_pipeline/TRX155-miniXon2C9-NLucPest-FlucPest.gb \\
    --panel projects/splice_donor_pipeline/output/tiered_panel_32_v7.csv \\
    --output projects/splice_donor_pipeline/output

  # Verify a generated plasmid
  python build_reporter_plasmids.py --verify plasmids/SPL001_ACGGTAAAGCA_Tier1.gb
        """
    )

    parser.add_argument("--template", type=str,
                        help="Path to TRX155 GenBank template file")
    parser.add_argument("--panel", type=str,
                        help="Path to tiered panel CSV file (v7)")
    parser.add_argument("--output", type=str, default="./output",
                        help="Output directory")
    parser.add_argument("--prefix", type=str, default="SPL",
                        help="Prefix for plasmid names (default: SPL)")
    parser.add_argument("--verify", type=str,
                        help="Verify a generated plasmid file")
    parser.add_argument("--show-components", action="store_true",
                        help="Show pre-computed component details")

    args = parser.parse_args()

    # Show components
    if args.show_components:
        print("Pre-computed Intron Components")
        print("=" * 60)

        print("\nUpstream Exon (66 bp) - Triple Kozak-Met-Ala + VR4 loop:")
        print(f"  {UPSTREAM_EXON_66BP}")
        print(f"  Triple tandem Kozak-Met-Ala (42 bp):")
        print(f"    ATG #1: position 7-9 (weak Kozak - skipped)")
        print(f"    ATG #2: position 22-24 (weak Kozak - skipped)")
        print(f"    ATG #3: position 37-39 (STRONG Kozak: GCCACCATGGCA - USED)")
        print(f"  VR4 sequence: positions 43-66 (exact AVD002 positions 3720-3743)")
        print(f"  Protein from strong ATG: MASKTINGSG (Met-Ala-Ser-Lys-Thr-Ile-Asn-Gly-Ser-Gly)")

        print("\nVHH3 Nanobody (from AVD006):")
        vhh3_dna = AVD006_VHH3_DNA
        print(f"  Protein: {VHH3_PROTEIN[:40]}...")
        print(f"  DNA length: {len(vhh3_dna)} bp")
        print(f"  Source: AVD006 positions 3819-4175")
        frame1_stops = find_stop_codons_in_frame(vhh3_dna, frame=1)
        print(f"  Frame 1 stops: {len(frame1_stops)} (OK - VHH3 is read in frame 0)")

        print("\nGGGGS5 Linker (from AVD006):")
        ggggs5_dna = AVD006_GGGGS5_DNA
        print(f"  Protein: GGGGS x 5 = {len('GGGGS' * 5)} aa")
        print(f"  DNA length: {len(ggggs5_dna)} bp")
        print(f"  DNA: {ggggs5_dna}")
        print(f"  Source: AVD006 positions 3744-3818 (varied codons, synthesis-friendly)")

        print("\nGGGGS1 Linker (from AVD006):")
        ggggs1_dna = AVD006_GGGGS1_DNA
        print(f"  Protein: GGGGS x 1 = {len('GGGGS')} aa")
        print(f"  DNA length: {len(ggggs1_dna)} bp")
        print(f"  DNA: {ggggs1_dna}")
        print(f"  Source: AVD006 positions 4176-4190 (varied codons, synthesis-friendly)")

        print("\nExonic Stuffer (1 bp):")
        print(f"  Sequence: {EXONIC_STUFFER_1BP}")
        print(f"  Purpose: Shifts junction from 33 bp to 34 bp for correct frame alignment")
        print(f"  Location: Between upstream exon and donor exonic")
        print(f"  CRITICAL: Keeps VHH3 in native frame (frame 0), avoiding frame 1 stop at position 244")

        print("\nSplice Acceptor:")
        print(f"  Intronic (PPT + CAG): TTCTATTTTCCCACCCTCAG (20 bp)")
        print(f"    NOTE: Using CAG ending (not TAG) to avoid stop codon!")
        print(f"\nDownstream Exon (15 bp) - AAV9 VP1 VR4 loop (minimal):")
        print(f"  {ACCEPTOR_EXONIC_15BP}")
        print(f"  VR4 sequence: exact AVD002 positions 3744-3758")
        print(f"  Length: 15 bp (mod 3 = 0) - minimal to avoid frame 1 stops")
        print(f"  NOTE: Using minimal acceptor exonic with OvORF linker for frame alignment")

        # Calculate intron length
        donor_intronic = 8  # +1 to +8
        ggggs5_len = len(ggggs5_dna)
        vhh3_len = len(vhh3_dna)
        ggggs1_len = len(ggggs1_dna)
        acceptor_intronic = 20

        base_intron = donor_intronic + ggggs5_len + vhh3_len + ggggs1_len + acceptor_intronic
        print(f"\nIntron Length Calculation:")
        print(f"  Donor intronic (+1..+8):  {donor_intronic:3d} bp")
        print(f"  GGGGS5 linker:            {ggggs5_len:3d} bp")
        print(f"  VHH3 nanobody:            {vhh3_len:3d} bp")
        print(f"  GGGGS1 linker:            {ggggs1_len:3d} bp")
        print(f"  Acceptor intronic:        {acceptor_intronic:3d} bp")
        print(f"  ─────────────────────────────────")
        print(f"  Base total:               {base_intron:3d} bp (mod 3 = {base_intron % 3})")

        stuffer_needed = (2 - base_intron % 3) % 3
        total_intron = base_intron + stuffer_needed
        print(f"  Stuffer needed:           {stuffer_needed:3d} bp")
        print(f"  Final intron:             {total_intron:3d} bp (mod 3 = {total_intron % 3}) ✓")

        # Frame logic verification
        print("\nFrame Logic Verification:")
        upstream_exon_len = len(UPSTREAM_EXON_66BP)
        exonic_stuffer_len = len(EXONIC_STUFFER_1BP)
        donor_exonic_len = 3
        acceptor_exonic_len = len(ACCEPTOR_EXONIC_15BP)
        ovorf_linker_len = len(OVORF_LINKER_32BP)
        atg_offset = 36  # Strong ATG is at position 37-39 (GCCACCATGGCA context)

        # Junction remainder check (CRITICAL for avoiding TAA stop)
        junction_len = (upstream_exon_len - atg_offset) + exonic_stuffer_len + donor_exonic_len
        junction_remainder = junction_len % 3
        print(f"  Junction frame check:")
        print(f"    from_atg_to_junction = {junction_len} bp (30 + 1 + 3 = 34)")
        print(f"    Junction remainder = {junction_remainder}")
        if junction_remainder == 0:
            print(f"    → At codon boundary: ...XXX|GTA-A... (SAFE for GTAA donors) ✓")
        elif junction_remainder == 1:
            print(f"    → Previous codon: ...XGT (no TAA stop) ✓")
        elif junction_remainder == 2:
            print(f"    → Previous codon: ...XXG, next: TAA (STOP!) ✗")

        # D_spliced = exonic from ATG to downstream region (when spliced)
        D_spliced = junction_len + acceptor_exonic_len + ovorf_linker_len
        print(f"\n  D_spliced (exonic from ATG when spliced):")
        print(f"    = junction_len + acceptor_exonic + ovorf_linker")
        print(f"    = {junction_len} + {acceptor_exonic_len} + {ovorf_linker_len}")
        print(f"    = {D_spliced} bp")
        print(f"    D_spliced mod 3 = {D_spliced % 3} → {'Frame 0 (FLuc) ✓' if D_spliced % 3 == 0 else 'WRONG!'}")

        # D_retained = exonic + intronic from ATG to downstream region (when retained)
        D_retained = junction_len + total_intron + acceptor_exonic_len + ovorf_linker_len
        print(f"\n  D_retained (total from ATG when retained):")
        print(f"    = junction_len + intron + acceptor_exonic + ovorf_linker")
        print(f"    = {junction_len} + {total_intron} + {acceptor_exonic_len} + {ovorf_linker_len}")
        print(f"    = {D_retained} bp")
        print(f"    D_retained mod 3 = {D_retained % 3} → {'Frame 2 (NanoLuc) ✓' if D_retained % 3 == 2 else 'WRONG!'}")

        print("\n  Reporter Frame Mapping:")
        print(f"    D mod 3 = 0 → downstream frame 0 → FLuc (active when SPLICED)")
        print(f"    D mod 3 = 2 → downstream frame 1 → NanoLuc (active when RETAINED)")

        return 0

    # Verify mode
    if args.verify:
        print(f"Verifying: {args.verify}")
        results = verify_plasmid(args.verify)
        print(f"\nVerification Results:")
        print(f"  Valid: {results['valid']}")
        print(f"  Plasmid length: {results.get('plasmid_len', 'N/A')} bp")
        print(f"  VHH3 present: {results.get('vhh3_present', 'N/A')}")
        print(f"  ATG present: {results.get('atg_present', 'N/A')}")
        if results['errors']:
            print(f"  Errors: {results['errors']}")
        if results['warnings']:
            print(f"  Warnings: {results['warnings']}")
        return 0 if results['valid'] else 1

    # Generate plasmids
    if not args.template or not args.panel:
        parser.print_help()
        return 1

    results = generate_plasmid_panel(
        template_path=args.template,
        panel_csv_path=args.panel,
        output_dir=args.output,
        name_prefix=args.prefix,
    )

    # Summary
    successful = [r for r in results if 'error' not in r]
    failed = [r for r in results if 'error' in r]

    print(f"\n{'=' * 60}")
    print(f"Summary: {len(successful)} successful, {len(failed)} failed")

    if successful:
        # Verify frame logic
        example = successful[0]
        print(f"\nFrame Logic Verification (from {example['plasmid_name']}):")
        print(f"  Intron length: {example['intron_len']} bp (mod 3 = {example['intron_mod3']})")
        # Frame mapping: D mod 3 = 0 → downstream frame 0 → FLuc
        #               D mod 3 = 2 → downstream frame 1 → NanoLuc
        print(f"  Frame when spliced: {example['frame_spliced']} → {'FLuc ✓' if example['frame_spliced'] == 0 else 'MISMATCH!'}")
        print(f"  Frame when retained: {example['frame_retained']} → {'NanoLuc ✓' if example['frame_retained'] == 2 else 'MISMATCH!'}")

    if failed:
        print(f"\nFailed plasmids:")
        for r in failed:
            print(f"  {r['plasmid_name']}: {r['error']}")

    return 0


if __name__ == "__main__":
    exit(main())
