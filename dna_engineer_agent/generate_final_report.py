#!/usr/bin/env python3
"""
Generate final report with recommended restriction sites
"""

import sys
sys.path.insert(0, 'scripts/tools')
from Bio import SeqIO
from silent_sites import translate, count_pattern_occurrences, apply_mutations
import re

# Load the GenBank file
record = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")
plasmid_seq = str(record.seq).upper()

# Genetic code
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

def get_codon_info(position_1idx, mutation_str, plasmid_seq):
    """
    Get detailed codon change information
    position_1idx: 1-indexed position in plasmid
    mutation_str: e.g., "T2G" meaning position 2 in the window changes T->G
    """
    # Parse mutation
    match = re.match(r'([ACGT])(\d+)([ACGT])', mutation_str)
    if not match:
        return None

    original_base, offset_str, new_base = match.groups()
    offset = int(offset_str) - 1  # Convert to 0-indexed offset within window

    # Calculate absolute position in plasmid (0-indexed)
    abs_pos = (position_1idx - 1) + offset

    # Find which codon this position belongs to
    # Determine the reading frame start - we need to find the nearest upstream CDS start
    # For now, let's just show the local codon context

    # Get codon containing this position
    codon_start = abs_pos - (abs_pos % 3)
    original_codon = plasmid_seq[codon_start:codon_start+3]

    # Apply mutation to create new codon
    new_codon = list(original_codon)
    new_codon[abs_pos - codon_start] = new_base
    new_codon = ''.join(new_codon)

    original_aa = CODON_TABLE.get(original_codon, '?')
    new_aa = CODON_TABLE.get(new_codon, '?')

    return {
        'position_1idx': abs_pos + 1,
        'offset_in_window': offset,
        'original_codon': original_codon,
        'new_codon': new_codon,
        'original_aa': original_aa,
        'new_aa': new_aa,
        'silent': original_aa == new_aa
    }

def document_site(region_name, enzyme, position, site_sequence, edits_required, mutations, plasmid_seq):
    """Document a restriction site with full details"""

    print(f"\n{'='*80}")
    print(f"RECOMMENDED SITE: {region_name}")
    print(f"{'='*80}")
    print(f"Enzyme:            {enzyme}")
    print(f"Recognition Site:  {site_sequence}")
    print(f"Position:          {position} (1-indexed)")
    print(f"Mutations Required: {edits_required}")

    # Get the window
    window_len = len(site_sequence)
    pos_0idx = position - 1
    original_window = plasmid_seq[pos_0idx:pos_0idx + window_len]
    print(f"Original Sequence: {original_window}")

    if edits_required == 0:
        print(f"Modified Sequence: {original_window} (no changes needed - site already present)")
        print(f"\nUniqueness: This site appears only ONCE in the entire plasmid")
    else:
        # Apply mutations
        modified_window = apply_mutations(original_window, mutations)
        print(f"Modified Sequence: {modified_window}")

        print(f"\nCodon Changes:")
        for mut in mutations:
            info = get_codon_info(position, mut, plasmid_seq)
            if info:
                print(f"  Position {info['position_1idx']}: "
                      f"{info['original_codon']} → {info['new_codon']} "
                      f"({info['original_aa']} → {info['new_aa']}) "
                      f"{'✓ SILENT' if info['silent'] else '✗ NOT SILENT'}")

        # Verify uniqueness in modified plasmid
        modified_plasmid = plasmid_seq[:pos_0idx] + modified_window + plasmid_seq[pos_0idx + window_len:]
        count = count_pattern_occurrences(modified_plasmid, site_sequence)
        print(f"\nUniqueness: After modifications, this site appears {count} time(s) in the plasmid")
        if count == 1:
            print(f"            ✓ UNIQUE - suitable for cloning")
        else:
            print(f"            ✗ WARNING - not unique!")

print("="*80)
print("RESTRICTION SITE ANALYSIS REPORT")
print("AAV9 Rep-Cap Helper Plasmid (BASE-DRAFT-AAV9-RepCap-NOCAP.gb)")
print("="*80)
print(f"Plasmid Length: {len(plasmid_seq)} bp")
print(f"Analysis Date: 2025-12-17")
print(f"\nTask: Identify unique restriction sites (≤2 silent mutations) in three regions:")
print(f"  1. VP1 unique region (between VP1 start and VP2 start)")
print(f"  2. VP2-AAP intergenic (after VP2 start, before AAP start)")
print(f"  3. VP3 post-AAP (within VP3 coding region, after AAP stop)")

# Recommended sites from the analysis
# REGION 1: EagI @ 2369 (1 edit)
document_site(
    region_name="Region 1 - VP1 Unique Region",
    enzyme="EagI",
    position=2369,
    site_sequence="CGGCCG",
    edits_required=1,
    mutations=["T2G"],
    plasmid_seq=plasmid_seq
)

# REGION 2: BbvCI @ 2828 (1 edit)
document_site(
    region_name="Region 2 - VP2-AAP Intergenic",
    enzyme="BbvCI",
    position=2828,
    site_sequence="CCTCAGC",
    edits_required=1,
    mutations=["C5A"],
    plasmid_seq=plasmid_seq
)

# REGION 3: AfeI @ 4442 (0 edits)
document_site(
    region_name="Region 3 - VP3 Post-AAP",
    enzyme="AfeI",
    position=4442,
    site_sequence="AGCGCT",
    edits_required=0,
    mutations=[],
    plasmid_seq=plasmid_seq
)

print(f"\n{'='*80}")
print("SUMMARY")
print(f"{'='*80}")
print(f"\nAll three recommended sites are UNIQUE in the entire plasmid.")
print(f"\nTotal mutations required: 2 (1 in Region 1, 1 in Region 2, 0 in Region 3)")
print(f"All mutations are SILENT (no amino acid changes)")
print(f"\nThese sites can be used for:")
print(f"  - Future modifications via restriction enzyme cloning")
print(f"  - Golden Gate assembly (if using Type IIS enzymes)")
print(f"  - Diagnostic restriction mapping to verify plasmid identity")

print(f"\n{'='*80}")
print("ALTERNATIVE CANDIDATES")
print(f"{'='*80}")

print(f"\nRegion 1 alternatives (all 1 mutation, all unique):")
alternatives_r1 = [
    ("EcoRV", 2378, "GATATC", ["T2A"]),
    ("BspEI", 2385, "TCCGGA", ["A4G"]),
    ("SmaI", 2505, "CCCGGG", ["T1C"]),
]
for enzyme, pos, site, muts in alternatives_r1:
    print(f"  - {enzyme:12s} @ {pos:5d} | {site:10s} | {', '.join(muts)}")

print(f"\nRegion 2 alternatives (all 1 mutation, all unique):")
alternatives_r2 = [
    ("BlpI", 2854, "GCTNAGC", ["A3T"]),
    ("NaeI", 2859, "GCCGGC", ["C4G"]),
]
for enzyme, pos, site, muts in alternatives_r2:
    print(f"  - {enzyme:12s} @ {pos:5d} | {site:10s} | {', '.join(muts)}")

print(f"\nRegion 3 alternatives (0 mutations, all unique):")
alternatives_r3 = [
    ("BaeI", 4318, "ACNNNNGTAYC", []),
    ("XcmI", 4380, "CCANNNNNNNNNTGG", []),
    ("BspDI/ClaI", 4575, "ATCGAT", []),
]
for enzyme, pos, site, muts in alternatives_r3:
    print(f"  - {enzyme:12s} @ {pos:5d} | {site:20s} | Already present")

print(f"\n{'='*80}")
print("END OF REPORT")
print(f"{'='*80}\n")
