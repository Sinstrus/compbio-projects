#!/usr/bin/env python3
"""
Debug: Check what protein sequence is actually in VP1 for region 2415-2775
and compare it to what we're generating
"""

import sys
sys.path.insert(0, 'scripts/tools')
from Bio import SeqIO
from Bio.Seq import Seq
from silent_sites import translate

# Load plasmid
record = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")
plasmid_seq = str(record.seq).upper()

print("="*80)
print("DEBUGGING PROTEIN EXTRACTION FOR REGION 2415-2775")
print("="*80)

# VP1 starts at position 2365 (1-indexed) and goes to 4575
vp1_start = 2365
vp1_end = 4575
vp1_dna = plasmid_seq[vp1_start-1:vp1_end]  # 0-indexed slicing

print(f"\nVP1 CDS: positions {vp1_start}-{vp1_end}")
print(f"VP1 DNA length: {len(vp1_dna)} bp")

# Translate VP1 (should be frame 0 since we extracted from the start)
vp1_protein = translate(vp1_dna, frame=0, stop_at_terminator=True)
print(f"VP1 protein length: {len(vp1_protein)} aa")

# Now extract the region 2415-2775
region1_start = 2415
region1_end = 2775
region1_dna = plasmid_seq[region1_start-1:region1_end]

print(f"\n{'='*80}")
print(f"REGION 1: positions {region1_start}-{region1_end}")
print(f"{'='*80}")
print(f"Region length: {len(region1_dna)} bp")

# Calculate which amino acids this corresponds to in VP1
# Position 2415 in plasmid is position (2415-2365) = 50 in VP1 DNA (0-indexed 49)
offset_in_vp1_dna = region1_start - vp1_start  # 50 bp from start of VP1
aa_start_in_vp1 = offset_in_vp1_dna // 3  # Which amino acid does this start at?
position_in_codon = offset_in_vp1_dna % 3  # Position within that codon

print(f"\nRegion starts at offset {offset_in_vp1_dna} in VP1 DNA")
print(f"This is amino acid {aa_start_in_vp1 + 1} of VP1 (1-indexed)")
print(f"Position within codon: {position_in_codon} (0=first, 1=second, 2=third base)")

# Extract the corresponding protein sequence from VP1
# If position_in_codon != 0, we start mid-codon, so we need to handle this carefully
if position_in_codon == 0:
    # Starts at beginning of a codon - simple case
    region1_protein_from_vp1 = vp1_protein[aa_start_in_vp1:]
    print(f"\nRegion starts at codon boundary")
else:
    # Starts mid-codon - we can't extract a clean protein sequence
    # We need to skip to the next codon boundary
    bases_to_skip = 3 - position_in_codon
    aa_start_in_vp1_adjusted = aa_start_in_vp1 + 1
    region1_protein_from_vp1 = vp1_protein[aa_start_in_vp1_adjusted:]
    print(f"\nRegion starts mid-codon, need to skip {bases_to_skip} bases to next codon boundary")
    print(f"Adjusted to start at amino acid {aa_start_in_vp1_adjusted + 1} of VP1")

# Calculate end position in VP1
offset_end_in_vp1_dna = region1_end - vp1_start  # Position in VP1 DNA
aa_end_in_vp1 = offset_end_in_vp1_dna // 3

region1_protein_from_vp1 = vp1_protein[aa_start_in_vp1:aa_end_in_vp1]
region1_protein_from_vp1_adjusted = vp1_protein[aa_start_in_vp1+1:aa_end_in_vp1] if position_in_codon != 0 else region1_protein_from_vp1

print(f"Region ends at offset {offset_end_in_vp1_dna} in VP1 DNA")
print(f"This is amino acid {aa_end_in_vp1 + 1} of VP1 (1-indexed)")

print(f"\n{'='*80}")
print("ACTUAL VP1 PROTEIN IN REGION 2415-2775")
print(f"{'='*80}")
if position_in_codon == 0:
    print(f"VP1 amino acids {aa_start_in_vp1+1} to {aa_end_in_vp1} ({len(region1_protein_from_vp1)} aa):")
    print(f"  {region1_protein_from_vp1[:50]}..." if len(region1_protein_from_vp1) > 50 else f"  {region1_protein_from_vp1}")
else:
    print(f"Starts mid-codon, so adjusted to VP1 amino acids {aa_start_in_vp1+2} to {aa_end_in_vp1} ({len(region1_protein_from_vp1_adjusted)} aa):")
    print(f"  {region1_protein_from_vp1_adjusted[:50]}..." if len(region1_protein_from_vp1_adjusted) > 50 else f"  {region1_protein_from_vp1_adjusted}")

# Now let's see what we get if we translate the region DNA with different frame offsets
print(f"\n{'='*80}")
print("TRANSLATING REGION DNA WITH DIFFERENT FRAME OFFSETS")
print(f"{'='*80}")

for frame_offset in [0, 1, 2]:
    region1_protein_translated = translate(region1_dna, frame=frame_offset, stop_at_terminator=False)
    matches = region1_protein_translated == region1_protein_from_vp1_adjusted if position_in_codon != 0 else region1_protein_translated == region1_protein_from_vp1
    match_str = "✓ MATCHES" if matches else "✗ DOESN'T MATCH"
    print(f"\nFrame offset {frame_offset}:")
    print(f"  Protein length: {len(region1_protein_translated)} aa")
    print(f"  First 50 aa: {region1_protein_translated[:50]}..." if len(region1_protein_translated) > 50 else f"  Protein: {region1_protein_translated}")
    print(f"  {match_str} VP1 protein in this region")

# Calculate what the correct frame offset should be
correct_frame_offset = (3 - position_in_codon) % 3 if position_in_codon != 0 else 0
print(f"\n{'='*80}")
print(f"CONCLUSION")
print(f"{'='*80}")
print(f"Region starts at position {position_in_codon} within a codon")
print(f"Correct frame offset should be: {correct_frame_offset}")
print(f"Current script uses: {(region1_start - vp1_start) % 3} = {offset_in_vp1_dna % 3}")
