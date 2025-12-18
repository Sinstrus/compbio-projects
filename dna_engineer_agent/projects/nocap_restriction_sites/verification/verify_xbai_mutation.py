#!/usr/bin/env python3
"""
Verify that the XbaI mutation at position 2478 is truly silent
"""

from Bio import SeqIO
from Bio.Seq import Seq

# Load plasmid
record = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")
plasmid_seq = str(record.seq).upper()

print("="*80)
print("VERIFYING XbaI MUTATION AT POSITION 2478")
print("="*80)

# VP1 starts at position 2365
vp1_start = 2365

# Check the sequence at position 2478
site_start = 2478
site_end = site_start + 6  # XbaI recognition is 6 bp (TCTAGA)

# Extract current sequence
current_seq = plasmid_seq[site_start-1:site_end-1]
print(f"\nCurrent sequence at positions {site_start}-{site_end-1}:")
print(f"  {current_seq}")

# XbaI recognition sequence
xbai_seq = "TCTAGA"
print(f"\nXbaI recognition sequence:")
print(f"  {xbai_seq}")

# Identify the mutation needed
print(f"\nMutation required:")
for i, (curr, target) in enumerate(zip(current_seq, xbai_seq), start=1):
    if curr != target:
        pos = site_start + i - 1
        print(f"  Position {pos}: {curr} → {target}")

# Now check if this mutation is silent in VP1
print(f"\n{'='*80}")
print("CHECKING IF MUTATION IS SILENT IN VP1")
print(f"{'='*80}")

# The mutation is at position 2480 (A→T based on "A3T" meaning position 3 of site)
mutation_pos = 2480  # site_start + 2 (position 3, 0-indexed 2)
offset_from_vp1 = mutation_pos - vp1_start
frame_position = offset_from_vp1 % 3

print(f"\nVP1 starts at position: {vp1_start}")
print(f"Mutation at position: {mutation_pos}")
print(f"Offset from VP1 start: {offset_from_vp1}")
print(f"Frame position (0=first, 1=second, 2=third): {frame_position}")

# Extract the codon containing this position
if frame_position == 0:
    codon_start = mutation_pos
elif frame_position == 1:
    codon_start = mutation_pos - 1
else:  # frame_position == 2
    codon_start = mutation_pos - 2

codon_end = codon_start + 3
codon_num = (codon_start - vp1_start) // 3 + 1

print(f"\nCodon {codon_num} spans positions {codon_start}-{codon_end-1}")

# Extract original codon
original_codon = plasmid_seq[codon_start-1:codon_end-1]
print(f"Original codon: {original_codon}")
print(f"  Translates to: {Seq(original_codon).translate()}")

# Create mutated sequence
mutated_plasmid = list(plasmid_seq)
mutated_plasmid[mutation_pos-1] = 'T'  # Change A to T
mutated_codon = ''.join(mutated_plasmid[codon_start-1:codon_end-1])
print(f"\nMutated codon: {mutated_codon}")
print(f"  Translates to: {Seq(mutated_codon).translate()}")

# Check if silent
if Seq(original_codon).translate() == Seq(mutated_codon).translate():
    print(f"\n✓ MUTATION IS SILENT")
    print(f"  Both translate to: {Seq(original_codon).translate()}")
else:
    print(f"\n✗ MUTATION IS NOT SILENT!")
    print(f"  This changes {Seq(original_codon).translate()} → {Seq(mutated_codon).translate()}")
