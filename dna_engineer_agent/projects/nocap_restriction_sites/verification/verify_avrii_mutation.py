#!/usr/bin/env python3
"""
Verify that the AvrII mutation at position 2459 is truly silent
"""

from Bio import SeqIO
from Bio.Seq import Seq

# Load plasmid
record = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")
plasmid_seq = str(record.seq).upper()

print("="*80)
print("VERIFYING AvrII MUTATION AT POSITION 2459")
print("="*80)

# VP1 starts at position 2365
vp1_start = 2365

# Check the sequence at position 2459
site_start = 2459
site_end = site_start + 6  # AvrII recognition is 6 bp (CCTAGG)

# Extract current sequence
current_seq = plasmid_seq[site_start-1:site_end-1]
print(f"\nCurrent sequence at positions {site_start}-{site_end-1}:")
print(f"  {current_seq}")

# AvrII recognition sequence
avrii_seq = "CCTAGG"
print(f"\nAvrII recognition sequence:")
print(f"  {avrii_seq}")

# Identify the mutation needed
print(f"\nMutation required:")
for i, (curr, target) in enumerate(zip(current_seq, avrii_seq), start=1):
    if curr != target:
        pos = site_start + i - 1
        print(f"  Position {pos}: {curr} → {target}")

# Now check if this mutation is silent in VP1
print(f"\n{'='*80}")
print("CHECKING IF MUTATION IS SILENT IN VP1")
print(f"{'='*80}")

# VP1 starts at 2365, so check what codon position 2461 is in
mutation_pos = 2461  # This is where A→T change occurs
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
else:
    print(f"\n✗ MUTATION IS NOT SILENT!")
    print(f"  This changes {Seq(original_codon).translate()} → {Seq(mutated_codon).translate()}")

# Let's also check the entire VP1 sequence context
print(f"\n{'='*80}")
print("VP1 SEQUENCE CONTEXT")
print(f"{'='*80}")

context_start = mutation_pos - 15
context_end = mutation_pos + 15

print(f"\nSequence context around position {mutation_pos}:")
print(f"Position {context_start}-{context_end}:")
print(plasmid_seq[context_start-1:context_end-1])
print(" " * 15 + "^")
print(" " * 15 + f"Position {mutation_pos}")

# Show the reading frame
print(f"\nReading frame from VP1 start (position {vp1_start}):")
vp1_seq = plasmid_seq[vp1_start-1:context_end-1]
offset_to_context = context_start - vp1_start
print(f"\nOffset to context start: {offset_to_context}")
print(f"Frame at context start: {offset_to_context % 3}")

context_for_display = vp1_seq[offset_to_context:offset_to_context + 30]
print(f"\nDNA: {context_for_display}")
print(f"AA:  ", end="")
for i in range(0, len(context_for_display), 3):
    if i + 3 <= len(context_for_display):
        codon = context_for_display[i:i+3]
        aa = Seq(codon).translate()
        print(f"{aa}", end="  ")
print()

# Show where the mutation occurs in this display
mutation_offset_in_context = mutation_pos - context_start
print(f"\n     ", end="")
print(" " * mutation_offset_in_context + "^")
print(f"     ", end="")
print(" " * mutation_offset_in_context + f"Mutation here (A→T)")
