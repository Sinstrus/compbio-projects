#!/usr/bin/env python3
"""
Verify alternative enzymes for Region 3 (instead of EcoRV)
"""

from Bio import SeqIO
from Bio.Seq import Seq

# Load plasmid
record = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")
plasmid_seq = str(record.seq).upper()

vp1_start = 2365

# Alternative candidates (excluding EcoRV):
# 2. PshAI @ 3533 (1 edit) - GACNNNNGTC
# 3. BsmBI @ 3566 (1 edit) - CGTCTC

alternatives = [
    ("PshAI", 3533, "GACNNNNGTC", "T2A"),
    ("BsmBI", 3566, "CGTCTC", "G5T"),
]

print("="*80)
print("VERIFYING REGION 3 ALTERNATIVE ENZYMES (instead of EcoRV)")
print("="*80)

for enzyme_name, site_start, recognition, mutation_desc in alternatives:
    print(f"\n{'='*80}")
    print(f"Checking: {enzyme_name} @ position {site_start}")
    print(f"Recognition: {recognition}")
    print(f"Mutation: {mutation_desc}")
    print(f"{'='*80}")

    # Extract current sequence
    site_len = len(recognition)
    site_seq = plasmid_seq[site_start-1:site_start-1+site_len]
    print(f"Current sequence: {site_seq}")

    # Parse mutation (format: "T2A" means position 2, T→A)
    mut_pos_in_site = int(''.join(filter(str.isdigit, mutation_desc)))
    mutation_pos = site_start + mut_pos_in_site - 1
    old_base = mutation_desc[0]
    new_base = mutation_desc[-1]

    print(f"\nMutation at position {mutation_pos}: {old_base} → {new_base}")

    # Check if silent
    offset_from_vp1 = mutation_pos - vp1_start
    frame_position = offset_from_vp1 % 3

    print(f"Offset from VP1 start: {offset_from_vp1}")
    print(f"Frame position: {frame_position} (0=first, 1=second, 2=third)")

    # Extract codon
    if frame_position == 0:
        codon_start = mutation_pos
    elif frame_position == 1:
        codon_start = mutation_pos - 1
    else:
        codon_start = mutation_pos - 2

    codon_end = codon_start + 3
    codon_num = (codon_start - vp1_start) // 3 + 1

    print(f"\nCodon {codon_num} spans positions {codon_start}-{codon_end-1}")

    # Original codon
    original_codon = plasmid_seq[codon_start-1:codon_end-1]
    print(f"Original codon: {original_codon} → {Seq(original_codon).translate()}")

    # Mutated codon
    mutated_plasmid = list(plasmid_seq)
    mutated_plasmid[mutation_pos-1] = new_base
    mutated_codon = ''.join(mutated_plasmid[codon_start-1:codon_end-1])
    print(f"Mutated codon:  {mutated_codon} → {Seq(mutated_codon).translate()}")

    # Check if silent
    if Seq(original_codon).translate() == Seq(mutated_codon).translate():
        print(f"✅ SILENT - Both encode {Seq(original_codon).translate()}")
    else:
        print(f"❌ NOT SILENT - Changes {Seq(original_codon).translate()} → {Seq(mutated_codon).translate()}")

print(f"\n{'='*80}")
print("RECOMMENDATION")
print(f"{'='*80}")
print("Based on common usage and reliability:")
print("  1st choice: PshAI @ 3533 (reliable, Type II restriction enzyme)")
print("  2nd choice: BsmBI @ 3566 (Type IIS, cuts outside recognition site)")
print("  Avoid: EcoRV (commonly used by gene synthesis companies)")
