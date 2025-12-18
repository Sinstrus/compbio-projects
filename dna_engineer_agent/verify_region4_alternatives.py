#!/usr/bin/env python3
"""
Verify alternative enzymes for Region 4 (instead of FseI)
"""

from Bio import SeqIO
from Bio.Seq import Seq

# Load plasmid
record = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")
plasmid_seq = str(record.seq).upper()

vp1_start = 2365

# Alternative candidates from the analysis:
# 2. NaeI @ 3760 (1 edit) - GCCGGC
# 3. NgoMIV @ 3760 (1 edit) - GCCGGC (same as NaeI)
# 4. RsrII @ 3762 (1 edit) - CGGWCCG
# 5. BsrGI @ 3780 (1 edit) - TGTACA

alternatives = [
    ("NaeI", 3760, "GCCGGC", "A6C"),
    ("RsrII", 3762, "CGGWCCG", "C7G"),
    ("BsrGI", 3780, "TGTACA", "C4A"),
]

print("="*80)
print("VERIFYING REGION 4 ALTERNATIVE ENZYMES (instead of FseI)")
print("="*80)

for enzyme_name, site_start, recognition, mutation_desc in alternatives:
    print(f"\n{'='*80}")
    print(f"Checking: {enzyme_name} @ position {site_start}")
    print(f"Recognition: {recognition}")
    print(f"Mutation: {mutation_desc}")
    print(f"{'='*80}")

    # Extract position from mutation description (e.g., "A6C" means position 6)
    if 'W' in recognition:
        # IUPAC code - need to handle differently
        print("Note: Contains IUPAC ambiguity code (W = A or T)")
        # For RsrII: CGGWCCG means CGG[A/T]CCG
        # Let's check what's actually there
        site_seq = plasmid_seq[site_start-1:site_start-1+len(recognition)]
        print(f"Current sequence: {site_seq}")

        # Parse mutation
        mut_pos_in_site = int(''.join(filter(str.isdigit, mutation_desc)))
        mutation_pos = site_start + mut_pos_in_site - 1
        old_base = mutation_desc[0]
        new_base = mutation_desc[-1]
    else:
        site_len = len(recognition)
        site_seq = plasmid_seq[site_start-1:site_start-1+site_len]
        print(f"Current sequence: {site_seq}")

        # Parse mutation (format: "A6C" means position 6, A→C)
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
print("Based on ease of use and reliability:")
print("  1st choice: NaeI @ 3760 (very common, reliable enzyme)")
print("  2nd choice: BsrGI @ 3780 (also reliable)")
print("  Avoid: FseI (fussy, requires specific conditions)")
