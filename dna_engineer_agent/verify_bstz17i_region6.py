#!/usr/bin/env python3
"""
Manual verification of BstZ17I @ 4163 in narrowed Region 6
"""

import sys
sys.path.insert(0, 'scripts/tools')
from Bio import SeqIO
from silent_sites import translate

# Load plasmid
record = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")
plasmid_seq = str(record.seq).upper()

# BstZ17I site at position 4163
site_pos = 4163  # 1-indexed
site_seq = "GTATAC"

# Extract context
context_start = site_pos - 1 - 12  # 0-indexed, 12 bases before
context_end = site_pos - 1 + len(site_seq) + 12  # 12 bases after

print("BstZ17I @ position 4163 (narrowed Region 6)")
print("="*80)

# Show original sequence
original_seq = plasmid_seq[context_start:context_end]
print(f"\nOriginal sequence (positions {context_start+1}-{context_end}):")
print(f"  {original_seq}")
print(f"  {' '*12}{'▼'*len(site_seq)} (site @ {site_pos})")

# Extract the actual site
actual_site = plasmid_seq[site_pos-1:site_pos-1+len(site_seq)]
print(f"\nActual sequence at {site_pos}: {actual_site}")
print(f"Target site (BstZ17I):        {site_seq}")

# Check if mutation is needed
if actual_site != site_seq:
    print(f"\nMutation needed:")
    for i, (a, t) in enumerate(zip(actual_site, site_seq), 1):
        if a != t:
            abs_pos = site_pos + i - 1
            print(f"  Position {i} of site (absolute {abs_pos}): {a} → {t}")

    # Apply mutation
    mutated_seq = list(original_seq)
    # Find position in context
    pos_in_context = (site_pos - 1) - context_start + 1  # Position of mutation in context
    print(f"  Position in context: {pos_in_context}")
    mutated_seq[pos_in_context] = site_seq[1]  # A→T at position 2
    mutated_seq_str = ''.join(mutated_seq)

    print(f"\nMutated sequence:")
    print(f"  {mutated_seq_str}")
else:
    print(f"\n✅ Site already exists (no mutation needed)")
    mutated_seq_str = original_seq

# Now check the codon context
# Region 6 starts at 4144 with frame offset 0
region_start = 4144
frame_offset = 0

# Calculate position of site relative to region start
site_offset = site_pos - region_start  # 4163 - 4144 = 19 (0-indexed)

# Extract a larger region for codon context (start from region beginning)
region_for_translation = plasmid_seq[region_start-1:context_end]

print(f"\n" + "="*80)
print("CODON ANALYSIS")
print("="*80)

# Find which codon contains the mutation site
codon_index = (site_offset + 1) // 3  # Which codon (0-indexed)
position_in_codon = (site_offset + 1) % 3  # Position within codon (0, 1, or 2)

print(f"Site position relative to region start: {site_offset + 1} (1-indexed)")
print(f"Codon index: {codon_index}")
print(f"Position in codon: {position_in_codon}")

# Extract the codon containing the mutation
codon_start = codon_index * 3
codon_original = region_for_translation[codon_start:codon_start+3]

print(f"\nOriginal codon: {codon_original} → {translate(codon_original, frame=0)}")

# Apply mutation to the codon
mutated_region = list(region_for_translation)
mutation_pos_in_region = site_offset + 1  # Position of A→T in region
mutated_region[mutation_pos_in_region] = 'T'
mutated_region_str = ''.join(mutated_region)

codon_mutated = mutated_region_str[codon_start:codon_start+3]
print(f"Mutated codon:  {codon_mutated} → {translate(codon_mutated, frame=0)}")

if translate(codon_original, frame=0) == translate(codon_mutated, frame=0):
    print(f"\n✅ SILENT MUTATION VERIFIED: {codon_original} → {codon_mutated} (both {translate(codon_original, frame=0)})")
else:
    print(f"\n❌ NOT SILENT: {codon_original} ({translate(codon_original, frame=0)}) → {codon_mutated} ({translate(codon_mutated, frame=0)})")

# Full protein check
print(f"\n" + "="*80)
print("FULL PROTEIN VERIFICATION")
print("="*80)

region_full = plasmid_seq[region_start-1:region_start+199]
original_protein = translate(region_full, frame=frame_offset, stop_at_terminator=False)

mutated_full = list(region_full)
mutated_full[site_offset + 1] = 'T'
mutated_full_str = ''.join(mutated_full)
mutated_protein = translate(mutated_full_str, frame=frame_offset, stop_at_terminator=False)

print(f"Original protein: {original_protein[:60]}...")
print(f"Mutated protein:  {mutated_protein[:60]}...")

if original_protein == mutated_protein:
    print(f"\n✅ FULL VERIFICATION PASSED: Proteins are identical")
else:
    print(f"\n❌ VERIFICATION FAILED: Proteins differ")
    for i, (o, m) in enumerate(zip(original_protein, mutated_protein)):
        if o != m:
            print(f"  Position {i}: {o} → {m}")
