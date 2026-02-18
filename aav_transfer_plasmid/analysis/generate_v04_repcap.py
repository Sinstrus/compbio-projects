#!/usr/bin/env python3
"""
Generate AAV9-RepCap-NOCAP-v04.gb with 6 unique restriction sites:
- Region 1: AvrII (replaces XbaI from v03) - position 2523
- Region 2: BspEI (kept from v03) - position 2850
- Region 3: BsmBI (kept from v03) - position 3574
- Region 4: BsrGI (kept from v03) - position 3782
- Region 5: BmtI (kept from v03) - position 3943
- Region 6: BstZ17I (kept from v03) - position 4167
"""

from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.Restriction import Analysis, RestrictionBatch, AvrII, BspEI, BsmBI, BsrGI, BmtI, BstZ17I
import sys

# Mutations to apply (1-indexed positions, verified from BASE-DRAFT and v03 comparison)
mutations = [
    # Region 1: AvrII (replaces XbaI from v03)
    (2523, 'T', 'A', 'AvrII', 'VP1 aa 53: CTT→CTA (L→L)'),

    # Region 2: BspEI (from v03)
    (2850, 'G', 'C', 'BspEI', 'VP1 aa 162: TCG→TCC (S→S)'),
    (2853, 'T', 'A', 'BspEI', 'VP1 aa 163: GGT→GGA (G→G)'),

    # Region 3: BsmBI (from v03)
    (3570, 'G', 'T', 'BsmBI', 'VP1 aa 402: TCG→TCT (S→S)'),

    # Region 4: BsrGI (from v03)
    (3783, 'C', 'A', 'BsrGI', 'VP1 aa 473: GTC→GTA (V→V)'),

    # Region 5: BmtI (from v03)
    (3939, 'C', 'T', 'BmtI', 'VP1 aa 525: GCC→GCT (A→A)'),

    # Region 6: BstZ17I (from v03)
    (4164, 'A', 'T', 'BstZ17I', 'VP1 aa 600: GGA→GGT (G→G)')
]

print("Loading BASE-DRAFT-AAV9-RepCap-NOCAP.gb...")
base_draft = SeqIO.read('test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb', 'genbank')

# Convert to mutable string for modifications
seq_str = str(base_draft.seq)

print("\nApplying mutations:")
print("=" * 70)

for pos, orig, new, enzyme, note in mutations:
    # Convert to 0-indexed
    idx = pos - 1

    # Verify original base
    if seq_str[idx] != orig:
        print(f"ERROR at position {pos}: Expected '{orig}', found '{seq_str[idx]}'")
        sys.exit(1)

    # Apply mutation
    seq_str = seq_str[:idx] + new + seq_str[idx+1:]
    print(f"{enzyme:10s} Position {pos:4d}: {orig}→{new} | {note}")

# Convert back to Seq
base_draft.seq = Seq(seq_str)

# Update description
base_draft.description = "AAV9-RepCap-NOCAP-v04 with 6 unique restriction sites"

print("\n" + "=" * 70)
print("Verifying restriction sites...")
print("=" * 70)

enzymes = [AvrII, BspEI, BsmBI, BsrGI, BmtI, BstZ17I]
batch = RestrictionBatch(enzymes)
analysis = Analysis(batch, base_draft.seq)
sites = analysis.full()

all_good = True
for enzyme in enzymes:
    enzyme_sites = sites.get(enzyme, [])
    count = len(enzyme_sites)
    positions = [p+1 for p in enzyme_sites]

    if count == 1:
        print(f"✅ {enzyme.__name__:10s} {count} site @ position {positions[0]}")
    else:
        print(f"❌ {enzyme.__name__:10s} {count} sites @ {positions}")
        all_good = False

if not all_good:
    print("\n❌ ERROR: Not all sites are unique!")
    sys.exit(1)

print("\n" + "=" * 70)
print("✅ All sites verified unique!")
print("=" * 70)

# Save
output_path = 'test_data/AAV9-RepCap-NOCAP-v04.gb'
SeqIO.write(base_draft, output_path, 'genbank')
print(f"\n✅ Saved: {output_path}")
print(f"   Size: {len(base_draft.seq)} bp")
