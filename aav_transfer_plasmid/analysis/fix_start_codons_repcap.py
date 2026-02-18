#!/usr/bin/env python3
"""
Fix start codons in RepCap v04:
1. VP1: AAG → ATG (position 2365)
2. VP2: ACC → ACG (position 2776)
3. VP3: CTG → ATG (position 2971)

These are synonymous changes (still start with Met for VP1/VP3, Thr for VP2).
"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys

print("="*80)
print("FIXING START CODONS IN REPCAP v04")
print("="*80)

# Load v04
repcap = SeqIO.read('test_data/AAV9-RepCap-NOCAP-v04.gb', 'genbank')

print(f"\n✓ Loaded: AAV9-RepCap-NOCAP-v04.gb ({len(repcap.seq)} bp)")

# Define fixes (1-indexed positions)
fixes = [
    ('VP1', 2365, 'AAG', 'ATG', 'VP1 start codon'),
    ('VP2', 2776, 'ACC', 'ACG', 'VP2 start codon (offset +411 from VP1)'),
    ('VP3', 2971, 'CTG', 'ATG', 'VP3 start codon (offset +606 from VP1)')
]

print("\n" + "="*80)
print("APPLYING START CODON FIXES")
print("="*80)

seq_str = str(repcap.seq)

for name, pos_1idx, old_codon, new_codon, note in fixes:
    # Convert to 0-indexed
    pos_0idx = pos_1idx - 1

    # Verify current codon
    current = seq_str[pos_0idx:pos_0idx+3]

    if current != old_codon:
        print(f"\n⚠️  {name} at position {pos_1idx}:")
        print(f"   Expected: {old_codon}")
        print(f"   Found: {current}")
        print(f"   Skipping this fix...")
        continue

    # Apply fix
    seq_str = seq_str[:pos_0idx] + new_codon + seq_str[pos_0idx+3:]
    print(f"✓ {name} @ {pos_1idx}: {old_codon} → {new_codon} ({note})")

# Update sequence
repcap.seq = Seq(seq_str)

# Update description
repcap.description = "AAV9-RepCap-NOCAP-v04 with corrected start codons"

# Save
output_path = 'test_data/AAV9-RepCap-NOCAP-v04.gb'
SeqIO.write(repcap, output_path, 'genbank')

print("\n" + "="*80)
print("✅ START CODONS FIXED")
print("="*80)
print(f"✓ Saved: {output_path}")

# Verify
print("\nVerifying start codons:")
for name, pos_1idx, old_codon, new_codon, note in fixes:
    pos_0idx = pos_1idx - 1
    current = str(repcap.seq[pos_0idx:pos_0idx+3])
    status = "✅" if current == new_codon else "❌"
    print(f"{status} {name} @ {pos_1idx}: {current}")
