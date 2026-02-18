#!/usr/bin/env python3
"""
Verify that all feature annotations are correctly preserved in the mutated plasmid
by checking that features still start and end with the same sequences
"""

from Bio import SeqIO

print("="*80)
print("VERIFYING FEATURE SEQUENCE PRESERVATION")
print("="*80)

# Load both files
original = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")
mutated = SeqIO.read("test_data/AAV9-RepCap-NOCAP-WITH-RESTRICTION-SITES.gb", "genbank")

original_seq = str(original.seq).upper()
mutated_seq = str(mutated.seq).upper()

print(f"\nOriginal plasmid: {len(original_seq)} bp, {len(original.features)} features")
print(f"Mutated plasmid:  {len(mutated_seq)} bp, {len(mutated.features)} features")

# Check that lengths are the same
if len(original_seq) != len(mutated_seq):
    print(f"\n❌ ERROR: Plasmid lengths differ!")
    exit(1)

# Count differences
differences = []
for i, (o, m) in enumerate(zip(original_seq, mutated_seq), 1):
    if o != m:
        differences.append((i, o, m))

print(f"\nTotal nucleotide differences: {len(differences)}")
print(f"Expected: 6 (one per restriction site)")

if len(differences) != 6:
    print(f"❌ ERROR: Expected 6 differences, found {len(differences)}")
else:
    print(f"✅ Correct number of mutations")

print(f"\nMutations:")
for pos, orig, mut in differences:
    print(f"  Position {pos}: {orig} → {mut}")

# Verify key features have correct sequences
print("\n" + "="*80)
print("VERIFYING KEY FEATURES")
print("="*80)

# Find matching features in both files
key_features = ["Rep68", "VP1", "VP2", "VP3", "AAV9 AAP"]

all_valid = True

for feat_name in key_features:
    # Find feature in original
    orig_feature = None
    mut_feature = None

    for f in original.features:
        if f.type == "CDS":
            label = f.qualifiers.get("label", [""])[0]
            if feat_name in label:
                orig_feature = f
                break

    for f in mutated.features:
        if f.type == "CDS":
            label = f.qualifiers.get("label", [""])[0]
            if feat_name in label:
                mut_feature = f
                break

    if not orig_feature or not mut_feature:
        print(f"\n{feat_name}: Not found in one or both files")
        continue

    # Extract sequences
    orig_start = int(orig_feature.location.start)
    orig_end = int(orig_feature.location.end)
    mut_start = int(mut_feature.location.start)
    mut_end = int(mut_feature.location.end)

    # Check coordinates are the same
    if orig_start != mut_start or orig_end != mut_end:
        print(f"\n{feat_name}:")
        print(f"  ❌ ERROR: Coordinates changed!")
        print(f"  Original: {orig_start+1}-{orig_end}")
        print(f"  Mutated:  {mut_start+1}-{mut_end}")
        all_valid = False
        continue

    # Extract sequences
    orig_seq_feat = str(original.seq[orig_start:orig_end])
    mut_seq_feat = str(mutated.seq[mut_start:mut_end])

    # Check first 20 and last 20 bases
    orig_first = orig_seq_feat[:20]
    mut_first = mut_seq_feat[:20]
    orig_last = orig_seq_feat[-20:]
    mut_last = mut_seq_feat[-20:]

    print(f"\n{feat_name} ({orig_start+1}-{orig_end}, {len(orig_seq_feat)} bp):")

    if orig_first != mut_first:
        print(f"  ❌ ERROR: First 20 bp differ!")
        print(f"    Original: {orig_first}")
        print(f"    Mutated:  {mut_first}")
        all_valid = False
    else:
        print(f"  ✅ First 20 bp: {orig_first}")

    if orig_last != mut_last:
        print(f"  ❌ ERROR: Last 20 bp differ!")
        print(f"    Original: {orig_last}")
        print(f"    Mutated:  {mut_last}")
        all_valid = False
    else:
        print(f"  ✅ Last 20 bp:  {orig_last}")

    # Count differences within feature
    feat_diffs = 0
    for i, (o, m) in enumerate(zip(orig_seq_feat, mut_seq_feat), orig_start+1):
        if o != m:
            feat_diffs += 1

    if feat_diffs > 0:
        print(f"  Changes within feature: {feat_diffs} nucleotide(s)")
    else:
        print(f"  No changes within this feature")

# Check ITRs (very important!)
print("\n" + "="*80)
print("VERIFYING ITRs (CRITICAL)")
print("="*80)

for f in original.features:
    if f.type == "misc_feature" or f.type == "repeat_region":
        label = f.qualifiers.get("label", [""])[0]
        if "ITR" in label:
            start = int(f.location.start)
            end = int(f.location.end)

            orig_itr = str(original.seq[start:end])
            mut_itr = str(mutated.seq[start:end])

            print(f"\n{label} ({start+1}-{end}):")
            if orig_itr == mut_itr:
                print(f"  ✅ ITR sequence unchanged")
                print(f"  First 20: {orig_itr[:20]}")
                print(f"  Last 20:  {orig_itr[-20:]}")
            else:
                print(f"  ❌ ERROR: ITR sequence changed!")
                all_valid = False

print("\n" + "="*80)
if all_valid:
    print("✅ ALL FEATURES VERIFIED - SEQUENCES CORRECTLY PRESERVED")
else:
    print("❌ SOME FEATURES FAILED VERIFICATION")
print("="*80)
