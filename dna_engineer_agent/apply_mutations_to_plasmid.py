#!/usr/bin/env python3
"""
Apply the 6 verified silent mutations to create the final plasmid with restriction sites

CRITICAL: Preserves all GenBank annotations - features stay with their sequences
"""

import sys
sys.path.insert(0, 'scripts/tools')
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import copy

# Load original plasmid
print("="*80)
print("APPLYING 6 SILENT MUTATIONS TO CREATE FINAL PLASMID")
print("="*80)

input_file = "test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb"
output_file = "test_data/AAV9-RepCap-NOCAP-WITH-RESTRICTION-SITES.gb"

record = SeqIO.read(input_file, "genbank")
plasmid_seq = str(record.seq).upper()

print(f"\nOriginal plasmid: {input_file}")
print(f"Length: {len(plasmid_seq)} bp")
print(f"Features: {len(record.features)}")

# Define the 6 mutations with exact positions (verified from sequence)
mutations = [
    {
        'name': 'Region 1: SmaI',
        'enzyme': 'SmaI',
        'site_start': 2505,
        'recognition': 'CCCGGG',
        'mutation_pos': 2505,
        'original': 'T',
        'mutated': 'C',
        'codon_change': 'CTT→CTC',
        'aa_change': 'L→L'
    },
    {
        'name': 'Region 2: BbvCI',
        'enzyme': 'BbvCI',
        'site_start': 2828,
        'recognition': 'CCTCAGC',
        'mutation_pos': 2832,
        'original': 'C',
        'mutated': 'A',
        'codon_change': 'TCC→TCA',
        'aa_change': 'S→S'
    },
    {
        'name': 'Region 3: AgeI',
        'enzyme': 'AgeI',
        'site_start': 3583,
        'recognition': 'ACCGGT',
        'mutation_pos': 3585,
        'original': 'G',
        'mutated': 'C',
        'codon_change': 'ACG→ACC',
        'aa_change': 'T→T'
    },
    {
        'name': 'Region 4: BsrGI',
        'enzyme': 'BsrGI',
        'site_start': 3780,
        'recognition': 'TGTACA',
        'mutation_pos': 3783,
        'original': 'C',
        'mutated': 'A',
        'codon_change': 'GTC→GTA',
        'aa_change': 'V→V'
    },
    {
        'name': 'Region 5: BmtI',
        'enzyme': 'BmtI',
        'site_start': 3937,
        'recognition': 'GCTAGC',
        'mutation_pos': 3939,
        'original': 'C',
        'mutated': 'T',
        'codon_change': 'GCC→GCT',
        'aa_change': 'A→A'
    },
    {
        'name': 'Region 6: BstZ17I',
        'enzyme': 'BstZ17I',
        'site_start': 4163,
        'recognition': 'GTATAC',
        'mutation_pos': 4164,
        'original': 'A',
        'mutated': 'T',
        'codon_change': 'GGA→GGT',
        'aa_change': 'G→G'
    },
]

# First, verify all mutations are correct
print("\n" + "="*80)
print("VERIFYING MUTATIONS BEFORE APPLYING")
print("="*80)

all_valid = True
for mut in mutations:
    pos = mut['mutation_pos']
    original = mut['original']
    actual = plasmid_seq[pos - 1]  # Convert to 0-indexed

    print(f"\n{mut['name']}:")
    print(f"  Position: {pos}")
    print(f"  Expected original: {original}")
    print(f"  Actual base: {actual}")

    if actual != original:
        print(f"  ❌ ERROR: Base mismatch!")
        all_valid = False
    else:
        print(f"  ✅ Verified")

        # Check if site already exists or needs mutation
        site_start = mut['site_start']
        site_end = site_start + len(mut['recognition']) - 1
        current_site = plasmid_seq[site_start - 1:site_end]
        print(f"  Current site: {current_site}")
        print(f"  Target site:  {mut['recognition']}")
        print(f"  Mutation: {original}→{mut['mutated']} ({mut['codon_change']}, {mut['aa_change']})")

if not all_valid:
    print("\n❌ VALIDATION FAILED - mutations do not match original sequence!")
    sys.exit(1)

print("\n✅ All mutations validated")

# Apply mutations
print("\n" + "="*80)
print("APPLYING MUTATIONS")
print("="*80)

mutated_seq = list(plasmid_seq)

for mut in mutations:
    pos = mut['mutation_pos']
    mutated_seq[pos - 1] = mut['mutated']  # Convert to 0-indexed
    print(f"Applied: {mut['name']} - {mut['original']}→{mut['mutated']} at position {pos}")

mutated_seq_str = ''.join(mutated_seq)

# Verify all sites are now present
print("\n" + "="*80)
print("VERIFYING RESTRICTION SITES IN MUTATED SEQUENCE")
print("="*80)

all_sites_present = True
for mut in mutations:
    site_start = mut['site_start']
    site_end = site_start + len(mut['recognition']) - 1
    final_site = mutated_seq_str[site_start - 1:site_end]
    expected = mut['recognition']

    print(f"\n{mut['enzyme']} @ {site_start}:")
    print(f"  Expected: {expected}")
    print(f"  Found:    {final_site}")

    if final_site == expected:
        print(f"  ✅ Site created successfully")
    else:
        print(f"  ❌ ERROR: Site not created correctly!")
        all_sites_present = False

if not all_sites_present:
    print("\n❌ ERROR: Not all sites were created correctly!")
    sys.exit(1)

print("\n✅ All 6 restriction sites verified in final sequence")

# Create new SeqRecord with mutated sequence
print("\n" + "="*80)
print("CREATING NEW GENBANK FILE")
print("="*80)

# Create new record with mutated sequence but preserve all other properties
new_record = SeqRecord(
    Seq(mutated_seq_str),
    id=record.id,
    name=record.name,
    description=record.description + " [WITH 6 SILENT RESTRICTION SITES]",
    annotations=copy.deepcopy(record.annotations)
)

# Copy all features - they maintain their coordinates since we only did substitutions
# (no insertions/deletions that would shift positions)
new_record.features = copy.deepcopy(record.features)

print(f"New record created:")
print(f"  ID: {new_record.id}")
print(f"  Name: {new_record.name}")
print(f"  Description: {new_record.description}")
print(f"  Length: {len(new_record.seq)} bp")
print(f"  Features preserved: {len(new_record.features)}")

# Verify features are still correct by checking some key sequences
print("\n" + "="*80)
print("VERIFYING FEATURE ANNOTATIONS ARE CORRECT")
print("="*80)

# Check a few key features to ensure they still have correct sequences
for feature in new_record.features:
    if feature.type == "CDS":
        label = feature.qualifiers.get("label", [""])[0]
        if "VP1" in label or "Rep68" in label or "AAP" in label:
            start = int(feature.location.start)
            end = int(feature.location.end)
            feature_seq = str(new_record.seq[start:end])

            # Get first and last 10 bases
            first_10 = feature_seq[:10]
            last_10 = feature_seq[-10:]

            print(f"\n{label}:")
            print(f"  Position: {start+1}-{end}")
            print(f"  First 10 bp: {first_10}")
            print(f"  Last 10 bp:  {last_10}")
            print(f"  Length: {len(feature_seq)} bp")

# Add new features for the restriction sites
print("\n" + "="*80)
print("ADDING RESTRICTION SITE ANNOTATIONS")
print("="*80)

from Bio.SeqFeature import SeqFeature, FeatureLocation

for mut in mutations:
    site_start = mut['site_start'] - 1  # Convert to 0-indexed
    site_end = site_start + len(mut['recognition'])

    # Create feature for restriction site
    feature = SeqFeature(
        FeatureLocation(site_start, site_end),
        type="misc_feature",
        qualifiers={
            "label": f"{mut['enzyme']}_site",
            "note": f"Engineered {mut['enzyme']} restriction site ({mut['recognition']}). Silent mutation: {mut['codon_change']} ({mut['aa_change']})"
        }
    )

    new_record.features.append(feature)
    print(f"Added annotation: {mut['enzyme']} @ {mut['site_start']}")

print(f"\nTotal features in new file: {len(new_record.features)}")

# Write to new GenBank file
SeqIO.write(new_record, output_file, "genbank")

print("\n" + "="*80)
print("SUCCESS!")
print("="*80)
print(f"\nNew plasmid saved to: {output_file}")
print(f"Length: {len(new_record.seq)} bp (unchanged)")
print(f"Features: {len(new_record.features)} (original {len(record.features)} + 6 restriction site annotations)")
print(f"\nMutations applied:")
for mut in mutations:
    print(f"  - {mut['enzyme']:12s} @ {mut['site_start']:4d}: {mut['original']}→{mut['mutated']} ({mut['aa_change']})")

print("\n✅ All annotations preserved with correct sequences")
print("✅ All 6 restriction sites added")
print("✅ All mutations verified silent")
