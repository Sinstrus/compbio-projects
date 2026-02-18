#!/usr/bin/env python3
"""
Verify exact mutation positions for all 6 restriction sites
"""

import sys
sys.path.insert(0, 'scripts/tools')
from Bio import SeqIO

# Load plasmid
record = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")
plasmid_seq = str(record.seq).upper()

# Define the 6 restriction sites
sites = [
    {
        'name': 'Region 1: SmaI',
        'enzyme': 'SmaI',
        'position': 2505,
        'recognition': 'CCCGGG',
    },
    {
        'name': 'Region 2: BbvCI',
        'enzyme': 'BbvCI',
        'position': 2828,
        'recognition': 'CCTCAGC',
    },
    {
        'name': 'Region 3: AgeI',
        'enzyme': 'AgeI',
        'position': 3583,
        'recognition': 'ACCGGT',
    },
    {
        'name': 'Region 4: BsrGI',
        'enzyme': 'BsrGI',
        'position': 3780,
        'recognition': 'TGTACA',
    },
    {
        'name': 'Region 5: BmtI/NheI',
        'enzyme': 'BmtI',
        'position': 3937,
        'recognition': 'GCTAGC',
    },
    {
        'name': 'Region 6: BstZ17I',
        'enzyme': 'BstZ17I',
        'position': 4163,
        'recognition': 'GTATAC',
    },
]

print("="*80)
print("VERIFYING EXACT MUTATION POSITIONS FOR ALL 6 SITES")
print("="*80)

mutations = []

for site in sites:
    pos = site['position']
    target = site['recognition']
    site_len = len(target)

    # Get current sequence at this position
    current = plasmid_seq[pos - 1:pos - 1 + site_len]

    print(f"\n{site['name']}:")
    print(f"  Position: {pos}")
    print(f"  Target:   {target}")
    print(f"  Current:  {current}")

    # Find differences
    if current == target:
        print(f"  ✅ Site already exists (no mutation needed)")
        mutations.append({
            'enzyme': site['enzyme'],
            'position': pos,
            'recognition': target,
            'mutations': []
        })
    else:
        print(f"  Mutations needed:")
        site_mutations = []
        for i, (c, t) in enumerate(zip(current, target)):
            if c != t:
                abs_pos = pos + i
                print(f"    Position {abs_pos} (site position {i+1}): {c} → {t}")
                site_mutations.append({
                    'abs_pos': abs_pos,
                    'site_rel_pos': i + 1,
                    'original': c,
                    'mutated': t
                })

        mutations.append({
            'enzyme': site['enzyme'],
            'position': pos,
            'recognition': target,
            'mutations': site_mutations
        })

# Print summary
print("\n" + "="*80)
print("SUMMARY OF ALL MUTATIONS")
print("="*80)

for mut_info in mutations:
    print(f"\n{mut_info['enzyme']} @ {mut_info['position']}:")
    if len(mut_info['mutations']) == 0:
        print(f"  No mutations needed (site exists)")
    else:
        for mut in mut_info['mutations']:
            print(f"  Position {mut['abs_pos']}: {mut['original']} → {mut['mutated']}")

# Export for use in other script
import json
with open('/tmp/mutation_positions.json', 'w') as f:
    json.dump(mutations, f, indent=2)

print(f"\n✅ Mutation details saved to /tmp/mutation_positions.json")
