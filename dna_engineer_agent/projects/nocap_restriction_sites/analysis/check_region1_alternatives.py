#!/usr/bin/env python3
"""
Check top alternatives for Region 1 that don't conflict with BbvCI recognition sequence
"""

import sys
sys.path.insert(0, 'scripts/tools')
import pickle

# Load results
with open('/tmp/region1_revised_results.pkl', 'rb') as f:
    results = pickle.load(f)

candidates = results['all_candidates']

# BbvCI family enzymes all recognize CCTCAGC
bbvci_family = {'BbvCI', 'Nb.BbvCI', 'Nt.BbvCI'}

# Filter out BbvCI family
non_bbvci = [c for c in candidates if c.enzyme not in bbvci_family]

print(f"{'='*80}")
print("TOP 20 CANDIDATES (excluding BbvCI family)")
print(f"{'='*80}")
print(f"Total candidates: {len(non_bbvci)}")
print()

for i, c in enumerate(non_bbvci[:20], 1):
    print(f"{i:2d}. {c.enzyme:15s} @ {c.position:5d} | "
          f"{c.edits_required} edit | {c.site_sequence:20s} | "
          f"Mutations: {', '.join(c.mutations) if c.mutations else 'None'}")

# Recommend AvrII
if len(non_bbvci) > 0:
    recommended = non_bbvci[0]

    print(f"\n{'='*80}")
    print("REVISED RECOMMENDATION FOR REGION 1:")
    print(f"{'='*80}")
    print(f"Enzyme:      {recommended.enzyme}")
    print(f"Position:    {recommended.position}")
    print(f"Recognition: {recommended.site_sequence}")
    print(f"Mutations:   {recommended.edits_required} ({', '.join(recommended.mutations)})")
    print(f"Uniqueness:  {recommended.uniqueness_dna}")

    print(f"\n{'='*80}")
    print("FINAL UPDATED SET OF 6 RESTRICTION SITES:")
    print(f"{'='*80}")
    print(f"Region 1 (Rep68-stop to VP2-start): {recommended.enzyme:10s} @ {recommended.position:5d} | {recommended.site_sequence}")
    print(f"Region 2 (VP2-AAP Intergenic):      BbvCI      @  2828 | CCTCAGC")
    print(f"Region 3 (AAP-stop to VR4):         EcoNI      @  3487 | CCTCAGTAAGG")
    print(f"Region 4 (VR4 to VR5):              FseI       @  3759 | GGCCGGCC")
    print(f"Region 5 (VR5 to VR8):              BmtI       @  3937 | GCTAGC")
    print(f"Region 6 (Post-VR8):                BaeI       @  4318 | ACACCTGTACC")
