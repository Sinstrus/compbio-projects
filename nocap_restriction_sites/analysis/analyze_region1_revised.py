#!/usr/bin/env python3
"""
Re-analyze Region 1 with narrowed boundaries:
- After Rep68 stop codon (position 2415)
- Before VP2 start codon (position 2776)
"""

import sys
sys.path.insert(0, 'scripts/tools')
from Bio import SeqIO
from silent_sites import find_candidates, translate
from dataclasses import dataclass
from typing import List, Optional

@dataclass
class Boundary:
    """Define a critical boundary that sites cannot overlap"""
    name: str
    position: int
    type: str

# Load plasmid
print("Loading GenBank file...")
record = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")
plasmid_seq = str(record.seq).upper()
print(f"Plasmid length: {len(plasmid_seq)} bp\n")

# Extract VR coordinates from GenBank annotations
vr_coords = {}
for feature in record.features:
    if feature.type == "misc_feature":
        label = feature.qualifiers.get("label", [""])[0]
        if "AAV9_VR-" in label:
            vr_name = label.replace("AAV9_", "").replace("-", "")
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            vr_coords[vr_name] = (start, end)

# Define critical boundaries
boundaries = [
    Boundary("VP1_start", 2365, "start"),
    Boundary("VP2_start", 2776, "start"),
    Boundary("AAP_start", 2891, "start"),
    Boundary("VP3_start", 2971, "start"),
    Boundary("AAP_stop", 3484, "end"),
    Boundary("VR4_start", vr_coords["VRIV"][0], "both"),
    Boundary("VR4_end", vr_coords["VRIV"][1], "both"),
    Boundary("VR5_start", vr_coords["VRV"][0], "both"),
    Boundary("VR5_end", vr_coords["VRV"][1], "both"),
    Boundary("VR8_start", vr_coords["VRVIII"][0], "both"),
    Boundary("VR8_end", vr_coords["VRVIII"][1], "both"),
]

# Previously selected enzymes for regions 2-6
# We need to ensure Region 1 doesn't conflict with these
previously_selected = {
    "Region 2": ("BbvCI", 2828, "CCTCAGC"),
    "Region 3": ("EcoNI", 3487, "CCTCAGTAAGG"),
    "Region 4": ("FseI", 3759, "GGCCGGCC"),
    "Region 5": ("BmtI", 3937, "GCTAGC"),
    "Region 6": ("BaeI", 4318, "ACACCTGTACC"),
}

print("Previously selected restriction sites:")
for region, (enzyme, pos, seq) in previously_selected.items():
    print(f"  {region}: {enzyme:10s} @ {pos:5d} ({seq})")

# NEW Region 1: Rep68 stop to VP2 start
region1_start = 2415  # After Rep68 stop at 2414
region1_end = 2775    # Before VP2 start at 2776
region1_length = region1_end - region1_start + 1

print(f"\n{'='*80}")
print(f"REVISED REGION 1 ANALYSIS")
print(f"{'='*80}")
print(f"New boundaries: {region1_start}-{region1_end} ({region1_length} bp)")
print(f"Previous boundaries: 2365-2775 (411 bp)")
print(f"Difference: Narrowed by {411 - region1_length} bp")

# Extract sequence
region1_dna = plasmid_seq[region1_start - 1:region1_end]

# Calculate frame offset from VP1 start (2365)
# CORRECTED: If region starts mid-codon, skip to next codon boundary
position_in_codon = (region1_start - 2365) % 3
frame_offset = (3 - position_in_codon) % 3 if position_in_codon != 0 else 0
region1_protein = translate(region1_dna, frame=frame_offset, stop_at_terminator=False)

print(f"\nRegion 1 protein length: {len(region1_protein)} aa")
print(f"Frame offset from VP1 start: {frame_offset}")

def check_boundary_overlap(site_position: int, site_length: int, boundaries: List[Boundary]) -> Optional[str]:
    """Check if a restriction site overlaps any critical boundary"""
    site_start = site_position
    site_end = site_position + site_length - 1

    for boundary in boundaries:
        if boundary.position >= site_start and boundary.position <= site_end:
            return f"Overlaps {boundary.name} at position {boundary.position}"

    return None

# Run silent_sites.py
print("\nSearching for restriction sites...")
candidates = find_candidates(
    dna_seq=plasmid_seq,
    protein_seq=region1_protein,
    max_mutations=2,
    min_length=6,
    roi_seq=region1_dna
)

print(f"Found {len(candidates)} total candidates")

# Filter for silent mutations
silent = [c for c in candidates if c.mutation_type == "Silent"]
print(f"Silent mutations only: {len(silent)} candidates")

# Filter for unique in full plasmid
unique = [c for c in silent if c.uniqueness_dna == "Unique"]
print(f"Unique in full plasmid: {len(unique)} candidates")

# Check for boundary overlaps
no_overlap = []
for c in unique:
    overlap = check_boundary_overlap(c.position, len(c.site_sequence), boundaries)
    if overlap is None:
        no_overlap.append(c)

print(f"No boundary overlaps: {len(no_overlap)} candidates")

# Additional filter: Check if enzyme is already used in other regions
# We want to avoid using the same enzyme name if possible
previously_used_enzymes = {enzyme for enzyme, _, _ in previously_selected.values()}

non_conflicting = []
for c in no_overlap:
    if c.enzyme not in previously_used_enzymes:
        non_conflicting.append(c)

print(f"Not conflicting with enzymes in Regions 2-6: {len(non_conflicting)} candidates")

if len(non_conflicting) == 0:
    print(f"\n⚠️  WARNING: All unique sites conflict with previously selected enzymes")
    print(f"   Will show all candidates including those with same enzyme names")
    final_candidates = no_overlap
else:
    final_candidates = non_conflicting

# Sort by edits required
final_candidates.sort(key=lambda c: c.edits_required)

# Show top 10
print(f"\n{'='*80}")
print("TOP 10 CANDIDATES (sorted by fewest mutations):")
print(f"{'='*80}")
for i, c in enumerate(final_candidates[:10], 1):
    conflict_marker = " ⚠️ CONFLICT" if c.enzyme in previously_used_enzymes else ""
    print(f"{i:2d}. {c.enzyme:15s} @ position {c.position:5d} | "
          f"{c.edits_required} edit(s) | {c.site_sequence:20s} | "
          f"Mutations: {', '.join(c.mutations) if c.mutations else 'None'}{conflict_marker}")

# Recommended site
if len(final_candidates) > 0:
    recommended = final_candidates[0]
    print(f"\n{'='*80}")
    print("RECOMMENDED SITE FOR REVISED REGION 1:")
    print(f"{'='*80}")
    print(f"Enzyme:      {recommended.enzyme}")
    print(f"Position:    {recommended.position}")
    print(f"Recognition: {recommended.site_sequence}")
    print(f"Mutations:   {recommended.edits_required} ({', '.join(recommended.mutations) if recommended.mutations else 'None'})")
    print(f"Uniqueness:  {recommended.uniqueness_dna}")

    # Show sequence context
    context_start = recommended.position - 30
    context_end = recommended.position + len(recommended.site_sequence) + 30
    context_seq = plasmid_seq[context_start - 1:context_end]

    print(f"\nSequence context (position {context_start}-{context_end}):")
    print(context_seq)
    print(" " * 30 + "^" * len(recommended.site_sequence))
    print(f" " * 30 + f"{recommended.enzyme} site @ {recommended.position}")

    # Check if this changes our complete set
    print(f"\n{'='*80}")
    print("UPDATED COMPLETE SET OF 6 RESTRICTION SITES:")
    print(f"{'='*80}")
    print(f"Region 1 (REVISED): {recommended.enzyme:10s} @ {recommended.position:5d} | {recommended.site_sequence}")
    for region, (enzyme, pos, seq) in previously_selected.items():
        print(f"{region:15s}: {enzyme:10s} @ {pos:5d} | {seq}")

    # Save results
    import pickle
    results = {
        'recommended': recommended,
        'all_candidates': final_candidates,
        'region_start': region1_start,
        'region_end': region1_end,
    }

    with open('/tmp/region1_revised_results.pkl', 'wb') as f:
        pickle.dump(results, f)

    print(f"\n{'='*80}")
    print("Analysis complete. Results saved to /tmp/region1_revised_results.pkl")
    print(f"{'='*80}")
else:
    print("\n❌ ERROR: No valid sites found for revised Region 1")
