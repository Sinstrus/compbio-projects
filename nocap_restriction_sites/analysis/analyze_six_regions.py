#!/usr/bin/env python3
"""
Comprehensive restriction site analysis for 6 regions in AAV9 Rep-Cap plasmid
with strict boundary overlap checking
"""

import sys
sys.path.insert(0, 'scripts/tools')
from Bio import SeqIO
from silent_sites import find_candidates, translate
from dataclasses import dataclass
from typing import List, Optional, Tuple

@dataclass
class Region:
    """Define a region with clear boundaries"""
    name: str
    start: int  # 1-indexed, inclusive
    end: int    # 1-indexed, inclusive
    description: str
    dna_seq: str = ""
    protein_seq: str = ""
    frame_offset: int = 0  # offset from VP1 start for correct translation

@dataclass
class Boundary:
    """Define a critical boundary that sites cannot overlap"""
    name: str
    position: int  # 1-indexed
    type: str  # "start", "end", or "both"

# Load plasmid
print("Loading GenBank file...")
record = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")
plasmid_seq = str(record.seq).upper()
print(f"Plasmid length: {len(plasmid_seq)} bp\n")

# Extract VR coordinates from GenBank annotations
print("Extracting Variable Region coordinates from GenBank annotations...")
vr_coords = {}
for feature in record.features:
    if feature.type == "misc_feature":
        label = feature.qualifiers.get("label", [""])[0]
        if "AAV9_VR-" in label:
            vr_name = label.replace("AAV9_", "").replace("-", "")
            start = int(feature.location.start) + 1  # Convert to 1-indexed
            end = int(feature.location.end)
            vr_coords[vr_name] = (start, end)
            print(f"  {vr_name:8s}: {start}-{end} ({end-start+1} bp)")

# Define all critical boundaries
print("\nDefining critical boundaries (sites cannot overlap these)...")
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

for b in boundaries:
    print(f"  {b.name:15s}: position {b.position}")

# Define the 6 regions
print("\nDefining target regions...")
regions = [
    Region(
        name="Region 1: VP1 Unique",
        start=2365,
        end=2775,
        description="VP1 unique region (VP1 start to VP2 start)",
        frame_offset=0
    ),
    Region(
        name="Region 2: VP2-AAP Intergenic",
        start=2779,
        end=2890,
        description="After VP2 start, before AAP start",
        frame_offset=(2779 - 2365) % 3
    ),
    Region(
        name="Region 3: AAP-stop to VR4",
        start=3485,
        end=vr_coords["VRIV"][0] - 1,  # 3717
        description="Between AAP stop and VR4 start",
        frame_offset=(3485 - 2365) % 3
    ),
    Region(
        name="Region 4: VR4 to VR5",
        start=vr_coords["VRIV"][1] + 1,  # 3745
        end=vr_coords["VRV"][0] - 1,     # 3825
        description="Between VR4 end and VR5 start",
        frame_offset=(vr_coords["VRIV"][1] + 1 - 2365) % 3
    ),
    Region(
        name="Region 5: VR5 to VR8",
        start=vr_coords["VRV"][1] + 1,    # 3880
        end=vr_coords["VRVIII"][0] - 1,   # 4104
        description="Between VR5 end and VR8 start",
        frame_offset=(vr_coords["VRV"][1] + 1 - 2365) % 3
    ),
    Region(
        name="Region 6: Post-VR8",
        start=vr_coords["VRVIII"][1] + 1,  # 4144
        end=4575,
        description="After VR8 end (as close to VR8 as possible)",
        frame_offset=(vr_coords["VRVIII"][1] + 1 - 2365) % 3
    ),
]

# Extract sequences for each region
print("\nExtracting sequences for each region...")
for region in regions:
    # Extract DNA (0-indexed slicing)
    region.dna_seq = plasmid_seq[region.start - 1:region.end]

    # Translate in correct frame
    # All regions are within VP1/VP2/VP3 coding sequence
    # VP1 starts at 2365, so we need to maintain frame relative to that
    region.protein_seq = translate(region.dna_seq, frame=region.frame_offset, stop_at_terminator=False)

    print(f"\n{region.name}")
    print(f"  Coordinates: {region.start}-{region.end} ({len(region.dna_seq)} bp, {len(region.protein_seq)} aa)")
    print(f"  Frame offset from VP1: {region.frame_offset}")
    print(f"  Description: {region.description}")

def check_boundary_overlap(site_position: int, site_length: int, boundaries: List[Boundary]) -> Optional[str]:
    """
    Check if a restriction site overlaps any critical boundary

    Args:
        site_position: 1-indexed start position of site
        site_length: length of recognition sequence
        boundaries: list of Boundary objects

    Returns:
        None if no overlap, or a string describing the overlap
    """
    site_start = site_position
    site_end = site_position + site_length - 1

    for boundary in boundaries:
        if boundary.position >= site_start and boundary.position <= site_end:
            return f"Overlaps {boundary.name} at position {boundary.position}"

    return None

def analyze_region(region: Region, plasmid_seq: str, boundaries: List[Boundary]):
    """Analyze a single region for restriction sites"""

    print(f"\n{'='*80}")
    print(f"ANALYZING: {region.name}")
    print(f"{'='*80}")

    # Run silent_sites.py
    candidates = find_candidates(
        dna_seq=plasmid_seq,
        protein_seq=region.protein_seq,
        max_mutations=2,
        min_length=6,
        roi_seq=region.dna_seq
    )

    print(f"\nFound {len(candidates)} total candidates")

    # Filter for silent mutations
    silent = [c for c in candidates if c.mutation_type == "Silent"]
    print(f"Silent mutations only: {len(silent)} candidates")

    # Filter for unique in full plasmid
    unique = [c for c in silent if c.uniqueness_dna == "Unique"]
    print(f"Unique in full plasmid: {len(unique)} candidates")

    # Check for boundary overlaps
    valid = []
    for c in unique:
        overlap = check_boundary_overlap(c.position, len(c.site_sequence), boundaries)
        if overlap is None:
            valid.append(c)

    print(f"No boundary overlaps: {len(valid)} candidates")

    if len(valid) == 0:
        print(f"\n⚠️  WARNING: No valid sites found for {region.name}")
        print(f"   All candidates either overlap boundaries or are non-unique")
        return None

    # Sort by edits required
    valid.sort(key=lambda c: c.edits_required)

    # Show top 5
    print(f"\nTop candidates (sorted by fewest mutations):")
    for i, c in enumerate(valid[:5], 1):
        print(f"{i}. {c.enzyme:15s} @ position {c.position:5d} | "
              f"{c.edits_required} edit(s) | {c.site_sequence:20s} | "
              f"Mutations: {', '.join(c.mutations) if c.mutations else 'None'}")

    return valid

# Analyze all regions
print(f"\n{'='*80}")
print("RUNNING ANALYSIS FOR ALL 6 REGIONS")
print(f"{'='*80}")

results = {}
for region in regions:
    valid_sites = analyze_region(region, plasmid_seq, boundaries)
    if valid_sites:
        results[region.name] = {
            'region': region,
            'candidates': valid_sites,
            'recommended': valid_sites[0]  # Top candidate
        }
    else:
        results[region.name] = {
            'region': region,
            'candidates': [],
            'recommended': None
        }

# Generate summary
print(f"\n{'='*80}")
print("SUMMARY OF RECOMMENDED SITES")
print(f"{'='*80}")

for region_name, data in results.items():
    region = data['region']
    rec = data['recommended']

    print(f"\n{region_name}")
    print(f"  Coordinates: {region.start}-{region.end}")
    if rec:
        print(f"  ✓ Recommended: {rec.enzyme} @ position {rec.position}")
        print(f"    Recognition: {rec.site_sequence}")
        print(f"    Mutations: {rec.edits_required} ({', '.join(rec.mutations) if rec.mutations else 'None'})")
    else:
        print(f"  ✗ No valid sites found")

# Save results for detailed report generation
import pickle
with open('/tmp/six_region_results.pkl', 'wb') as f:
    pickle.dump(results, f)

print(f"\n{'='*80}")
print("Analysis complete. Results saved to /tmp/six_region_results.pkl")
print("Run generate_detailed_report.py to create the final markdown report.")
print(f"{'='*80}")
