#!/usr/bin/env python3
"""
Re-analyze Region 6 with NARROWED boundaries: 200 bases after VR8
VR8 ends at 4143, so Region 6 is now 4144-4344 (200 bp instead of 432 bp)
"""

import sys
sys.path.insert(0, 'scripts/tools')
from Bio import SeqIO
from silent_sites import translate
from silent_sites_curated import find_candidates_curated as find_candidates, RESTRICTION_ENZYMES
from dataclasses import dataclass
from typing import List, Optional

def calc_frame_offset(region_start: int, vp1_start: int = 2365) -> int:
    """Calculate correct frame offset to maintain VP1 reading frame"""
    position_in_codon = (region_start - vp1_start) % 3
    return (3 - position_in_codon) % 3 if position_in_codon != 0 else 0

@dataclass
class Boundary:
    """Define a critical boundary that sites cannot overlap"""
    name: str
    position: int
    type: str

# Load plasmid
print("="*80)
print("RE-ANALYZING REGION 6 WITH NARROWED BOUNDARIES")
print("="*80)
print("New Region 6: 200 bases after VR8 (4144-4344)")
print("="*80)

record = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")
plasmid_seq = str(record.seq).upper()

# Extract VR8 coordinates
vr8_coords = None
for feature in record.features:
    if feature.type == "misc_feature":
        label = feature.qualifiers.get("label", [""])[0]
        if "AAV9_VR-VIII" in label:
            vr8_start = int(feature.location.start) + 1
            vr8_end = int(feature.location.end)
            vr8_coords = (vr8_start, vr8_end)
            print(f"VR8: {vr8_start}-{vr8_end}")
            break

# Define Region 6 with new boundaries
region_start = vr8_coords[1] + 1  # 4144
region_end = region_start + 199    # 4344 (200 bases total)

print(f"Region 6: {region_start}-{region_end} ({region_end - region_start + 1} bp)")

# Extract sequence
region_dna = plasmid_seq[region_start - 1:region_end]
frame_offset = calc_frame_offset(region_start)
region_protein = translate(region_dna, frame=frame_offset, stop_at_terminator=False)

print(f"Frame offset: {frame_offset}")
print(f"DNA length: {len(region_dna)} bp")
print(f"Protein length: {len(region_protein)} aa")

# Define critical boundaries (none in this region, but keep for consistency)
boundaries = [
    Boundary("VR8_end", vr8_coords[1], "both"),
]

def check_boundary_overlap(site_position: int, site_length: int, boundaries: List[Boundary]) -> Optional[str]:
    """Check if a restriction site overlaps any critical boundary"""
    site_start = site_position
    site_end = site_position + site_length - 1

    for boundary in boundaries:
        if boundary.position >= site_start and boundary.position <= site_end:
            return f"Overlaps {boundary.name} at position {boundary.position}"

    return None

# Run analysis
print("\n" + "="*80)
print("SEARCHING FOR RESTRICTION SITES")
print("="*80)

candidates = find_candidates(
    dna_seq=plasmid_seq,
    protein_seq=region_protein,
    max_mutations=2,
    min_length=6,
    roi_seq=region_dna
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
    print("\n⚠️  WARNING: No valid sites found in narrowed Region 6!")
    print("    May need to expand region or relax constraints")
else:
    # Sort by edits required
    valid.sort(key=lambda c: c.edits_required)

    print(f"\n" + "="*80)
    print(f"TOP CANDIDATES (sorted by fewest mutations)")
    print("="*80)
    for i, c in enumerate(valid[:10], 1):
        print(f"{i:2d}. {c.enzyme:15s} @ position {c.position:5d} | "
              f"{c.edits_required} edit(s) | {c.site_sequence:20s} | "
              f"Mutations: {', '.join(c.mutations) if c.mutations else 'None (natural site)'}")

    # Show recommended site
    rec = valid[0]
    print(f"\n" + "="*80)
    print("RECOMMENDED SITE FOR NARROWED REGION 6")
    print("="*80)
    print(f"Enzyme: {rec.enzyme}")
    print(f"Position: {rec.position}")
    print(f"Recognition: {rec.site_sequence}")
    print(f"Mutations needed: {rec.edits_required}")
    if rec.mutations:
        print(f"Mutations: {', '.join(rec.mutations)}")
    else:
        print(f"Mutations: None (natural site)")

    # Manual verification
    print(f"\n" + "="*80)
    print("MANUAL VERIFICATION")
    print("="*80)

    # Get the exact position in region
    pos_in_region = rec.position - region_start
    site_len = len(rec.site_sequence)

    # Show original context
    context_start = max(0, pos_in_region - 15)
    context_end = min(len(region_dna), pos_in_region + site_len + 15)
    context = region_dna[context_start:context_end]

    print(f"Original context (position {region_start + context_start}-{region_start + context_end - 1}):")
    print(f"  {context}")
    print(f"  {' '*15}{'▼'*site_len}")
    print(f"  Site at position {rec.position}: {region_dna[pos_in_region:pos_in_region+site_len]}")

    if rec.edits_required > 0:
        # Apply mutations
        mutated_region = list(region_dna)
        for mut in rec.mutations:
            # Parse mutation like "A2T" (position relative to site, 1-indexed)
            old_base = mut[0]
            site_relative_pos = int(mut[1:-1])  # 1-indexed position in site
            new_base = mut[-1]

            # Convert to position in region (0-indexed)
            pos_in_region_mut = pos_in_region + (site_relative_pos - 1)

            # Absolute position in plasmid
            abs_pos = region_start + pos_in_region_mut

            print(f"\nMutation: {mut}")
            print(f"  Position in site (1-indexed): {site_relative_pos}")
            print(f"  Absolute position in plasmid: {abs_pos}")
            print(f"  Position in region (0-indexed): {pos_in_region_mut}")
            print(f"  Original base: {mutated_region[pos_in_region_mut]} (expected {old_base})")

            mutated_region[pos_in_region_mut] = new_base

        mutated_dna = ''.join(mutated_region)

        # Translate original and mutated
        original_protein = translate(region_dna, frame=frame_offset, stop_at_terminator=False)
        mutated_protein = translate(mutated_dna, frame=frame_offset, stop_at_terminator=False)

        print(f"\nProtein translation check:")
        print(f"  Original: {original_protein[:50]}..." if len(original_protein) > 50 else f"  Original: {original_protein}")
        print(f"  Mutated:  {mutated_protein[:50]}..." if len(mutated_protein) > 50 else f"  Mutated:  {mutated_protein}")

        if original_protein == mutated_protein:
            print(f"\n✅ VERIFIED SILENT: Protein sequence unchanged")
        else:
            print(f"\n❌ NOT SILENT: Protein sequence changed!")
            # Find differences
            for i, (o, m) in enumerate(zip(original_protein, mutated_protein)):
                if o != m:
                    print(f"   Position {i}: {o} → {m}")
    else:
        print(f"\n✅ Natural site - no mutations needed")
