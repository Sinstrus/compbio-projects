#!/usr/bin/env python3
"""
Analyze restriction sites in AAV9 Rep-Cap plasmid
Identify candidate sites in three specific regions:
1. VP1 unique region (VP1 start to VP2 start)
2. VP2-AAP intergenic (VP2 start to AAP start)
3. VP3 post-AAP (VP3 region after AAP stop)
"""

import sys
sys.path.insert(0, 'scripts/tools')
from Bio import SeqIO
from silent_sites import find_candidates, translate

# Load the GenBank file
print("Loading GenBank file...")
record = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")
plasmid_seq = str(record.seq)

print(f"Plasmid length: {len(plasmid_seq)} bp\n")

# Define regions based on GenBank annotations
# VP1: 2365-4575
# VP2: join(2776..2778, 2779..4575)
# VP3: 2971-4575
# AAP: 2891-3484

print("="*80)
print("REGION 1: VP1 UNIQUE REGION")
print("="*80)
print("Definition: Between VP1 start and VP2 start")
print("VP1 starts at position 2365 (1-indexed)")
print("VP2 starts at position 2776 (1-indexed)")

# Convert to 0-indexed for Python slicing
vp1_start_0idx = 2365 - 1  # 2364
vp2_start_0idx = 2776 - 1  # 2775

region1_dna = plasmid_seq[vp1_start_0idx:vp2_start_0idx]
region1_protein = translate(region1_dna, frame=0, stop_at_terminator=False)

print(f"Region coordinates (1-indexed): {vp1_start_0idx + 1} to {vp2_start_0idx}")
print(f"Region length: {len(region1_dna)} bp ({len(region1_protein)} aa)")
print(f"DNA sequence: {region1_dna[:60]}...{region1_dna[-60:]}")
print(f"Protein sequence: {region1_protein[:30]}...{region1_protein[-30:]}")

print("\n" + "="*80)
print("REGION 2: VP2-AAP INTERGENIC")
print("="*80)
print("Definition: After VP2 start but before AAP start codon")
print("VP2 starts at position 2776-2778 (ACG codon), then continues at 2779")
print("AAP starts at position 2891 (1-indexed)")

# VP2 frame continues after the join at 2779
# We want from 2779 to just before AAP start (2891)
region2_start_0idx = 2779 - 1  # 2778
aap_start_0idx = 2891 - 1      # 2890

region2_dna = plasmid_seq[region2_start_0idx:aap_start_0idx]
# VP2 is in the same frame as VP1, starting from position 2776
# To get the correct protein frame, we need to consider the offset from VP1 start
# VP1 starts at 2365, so VP2 frame offset is (2779-2365) % 3 = 414 % 3 = 0
region2_protein = translate(region2_dna, frame=0, stop_at_terminator=False)

print(f"Region coordinates (1-indexed): {region2_start_0idx + 1} to {aap_start_0idx}")
print(f"Region length: {len(region2_dna)} bp ({len(region2_protein)} aa)")
print(f"DNA sequence: {region2_dna}")
print(f"Protein sequence: {region2_protein}")

print("\n" + "="*80)
print("REGION 3: VP3 POST-AAP")
print("="*80)
print("Definition: Within VP3 coding region but after AAP stop codon")
print("VP3 starts at position 2971 (1-indexed)")
print("AAP stops at position 3484 (1-indexed, includes stop codon)")
print("VP3 ends at position 4575 (1-indexed)")

vp3_start_0idx = 2971 - 1  # 2970
aap_stop_0idx = 3484       # 3483 (0-indexed, stop codon ends here)
vp3_end_0idx = 4575        # 4574 (0-indexed)

# Region starts after AAP stop
region3_start_0idx = aap_stop_0idx  # Position after AAP stop
region3_dna = plasmid_seq[region3_start_0idx:vp3_end_0idx]

# VP3 frame: VP1 starts at 2365, VP3 starts at 2971
# Frame offset: (2971 - 2365) % 3 = 606 % 3 = 0
# Region starts at 3485, offset from VP1: (3485 - 2365) % 3 = 1120 % 3 = 2
# So we need frame=2 for this region relative to position 3485
# Actually, let's calculate relative to VP3 start:
# VP3 starts at 2971 (0-indexed: 2970)
# Region starts at 3485 (0-indexed: 3484)
# Offset: 3484 - 2970 = 514, 514 % 3 = 1
region3_protein = translate(region3_dna, frame=1, stop_at_terminator=False)

print(f"Region coordinates (1-indexed): {region3_start_0idx + 1} to {vp3_end_0idx}")
print(f"Region length: {len(region3_dna)} bp ({len(region3_protein)} aa)")
print(f"DNA sequence: {region3_dna[:60]}...{region3_dna[-60:]}")
print(f"Protein sequence: {region3_protein[:30]}...{region3_protein[-30:]}")

print("\n" + "="*80)
print("RUNNING SILENT_SITES.PY FOR EACH REGION")
print("="*80)

# Save regions to files for analysis
with open("/tmp/region1_vp1_unique.txt", "w") as f:
    f.write(f"# VP1 Unique Region (positions {vp1_start_0idx + 1}-{vp2_start_0idx})\n")
    f.write(f"DNA: {region1_dna}\n")
    f.write(f"Protein: {region1_protein}\n")

with open("/tmp/region2_vp2_aap.txt", "w") as f:
    f.write(f"# VP2-AAP Intergenic (positions {region2_start_0idx + 1}-{aap_start_0idx})\n")
    f.write(f"DNA: {region2_dna}\n")
    f.write(f"Protein: {region2_protein}\n")

with open("/tmp/region3_vp3_post_aap.txt", "w") as f:
    f.write(f"# VP3 Post-AAP (positions {region3_start_0idx + 1}-{vp3_end_0idx})\n")
    f.write(f"DNA: {region3_dna}\n")
    f.write(f"Protein: {region3_protein}\n")

print("\nRegion information saved to /tmp/region*.txt files")
print("\nNow running silent_sites.py for each region...")

# Region 1: VP1 Unique
print("\n" + "="*80)
print("REGION 1 CANDIDATES: VP1 UNIQUE")
print("="*80)
candidates_r1 = find_candidates(
    dna_seq=plasmid_seq,
    protein_seq=region1_protein,
    max_mutations=2,
    min_length=6,
    roi_seq=region1_dna
)
print(f"\nFound {len(candidates_r1)} total candidates in Region 1")

# Filter for silent mutations only
silent_r1 = [c for c in candidates_r1 if c.mutation_type == "Silent"]
print(f"Silent mutations only: {len(silent_r1)} candidates")

# Filter for unique sites in full plasmid
unique_r1 = [c for c in silent_r1 if c.uniqueness_dna == "Unique"]
print(f"Unique in full plasmid: {len(unique_r1)} candidates")

# Sort by edits required (0 > 1 > 2)
unique_r1.sort(key=lambda c: c.edits_required)

print(f"\nTop candidates (sorted by fewest mutations):")
for i, c in enumerate(unique_r1[:10], 1):
    print(f"{i}. {c.enzyme:12s} @ position {c.position:5d} | "
          f"{c.edits_required} edit(s) | {c.site_sequence} | "
          f"Mutations: {', '.join(c.mutations) if c.mutations else 'None'}")

# Region 2: VP2-AAP Intergenic
print("\n" + "="*80)
print("REGION 2 CANDIDATES: VP2-AAP INTERGENIC")
print("="*80)
candidates_r2 = find_candidates(
    dna_seq=plasmid_seq,
    protein_seq=region2_protein,
    max_mutations=2,
    min_length=6,
    roi_seq=region2_dna
)
print(f"\nFound {len(candidates_r2)} total candidates in Region 2")

silent_r2 = [c for c in candidates_r2 if c.mutation_type == "Silent"]
print(f"Silent mutations only: {len(silent_r2)} candidates")

unique_r2 = [c for c in silent_r2 if c.uniqueness_dna == "Unique"]
print(f"Unique in full plasmid: {len(unique_r2)} candidates")

unique_r2.sort(key=lambda c: c.edits_required)

print(f"\nTop candidates (sorted by fewest mutations):")
for i, c in enumerate(unique_r2[:10], 1):
    print(f"{i}. {c.enzyme:12s} @ position {c.position:5d} | "
          f"{c.edits_required} edit(s) | {c.site_sequence} | "
          f"Mutations: {', '.join(c.mutations) if c.mutations else 'None'}")

# Region 3: VP3 Post-AAP
print("\n" + "="*80)
print("REGION 3 CANDIDATES: VP3 POST-AAP")
print("="*80)

# Wait, I need to be more careful about the frame here. Let me recalculate.
# The region3_dna needs to be translated in the correct frame.
# VP3 starts at 2971 (1-indexed), which is the ATG start codon
# Region starts at 3485 (1-indexed)
# From VP3 start to region start: 3485 - 2971 = 514 bp
# 514 % 3 = 1, so frame = 1

# But the protein sequence needs to match! Let me extract the actual VP3 protein
# from the region to verify
vp3_full_dna = plasmid_seq[vp3_start_0idx:vp3_end_0idx]
vp3_full_protein = translate(vp3_full_dna, frame=0, stop_at_terminator=False)

# Now get the portion after AAP
# AAP ends at position 3484 (1-indexed), which is position 3483 (0-indexed) + 1 = 3484
# Offset from VP3 start: 3484 - 2970 = 514 bp, which is 514/3 = 171 aa + 1 bp
# So the protein starts at position 172 (1-indexed) in VP3

offset_aa = (aap_stop_0idx - vp3_start_0idx) // 3
region3_protein_correct = vp3_full_protein[offset_aa:]

print(f"VP3 full protein length: {len(vp3_full_protein)} aa")
print(f"Region 3 starts at AA position {offset_aa + 1} of VP3")
print(f"Region 3 protein length: {len(region3_protein_correct)} aa")

candidates_r3 = find_candidates(
    dna_seq=plasmid_seq,
    protein_seq=region3_protein_correct,
    max_mutations=2,
    min_length=6,
    roi_seq=region3_dna
)
print(f"\nFound {len(candidates_r3)} total candidates in Region 3")

silent_r3 = [c for c in candidates_r3 if c.mutation_type == "Silent"]
print(f"Silent mutations only: {len(silent_r3)} candidates")

unique_r3 = [c for c in silent_r3 if c.uniqueness_dna == "Unique"]
print(f"Unique in full plasmid: {len(unique_r3)} candidates")

unique_r3.sort(key=lambda c: c.edits_required)

print(f"\nTop candidates (sorted by fewest mutations):")
for i, c in enumerate(unique_r3[:10], 1):
    print(f"{i}. {c.enzyme:12s} @ position {c.position:5d} | "
          f"{c.edits_required} edit(s) | {c.site_sequence} | "
          f"Mutations: {', '.join(c.mutations) if c.mutations else 'None'}")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
