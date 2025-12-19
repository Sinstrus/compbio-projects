#!/usr/bin/env python3
"""
Debug script to investigate AAP frame issues
"""

from Bio import SeqIO
from Bio.Seq import Seq

INPUT_FILE = "/home/cnguy/projects/dna_engineer_agent/test_data/TTRC004_rep2mut02-cap9-p5.gb"
record = SeqIO.read(INPUT_FILE, "genbank")
seq = record.seq

# Find VP1 and AAP annotations
vp1_feature = None
aap_feature = None

for feat in record.features:
    if feat.type == "CDS":
        label = feat.qualifiers.get('label', [''])[0]
        if label == 'VP1':
            vp1_feature = feat
        elif label == 'AAP':
            aap_feature = feat

print("=" * 80)
print("AAP FRAME DEBUG")
print("=" * 80)

if vp1_feature:
    vp1_start = int(vp1_feature.location.start)  # 0-based
    vp1_end = int(vp1_feature.location.end)
    print(f"VP1: {vp1_start+1}-{vp1_end} bp (0-based: {vp1_start}-{vp1_end})")

    # Get VP1 DNA sequence
    vp1_dna = seq[vp1_start:vp1_end]
    print(f"VP1 length: {len(vp1_dna)} bp ({len(vp1_dna)//3} codons)")

    # Translate in all 3 frames
    print("\n--- VP1 Frame 0 (VP frame) ---")
    vp_protein = vp1_dna.translate(to_stop=False)
    print(f"First 50aa: {vp_protein[:50]}")

    print("\n--- VP1 Frame +1 (AAP frame) ---")
    aap_frame_dna = vp1_dna[1:]
    aap_frame_protein = aap_frame_dna.translate(to_stop=False)
    print(f"First 50aa: {aap_frame_protein[:50]}")

    print("\n--- VP1 Frame +2 ---")
    frame2_dna = vp1_dna[2:]
    frame2_protein = frame2_dna.translate(to_stop=False)
    print(f"First 50aa: {frame2_protein[:50]}")

if aap_feature:
    aap_start = int(aap_feature.location.start)  # 0-based
    aap_end = int(aap_feature.location.end)
    print(f"\n\nAAP annotation: {aap_start+1}-{aap_end} bp (0-based: {aap_start}-{aap_end})")
    print(f"AAP annotation translation from file: {aap_feature.qualifiers.get('translation', [''])[0][:50]}")

    # Extract AAP DNA and translate
    aap_dna = seq[aap_start:aap_end]
    aap_protein = aap_dna.translate(to_stop=False)
    print(f"AAP direct translation: {aap_protein[:50]}")

    # Check what frame AAP is relative to VP1
    offset_from_vp1 = aap_start - vp1_start
    frame = offset_from_vp1 % 3
    print(f"\nAAP offset from VP1 start: {offset_from_vp1} bp")
    print(f"AAP frame relative to VP1: +{frame}")

    # Now translate from VP1 in that frame
    if frame == 0:
        print("ERROR: AAP should not be in frame 0 (that's the VP frame)")
    elif frame == 1:
        print("\n--- Verifying AAP is in VP1 +1 frame ---")
        # Calculate position in the +1 frame
        aa_position_in_frame = offset_from_vp1 // 3
        print(f"AAP starts at codon position {aa_position_in_frame} in the +1 frame")

        # Translate +1 frame starting from VP1
        aap_frame_dna = vp1_dna[1:]
        aap_frame_protein = aap_frame_dna.translate(to_stop=False)

        # Extract the AAP portion
        aap_length_aa = (aap_end - aap_start) // 3
        extracted_aap = str(aap_frame_protein)[aa_position_in_frame:aa_position_in_frame + aap_length_aa]
        print(f"Extracted AAP from +1 frame: {extracted_aap[:50]}...")

        # Compare with annotation
        annotated_aap = aap_feature.qualifiers.get('translation', [''])[0]
        if extracted_aap == annotated_aap:
            print("✅ Match! Annotation is correct.")
        else:
            print("❌ Mismatch! Annotation doesn't match +1 frame translation.")
            print(f"Annotated: {annotated_aap[:50]}")
            print(f"Extracted: {extracted_aap[:50]}")
    elif frame == 2:
        print("AAP is in +2 frame relative to VP1")
        aa_position_in_frame = offset_from_vp1 // 3
        aap_frame_dna = vp1_dna[2:]
        aap_frame_protein = aap_frame_dna.translate(to_stop=False)
        extracted_aap = str(aap_frame_protein)[aa_position_in_frame:aa_position_in_frame + 50]
        print(f"Extracted AAP from +2 frame: {extracted_aap}")

# Search for known AAP motifs
print("\n\n" + "=" * 80)
print("SEARCHING FOR AAP START MOTIFS")
print("=" * 80)

motifs_to_search = [
    "MPGFYEIVIKVPSD",
    "MPGFYEIVIKV",
    "MATQSQSQTL",  # From the annotation in the file
    "MATQ",
]

if vp1_feature:
    vp1_dna = seq[vp1_start:vp1_end]

    for frame in [0, 1, 2]:
        frame_dna = vp1_dna[frame:]
        frame_protein = str(frame_dna.translate(to_stop=False))

        print(f"\n--- Frame +{frame} ---")
        for motif in motifs_to_search:
            pos = frame_protein.find(motif)
            if pos != -1:
                print(f"✅ Found '{motif}' at position {pos} aa (context: ...{frame_protein[max(0,pos-5):pos+len(motif)+5]}...)")
                # Calculate bp position in genome
                bp_in_genome = vp1_start + frame + (pos * 3) + 1  # +1 for 1-based
                print(f"   Genome position: {bp_in_genome} bp")
            else:
                print(f"❌ Not found: '{motif}'")

print("\n" + "=" * 80)
