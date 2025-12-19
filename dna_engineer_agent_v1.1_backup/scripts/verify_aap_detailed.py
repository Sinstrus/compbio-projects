#!/usr/bin/env python3
"""
Detailed AAP verification - check actual start codon
"""

from Bio import SeqIO
from Bio.Seq import Seq

INPUT_FILE = "/home/cnguy/projects/dna_engineer_agent/test_data/TTRC004_rep2mut02-cap9-p5.gb"
record = SeqIO.read(INPUT_FILE, "genbank")
seq = record.seq

# Find annotations
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
print("DETAILED AAP START ANALYSIS")
print("=" * 80)

if aap_feature and vp1_feature:
    # AAP annotation coordinates
    aap_start_0based = int(aap_feature.location.start)
    aap_end_0based = int(aap_feature.location.end)
    aap_start_1based = aap_start_0based + 1

    print(f"AAP annotation: {aap_start_1based}-{aap_end_0based} bp (GenBank 1-based)")
    print(f"AAP annotation: {aap_start_0based}-{aap_end_0based} bp (Python 0-based)")

    # Extract AAP DNA from annotation
    aap_dna = seq[aap_start_0based:aap_end_0based]
    print(f"\nAAP DNA (first 30 bp): {aap_dna[:30]}")
    print(f"First codon: {aap_dna[:3]}")
    print(f"First codon translates to: {aap_dna[:3].translate()}")

    # Translate AAP
    aap_protein_direct = aap_dna.translate(to_stop=False)
    print(f"AAP direct translation (first 30aa): {aap_protein_direct[:30]}")

    # What does the annotation say?
    aap_annotated_protein = aap_feature.qualifiers.get('translation', [''])[0]
    print(f"AAP annotation protein (first 30aa): {aap_annotated_protein[:30]}")

    # Check if annotation has extra M
    if aap_annotated_protein[0] == 'M' and str(aap_protein_direct)[0] == 'L':
        print("\n⚠️  ANNOTATION ISSUE: Annotation starts with M, but actual DNA translates to L")
        print("   This suggests the annotation includes a non-existent start codon")

        # Check if there's an ATG upstream
        upstream_seq = seq[aap_start_0based-10:aap_start_0based+10]
        print(f"   Context (10bp before and after): {upstream_seq}")

        # Look for ATG in the region
        search_region = seq[aap_start_0based-30:aap_start_0based]
        atg_positions = []
        for i in range(len(search_region) - 2):
            if str(search_region[i:i+3]) == 'ATG':
                atg_positions.append(i)

        if atg_positions:
            print(f"   Found ATG codons at positions (relative to AAP start): {[-30+i for i in atg_positions]}")

    # VP1 information
    vp1_start_0based = int(vp1_feature.location.start)
    vp1_end_0based = int(vp1_feature.location.end)

    # AAP offset from VP1
    offset = aap_start_0based - vp1_start_0based
    frame = offset % 3

    print(f"\n--- Relationship to VP1 ---")
    print(f"VP1 starts at: {vp1_start_0based} (0-based)")
    print(f"AAP starts at: {aap_start_0based} (0-based)")
    print(f"Offset: {offset} bp")
    print(f"Frame: +{frame}")

    # Extract from VP1 +1 frame
    vp1_dna = seq[vp1_start_0based:vp1_end_0based]
    vp1_plus1_dna = vp1_dna[1:]
    vp1_plus1_protein = vp1_plus1_dna.translate(to_stop=False)

    # AAP position in +1 frame
    aap_aa_position = offset // 3
    print(f"AAP starts at codon {aap_aa_position} in VP1 +1 frame")

    # Extract AAP from +1 frame
    aap_from_plus1 = str(vp1_plus1_protein)[aap_aa_position:aap_aa_position+50]
    print(f"AAP from VP1 +1 frame: {aap_from_plus1}")

    # Find the actual start by looking for stop codon before AAP
    plus1_before_aap = str(vp1_plus1_protein)[:aap_aa_position+10]
    print(f"\nVP1 +1 frame before AAP (last 50aa): ...{plus1_before_aap[-50:]}")

    # Search for last stop codon before AAP
    last_stop = -1
    for i in range(aap_aa_position - 1, -1, -1):
        if vp1_plus1_protein[i] == '*':
            last_stop = i
            break

    if last_stop != -1:
        print(f"Last stop codon in +1 frame before AAP: position {last_stop}")
        print(f"Distance from stop to AAP: {aap_aa_position - last_stop - 1} aa")

        # The true AAP start should be right after this stop
        true_aap_start_aa = last_stop + 1
        true_aap_start_bp_0based = vp1_start_0based + 1 + (true_aap_start_aa * 3)
        true_aap_start_bp_1based = true_aap_start_bp_0based + 1

        print(f"\n--- VERIFIED AAP START ---")
        print(f"True AAP start (after stop codon): position {true_aap_start_aa} in +1 frame")
        print(f"True AAP start (bp, 1-based): {true_aap_start_bp_1based}")
        print(f"Annotated AAP start (bp, 1-based): {aap_start_1based}")
        print(f"Difference: {true_aap_start_bp_1based - aap_start_1based} bp")

        # Extract true AAP sequence
        true_aap_protein = str(vp1_plus1_protein)[true_aap_start_aa:]

        # Find stop codon
        aap_stop_pos = true_aap_protein.find('*')
        if aap_stop_pos != -1:
            true_aap_protein = true_aap_protein[:aap_stop_pos]
            true_aap_end_aa = true_aap_start_aa + aap_stop_pos
            true_aap_end_bp_0based = vp1_start_0based + 1 + (true_aap_end_aa * 3) + 3  # include stop codon
            true_aap_end_bp_1based = true_aap_end_bp_0based

            print(f"True AAP end (bp, 1-based): {true_aap_end_bp_1based}")
            print(f"True AAP length: {len(true_aap_protein)} aa")
            print(f"True AAP sequence (first 50aa): {true_aap_protein[:50]}")
        else:
            print("WARNING: No stop codon found for AAP")
    else:
        print("No stop codon found before AAP in +1 frame")

print("\n" + "=" * 80)
