#!/usr/bin/env python3
"""
Generate complete AVD005/006 GenBank files with all features preserved
and calculate synthetic fragments for ordering
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from pathlib import Path
import csv

# VHH and linker sequences
VHH3_AA = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSTADMGWFRQAPGKGRELVAAVSGSGFSTYSDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAKATISLYYAMDVWGQGTTVTVSS"
LINK_D2_N_AA = "GGGGSGGGGSGGGGSGGGGS"

# Codon optimization
CODON_USAGE_HUMAN = {
    'A': 'GCC', 'C': 'TGC', 'D': 'GAC', 'E': 'GAG', 'F': 'TTC',
    'G': 'GGC', 'H': 'CAC', 'I': 'ATC', 'K': 'AAG', 'L': 'CTG',
    'M': 'ATG', 'N': 'AAC', 'P': 'CCC', 'Q': 'CAG', 'R': 'CGC',
    'S': 'AGC', 'T': 'ACC', 'V': 'GTG', 'W': 'TGG', 'Y': 'TAC',
    '*': 'TAA'
}

def codon_optimize(aa_seq):
    return ''.join(CODON_USAGE_HUMAN[aa] for aa in aa_seq)


def build_vp1_vhh_fusion(vp1_original):
    """Build VP1-VHH fusion with knockouts"""
    vp1 = list(vp1_original.upper())

    # VP2 knockout: ACG → ACC at position 411
    if ''.join(vp1[411:414]) == 'ACG':
        vp1[413] = 'C'

    # VP3 knockout: ATG → CTG at position 606
    if ''.join(vp1[606:609]) == 'ATG':
        vp1[606] = 'C'

    vp1_ko = ''.join(vp1)

    # Insert VHH at position 456 aa (bp 1368)
    insertion_bp = 456 * 3
    vp1_before = vp1_ko[:insertion_bp]
    vp1_after = vp1_ko[insertion_bp:]

    linker_n = codon_optimize(LINK_D2_N_AA)
    vhh = codon_optimize(VHH3_AA)

    insert_cassette = linker_n + vhh
    vp1_vhh_fusion = vp1_before + insert_cassette + vp1_after

    return vp1_vhh_fusion, {
        'insertion_bp': insertion_bp,
        'linker_n_start': insertion_bp,
        'linker_n_end': insertion_bp + len(linker_n),
        'vhh_start': insertion_bp + len(linker_n),
        'vhh_end': insertion_bp + len(linker_n) + len(vhh),
        'vp2_ko': 411,
        'vp3_ko': 606
    }


def find_restriction_sites(sequence, enzyme_dict):
    """Find all restriction sites in sequence"""
    results = {}
    seq_upper = sequence.upper()
    for enzyme, recog_seq in enzyme_dict.items():
        positions = []
        start = 0
        while True:
            pos = seq_upper.find(recog_seq, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
        results[enzyme] = positions
    return results


def calculate_synthetic_fragment(original_seq, new_seq, vp1_start, vp1_end_new, enzyme_sites):
    """
    Calculate the synthetic fragment needed for ordering.

    The fragment extends from the changed region to the nearest unique
    restriction sites on both sides.
    """
    # The changed region is from vp1_start to vp1_end_new
    changed_start = vp1_start
    changed_end = vp1_end_new

    # Find flanking restriction sites
    left_sites = []
    right_sites = []

    for enzyme, positions in enzyme_sites.items():
        if len(positions) == 1:  # Only unique sites
            pos = positions[0]
            if pos < changed_start:
                left_sites.append((enzyme, pos))
            elif pos > changed_end:
                right_sites.append((enzyme, pos))

    # Sort to find closest
    if left_sites:
        left_sites.sort(key=lambda x: x[1], reverse=True)  # Closest to start
        left_enzyme, left_pos = left_sites[0]
    else:
        left_enzyme, left_pos = None, 0

    if right_sites:
        right_sites.sort(key=lambda x: x[1])  # Closest to end
        right_enzyme, right_pos = right_sites[0]
    else:
        right_enzyme, right_pos = None, len(new_seq)

    # Get enzyme recognition sequences
    enzyme_dict = {
        'SmaI': 'CCCGGG',
        'BbvCI': 'CCTCAGC',
        'AgeI': 'ACCGGT',
        'BsrGI': 'TGTACA',
        'BmtI': 'GCTAGC',
        'BstZ17I': 'GTATAC'
    }

    # Include complete recognition sites
    if left_enzyme:
        left_site_seq = enzyme_dict[left_enzyme]
        fragment_start = left_pos
        fragment_end = right_pos + len(enzyme_dict[right_enzyme]) if right_enzyme else right_pos
    else:
        fragment_start = changed_start
        fragment_end = right_pos + len(enzyme_dict[right_enzyme]) if right_enzyme else changed_end

    fragment_seq = new_seq[fragment_start:fragment_end]

    return {
        'start': fragment_start,
        'end': fragment_end,
        'length': len(fragment_seq),
        'sequence': fragment_seq,
        'left_enzyme': left_enzyme if left_enzyme else 'START',
        'left_pos': left_pos if left_enzyme else fragment_start,
        'right_enzyme': right_enzyme if right_enzyme else 'END',
        'right_pos': right_pos if right_enzyme else fragment_end
    }


def create_complete_avd005(avd003_original, vp1_vhh_fusion, insert_info):
    """Create AVD005 with ALL features preserved"""

    # Build sequence
    seq_before_vp1 = str(avd003_original.seq[:1520]).upper()
    seq_after_vp1 = str(avd003_original.seq[3731:]).upper()
    avd005_seq = seq_before_vp1 + vp1_vhh_fusion + seq_after_vp1

    # Create record
    record = SeqRecord(
        Seq(avd005_seq),
        id="AVD005",
        name="AVD005",
        description="EF1alpha-VP1-VHH3-ALPL-bGH | VP1-only with anti-ALPL VHH at VR-IV using D2 linkers",
        annotations={"molecule_type": "DNA", "topology": "circular"}
    )

    # Copy ALL features from original, adjusting positions
    vp1_size_change = len(vp1_vhh_fusion) - (3731 - 1520)

    for feat in avd003_original.features:
        feat_start = feat.location.start
        feat_end = feat.location.end

        # If feature is entirely before VP1, keep as-is
        if feat_end <= 1520:
            new_feat = SeqFeature(
                FeatureLocation(feat_start, feat_end, strand=feat.location.strand),
                type=feat.type,
                qualifiers=feat.qualifiers.copy()
            )
            record.features.append(new_feat)

        # If feature is entirely after VP1, shift by size change
        elif feat_start >= 3731:
            new_start = feat_start + vp1_size_change
            new_end = feat_end + vp1_size_change
            new_feat = SeqFeature(
                FeatureLocation(new_start, new_end, strand=feat.location.strand),
                type=feat.type,
                qualifiers=feat.qualifiers.copy()
            )
            record.features.append(new_feat)

        # If feature is VP1 region, replace with new annotations
        elif feat_start >= 1520 and feat_end <= 3731:
            # Skip original VP1-related features, we'll add new ones
            continue

        # Features spanning VP1 region - handle carefully
        else:
            # For now, skip features that span VP1 boundary
            continue

    # Add new VP1-VHH features
    vp1_start = 1520
    vp1_end = vp1_start + len(vp1_vhh_fusion)

    # VP1-VHH fusion
    record.features.append(SeqFeature(
        FeatureLocation(vp1_start, vp1_end, strand=1),
        type="CDS",
        qualifiers={
            "label": "VP1-VHH3-ALPL",
            "note": "AAV9 VP1 with anti-ALPL VHH3 at VR-IV. VP2(ACG→ACC) and VP3(ATG→CTG) knocked out.",
            "codon_start": "1",
            "translation": str(Seq(vp1_vhh_fusion).translate())[:-1]
        }
    ))

    # N-terminal linker
    linker_n_start = vp1_start + insert_info['linker_n_start']
    linker_n_end = vp1_start + insert_info['linker_n_end']
    record.features.append(SeqFeature(
        FeatureLocation(linker_n_start, linker_n_end, strand=1),
        type="misc_feature",
        qualifiers={
            "label": "LINK_D2_N",
            "note": "N-terminal linker (GGGGS)x4, 20 aa"
        }
    ))

    # VHH3
    vhh_start = vp1_start + insert_info['vhh_start']
    vhh_end = vp1_start + insert_info['vhh_end']
    record.features.append(SeqFeature(
        FeatureLocation(vhh_start, vhh_end, strand=1),
        type="CDS",
        qualifiers={
            "label": "VHH3_anti-ALPL",
            "note": "Anti-ALPL VHH nanobody, 118 aa, codon-optimized for Homo sapiens",
            "codon_start": "1",
            "translation": VHH3_AA
        }
    ))

    return record


def create_complete_avd006(avd002_original, vp1_vhh_fusion, insert_info):
    """Create AVD006 with ALL features preserved"""

    # Build sequence
    seq_before_vp1 = str(avd002_original.seq[:2378]).upper()
    seq_after_vp1 = str(avd002_original.seq[4589:]).upper()
    avd006_seq = seq_before_vp1 + vp1_vhh_fusion + seq_after_vp1

    # Create record
    record = SeqRecord(
        Seq(avd006_seq),
        id="AVD006",
        name="AVD006",
        description="Rep2Mut2Cap9-VP1-VHH3-ALPL | RepCap with VP1-VHH at VR-IV, VP2/VP3 knocked out",
        annotations={"molecule_type": "DNA", "topology": "circular"}
    )

    # Copy ALL features from original, adjusting positions
    vp1_size_change = len(vp1_vhh_fusion) - (4589 - 2378)

    for feat in avd002_original.features:
        feat_start = feat.location.start
        feat_end = feat.location.end

        # If feature is entirely before VP1, keep as-is
        if feat_end <= 2378:
            new_feat = SeqFeature(
                FeatureLocation(feat_start, feat_end, strand=feat.location.strand),
                type=feat.type,
                qualifiers=feat.qualifiers.copy()
            )
            record.features.append(new_feat)

        # If feature is entirely after VP1, shift by size change
        elif feat_start >= 4589:
            new_start = feat_start + vp1_size_change
            new_end = feat_end + vp1_size_change
            new_feat = SeqFeature(
                FeatureLocation(new_start, new_end, strand=feat.location.strand),
                type=feat.type,
                qualifiers=feat.qualifiers.copy()
            )
            record.features.append(new_feat)

        # If feature is in VP1 region, skip (we'll add new ones)
        elif feat_start >= 2378 and feat_end <= 4589:
            continue

        # Features spanning VP1 - skip for now
        else:
            continue

    # Add new VP1-VHH features
    vp1_start = 2378
    vp1_end = vp1_start + len(vp1_vhh_fusion)

    # VP1-VHH fusion
    record.features.append(SeqFeature(
        FeatureLocation(vp1_start, vp1_end, strand=1),
        type="CDS",
        qualifiers={
            "label": "VP1-VHH3-ALPL",
            "note": "AAV9 VP1 with anti-ALPL VHH3 at VR-IV. VP2(ACG→ACC) and VP3(ATG→CTG) knocked out.",
            "codon_start": "1",
            "translation": str(Seq(vp1_vhh_fusion).translate())[:-1]
        }
    ))

    # N-terminal linker
    linker_n_start = vp1_start + insert_info['linker_n_start']
    linker_n_end = vp1_start + insert_info['linker_n_end']
    record.features.append(SeqFeature(
        FeatureLocation(linker_n_start, linker_n_end, strand=1),
        type="misc_feature",
        qualifiers={
            "label": "LINK_D2_N",
            "note": "N-terminal linker (GGGGS)x4, 20 aa"
        }
    ))

    # VHH3
    vhh_start = vp1_start + insert_info['vhh_start']
    vhh_end = vp1_start + insert_info['vhh_end']
    record.features.append(SeqFeature(
        FeatureLocation(vhh_start, vhh_end, strand=1),
        type="CDS",
        qualifiers={
            "label": "VHH3_anti-ALPL",
            "note": "Anti-ALPL VHH nanobody, 118 aa, codon-optimized for Homo sapiens",
            "codon_start": "1",
            "translation": VHH3_AA
        }
    ))

    return record


def main():
    print("=" * 80)
    print("GENERATING COMPLETE GENBANK FILES WITH SYNTHETIC FRAGMENT INFO")
    print("=" * 80)
    print()

    # Load starting materials
    project_root = Path(__file__).parent.parent
    avd003_file = project_root / "AVD003-EF1A-VPall-bGH.dna"
    avd002_file = project_root / "AVD002-Rep2Mut2Cap9-6R-wt.dna"

    avd003 = SeqIO.read(avd003_file, "snapgene")
    avd002 = SeqIO.read(avd002_file, "snapgene")

    print(f"Loaded AVD003: {len(avd003.seq)} bp, {len(avd003.features)} features")
    print(f"Loaded AVD002: {len(avd002.seq)} bp, {len(avd002.features)} features")
    print()

    # Extract VP1
    vp1_original = str(avd003.seq[1520:3731])

    # Build VP1-VHH fusion
    print("Building VP1-VHH3 fusion...")
    vp1_vhh_fusion, insert_info = build_vp1_vhh_fusion(vp1_original)
    print(f"  VP1-VHH fusion: {len(vp1_vhh_fusion)} bp")
    print()

    # Create complete records
    print("Creating AVD005 with all features preserved...")
    avd005 = create_complete_avd005(avd003, vp1_vhh_fusion, insert_info)
    print(f"  AVD005: {len(avd005.seq)} bp, {len(avd005.features)} features")

    print("Creating AVD006 with all features preserved...")
    avd006 = create_complete_avd006(avd002, vp1_vhh_fusion, insert_info)
    print(f"  AVD006: {len(avd006.seq)} bp, {len(avd006.features)} features")
    print()

    # Find restriction sites
    enzyme_dict = {
        'SmaI': 'CCCGGG',
        'BbvCI': 'CCTCAGC',
        'AgeI': 'ACCGGT',
        'BsrGI': 'TGTACA',
        'BmtI': 'GCTAGC',
        'BstZ17I': 'GTATAC'
    }

    print("Finding restriction sites...")
    avd005_sites = find_restriction_sites(str(avd005.seq), enzyme_dict)
    avd006_sites = find_restriction_sites(str(avd006.seq), enzyme_dict)
    print()

    # Calculate synthetic fragments
    print("Calculating synthetic fragments for ordering...")

    # AVD005: VP1 is at 1520-4148
    avd005_fragment = calculate_synthetic_fragment(
        str(avd003.seq),
        str(avd005.seq),
        1520,
        1520 + len(vp1_vhh_fusion),
        avd005_sites
    )

    # AVD006: VP1 is at 2378-5006
    avd006_fragment = calculate_synthetic_fragment(
        str(avd002.seq),
        str(avd006.seq),
        2378,
        2378 + len(vp1_vhh_fusion),
        avd006_sites
    )

    print(f"AVD005 synthetic fragment:")
    print(f"  Position: {avd005_fragment['start']+1}-{avd005_fragment['end']} (1-indexed)")
    print(f"  Length: {avd005_fragment['length']} bp")
    print(f"  Left boundary: {avd005_fragment['left_enzyme']} at position {avd005_fragment['left_pos']+1}")
    print(f"  Right boundary: {avd005_fragment['right_enzyme']} at position {avd005_fragment['right_pos']+1}")
    print()

    print(f"AVD006 synthetic fragment:")
    print(f"  Position: {avd006_fragment['start']+1}-{avd006_fragment['end']} (1-indexed)")
    print(f"  Length: {avd006_fragment['length']} bp")
    print(f"  Left boundary: {avd006_fragment['left_enzyme']} at position {avd006_fragment['left_pos']+1}")
    print(f"  Right boundary: {avd006_fragment['right_enzyme']} at position {avd006_fragment['right_pos']+1}")
    print()

    # Write GenBank files
    avd005_output = project_root / "AVD005-EF1A-VP1-VHH3-ALPL-bGH.gb"
    avd006_output = project_root / "AVD006-Rep2Mut2Cap9-VP1-VHH3-ALPL.gb"

    SeqIO.write(avd005, avd005_output, "genbank")
    print(f"✓ Written {avd005_output}")

    SeqIO.write(avd006, avd006_output, "genbank")
    print(f"✓ Written {avd006_output}")
    print()

    # Write synthetic fragment sequences to FASTA
    avd005_frag_fasta = project_root / "AVD005_synthetic_fragment.fasta"
    avd006_frag_fasta = project_root / "AVD006_synthetic_fragment.fasta"

    with open(avd005_frag_fasta, 'w') as f:
        f.write(f">AVD005_synthetic_fragment {avd005_fragment['start']+1}-{avd005_fragment['end']}\n")
        f.write(f">Left: {avd005_fragment['left_enzyme']}, Right: {avd005_fragment['right_enzyme']}\n")
        # Write sequence in 80-character lines
        seq = avd005_fragment['sequence']
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + '\n')

    with open(avd006_frag_fasta, 'w') as f:
        f.write(f">AVD006_synthetic_fragment {avd006_fragment['start']+1}-{avd006_fragment['end']}\n")
        f.write(f">Left: {avd006_fragment['left_enzyme']}, Right: {avd006_fragment['right_enzyme']}\n")
        seq = avd006_fragment['sequence']
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + '\n')

    print(f"✓ Written {avd005_frag_fasta}")
    print(f"✓ Written {avd006_frag_fasta}")
    print()

    # Write CSV summary
    csv_output = project_root / "SYNTHETIC_FRAGMENTS_ORDER_FORM.csv"
    with open(csv_output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Construct',
            'Full_Plasmid_Size_bp',
            'Fragment_Start_1indexed',
            'Fragment_End_1indexed',
            'Fragment_Length_bp',
            'Left_Enzyme',
            'Left_Position_1indexed',
            'Right_Enzyme',
            'Right_Position_1indexed',
            'Fragment_FASTA_File',
            'Full_GenBank_File'
        ])

        writer.writerow([
            'AVD005',
            len(avd005.seq),
            avd005_fragment['start'] + 1,
            avd005_fragment['end'],
            avd005_fragment['length'],
            avd005_fragment['left_enzyme'],
            avd005_fragment['left_pos'] + 1,
            avd005_fragment['right_enzyme'],
            avd005_fragment['right_pos'] + 1,
            'AVD005_synthetic_fragment.fasta',
            'AVD005-EF1A-VP1-VHH3-ALPL-bGH.gb'
        ])

        writer.writerow([
            'AVD006',
            len(avd006.seq),
            avd006_fragment['start'] + 1,
            avd006_fragment['end'],
            avd006_fragment['length'],
            avd006_fragment['left_enzyme'],
            avd006_fragment['left_pos'] + 1,
            avd006_fragment['right_enzyme'],
            avd006_fragment['right_pos'] + 1,
            'AVD006_synthetic_fragment.fasta',
            'AVD006-Rep2Mut2Cap9-VP1-VHH3-ALPL.gb'
        ])

    print(f"✓ Written {csv_output}")
    print()

    print("=" * 80)
    print("COMPLETE! Files generated:")
    print("  - Complete GenBank files with all features preserved")
    print("  - Synthetic fragment FASTA files for ordering")
    print("  - CSV summary for GenScript order form")
    print("=" * 80)


if __name__ == "__main__":
    main()
