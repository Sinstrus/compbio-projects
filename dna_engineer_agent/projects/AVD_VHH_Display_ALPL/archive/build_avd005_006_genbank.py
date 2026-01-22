#!/usr/bin/env python3
"""
Build AVD005 and AVD006 complete constructs and generate GenBank files

v2.0 - Uses dnachisel-optimized insert sequence to avoid:
  - High GC content (original: 93% in linker, 70% in VHH)
  - Repetitive sequences causing synthesis problems
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from datetime import datetime
import sys
from pathlib import Path

# Import design functions from the design script
sys.path.insert(0, str(Path(__file__).parent))

# VHH3 anti-ALPL sequence (118 amino acids)
VHH3_AA = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSTADMGWFRQAPGKGRELVAAVSGSGFSTYSDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAKATISLYYAMDVWGQGTTVTVSS"

# D2 Linkers (Biogen patent specification)
LINK_D2_N_AA = "GGGGSGGGGSGGGGSGGGGS"  # (GGGGS)x4 = 20 aa
LINK_D2_C_AA = "GGGGS"  # (GGGGS)x1 = 5 aa (ASYMMETRIC, not direct fusion)

# DNACHISEL-OPTIMIZED INSERT SEQUENCE (N-linker + VHH + C-linker)
# Optimized to:
#   - GC content: 60.2% (was 73.2%)
#   - Repeats >= 9bp: minimal (was 718)
#   - N-linker: varied codons (GGT, GGA, GGC for Gly; TCT, TCA, AGT for Ser)
#   - C-linker: (GGGGS)x1 with varied codons
#   - VHH optimized with dnachisel EnforceGCContent(35-65%) and AvoidHairpins
# Total: 429 bp encoding 143 aa (20 aa N-linker + 118 aa VHH + 5 aa C-linker)
OPTIMIZED_INSERT = (
    # N-terminal linker (60 bp, 20 aa) - manually designed with varied codons
    "GGTGGAGGCGGATCTGGAGGCGGTGGTTCAGGCGGTGGAGGAAGTGGTGGCGGAGGTTCT"
    # VHH (354 bp, 118 aa) - dnachisel optimized
    "GAGGTGCAACTGGTTGAAAGCGGCGGAGGACTTGTTCAACCCGGCGGCAGCCTTAGGCTTTCTTGCGCTGCCAGCGGCTTCACCTTTAGCACCGCCGACATGGGCTGGTTTAGGCAAGCTCCCGGAAAAGGCAGGGAACTTGTTGCCGCTGTGAGCGGCAGCGGCTTCAGCACCTACTCTGATAGCGTTGAGGGCAGGTTCACCATCAGCAGGGACAACGCCAAGAGGATGGTGTACCTGCAGATGAACAGCTTGAGGGCCGAGGACACCGCCGTGTACTACTGCGCCAAGGCCACAATTAGCCTGTACTACGCCATGGATGTGTGGGGACAGGGCACCACCGTGACCGTGAGCAGC"
    # C-terminal linker (15 bp, 5 aa) - varied codons (GGGGS = GGT GGA GGC GGT TCT)
    "GGTGGAGGCGGTTCT"
)


def build_vp1_vhh_fusion(vp1_original, insertion_aa_pos=456):
    """
    Build VP1-VHH fusion with D2 linkers.

    Steps:
    1. Split VP1 at insertion position
    2. Insert: N-linker + VHH (using OPTIMIZED_INSERT for synthesis)
    3. Knock out VP2 (ACG→ACC at codon 137) and VP3 (ATG→CTG at codon 202)
    """
    # Convert to uppercase and list for mutation
    vp1 = list(vp1_original.upper())

    # First, knock out VP2 and VP3 BEFORE insertion
    # VP2: ACG at codon 137 (bp 411) → ACC
    if ''.join(vp1[411:414]) == 'ACG':
        vp1[413] = 'C'  # ACG → ACC
        print("  ✓ VP2 knockout: ACG → ACC at codon 137 (bp 411)")
    else:
        print(f"  ⚠ VP2 codon at 411 is {''.join(vp1[411:414])}, not ACG")

    # VP3: ATG at codon 202 (bp 606) → CTG
    if ''.join(vp1[606:609]) == 'ATG':
        vp1[606] = 'C'  # ATG → CTG
        print("  ✓ VP3 knockout: ATG → CTG at codon 202 (bp 606)")
    else:
        print(f"  ⚠ VP3 codon at 606 is {''.join(vp1[606:609])}, not ATG")

    # Convert back to string
    vp1_ko = ''.join(vp1)

    # Now insert VHH at position 456 aa (bp 1368)
    insertion_bp = insertion_aa_pos * 3
    vp1_before = vp1_ko[:insertion_bp]
    vp1_after = vp1_ko[insertion_bp:]

    # Use dnachisel-optimized insert (N-linker + VHH + C-linker)
    insert_cassette = OPTIMIZED_INSERT
    linker_n_len = 60  # 20 aa * 3 bp
    vhh_len = 354      # 118 aa * 3 bp
    linker_c_len = 15  # 5 aa * 3 bp

    # Verify translation matches expected amino acids
    insert_aa = str(Seq(insert_cassette).translate())
    expected_aa = LINK_D2_N_AA + VHH3_AA + LINK_D2_C_AA
    if insert_aa != expected_aa:
        print(f"  ⚠ WARNING: Optimized insert translates to different AA sequence!")
        print(f"    Expected: {expected_aa}")
        print(f"    Got:      {insert_aa}")
    else:
        print(f"  ✓ Optimized insert translation verified (143 aa)")

    vp1_vhh_fusion = vp1_before + insert_cassette + vp1_after

    print(f"  ✓ VHH insertion at aa {insertion_aa_pos}:")
    print(f"    N-linker: {linker_n_len} bp ({len(LINK_D2_N_AA)} aa) - asymmetric long")
    print(f"    VHH3: {vhh_len} bp ({len(VHH3_AA)} aa)")
    print(f"    C-linker: {linker_c_len} bp ({len(LINK_D2_C_AA)} aa) - asymmetric short")
    print(f"    Total insert: {len(insert_cassette)} bp (dnachisel-optimized)")
    print(f"    Insert GC content: {100*sum(1 for c in insert_cassette if c in 'GC')/len(insert_cassette):.1f}%")
    print(f"    Final VP1-VHH: {len(vp1_vhh_fusion)} bp ({len(vp1_vhh_fusion)//3} aa)")

    return vp1_vhh_fusion, {
        'insertion_bp': insertion_bp,
        'linker_n_start': insertion_bp,
        'linker_n_end': insertion_bp + linker_n_len,
        'vhh_start': insertion_bp + linker_n_len,
        'vhh_end': insertion_bp + linker_n_len + vhh_len,
        'linker_c_start': insertion_bp + linker_n_len + vhh_len,
        'linker_c_end': insertion_bp + len(insert_cassette),
        'vp2_knockout_bp': 411,
        'vp3_knockout_bp': 606,
    }


def check_restriction_sites_full(full_seq, sites):
    """Check restriction sites in full construct"""
    results = {}
    for enzyme, recog_seq in sites.items():
        positions = []
        seq_upper = full_seq.upper()
        start = 0
        while True:
            pos = seq_upper.find(recog_seq, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1

        status = "✓ UNIQUE" if len(positions) == 1 else f"⚠ {len(positions)} sites" if len(positions) > 1 else "○ ABSENT"
        results[enzyme] = {
            'count': len(positions),
            'positions': positions,
            'status': status
        }
    return results


def create_avd005(avd003_original, vp1_vhh_fusion, insert_info):
    """
    Create AVD005: EF1α-VP1-VHH3-ALPL-bGH
    Based on AVD003 with VP1 replaced by VP1-VHH fusion
    """
    print("\nBuilding AVD005 (EF1α-VP1-VHH3-ALPL-bGH)...")

    # Extract original sequence parts
    seq_before_vp1 = str(avd003_original.seq[:1520]).upper()  # Up to VP1 start
    seq_after_vp1 = str(avd003_original.seq[3731:]).upper()   # After VP1 end

    # Assemble new construct
    avd005_seq = seq_before_vp1 + vp1_vhh_fusion + seq_after_vp1

    print(f"  Backbone before VP1: {len(seq_before_vp1)} bp")
    print(f"  VP1-VHH fusion: {len(vp1_vhh_fusion)} bp")
    print(f"  Backbone after VP1: {len(seq_after_vp1)} bp")
    print(f"  Total AVD005: {len(avd005_seq)} bp")

    # Create SeqRecord
    record = SeqRecord(
        Seq(avd005_seq),
        id="AVD005",
        name="AVD005",
        description="EF1alpha-VP1-VHH3-ALPL-bGH | VP1-only with anti-ALPL VHH at VR-IV using D2 linkers",
        annotations={"molecule_type": "DNA", "topology": "circular"}
    )

    # Add features (simplified - main components only)
    vp1_start_in_construct = 1520

    features = []

    # EF1α promoter (original position)
    features.append(SeqFeature(
        FeatureLocation(312, 1506),
        type="promoter",
        qualifiers={"label": "EF1α promoter",
                   "note": "Human EEF1A1 promoter with first intron (1194 bp)"}
    ))

    # VP1-VHH fusion CDS
    vp1_end_in_construct = vp1_start_in_construct + len(vp1_vhh_fusion)
    features.append(SeqFeature(
        FeatureLocation(vp1_start_in_construct, vp1_end_in_construct, strand=+1),
        type="CDS",
        qualifiers={"label": "VP1-VHH3-ALPL fusion",
                   "note": "AAV9 VP1 with anti-ALPL VHH inserted at VR-IV using D2 linkers",
                   "translation": str(Seq(vp1_vhh_fusion).translate())[:-1]}  # Remove stop
    ))

    # N-terminal linker
    linker_n_start = vp1_start_in_construct + insert_info['linker_n_start']
    linker_n_end = vp1_start_in_construct + insert_info['linker_n_end']
    features.append(SeqFeature(
        FeatureLocation(linker_n_start, linker_n_end),
        type="misc_feature",
        qualifiers={"label": "LINK_D2_N",
                   "note": f"N-terminal linker (GGGGS)x4, {len(LINK_D2_N_AA)} aa"}
    ))

    # VHH3
    vhh_start = vp1_start_in_construct + insert_info['vhh_start']
    vhh_end = vp1_start_in_construct + insert_info['vhh_end']
    features.append(SeqFeature(
        FeatureLocation(vhh_start, vhh_end, strand=+1),
        type="CDS",
        qualifiers={"label": "VHH3_anti-ALPL",
                   "note": f"Anti-ALPL VHH nanobody, {len(VHH3_AA)} aa, codon-optimized for Homo sapiens",
                   "translation": VHH3_AA}
    ))

    # bGH polyA (original position adjusted)
    polya_start = vp1_end_in_construct
    polya_end = polya_start + 132  # Original bGH polyA is 132 bp
    features.append(SeqFeature(
        FeatureLocation(polya_start, polya_end),
        type="polyA_signal",
        qualifiers={"label": "bGH polyA"}
    ))

    record.features = features

    return record


def create_avd006(avd002_original, vp1_vhh_fusion, insert_info):
    """
    Create AVD006: Rep2Mut2Cap9-VP1-VHH3-ALPL
    Based on AVD002 with VP1 replaced by VP1-VHH fusion
    """
    print("\nBuilding AVD006 (Rep2Mut2Cap9-VP1-VHH3-ALPL)...")

    # Extract original sequence parts
    seq_before_vp1 = str(avd002_original.seq[:2378]).upper()  # Up to VP1 start
    seq_after_vp1 = str(avd002_original.seq[4589:]).upper()   # After VP1 end

    # Assemble new construct
    avd006_seq = seq_before_vp1 + vp1_vhh_fusion + seq_after_vp1

    print(f"  Backbone before VP1: {len(seq_before_vp1)} bp")
    print(f"  VP1-VHH fusion: {len(vp1_vhh_fusion)} bp")
    print(f"  Backbone after VP1: {len(seq_after_vp1)} bp")
    print(f"  Total AVD006: {len(avd006_seq)} bp")

    # Create SeqRecord
    record = SeqRecord(
        Seq(avd006_seq),
        id="AVD006",
        name="AVD006",
        description="Rep2Mut2Cap9-VP1-VHH3-ALPL | RepCap helper with VP1-VHH fusion at VR-IV using D2 linkers",
        annotations={"molecule_type": "DNA", "topology": "circular"}
    )

    # Add features (simplified)
    vp1_start_in_construct = 2378
    vp1_end_in_construct = vp1_start_in_construct + len(vp1_vhh_fusion)

    features = []

    # Rep78 (original position)
    features.append(SeqFeature(
        FeatureLocation(496, 2362, strand=+1),
        type="CDS",
        qualifiers={"label": "Rep78",
                   "note": "AAV2 Rep78 replication protein"}
    ))

    # Rep52 (original position)
    features.append(SeqFeature(
        FeatureLocation(1168, 2362, strand=+1),
        type="CDS",
        qualifiers={"label": "Rep52",
                   "note": "AAV2 Rep52 packaging protein"}
    ))

    # p40 promoter (original position)
    features.append(SeqFeature(
        FeatureLocation(1875, 2028),
        type="promoter",
        qualifiers={"label": "p40 promoter",
                   "note": "AAV p40 promoter drives VP1/2/3 and AAP expression"}
    ))

    # VP1-VHH fusion CDS
    features.append(SeqFeature(
        FeatureLocation(vp1_start_in_construct, vp1_end_in_construct, strand=+1),
        type="CDS",
        qualifiers={"label": "VP1-VHH3-ALPL fusion",
                   "note": "AAV9 VP1 with anti-ALPL VHH inserted at VR-IV using D2 linkers. VP2 and VP3 knocked out.",
                   "translation": str(Seq(vp1_vhh_fusion).translate())[:-1]}
    ))

    # N-terminal linker
    linker_n_start = vp1_start_in_construct + insert_info['linker_n_start']
    linker_n_end = vp1_start_in_construct + insert_info['linker_n_end']
    features.append(SeqFeature(
        FeatureLocation(linker_n_start, linker_n_end),
        type="misc_feature",
        qualifiers={"label": "LINK_D2_N",
                   "note": f"N-terminal linker (GGGGS)x4, {len(LINK_D2_N_AA)} aa"}
    ))

    # VHH3
    vhh_start = vp1_start_in_construct + insert_info['vhh_start']
    vhh_end = vp1_start_in_construct + insert_info['vhh_end']
    features.append(SeqFeature(
        FeatureLocation(vhh_start, vhh_end, strand=+1),
        type="CDS",
        qualifiers={"label": "VHH3_anti-ALPL",
                   "note": f"Anti-ALPL VHH nanobody, {len(VHH3_AA)} aa, codon-optimized for Homo sapiens",
                   "translation": VHH3_AA}
    ))

    # AAP (in +1 frame, adjust position)
    aap_start = vp1_start_in_construct + 526  # ~codon 50 of VP1 in +1 frame = (50*3)+1 = 151
    aap_end = aap_start + 594  # AAP is typically ~198 aa = 594 bp
    features.append(SeqFeature(
        FeatureLocation(aap_start, aap_end, strand=+1),
        type="CDS",
        qualifiers={"label": "AAP",
                   "note": "Assembly-Activating Protein, +1 frame relative to VP1"}
    ))

    record.features = features

    return record


def main():
    print("=" * 80)
    print("BUILDING AVD005 & AVD006 GENBANK FILES")
    print("=" * 80)
    print()

    # Load starting materials
    project_root = Path(__file__).parent.parent
    plasmids_dir = project_root / "plasmids"
    avd003_file = plasmids_dir / "AVD003-EF1A-VPall-bGH.dna"
    avd002_file = plasmids_dir / "AVD002-Rep2Mut2Cap9-6R-wt.dna"

    avd003 = SeqIO.read(avd003_file, "snapgene")
    avd002 = SeqIO.read(avd002_file, "snapgene")

    print(f"Loaded AVD003: {len(avd003.seq)} bp")
    print(f"Loaded AVD002: {len(avd002.seq)} bp")
    print()

    # Extract VP1 from EACH parent plasmid (they have different sequences!)
    vp1_from_avd003 = str(avd003.seq[1520:3731])
    vp1_from_avd002 = str(avd002.seq[2378:4589])
    print(f"VP1 from AVD003: {len(vp1_from_avd003)} bp ({len(vp1_from_avd003)//3} aa)")
    print(f"VP1 from AVD002: {len(vp1_from_avd002)} bp ({len(vp1_from_avd002)//3} aa)")
    print()

    # Build VP1-VHH fusion for AVD005 (using AVD003's VP1)
    print("Creating VP1-VHH3 fusion for AVD005 (from AVD003 VP1)...")
    vp1_vhh_for_avd005, insert_info_005 = build_vp1_vhh_fusion(vp1_from_avd003, insertion_aa_pos=456)
    print()

    # Build VP1-VHH fusion for AVD006 (using AVD002's VP1)
    print("Creating VP1-VHH3 fusion for AVD006 (from AVD002 VP1)...")
    vp1_vhh_for_avd006, insert_info_006 = build_vp1_vhh_fusion(vp1_from_avd002, insertion_aa_pos=456)
    print()

    # Build AVD005 using AVD003's VP1-VHH
    avd005 = create_avd005(avd003, vp1_vhh_for_avd005, insert_info_005)

    # Build AVD006 using AVD002's VP1-VHH
    avd006 = create_avd006(avd002, vp1_vhh_for_avd006, insert_info_006)
    print()

    # Check 6R sites in both constructs
    restriction_sites_6r = {
        'SmaI': 'CCCGGG',
        'BbvCI': 'CCTCAGC',
        'AgeI': 'ACCGGT',
        'BsrGI': 'TGTACA',
        'BmtI': 'GCTAGC',
        'BstZ17I': 'GTATAC'
    }

    print("Checking 6R restriction sites in AVD005...")
    avd005_sites = check_restriction_sites_full(str(avd005.seq), restriction_sites_6r)
    for enzyme, info in avd005_sites.items():
        print(f"  {enzyme:10s} ({restriction_sites_6r[enzyme]:7s}): {info['status']}")
    print()

    print("Checking 6R restriction sites in AVD006...")
    avd006_sites = check_restriction_sites_full(str(avd006.seq), restriction_sites_6r)
    for enzyme, info in avd006_sites.items():
        print(f"  {enzyme:10s} ({restriction_sites_6r[enzyme]:7s}): {info['status']}")
    print()

    # Write GenBank files
    output_dir = plasmids_dir
    avd005_output = output_dir / "AVD005-EF1A-VP1-VHH3-ALPL-bGH.gb"
    avd006_output = output_dir / "AVD006-Rep2Mut2Cap9-VP1-VHH3-ALPL.gb"

    SeqIO.write(avd005, avd005_output, "genbank")
    print(f"✓ Written {avd005_output}")
    print(f"  Size: {len(avd005.seq)} bp")

    SeqIO.write(avd006, avd006_output, "genbank")
    print(f"✓ Written {avd006_output}")
    print(f"  Size: {len(avd006.seq)} bp")
    print()

    print("=" * 80)
    print("SUCCESS: AVD005 and AVD006 GenBank files generated!")
    print("=" * 80)


if __name__ == "__main__":
    main()
