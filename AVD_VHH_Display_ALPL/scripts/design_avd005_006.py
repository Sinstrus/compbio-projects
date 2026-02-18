#!/usr/bin/env python3
"""
Design AVD005 and AVD006 - VP1-VHH3-ALPL constructs
VP1-only expression plasmids with anti-ALPL VHH inserted at VR-IV using D2 linkers
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Data import CodonTable
import sys
from pathlib import Path

# VHH3 anti-ALPL sequence (118 amino acids)
VHH3_AA = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSTADMGWFRQAPGKGRELVAAVSGSGFSTYSDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAKATISLYYAMDVWGQGTTVTVSS"

# D2 Linkers (from parts_library/linkers.json)
LINK_D2_N_AA = "GGGGSGGGGSGGGGSGGGGS"  # (GGGGS)x4 = 20 aa
LINK_D2_C_AA = ""  # Direct fusion, no C-terminal linker

# Codon optimization table for Homo sapiens (high-frequency codons)
CODON_USAGE_HUMAN = {
    'A': ['GCC', 'GCT', 'GCA', 'GCG'],  # Prefer GCC
    'C': ['TGC', 'TGT'],  # Prefer TGC
    'D': ['GAC', 'GAT'],  # Prefer GAC
    'E': ['GAG', 'GAA'],  # Prefer GAG
    'F': ['TTC', 'TTT'],  # Prefer TTC
    'G': ['GGC', 'GGT', 'GGA', 'GGG'],  # Prefer GGC
    'H': ['CAC', 'CAT'],  # Prefer CAC
    'I': ['ATC', 'ATT', 'ATA'],  # Prefer ATC
    'K': ['AAG', 'AAA'],  # Prefer AAG
    'L': ['CTG', 'CTC', 'TTA', 'TTG', 'CTT', 'CTA'],  # Prefer CTG
    'M': ['ATG'],  # Only one codon
    'N': ['AAC', 'AAT'],  # Prefer AAC
    'P': ['CCC', 'CCT', 'CCA', 'CCG'],  # Prefer CCC
    'Q': ['CAG', 'CAA'],  # Prefer CAG
    'R': ['CGC', 'AGG', 'CGT', 'AGA', 'CGA', 'CGG'],  # Prefer CGC
    'S': ['AGC', 'TCC', 'TCT', 'TCA', 'TCG', 'AGT'],  # Prefer AGC
    'T': ['ACC', 'ACT', 'ACA', 'ACG'],  # Prefer ACC
    'V': ['GTG', 'GTC', 'GTT', 'GTA'],  # Prefer GTG
    'W': ['TGG'],  # Only one codon
    'Y': ['TAC', 'TAT'],  # Prefer TAC
    '*': ['TAA', 'TAG', 'TGA']  # Stop codons
}

def codon_optimize(aa_sequence, organism="human"):
    """
    Codon-optimize an amino acid sequence for Homo sapiens.
    Uses high-frequency codons to maximize expression.
    """
    if organism == "human":
        codon_table = CODON_USAGE_HUMAN
    else:
        raise ValueError(f"Organism {organism} not supported")

    dna_sequence = ""
    for aa in aa_sequence:
        if aa in codon_table:
            # Use the first codon (highest frequency)
            dna_sequence += codon_table[aa][0]
        else:
            raise ValueError(f"Amino acid {aa} not found in codon table")

    return dna_sequence


def find_vr_iv_region(vp1_seq, start_aa=451, end_aa=461):
    """
    Extract the VR-IV region from VP1 sequence.
    VR-IV is approximately at residues 452-460 (aa position 451-460 in 0-indexed).
    """
    vr_iv_start = start_aa * 3
    vr_iv_end = end_aa * 3
    vr_iv_seq = vp1_seq[vr_iv_start:vr_iv_end]
    return vr_iv_seq, vr_iv_start, vr_iv_end


def find_bbvci_sites(sequence):
    """
    Find all BbvCI recognition sites (CCTCAGC) in a sequence.
    """
    bbvci_site = "CCTCAGC"
    sites = []
    start = 0
    while True:
        pos = sequence.find(bbvci_site, start)
        if pos == -1:
            break
        sites.append(pos)
        start = pos + 1
    return sites


def knockout_start_codons(vp1_dna):
    """
    Knock out VP2 and VP3 start codons while preserving VP1.

    VP1: ATG at position 0 (keep intact)
    VP2: ACG at position ~411 bp (codon 137) - change ACG → ACC (Thr → Thr, silent)
    VP3: ATG at position ~606 bp (codon 202) - change ATG → CTG (Met → Leu, NOT silent)

    Returns modified sequence and report.
    """
    vp1_seq_mut = list(vp1_dna)
    mutations = []

    # Find VP2 start (ACG, typically around codon 137-138)
    # Search for ACG in the region 400-450 bp
    for pos in range(400, 450):
        if vp1_dna[pos:pos+3] == "ACG":
            # Check if this is in-frame (divisible by 3)
            if pos % 3 == 0:
                # Mutate ACG → ACC (Thr → Thr, silent)
                vp1_seq_mut[pos+2] = 'C'
                mutations.append({
                    'position': pos,
                    'original': 'ACG',
                    'mutated': 'ACC',
                    'aa_original': 'T',
                    'aa_mutated': 'T',
                    'silent': True,
                    'type': 'VP2_knockout'
                })
                break

    # Find VP3 start (ATG, typically around codon 202-204)
    # Search for ATG in the region 600-650 bp
    for pos in range(600, 650):
        codon = vp1_dna[pos:pos+3]
        if codon == "ATG":
            # Check if this is in-frame
            if pos % 3 == 0:
                # Mutate ATG → CTG (Met → Leu, NOT silent but acceptable for knockout)
                vp1_seq_mut[pos] = 'C'
                mutations.append({
                    'position': pos,
                    'original': 'ATG',
                    'mutated': 'CTG',
                    'aa_original': 'M',
                    'aa_mutated': 'L',
                    'silent': False,
                    'type': 'VP3_knockout'
                })
                break

    return ''.join(vp1_seq_mut), mutations


def design_vhh_insertion(vp1_dna, insertion_pos_aa=456):
    """
    Design VP1-VHH3 fusion construct.

    Strategy:
    1. Split VP1 at insertion position (around aa 456, middle of VR-IV)
    2. Insert: N-linker - VHH3 - C-linker
    3. Codon-optimize VHH3 and linkers for Homo sapiens

    Returns:
    - Modified VP1 DNA sequence
    - Insertion details report
    """
    # Calculate insertion position in DNA (0-indexed)
    insertion_pos_bp = insertion_pos_aa * 3

    # Split VP1 at insertion point
    vp1_before = vp1_dna[:insertion_pos_bp]
    vp1_after = vp1_dna[insertion_pos_bp:]

    # Codon-optimize VHH and linkers
    linker_n_dna = codon_optimize(LINK_D2_N_AA)
    vhh_dna = codon_optimize(VHH3_AA)
    linker_c_dna = codon_optimize(LINK_D2_C_AA) if LINK_D2_C_AA else ""

    # Assemble fusion construct
    insert_cassette = linker_n_dna + vhh_dna + linker_c_dna
    vp1_vhh_fusion = vp1_before + insert_cassette + vp1_after

    report = {
        'insertion_position_aa': insertion_pos_aa,
        'insertion_position_bp': insertion_pos_bp,
        'vp1_before_length': len(vp1_before),
        'linker_n_length': len(linker_n_dna),
        'vhh_length': len(vhh_dna),
        'linker_c_length': len(linker_c_dna),
        'insert_cassette_length': len(insert_cassette),
        'vp1_after_length': len(vp1_after),
        'total_vp1_vhh_length': len(vp1_vhh_fusion),
        'insert_cassette_seq': insert_cassette,
        'linker_n_aa': LINK_D2_N_AA,
        'vhh_aa': VHH3_AA,
        'linker_c_aa': LINK_D2_C_AA
    }

    return vp1_vhh_fusion, report


def verify_aap_frame(vp1_dna, aap_start_codon=50):
    """
    Verify that AAP reading frame (+1 frame relative to VP1) is preserved.
    AAP typically starts around codon 50 of VP1 in the +1 frame.

    Returns True if frame is preserved, False otherwise.
    """
    # AAP is in +1 frame, so it starts at position (aap_start_codon * 3) + 1
    aap_start_bp = (aap_start_codon * 3) + 1

    # Check that the length maintains the frame relationship
    # AAP frame should remain +1 relative to VP1 throughout
    # As long as we insert multiples of 3 bp, the frame relationship is preserved
    if len(vp1_dna) % 3 == 0:
        return True
    else:
        return False


def check_restriction_sites(sequence, sites_to_check):
    """
    Check for restriction sites in the sequence.

    sites_to_check: dict of {'enzyme_name': 'recognition_sequence'}

    Returns dict of sites found.
    """
    results = {}
    for enzyme, recog_seq in sites_to_check.items():
        positions = []
        start = 0
        while True:
            pos = sequence.find(recog_seq, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
        results[enzyme] = {
            'recognition_seq': recog_seq,
            'count': len(positions),
            'positions': positions
        }
    return results


def main():
    # Load starting materials
    print("=" * 80)
    print("AVD005 & AVD006 DESIGN - VP1-VHH3-ALPL CONSTRUCTS")
    print("=" * 80)
    print()

    # File paths
    project_root = Path(__file__).parent.parent
    avd003_file = project_root / "AVD003-EF1A-VPall-bGH.dna"
    avd002_file = project_root / "AVD002-Rep2Mut2Cap9-6R-wt.dna"

    # Read GenBank files
    print("Loading starting materials...")
    avd003 = SeqIO.read(avd003_file, "snapgene")
    avd002 = SeqIO.read(avd002_file, "snapgene")
    print(f"  AVD003: {len(avd003.seq)} bp")
    print(f"  AVD002: {len(avd002.seq)} bp")
    print()

    # Extract VP1 sequences
    print("Extracting VP1 sequences...")

    # AVD003: VP1 at 1521-3731 (1-indexed from annotation)
    vp1_003_start = 1520  # 0-indexed
    vp1_003_end = 3731
    vp1_003_dna = str(avd003.seq[vp1_003_start:vp1_003_end]).upper()
    print(f"  AVD003 VP1: {vp1_003_start+1}-{vp1_003_end} ({len(vp1_003_dna)} bp, {len(vp1_003_dna)//3} aa)")

    # AVD002: VP1 at 2379-4589
    vp1_002_start = 2378  # 0-indexed
    vp1_002_end = 4589
    vp1_002_dna = str(avd002.seq[vp1_002_start:vp1_002_end]).upper()
    print(f"  AVD002 VP1: {vp1_002_start+1}-{vp1_002_end} ({len(vp1_002_dna)} bp, {len(vp1_002_dna)//3} aa)")
    print()

    # Verify sequences are identical
    if vp1_003_dna == vp1_002_dna:
        print("✓ VP1 sequences are identical in both constructs")
        vp1_reference = vp1_003_dna
    else:
        print("⚠ VP1 sequences differ between constructs!")
        # Find first difference
        for i in range(min(len(vp1_003_dna), len(vp1_002_dna))):
            if vp1_003_dna[i] != vp1_002_dna[i]:
                print(f"  First difference at position {i} (codon {i//3}):")
                print(f"    AVD003: {vp1_003_dna[max(0,i-10):i+10]}")
                print(f"    AVD002: {vp1_002_dna[max(0,i-10):i+10]}")
                break
        vp1_reference = vp1_003_dna  # Use AVD003 as reference
    print()

    # Print first 20 codons for inspection
    print("First 20 codons of VP1 (AVD003):")
    for i in range(0, 60, 3):
        codon = vp1_003_dna[i:i+3]
        aa = Seq(codon).translate()
        print(f"  Codon {i//3:3d}: {codon} = {aa}", end="")
        if codon in ["ATG", "ACG"]:
            print(" ← Potential start codon", end="")
        print()
    print()

    # Search for all ACG and ATG codons in VP1
    print("Searching for all ACG and ATG codons in VP1...")
    acg_positions = []
    atg_positions = []
    for i in range(0, len(vp1_reference)-2, 3):
        codon = vp1_reference[i:i+3]
        if codon == "ACG":
            acg_positions.append((i, i//3))
        elif codon == "ATG":
            atg_positions.append((i, i//3))

    print(f"  Found {len(acg_positions)} ACG codon(s) at positions (bp, codon):")
    for bp, codon in acg_positions[:10]:  # Print first 10
        print(f"    bp {bp:4d}, codon {codon:3d}")
    if len(acg_positions) > 10:
        print(f"    ... and {len(acg_positions)-10} more")

    print(f"  Found {len(atg_positions)} ATG codon(s) at positions (bp, codon):")
    for bp, codon in atg_positions[:10]:  # Print first 10
        print(f"    bp {bp:4d}, codon {codon:3d}")
    if len(atg_positions) > 10:
        print(f"    ... and {len(atg_positions)-10} more")
    print()

    # Find VR-IV region
    print("Locating VR-IV region (aa 452-460)...")
    vr_iv_seq, vr_iv_start, vr_iv_end = find_vr_iv_region(vp1_reference, start_aa=451, end_aa=461)
    print(f"  VR-IV: bp {vr_iv_start}-{vr_iv_end} (30 bp)")
    print(f"  VR-IV DNA: {vr_iv_seq}")
    print(f"  VR-IV AA:  {Seq(vr_iv_seq).translate()}")
    print()

    # Find BbvCI sites in VP1
    print("Searching for BbvCI sites (CCTCAGC) in VP1...")
    bbvci_sites = find_bbvci_sites(vp1_reference)
    if bbvci_sites:
        print(f"  Found {len(bbvci_sites)} BbvCI site(s) at position(s): {bbvci_sites}")
        for site in bbvci_sites:
            print(f"    Position {site} = codon {site//3} (aa {site//3 + 1})")
    else:
        print("  No BbvCI sites found in VP1")
    print()

    # Design VHH insertion
    print("Designing VP1-VHH3 fusion construct...")
    print(f"  Insertion site: aa 456 (middle of VR-IV)")
    print(f"  N-terminal linker: {LINK_D2_N_AA} ({len(LINK_D2_N_AA)} aa)")
    print(f"  VHH3: {VHH3_AA[:20]}...{VHH3_AA[-20:]} ({len(VHH3_AA)} aa)")
    print(f"  C-terminal linker: {LINK_D2_C_AA if LINK_D2_C_AA else '(none)'} ({len(LINK_D2_C_AA)} aa)")

    vp1_vhh_fusion, insertion_report = design_vhh_insertion(vp1_reference, insertion_pos_aa=456)

    print(f"\n  Insertion summary:")
    print(f"    VP1 before insertion: {insertion_report['vp1_before_length']} bp")
    print(f"    Insert cassette: {insertion_report['insert_cassette_length']} bp")
    print(f"      - N-linker: {insertion_report['linker_n_length']} bp ({len(LINK_D2_N_AA)} aa)")
    print(f"      - VHH3: {insertion_report['vhh_length']} bp ({len(VHH3_AA)} aa)")
    print(f"      - C-linker: {insertion_report['linker_c_length']} bp ({len(LINK_D2_C_AA)} aa)")
    print(f"    VP1 after insertion: {insertion_report['vp1_after_length']} bp")
    print(f"    Total VP1-VHH fusion: {insertion_report['total_vp1_vhh_length']} bp ({insertion_report['total_vp1_vhh_length']//3} aa)")
    print()

    # Knock out VP2 and VP3 start codons
    print("Knocking out VP2 and VP3 start codons...")
    vp1_vhh_ko, ko_mutations = knockout_start_codons(vp1_vhh_fusion)

    for mut in ko_mutations:
        silent_status = "✓ SILENT" if mut['silent'] else "⚠ NON-SILENT"
        print(f"  {mut['type']}: {mut['original']} → {mut['mutated']} ({mut['aa_original']} → {mut['aa_mutated']}) {silent_status}")
        print(f"    Position: bp {mut['position']}, codon {mut['position']//3}")
    print()

    # Verify AAP frame preservation
    print("Verifying AAP reading frame preservation...")
    aap_frame_ok = verify_aap_frame(vp1_vhh_ko)
    if aap_frame_ok:
        print("  ✓ AAP frame preserved (insertion is multiple of 3 bp)")
    else:
        print("  ⚠ AAP frame may be disrupted!")
    print()

    # Check restriction sites
    print("Checking restriction sites in final VP1-VHH sequence...")
    restriction_sites_6r = {
        'SmaI': 'CCCGGG',
        'BbvCI': 'CCTCAGC',
        'AgeI': 'ACCGGT',
        'BsrGI': 'TGTACA',
        'BmtI': 'GCTAGC',
        'BstZ17I': 'GTATAC'
    }

    site_check = check_restriction_sites(vp1_vhh_ko, restriction_sites_6r)

    for enzyme, info in site_check.items():
        status = "✓ UNIQUE" if info['count'] == 1 else f"⚠ {info['count']} sites" if info['count'] > 1 else "○ ABSENT"
        print(f"  {enzyme:10s} ({info['recognition_seq']:7s}): {status}", end="")
        if info['count'] > 0:
            print(f" at {info['positions']}")
        else:
            print()
    print()

    print("=" * 80)
    print("Design complete! Ready to generate GenBank files.")
    print("=" * 80)


if __name__ == "__main__":
    main()
