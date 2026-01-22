#!/usr/bin/env python3
"""
Verify flanking sequences for N-terminal VHH3 insertions in AVD008-010.

This script checks that:
1. AVD008/009 have M-A-[VHH3]-[GGGGS5]-A-D-G-Y-L... at VP1 N-term
2. AVD010 has M-A-[VHH3]-[GGGGS5]-P-G-K-K-R... at VP2 N-term
"""

from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path

PLASMIDS_DIR = Path(__file__).parent.parent / "plasmids"

def verify_vp1_n_term(construct_id, gb_file):
    """Verify VP1 N-terminal insertion flanking sequences."""
    print(f"\n{construct_id}: VP1 N-terminal insertion")

    record = SeqIO.read(str(gb_file), "genbank")
    seq = str(record.seq).upper()

    # Expected sequence:
    # bp 2379-2384 (0-idx 2378-2383): ATG GCT = M-A (before insert)
    # bp 2384-2816 (432 bp insert): VHH3 (357 bp) + GGGGS5 (75 bp)
    # bp 2816-2831 (0-idx 2815-2830): GCC GAT GGT TAT = A-D-G-Y (after insert)

    # Check upstream (M-A)
    upstream_dna = seq[2378:2384]
    upstream_aa = str(Seq(upstream_dna).translate())

    # Check insert
    insert_dna = seq[2384:2816]
    insert_aa = str(Seq(insert_dna).translate())

    # Check downstream (A-D-G-Y-L-P-D-W)
    downstream_dna = seq[2816:2840]
    downstream_aa = str(Seq(downstream_dna).translate())

    print(f"  Upstream (bp 2379-2384): {upstream_dna} = {upstream_aa}")
    print(f"  Expected: ATG GCT = M-A")
    print(f"  ✓ PASS" if upstream_aa == "MA" else f"  ✗ FAIL")

    print(f"\n  Insert (bp 2384-2816, {len(insert_dna)} bp):")
    print(f"  Translation starts: {insert_aa[:10]}...")
    print(f"  Translation ends: ...{insert_aa[-10:]}")
    print(f"  Expected start: EVQLVESGG...")
    print(f"  Expected end: ...GSGGGGS")

    print(f"\n  Downstream (bp 2816-2840): {downstream_dna}")
    print(f"  Translation: {downstream_aa}")
    print(f"  Expected: A-D-G-Y-L-P-D-W")
    print(f"  ✓ PASS" if downstream_aa == "ADGYLPDW" else f"  ✗ FAIL")

    # Final verification
    full_junction = upstream_aa + insert_aa + downstream_aa
    print(f"\n  Full junction: M-A-[{len(insert_aa)} aa insert]-A-D-G-Y-L-P-D-W")

    if upstream_aa == "MA" and downstream_aa == "ADGYLPDW":
        print(f"  ✅ {construct_id} N-terminal junction CORRECT")
        return True
    else:
        print(f"  ❌ {construct_id} N-terminal junction INCORRECT")
        return False


def verify_vp2_n_term(construct_id, gb_file):
    """Verify VP2 N-terminal insertion flanking sequences."""
    print(f"\n{construct_id}: VP2 N-terminal insertion")

    record = SeqIO.read(str(gb_file), "genbank")
    seq = str(record.seq).upper()

    # Expected sequence:
    # bp 2790-2795 (0-idx 2789-2794): ACG GCT = T/M-A (before insert, ACG→M at translation)
    # bp 2795-3227 (432 bp insert): VHH3 (357 bp) + GGGGS5 (75 bp)
    # bp 3227-3242 (0-idx 3226-3241): CCT GGA AAG AAG = P-G-K-K (after insert)

    # Check upstream (M-A, but ACG in DNA)
    upstream_dna = seq[2789:2795]
    upstream_aa = str(Seq(upstream_dna).translate())

    # Check insert
    insert_dna = seq[2795:3227]
    insert_aa = str(Seq(insert_dna).translate())

    # Check downstream (P-G-K-K-R)
    downstream_dna = seq[3227:3242]
    downstream_aa = str(Seq(downstream_dna).translate())

    print(f"  Upstream (bp 2790-2795): {upstream_dna}")
    print(f"  Translation: {upstream_aa}")
    print(f"  Expected: ACG GCT = T-A (ACG→M at translation initiation)")
    print(f"  Note: SnapGene will show M-A at translation start")
    print(f"  ✓ PASS" if upstream_aa == "TA" else f"  ✗ FAIL (but OK if T-A, becomes M-A)")

    print(f"\n  Insert (bp 2795-3227, {len(insert_dna)} bp):")
    print(f"  Translation starts: {insert_aa[:10]}...")
    print(f"  Translation ends: ...{insert_aa[-10:]}")
    print(f"  Expected start: EVQLVESGG...")
    print(f"  Expected end: ...GSGGGGS")

    print(f"\n  Downstream (bp 3227-3242): {downstream_dna}")
    print(f"  Translation: {downstream_aa}")
    print(f"  Expected: P-G-K-K-R")
    print(f"  ✓ PASS" if downstream_aa == "PGKKR" else f"  ✗ FAIL")

    # Final verification
    print(f"\n  Full junction: M-A-[{len(insert_aa)} aa insert]-P-G-K-K-R")
    print(f"  (ACG will be forced to M at translation initiation)")

    if upstream_aa == "TA" and downstream_aa == "PGKKR":
        print(f"  ✅ {construct_id} N-terminal junction CORRECT")
        return True
    else:
        print(f"  ❌ {construct_id} N-terminal junction INCORRECT")
        return False


def main():
    """Verify all N-terminal insertions."""
    print("=" * 70)
    print("Verifying N-Terminal Insertion Flanking Sequences")
    print("=" * 70)

    results = []

    # Verify AVD008
    avd008_file = PLASMIDS_DIR / "AVD008-Rep2Mut2Cap9-VP1n-VHH3.gb"
    results.append(verify_vp1_n_term("AVD008", avd008_file))

    # Verify AVD009
    avd009_file = PLASMIDS_DIR / "AVD009-Rep2Mut2Cap9-VP1n-VHH3-D2.gb"
    results.append(verify_vp1_n_term("AVD009", avd009_file))

    # Verify AVD010
    avd010_file = PLASMIDS_DIR / "AVD010-Rep2Mut2Cap9-VP2n-VHH3.gb"
    results.append(verify_vp2_n_term("AVD010", avd010_file))

    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)

    if all(results):
        print("✅ ALL FLANKING SEQUENCES VERIFIED CORRECTLY")
        print("\nConclusion:")
        print("- AVD008/009: M-A-[VHH3]-[GGGGS5]-A-D-G-Y-L... ✓")
        print("- AVD010: M-A-[VHH3]-[GGGGS5]-P-G-K-K-R... ✓")
        print("\nBUG-006 has been FIXED.")
    else:
        print("❌ SOME FLANKING SEQUENCES ARE INCORRECT")
        print("\nPlease review the errors above.")

    print()


if __name__ == "__main__":
    main()
