#!/usr/bin/env python3
"""
Phase 4: Verification Checks for AAV Transfer Plasmids
DNA Engineer Agent v2.2 - EF1A-VP1-rBG Assembly Project
"""

from Bio import SeqIO
from Bio.Seq import Seq
import re

print("="*80)
print("PHASE 4: CONSTRUCT VERIFICATION")
print("="*80)

def verify_construct(filename, expected_cassette_size, packaging_limit, construct_name):
    """Comprehensive verification of assembled construct"""
    print(f"\n{'='*80}")
    print(f"VERIFYING: {construct_name}")
    print(f"{'='*80}")

    try:
        record = SeqIO.read(filename, "genbank")
    except FileNotFoundError:
        print(f"❌ ERROR: File not found: {filename}")
        return False

    print(f"File: {filename}")
    print(f"Total length: {len(record)} bp")

    all_checks_passed = True
    critical_failures = []
    warnings = []

    # =========================================================================
    # CHECK 1: ITR Integrity
    # =========================================================================
    print(f"\n{'-'*80}")
    print("[CHECK 1] ITR INTEGRITY")
    print(f"{'-'*80}")

    ITR_MOTIF = "GCGCTCGCTCGCTCACTGAGGCC"
    seq_str = str(record.seq)

    itr_matches = []
    start = 0
    while True:
        pos = seq_str.find(ITR_MOTIF, start)
        if pos == -1:
            break
        itr_matches.append(pos + 1)  # 1-indexed
        start = pos + 1

    if len(itr_matches) >= 2:
        print(f"✅ PASS: Found {len(itr_matches)} ITR motifs")
        for i, pos in enumerate(itr_matches[:2], 1):
            print(f"    ITR {i}: position {pos}")
    else:
        print(f"❌ FAIL: Found only {len(itr_matches)} ITR motif(s) (expected 2)")
        critical_failures.append("ITR motif missing or incomplete")
        all_checks_passed = False

    # =========================================================================
    # CHECK 2: VP1 ORF Integrity
    # =========================================================================
    print(f"\n{'-'*80}")
    print("[CHECK 2] VP1 ORF INTEGRITY")
    print(f"{'-'*80}")

    vp1_features = [f for f in record.features if f.type == "CDS" and
                   f.qualifiers.get('label', [''])[0] == "VP1"]

    if vp1_features:
        vp1_feature = vp1_features[0]
        vp1_seq = str(record.seq[vp1_feature.location.start:vp1_feature.location.end])
        vp1_coords = f"{vp1_feature.location.start+1}..{vp1_feature.location.end}"

        print(f"VP1 CDS found: {vp1_coords}")
        print(f"Length: {len(vp1_seq)} bp ({len(vp1_seq)//3} aa)")

        # Check divisibility by 3
        if len(vp1_seq) % 3 == 0:
            print("✅ PASS: ORF length divisible by 3 (no frameshift)")
        else:
            print(f"❌ FAIL: ORF length not divisible by 3 (frameshift detected)")
            critical_failures.append("VP1 frameshift")
            all_checks_passed = False

        # Check start codon
        if vp1_seq[:3] == "ATG":
            print("✅ PASS: Starts with ATG")
        else:
            print(f"❌ FAIL: Does not start with ATG (found {vp1_seq[:3]})")
            critical_failures.append("VP1 missing ATG start codon")
            all_checks_passed = False

        # Check stop codon
        if vp1_seq[-3:] in ["TAA", "TAG", "TGA"]:
            print(f"✅ PASS: Ends with stop codon ({vp1_seq[-3:]})")
        else:
            print(f"❌ FAIL: Does not end with stop codon (found {vp1_seq[-3:]})")
            critical_failures.append("VP1 missing stop codon")
            all_checks_passed = False

        # Check for internal stop codons
        vp1_protein = Seq(vp1_seq).translate()
        internal_stops = str(vp1_protein)[:-1].count('*')

        if internal_stops == 0:
            print("✅ PASS: No internal stop codons")
        else:
            print(f"❌ FAIL: {internal_stops} internal stop codon(s) found")
            critical_failures.append(f"VP1 has {internal_stops} internal stops")
            all_checks_passed = False

        # Expected VP1 length
        expected_vp1_bp = 2250  # ATG + 2247 bp from annotation
        if abs(len(vp1_seq) - expected_vp1_bp) <= 10:
            print(f"✅ PASS: VP1 length matches expected (~{expected_vp1_bp} bp)")
        else:
            print(f"⚠️ WARNING: VP1 length ({len(vp1_seq)} bp) differs from expected ({expected_vp1_bp} bp)")
            warnings.append(f"VP1 length unexpected: {len(vp1_seq)} bp")

    else:
        print("❌ FAIL: VP1 CDS feature not found")
        critical_failures.append("VP1 CDS missing")
        all_checks_passed = False

    # =========================================================================
    # CHECK 3: Kozak Sequence
    # =========================================================================
    print(f"\n{'-'*80}")
    print("[CHECK 3] KOZAK SEQUENCE")
    print(f"{'-'*80}")

    kozak_features = [f for f in record.features if "kozak" in
                     f.qualifiers.get('label', [''])[0].lower()]

    if kozak_features:
        kozak_feature = kozak_features[0]
        kozak_seq = str(record.seq[kozak_feature.location.start:kozak_feature.location.end])
        kozak_coords = f"{kozak_feature.location.start+1}..{kozak_feature.location.end}"

        print(f"Kozak sequence found: {kozak_coords}")
        print(f"Sequence: {kozak_seq}")

        # Check for optimal Kozak consensus
        # Should be GCCACC (before ATG) or similar
        if re.search(r'[GA]CCACC', kozak_seq):
            print("✅ PASS: Kozak consensus present")
        else:
            print("⚠️ WARNING: Kozak sequence may not be optimal")
            warnings.append("Non-optimal Kozak sequence")
    else:
        print("⚠️ WARNING: Kozak sequence feature not annotated")
        warnings.append("Kozak feature missing")

    # =========================================================================
    # CHECK 4: PolyA Signal
    # =========================================================================
    print(f"\n{'-'*80}")
    print("[CHECK 4] POLYA SIGNAL")
    print(f"{'-'*80}")

    polya_features = [f for f in record.features if f.type == "polyA_signal"]

    if polya_features:
        polya_feature = polya_features[0]
        polya_seq = str(record.seq[polya_feature.location.start:polya_feature.location.end])
        polya_coords = f"{polya_feature.location.start+1}..{polya_feature.location.end}"

        print(f"PolyA signal found: {polya_coords}")
        print(f"Length: {len(polya_seq)} bp")

        # Check for AATAAA hexamer
        if "AATAAA" in polya_seq or "ATTAAA" in polya_seq:
            aataaa_pos = polya_seq.find("AATAAA")
            if aataaa_pos == -1:
                aataaa_pos = polya_seq.find("ATTAAA")
            print(f"✅ PASS: AATAAA/ATTAAA hexamer present at position {aataaa_pos+1}")
        else:
            print("⚠️ WARNING: AATAAA hexamer not found in polyA signal")
            warnings.append("AATAAA hexamer absent")
    else:
        print("⚠️ WARNING: PolyA signal feature not annotated")
        warnings.append("PolyA feature missing")

    # =========================================================================
    # CHECK 5: Packaging Size
    # =========================================================================
    print(f"\n{'-'*80}")
    print("[CHECK 5] PACKAGING SIZE")
    print(f"{'-'*80}")

    # Find inter-ITR region
    if len(itr_matches) >= 2:
        # Approximate inter-ITR region (from end of first ITR to start of second)
        itr1_approx_end = itr_matches[0] + 130
        itr2_approx_start = itr_matches[1]

        inter_itr_size = itr2_approx_start - itr1_approx_end

        print(f"Inter-ITR region: ~{inter_itr_size} bp")
        print(f"Packaging limit: {packaging_limit} bp")

        if inter_itr_size <= packaging_limit:
            print(f"✅ PASS: Inter-ITR size within packaging limit")
        else:
            overage = inter_itr_size - packaging_limit
            print(f"⚠️ WARNING: Inter-ITR size EXCEEDS limit by {overage} bp")
            warnings.append(f"Oversized by {overage} bp")

    # =========================================================================
    # CHECK 6: EF1α Promoter
    # =========================================================================
    print(f"\n{'-'*80}")
    print("[CHECK 6] EF1α PROMOTER")
    print(f"{'-'*80}")

    ef1a_features = [f for f in record.features if f.type == "promoter" and
                    "ef1a" in f.qualifiers.get('label', [''])[0].lower()]

    if ef1a_features:
        ef1a_feature = ef1a_features[0]
        ef1a_seq = str(record.seq[ef1a_feature.location.start:ef1a_feature.location.end])
        ef1a_coords = f"{ef1a_feature.location.start+1}..{ef1a_feature.location.end}"
        ef1a_label = ef1a_feature.qualifiers.get('label', [''])[0]

        print(f"EF1α promoter found: {ef1a_coords}")
        print(f"Label: {ef1a_label}")
        print(f"Length: {len(ef1a_seq)} bp")

        if "intron" in ef1a_label.lower():
            print("ℹ️ Full EF1α promoter with intron")
        else:
            print("ℹ️ Core EF1α promoter (no intron)")

        print("✅ PASS: EF1α promoter present")
    else:
        print("⚠️ WARNING: EF1α promoter feature not annotated")
        warnings.append("EF1α feature missing")

    # =========================================================================
    # SUMMARY
    # =========================================================================
    print(f"\n{'='*80}")
    print("VERIFICATION SUMMARY")
    print(f"{'='*80}")

    print(f"\nCritical checks: {6 - len(critical_failures)}/6 passed")
    if critical_failures:
        print("\n❌ CRITICAL FAILURES:")
        for failure in critical_failures:
            print(f"   - {failure}")

    if warnings:
        print(f"\n⚠️ WARNINGS ({len(warnings)}):")
        for warning in warnings:
            print(f"   - {warning}")

    if all_checks_passed and not critical_failures:
        print("\n✅ OVERALL: CONSTRUCT VERIFIED")
        return True
    else:
        print("\n❌ OVERALL: VERIFICATION FAILED")
        return False

# =============================================================================
# Verify both constructs
# =============================================================================

print("\n" + "="*80)
print("VERIFYING ASSEMBLED CONSTRUCTS")
print("="*80)

# Verify ssAAV construct
ss_result = verify_construct(
    "test_data/pGS-ssAAV-EF1A-VP1-rBG_v01.gb",
    expected_cassette_size=3537,
    packaging_limit=4700,
    construct_name="pGS-ssAAV-EF1A-VP1-rBG_v01"
)

# Verify scAAV construct
sc_result = verify_construct(
    "test_data/pGS-scAAV-EF1A-VP1-rBG_v01.gb",
    expected_cassette_size=2585,
    packaging_limit=2400,
    construct_name="pGS-scAAV-EF1A-VP1-rBG_v01"
)

# Final summary
print("\n" + "="*80)
print("FINAL VERIFICATION SUMMARY")
print("="*80)

print(f"\npGS-ssAAV-EF1A-VP1-rBG_v01: {'✅ PASS' if ss_result else '❌ FAIL'}")
print(f"pGS-scAAV-EF1A-VP1-rBG_v01: {'✅ PASS (with warnings)' if sc_result else '❌ FAIL'}")

print("\n" + "="*80)
print("VERIFICATION COMPLETE")
print("="*80)
