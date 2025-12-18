#!/usr/bin/env python3
"""
Phase 1: Pre-Flight Analysis of pGS-scAAV and pGS-ssAAV Backbones
DNA Engineer Agent v2.2 - EF1A-VP1-rBG Assembly Project
"""

from Bio import SeqIO, Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re

print("="*80)
print("DNA ENGINEER AGENT v2.2 - PRE-FLIGHT ANALYSIS")
print("="*80)
print("\nProject: AAV Transfer Plasmid Assembly - EF1A-VP1-rBG Expression Cassette")
print("Date: 2025-12-17")
print("\n" + "="*80)
print("PHASE 1: BACKBONE PLASMID ANALYSIS")
print("="*80)

# Load backbone plasmids
sc_backbone = SeqIO.read("test_data/pGS-scAAV-ITR128-Amp-empty.gb", "genbank")
ss_backbone = SeqIO.read("test_data/pGS-ssAAV-ITR128-Amp-empty.gb", "genbank")

# Load VP1 source
vp1_source = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")

# ITR signature sequences for verification
ITR_LEFT_MOTIF = "GCGCTCGCTCGCTCACTGAGGCC"  # Diagnostic left ITR sequence

def analyze_backbone(record, name):
    """Comprehensive backbone analysis"""
    print(f"\n{'='*80}")
    print(f"ANALYZING: {name}")
    print(f"{'='*80}")
    print(f"\nFile: {record.name}")
    print(f"Length: {len(record)} bp")
    print(f"Topology: {record.annotations.get('topology', 'unknown')}")

    # Find ITRs
    print(f"\n{'-'*80}")
    print("ITR IDENTIFICATION")
    print(f"{'-'*80}")

    itr_features = []
    for feature in record.features:
        label = feature.qualifiers.get('label', [''])[0].lower()
        if 'itr' in label:
            itr_features.append(feature)
            itr_type = "LEFT" if '5' in label or 'left' in label else "RIGHT"
            coords = f"{feature.location.start+1}..{feature.location.end}"
            length = feature.location.end - feature.location.start
            print(f"  {itr_type} ITR: {coords} ({length} bp)")
            print(f"    Label: {feature.qualifiers.get('label', [''])[0]}")

            # Extract and show first/last 30bp of ITR
            itr_seq = str(record.seq[feature.location.start:feature.location.end])
            print(f"    First 30bp: {itr_seq[:30]}")
            print(f"    Last 30bp:  {itr_seq[-30:]}")

    if not itr_features:
        print("  WARNING: No ITR features annotated - searching for ITR motif...")
        # Search for ITR diagnostic sequence
        forward_match = record.seq.find(ITR_LEFT_MOTIF)
        if forward_match != -1:
            print(f"  Found ITR motif at position {forward_match+1}")

    # Determine inter-ITR region
    print(f"\n{'-'*80}")
    print("INTER-ITR REGION")
    print(f"{'-'*80}")

    if len(itr_features) >= 2:
        # Sort by start position
        itr_features_sorted = sorted(itr_features, key=lambda f: f.location.start)
        left_itr = itr_features_sorted[0]
        right_itr = itr_features_sorted[-1]

        inter_itr_start = left_itr.location.end
        inter_itr_end = right_itr.location.start
        inter_itr_length = inter_itr_end - inter_itr_start

        print(f"  Left ITR ends at: {inter_itr_start}")
        print(f"  Right ITR starts at: {inter_itr_end+1}")
        print(f"  Inter-ITR region: {inter_itr_start+1}..{inter_itr_end} ({inter_itr_length} bp)")
        print(f"  Current payload: {inter_itr_length} bp")

        # Check for MCS or stuffer sequence
        inter_itr_seq = str(record.seq[inter_itr_start:inter_itr_end])

        # Look for common restriction sites
        common_sites = {
            'EcoRI': 'GAATTC',
            'BamHI': 'GGATCC',
            'XbaI': 'TCTAGA',
            'NotI': 'GCGGCCGC',
            'HindIII': 'AAGCTT',
            'PstI': 'CTGCAG',
            'SalI': 'GTCGAC',
        }

        print(f"\n  Restriction sites in inter-ITR region:")
        for enzyme, site in common_sites.items():
            count = inter_itr_seq.count(site)
            if count > 0:
                positions = [m.start() + inter_itr_start + 1 for m in re.finditer(site, inter_itr_seq)]
                print(f"    {enzyme}: {count} site(s) at position(s) {positions}")

        return {
            'left_itr': left_itr,
            'right_itr': right_itr,
            'inter_itr_start': inter_itr_start,
            'inter_itr_end': inter_itr_end,
            'inter_itr_length': inter_itr_length
        }
    else:
        print("  ERROR: Could not identify both ITRs")
        return None

    # Check for bacterial backbone elements
    print(f"\n{'-'*80}")
    print("BACTERIAL BACKBONE ELEMENTS")
    print(f"{'-'*80}")

    for feature in record.features:
        if feature.type in ['rep_origin', 'CDS']:
            label = feature.qualifiers.get('label', ['unknown'])[0]
            if 'ori' in label.lower() or 'amp' in label.lower():
                coords = f"{feature.location.start+1}..{feature.location.end}"
                print(f"  {feature.type}: {label} at {coords}")

# Analyze both backbones
sc_info = analyze_backbone(sc_backbone, "pGS-scAAV-ITR128-Amp-empty.gb")
ss_info = analyze_backbone(ss_backbone, "pGS-ssAAV-ITR128-Amp-empty.gb")

# Analyze VP1 source
print(f"\n{'='*80}")
print("VP1 SOURCE ANALYSIS")
print(f"{'='*80}")

vp1_features = [f for f in vp1_source.features if f.type == "CDS" and
                f.qualifiers.get('label', [''])[0] == "VP1"]

if vp1_features:
    vp1_feature = vp1_features[0]
    vp1_coords = f"{vp1_feature.location.start+1}..{vp1_feature.location.end}"
    vp1_length = vp1_feature.location.end - vp1_feature.location.start
    vp1_translation = vp1_feature.qualifiers.get('translation', [''])[0]

    print(f"\nVP1 CDS found:")
    print(f"  Coordinates: {vp1_coords}")
    print(f"  Length: {vp1_length} bp ({len(vp1_translation)} aa)")
    print(f"  First 60 aa: {vp1_translation[:60]}")

    # Extract the nucleotide sequence
    vp1_nt_seq = str(vp1_source.seq[vp1_feature.location.start:vp1_feature.location.end])
    print(f"\n  Start codon context: {vp1_nt_seq[:12]} (should be ATG)")
    print(f"  Stop codon: {vp1_nt_seq[-3:]} (should be TAA/TAG/TGA)")

    # Check Kozak sequence
    upstream_context = str(vp1_source.seq[vp1_feature.location.start-6:vp1_feature.location.start+6])
    print(f"  Kozak context (-6 to +6): {upstream_context}")
    print(f"    Optimal Kozak: (G/A)CCATGG")

    # Check for internal stop codons
    test_translation = Seq.Seq(vp1_nt_seq).translate()
    internal_stops = test_translation[:-1].count('*')
    print(f"  Internal stop codons: {internal_stops}")
    if internal_stops > 0:
        print(f"    WARNING: VP1 ORF contains internal stops!")
else:
    print("\nERROR: VP1 CDS not found in source file")

# Size calculations
print(f"\n{'='*80}")
print("TRANSGENE CASSETTE SIZE ESTIMATION")
print(f"{'='*80}")

EF1A_FULL_LENGTH = 1172  # Full EF1α promoter with intron
EF1A_CORE_LENGTH = 212   # Core promoter only (no intron)
VP1_LENGTH = 2211
RBG_POLYA_LENGTH = 127
KOZAK_LENGTH = 9

cassette_full = EF1A_FULL_LENGTH + KOZAK_LENGTH + VP1_LENGTH + RBG_POLYA_LENGTH
cassette_core = EF1A_CORE_LENGTH + KOZAK_LENGTH + VP1_LENGTH + RBG_POLYA_LENGTH

print(f"\nComponent sizes:")
print(f"  EF1α promoter (full with intron): {EF1A_FULL_LENGTH} bp")
print(f"  EF1α promoter (core only): {EF1A_CORE_LENGTH} bp")
print(f"  Kozak sequence: {KOZAK_LENGTH} bp")
print(f"  VP1 ORF: {VP1_LENGTH} bp")
print(f"  rBG polyA: {RBG_POLYA_LENGTH} bp")

print(f"\nCassette size (full EF1α): {cassette_full} bp")
print(f"Cassette size (core EF1α): {cassette_core} bp")

print(f"\nPackaging capacity check:")
print(f"  scAAV limit: 2400 bp")
print(f"  ssAAV limit: 4700 bp")

if cassette_full > 2400:
    print(f"  ⚠️ FULL cassette ({cassette_full} bp) EXCEEDS scAAV limit by {cassette_full - 2400} bp")
    print(f"     Recommendation: Use core EF1α promoter for scAAV")

if cassette_core <= 2400:
    print(f"  ✅ CORE cassette ({cassette_core} bp) fits within scAAV limit")

if cassette_full <= 4700:
    print(f"  ✅ FULL cassette ({cassette_full} bp) fits within ssAAV limit")

print(f"\n{'='*80}")
print("PRE-FLIGHT ANALYSIS COMPLETE")
print(f"{'='*80}")
print("\nNext steps:")
print("  1. Retrieve EF1α promoter sequence")
print("  2. Retrieve rBG polyA sequence")
print("  3. Assemble cassette components")
print("  4. Write assembly script for both backbones")
print("  5. Generate constructs with verification")
