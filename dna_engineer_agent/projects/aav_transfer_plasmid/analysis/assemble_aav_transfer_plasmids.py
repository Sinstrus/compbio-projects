#!/usr/bin/env python3
"""
Phase 3: Assembly Script for AAV Transfer Plasmids
DNA Engineer Agent v2.2 - EF1A-VP1-rBG Assembly Project

Generates:
1. pGS-scAAV-EF1A-VP1-rBG_v01.gb (with core EF1α promoter - WARNING: oversized)
2. pGS-ssAAV-EF1A-VP1-rBG_v01.gb (with full EF1α promoter + intron)
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from datetime import datetime
import cassette_components as comp

print("="*80)
print("PHASE 3: AAV TRANSFER PLASMID ASSEMBLY")
print("="*80)
print(f"Date: {datetime.now().strftime('%Y-%m-%d')}")
print(f"Agent: DNA Engineer v2.2")

# Load backbone plasmids
sc_backbone = SeqIO.read("test_data/pGS-scAAV-ITR128-Amp-empty.gb", "genbank")
ss_backbone = SeqIO.read("test_data/pGS-ssAAV-ITR128-Amp-empty.gb", "genbank")

def find_itr_boundaries(record):
    """Find ITR coordinates in the backbone"""
    itr_features = []
    for feature in record.features:
        label = feature.qualifiers.get('label', [''])[0].lower()
        if 'itr' in label:
            itr_features.append(feature)

    if len(itr_features) >= 2:
        # Sort by start position
        itr_features_sorted = sorted(itr_features, key=lambda f: f.location.start)
        left_itr = itr_features_sorted[0]
        right_itr = itr_features_sorted[-1]

        inter_itr_start = left_itr.location.end
        inter_itr_end = right_itr.location.start

        return {
            'left_itr': left_itr,
            'right_itr': right_itr,
            'inter_itr_start': inter_itr_start,
            'inter_itr_end': inter_itr_end,
        }

    # If ITRs not annotated, search for the diagnostic sequence
    print("  ITRs not annotated, searching for ITR motif...")
    ITR_MOTIF = "GCGCTCGCTCGCTCACTGAGGCC"

    # Search for motif
    seq_str = str(record.seq)
    first_match = seq_str.find(ITR_MOTIF)

    if first_match != -1:
        # Look for second occurrence
        second_match = seq_str.find(ITR_MOTIF, first_match + 1)
        if second_match != -1:
            print(f"  Found ITR motifs at positions {first_match+1} and {second_match+1}")

            # Assume ITRs are ~130 bp and ~90 bp
            # Based on the scAAV backbone analysis
            left_itr_end = first_match + 130
            right_itr_start = second_match

            # Create pseudo-features
            left_itr = create_feature("misc_feature", first_match, left_itr_end, "5ITR", "Detected ITR")
            right_itr = create_feature("misc_feature", right_itr_start, second_match + 90, "3ITR", "Detected ITR")

            return {
                'left_itr': left_itr,
                'right_itr': right_itr,
                'inter_itr_start': left_itr_end,
                'inter_itr_end': right_itr_start,
            }

    return None

def create_feature(feature_type, start, end, label, note="", strand=1, **kwargs):
    """Create a SeqFeature with standard qualifiers"""
    location = FeatureLocation(start, end, strand=strand)
    qualifiers = {'label': [label]}
    if note:
        qualifiers['note'] = [note]
    for key, value in kwargs.items():
        qualifiers[key] = value if isinstance(value, list) else [value]
    return SeqFeature(location, type=feature_type, qualifiers=qualifiers)

def assemble_construct(backbone, cassette_seq, promoter_type, backbone_name):
    """Assemble the expression cassette into the backbone"""
    print(f"\n{'='*80}")
    print(f"ASSEMBLING: {backbone_name}")
    print(f"{'='*80}")

    # Find ITR boundaries
    itr_info = find_itr_boundaries(backbone)
    if not itr_info:
        print("ERROR: Could not find ITRs in backbone")
        return None

    left_itr = itr_info['left_itr']
    right_itr = itr_info['right_itr']
    inter_itr_start = itr_info['inter_itr_start']
    inter_itr_end = itr_info['inter_itr_end']

    print(f"Backbone: {len(backbone)} bp")
    print(f"Left ITR: {left_itr.location.start+1}..{left_itr.location.end}")
    print(f"Right ITR: {right_itr.location.start+1}..{right_itr.location.end}")
    print(f"Inter-ITR region: {inter_itr_start+1}..{inter_itr_end}")

    # Build new sequence
    # Keep: backbone up to end of left ITR
    # Insert: cassette
    # Keep: right ITR to end of backbone
    new_seq = (backbone.seq[:inter_itr_start] +
               Seq(cassette_seq) +
               backbone.seq[inter_itr_end:])

    print(f"\nNew construct length: {len(new_seq)} bp")
    print(f"Cassette inserted: {len(cassette_seq)} bp")

    # Create new record
    new_record = SeqRecord(
        new_seq,
        id=backbone.id,
        name=backbone.name,
        description=f"AAV transfer plasmid with EF1A-VP1-rBG expression cassette",
        annotations=backbone.annotations.copy()
    )
    new_record.annotations['topology'] = 'circular'

    # Add features
    # 1. Keep all original backbone features, adjusting positions after insertion
    cassette_length = len(cassette_seq)
    insertion_point = inter_itr_start

    for feature in backbone.features:
        # Keep features before insertion point as-is
        if feature.location.end <= insertion_point:
            new_record.features.append(feature)
        # Shift features after insertion point
        elif feature.location.start >= inter_itr_end:
            new_start = feature.location.start - (inter_itr_end - inter_itr_start) + cassette_length
            new_end = feature.location.end - (inter_itr_end - inter_itr_start) + cassette_length
            new_loc = FeatureLocation(new_start, new_end, strand=feature.location.strand)
            new_feature = SeqFeature(new_loc, type=feature.type, qualifiers=feature.qualifiers.copy())
            new_record.features.append(new_feature)

    # 2. Add cassette component features
    current_pos = inter_itr_start

    # EF1α Promoter
    if promoter_type == "FULL":
        promoter_seq = comp.EF1A_FULL
        promoter_label = "EF1A_promoter_with_intron"
        promoter_note = (f"Human EEF1A1 promoter with first intron ({len(promoter_seq)} bp). "
                        "[Action]: Inserted | [Reason]: Strong ubiquitous mammalian promoter "
                        "with intron-mediated enhancement | [Risk]: Low | "
                        f"[Date]: {datetime.now().strftime('%Y-%m-%d')}")
    else:  # CORE
        promoter_seq = comp.EF1A_CORE
        promoter_label = "EF1A_core_promoter"
        promoter_note = (f"Human EEF1A1 core promoter ({len(promoter_seq)} bp). "
                        "[Action]: Inserted | [Reason]: Compact strong promoter for scAAV | "
                        "[Risk]: Low | "
                        f"[Date]: {datetime.now().strftime('%Y-%m-%d')}")

    promoter_feature = create_feature(
        "promoter",
        current_pos,
        current_pos + len(promoter_seq),
        promoter_label,
        note=promoter_note
    )
    new_record.features.append(promoter_feature)
    current_pos += len(promoter_seq)

    # Kozak sequence (part of 5'UTR)
    kozak_full = comp.KOZAK + "ATG"  # GCCACCATG
    kozak_feature = create_feature(
        "misc_feature",
        current_pos,
        current_pos + len(comp.KOZAK),
        "Kozak_sequence",
        note="Kozak consensus for optimal translation initiation"
    )
    new_record.features.append(kozak_feature)
    current_pos += len(comp.KOZAK)

    # VP1 CDS (including the ATG from Kozak)
    vp1_cds_with_atg = "ATG" + comp.VP1_ORF
    vp1_protein = Seq(vp1_cds_with_atg).translate()

    vp1_feature = create_feature(
        "CDS",
        current_pos,  # Start at the ATG from Kozak
        current_pos + len(vp1_cds_with_atg),
        "VP1",
        note=(f"AAV9 VP1 capsid protein coding sequence. "
              f"[Action]: Inserted | [Reason]: VP1 expression cassette for protein production | "
              f"[Risk]: Low | [Date]: {datetime.now().strftime('%Y-%m-%d')}"),
        translation=str(vp1_protein)[:-1],  # Remove stop codon asterisk
        codon_start=1
    )
    new_record.features.append(vp1_feature)
    current_pos += len(vp1_cds_with_atg)

    # Rabbit β-globin polyA
    polya_feature = create_feature(
        "polyA_signal",
        current_pos,
        current_pos + len(comp.RBG_POLYA),
        "rBG_polyA",
        note=("Rabbit beta-globin polyadenylation signal. "
              "[Action]: Inserted | [Reason]: Efficient transcription termination "
              "and mRNA stabilization | [Risk]: Low | "
              f"[Date]: {datetime.now().strftime('%Y-%m-%d')}")
    )
    new_record.features.append(polya_feature)

    # Sort features by start position
    new_record.features.sort(key=lambda f: f.location.start)

    return new_record

# =============================================================================
# ASSEMBLE ssAAV CONSTRUCT (Full EF1α + Intron)
# =============================================================================
print("\n" + "="*80)
print("CONSTRUCT 1: ssAAV with FULL EF1α Promoter")
print("="*80)

cassette_full = comp.EF1A_FULL + comp.KOZAK + "ATG" + comp.VP1_ORF + comp.RBG_POLYA
print(f"Full cassette: {len(cassette_full)} bp")
print(f"ssAAV packaging limit: 4700 bp")

if len(cassette_full) <= 4700:
    print("✅ Cassette fits within ssAAV limit")
    ss_construct = assemble_construct(ss_backbone, cassette_full, "FULL", "ssAAV")

    if ss_construct:
        output_file_ss = "test_data/pGS-ssAAV-EF1A-VP1-rBG_v01.gb"
        ss_construct.name = "pGS-ssAAV-EF1A-VP1-rBG_v01"
        ss_construct.id = "pGS-ssAAV-EF1A-VP1-rBG_v01"

        # Add version history to comment
        comment = f"""
VERSION HISTORY:
v01 - {datetime.now().strftime('%Y-%m-%d')}
  - Initial assembly
  - Inserted EF1α-VP1-rBG expression cassette ({len(cassette_full)} bp)
  - EF1α promoter with first intron (1194 bp)
  - AAV9 VP1 ORF (2211 bp + ATG)
  - Rabbit β-globin polyA signal (123 bp)
  - Agent: DNA Engineer v2.2
"""
        ss_construct.annotations['comment'] = comment

        SeqIO.write(ss_construct, output_file_ss, "genbank")
        print(f"\n✅ Saved: {output_file_ss}")
        print(f"   Total size: {len(ss_construct)} bp")
        print(f"   Cassette size: {len(cassette_full)} bp")
else:
    print(f"❌ Cassette ({len(cassette_full)} bp) exceeds ssAAV limit")

# =============================================================================
# ASSEMBLE scAAV CONSTRUCT (Core EF1α - NO Intron)
# =============================================================================
print("\n" + "="*80)
print("CONSTRUCT 2: scAAV with CORE EF1α Promoter")
print("="*80)

cassette_core = comp.EF1A_CORE + comp.KOZAK + "ATG" + comp.VP1_ORF + comp.RBG_POLYA
print(f"Core cassette: {len(cassette_core)} bp")
print(f"scAAV packaging limit: 2400 bp")

if len(cassette_core) > 2400:
    print(f"⚠️ WARNING: Cassette ({len(cassette_core)} bp) EXCEEDS scAAV limit by {len(cassette_core) - 2400} bp")
    print("   This construct will be generated but may not package efficiently.")
    print("   Consider:")
    print("   - Using a shorter promoter (e.g., CAG, SV40)")
    print("   - Truncating VP1 (e.g., VP2 or VP3 region)")
    print("   - Using ssAAV backbone instead")

sc_construct = assemble_construct(sc_backbone, cassette_core, "CORE", "scAAV")

if sc_construct:
    output_file_sc = "test_data/pGS-scAAV-EF1A-VP1-rBG_v01.gb"
    sc_construct.name = "pGS-scAAV-EF1A-VP1-rBG_v01"
    sc_construct.id = "pGS-scAAV-EF1A-VP1-rBG_v01"

    # Add version history with WARNING
    comment = f"""
VERSION HISTORY:
v01 - {datetime.now().strftime('%Y-%m-%d')}
  - Initial assembly
  - Inserted EF1α-VP1-rBG expression cassette ({len(cassette_core)} bp)
  - EF1α core promoter (242 bp)
  - AAV9 VP1 ORF (2211 bp + ATG)
  - Rabbit β-globin polyA signal (123 bp)
  - Agent: DNA Engineer v2.2

⚠️ WARNING: OVERSIZED FOR scAAV PACKAGING
  - Cassette: {len(cassette_core)} bp
  - Limit: 2400 bp
  - Overage: {len(cassette_core) - 2400} bp
  - This construct may not package efficiently as self-complementary AAV
  - Recommendation: Use ssAAV backbone or reduce cassette size
"""
    sc_construct.annotations['comment'] = comment

    SeqIO.write(sc_construct, output_file_sc, "genbank")
    print(f"\n⚠️ Saved (with warning): {output_file_sc}")
    print(f"   Total size: {len(sc_construct)} bp")
    print(f"   Cassette size: {len(cassette_core)} bp")
    print(f"   OVERAGE: {len(cassette_core) - 2400} bp over scAAV limit")

print("\n" + "="*80)
print("ASSEMBLY COMPLETE")
print("="*80)
print("\nGenerated files:")
print(f"  1. test_data/pGS-ssAAV-EF1A-VP1-rBG_v01.gb")
print(f"  2. test_data/pGS-scAAV-EF1A-VP1-rBG_v01.gb (⚠️ OVERSIZED)")
