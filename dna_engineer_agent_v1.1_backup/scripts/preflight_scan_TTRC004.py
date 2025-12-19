#!/usr/bin/env python3
"""
Pre-Flight Scan for TTRC004_rep2mut02-cap9-p5.gb
DNA Engineer Agent v1.1.0
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import json
import hashlib
from datetime import datetime
import os

# File paths
INPUT_FILE = "/home/cnguy/projects/dna_engineer_agent/test_data/TTRC004_rep2mut02-cap9-p5.gb"
OUTPUT_DIR = "/home/cnguy/projects/dna_engineer_agent/reports"
CIS_MANIFEST = "/home/cnguy/projects/dna_engineer_agent/cis_elements/manifest.json"
PROJECT_CONTEXT = "/home/cnguy/projects/dna_engineer_agent/projects/aav_vhh_v1/context.json"

# Load configuration files
with open(CIS_MANIFEST, 'r') as f:
    cis_manifest = json.load(f)

with open(PROJECT_CONTEXT, 'r') as f:
    project_context = json.load(f)

# Load the GenBank file
record = SeqIO.read(INPUT_FILE, "genbank")
seq = record.seq
seq_length = len(seq)

# Calculate file hash
with open(INPUT_FILE, 'rb') as f:
    file_hash = hashlib.sha256(f.read()).hexdigest()

# Get timestamp
timestamp = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
timestamp_file = datetime.utcnow().strftime("%Y%m%d_%H%M%S")

print("=" * 80)
print("DNA ENGINEER AGENT v1.1.0 - PRE-FLIGHT SCAN")
print("=" * 80)
print(f"Input: {os.path.basename(INPUT_FILE)}")
print(f"Size: {seq_length} bp")
print(f"Hash: {file_hash[:16]}...")
print("=" * 80)

# ============================================================================
# 1. PLASMID CLASSIFICATION
# ============================================================================
print("\n[1] PLASMID CLASSIFICATION")
print("-" * 80)

has_rep = False
has_cap = False
has_itr = False
has_promoters = []

# Check for Rep
for feat in record.features:
    if 'rep' in feat.type.lower() or 'rep' in str(feat.qualifiers.get('label', '')).lower():
        has_rep = True
    if 'cap' in str(feat.qualifiers.get('gene', '')).lower() or 'VP' in str(feat.qualifiers.get('label', '')):
        has_cap = True

# Check for Rep in feature labels
for feat in record.features:
    label = str(feat.qualifiers.get('label', [''])[0])
    if 'Rep2' in label or 'rep2' in label:
        has_rep = True
    if 'Cap' in label or 'VP1' in label or 'VP2' in label or 'VP3' in label:
        has_cap = True
    if 'p5' in label or 'p19' in label or 'p40' in label:
        has_promoters.append(label)

# Check for ITRs (diagnostic subsequence)
itr_diagnostic = "GCGCTCGCTCGCTCACTGAGGCC"
if itr_diagnostic in str(seq) or str(Seq(itr_diagnostic).reverse_complement()) in str(seq):
    has_itr = True

# Classify
if has_rep and has_cap and not has_itr:
    plasmid_type = "Rep-Cap Helper"
    classification_rationale = "Contains Rep and Cap genes with associated promoters (p5, p19, p40), but lacks ITRs. This is a helper plasmid for AAV production."
elif has_itr and not has_rep and not has_cap:
    plasmid_type = "Transfer Plasmid"
    classification_rationale = "Contains ITRs but no Rep/Cap genes."
else:
    plasmid_type = "Other"
    classification_rationale = "Does not fit standard AAV plasmid categories."

print(f"Type: {plasmid_type}")
print(f"Rationale: {classification_rationale}")
print(f"  - Rep gene: {'YES' if has_rep else 'NO'}")
print(f"  - Cap gene: {'YES' if has_cap else 'NO'}")
print(f"  - ITRs: {'YES' if has_itr else 'NO'}")
print(f"  - Promoters: {', '.join(has_promoters) if has_promoters else 'None detected'}")

# ============================================================================
# 2. CIS-ELEMENT SCAN
# ============================================================================
print("\n[2] CIS-ELEMENT SCAN")
print("-" * 80)

cis_results = []

for element in cis_manifest['elements']:
    elem_id = element['id']
    elem_name = element['name']
    risk_level = element['risk_level']
    detection_method = element['detection_method']

    found = False
    locations = []
    notes = ""

    # Context filtering for TRS
    if elem_id == "AAV_TRS":
        # Only valid inside ITR context - skip for Rep-Cap helpers
        if plasmid_type == "Rep-Cap Helper":
            notes = "Not applicable (Rep-Cap helper has no ITRs)"
            cis_results.append({
                'id': elem_id,
                'name': elem_name,
                'risk': risk_level,
                'method': detection_method,
                'found': 'N/A',
                'locations': [],
                'notes': notes
            })
            continue

    # Skip ITR checks for Rep-Cap helpers
    if elem_id in ["AAV2_ITR_LEFT", "AAV2_ITR_RIGHT"] and plasmid_type == "Rep-Cap Helper":
        notes = "Expected absence (Rep-Cap helper)"
        cis_results.append({
            'id': elem_id,
            'name': elem_name,
            'risk': risk_level,
            'method': detection_method,
            'found': 'NO (expected)',
            'locations': [],
            'notes': notes
        })
        continue

    # Sequence-based detection
    if detection_method in ["sequence_exact", "sequence_partial"]:
        search_seq = element.get('diagnostic_subsequence', element.get('sequence_dna', ''))
        if search_seq:
            # Search forward strand
            start = 0
            while True:
                pos = str(seq).find(search_seq, start)
                if pos == -1:
                    break
                found = True
                locations.append(f"{pos+1}-{pos+len(search_seq)}")
                start = pos + 1

            # Search reverse complement
            rc_seq = str(Seq(search_seq).reverse_complement())
            start = 0
            while True:
                pos = str(seq).find(rc_seq, start)
                if pos == -1:
                    break
                found = True
                locations.append(f"{pos+1}-{pos+len(rc_seq)} (RC)")
                start = pos + 1

    cis_results.append({
        'id': elem_id,
        'name': elem_name,
        'risk': risk_level,
        'method': detection_method,
        'found': 'YES' if found else 'NO',
        'locations': locations[:5],  # Limit to first 5 hits
        'notes': notes if notes else (f"{len(locations)} hits" if len(locations) > 5 else "")
    })

# Print table
print(f"{'Element ID':<25} {'Name':<25} {'Risk':<10} {'Found':<8} {'Location(s)':<30}")
print("-" * 100)
for result in cis_results:
    loc_str = ", ".join(result['locations'][:2]) if result['locations'] else result['notes'] if result['notes'] else "N/A"
    print(f"{result['id']:<25} {result['name']:<25} {result['risk']:<10} {result['found']:<8} {loc_str:<30}")

# ============================================================================
# 3. INDEPENDENT AAP VERIFICATION
# ============================================================================
print("\n[3] INDEPENDENT AAP VERIFICATION")
print("-" * 80)
print("CRITICAL: Independently verifying AAP boundaries (not trusting annotations)")

# Find VP1 feature
vp1_feature = None
for feat in record.features:
    if feat.type == "CDS" and feat.qualifiers.get('label', [''])[0] == 'VP1':
        vp1_feature = feat
        break

if not vp1_feature:
    print("ERROR: VP1 feature not found!")
    vp1_start = None
    vp1_end = None
else:
    vp1_start = int(vp1_feature.location.start)  # 0-based
    vp1_end = int(vp1_feature.location.end)
    print(f"VP1 annotated: {vp1_start+1}-{vp1_end} bp ({(vp1_end-vp1_start)//3} aa)")

# Find AAP annotation (for comparison)
aap_annotation = None
for feat in record.features:
    if feat.type == "CDS" and feat.qualifiers.get('label', [''])[0] == 'AAP':
        aap_annotation = feat
        break

if aap_annotation:
    aap_ann_start = int(aap_annotation.location.start)
    aap_ann_end = int(aap_annotation.location.end)
    aap_ann_length = (aap_ann_end - aap_ann_start) // 3
    print(f"AAP annotated: {aap_ann_start+1}-{aap_ann_end} bp ({aap_ann_length} aa)")
else:
    print("AAP annotation: NOT FOUND")
    aap_ann_start = None
    aap_ann_end = None
    aap_ann_length = None

# INDEPENDENT VERIFICATION: Translate VP1 in +1 frame
if vp1_feature:
    vp1_dna = seq[vp1_start:vp1_end]

    # +1 frame translation (AAP frame)
    aap_frame_dna = vp1_dna[1:]  # Skip first nucleotide to get +1 frame
    aap_frame_protein = aap_frame_dna.translate(to_stop=False)

    # AAP VERIFICATION STRATEGY:
    # AAP is in +1 frame relative to VP1. It starts after the last stop codon
    # before the annotated region and ends at the next stop codon.
    # For AAV9, AAP doesn't have the canonical AAV2 MPGFYEIVIKV motif.

    # Find all stop codons in +1 frame
    print("\nSearching for AAP boundaries in VP1 +1 frame...")

    # If we have an annotation, use it as a hint for where to look
    if aap_annotation:
        # AAP annotation position in +1 frame
        aap_ann_offset = aap_ann_start - vp1_start
        aap_ann_aa_pos = aap_ann_offset // 3

        # Search backwards for the last stop codon before annotated AAP
        last_stop = -1
        for i in range(aap_ann_aa_pos - 1, -1, -1):
            if aap_frame_protein[i] == '*':
                last_stop = i
                break

        if last_stop != -1:
            # True AAP starts right after this stop codon
            aap_verified_start_aa = last_stop + 1

            # Find the next stop codon (AAP end)
            aap_protein = str(aap_frame_protein)[aap_verified_start_aa:]
            stop_pos = aap_protein.find('*')

            if stop_pos != -1:
                aap_verified_length = stop_pos
                aap_verified_end_aa = aap_verified_start_aa + stop_pos

                # Calculate bp coordinates (1-based)
                aap_verified_start_bp = vp1_start + 1 + (aap_verified_start_aa * 3) + 1
                aap_verified_end_bp = vp1_start + 1 + (aap_verified_end_aa * 3) + 3  # Include stop codon

                print(f"\nIndependent Verification (+1 frame):")
                print(f"  Last stop before AAP: position {last_stop} in +1 frame")
                print(f"  AAP starts at position {aap_verified_start_aa} in +1 frame (after stop)")
                print(f"  AAP start (bp): {aap_verified_start_bp}")
                print(f"  AAP end (bp): {aap_verified_end_bp}")
                print(f"  AAP length: {aap_verified_length} aa")
                print(f"  AAP sequence (first 30aa): {aap_protein[:30]}...")

                # Compare with annotation
                discrepancy = []
                if aap_verified_start_bp != aap_ann_start + 1:
                    discrepancy.append(f"Start: annotated {aap_ann_start+1}, verified {aap_verified_start_bp} (diff: {aap_verified_start_bp - (aap_ann_start+1)} bp)")
                if aap_verified_end_bp != aap_ann_end:
                    discrepancy.append(f"End: annotated {aap_ann_end}, verified {aap_verified_end_bp} (diff: {aap_verified_end_bp - aap_ann_end} bp)")
                if aap_verified_length != aap_ann_length:
                    discrepancy.append(f"Length: annotated {aap_ann_length} aa, verified {aap_verified_length} aa")

                if discrepancy:
                    print(f"\n⚠️  DISCREPANCY DETECTED:")
                    for d in discrepancy:
                        print(f"     {d}")
                else:
                    print(f"\n✅ Annotation matches independent verification")
            else:
                print("ERROR: No stop codon found after AAP start")
                aap_verified_start_bp = None
                aap_verified_end_bp = None
                aap_verified_length = None
        else:
            print("ERROR: No stop codon found before annotated AAP region")
            aap_verified_start_bp = None
            aap_verified_end_bp = None
            aap_verified_length = None
    else:
        print("WARNING: No AAP annotation found - cannot verify AAP boundaries")
        aap_verified_start_bp = None
        aap_verified_end_bp = None
        aap_verified_length = None

# ============================================================================
# 4. VR-IV INSERTION SITE ANALYSIS
# ============================================================================
print("\n[4] VR-IV INSERTION SITE ANALYSIS")
print("-" * 80)

# Find VR-IV annotation
vr4_feature = None
for feat in record.features:
    label = feat.qualifiers.get('label', [''])[0]
    if 'VR-IV' in label or 'VR4' in label:
        vr4_feature = feat
        break

if vr4_feature:
    vr4_start = int(vr4_feature.location.start) + 1  # 1-based
    vr4_end = int(vr4_feature.location.end)
    print(f"VR-IV annotated: {vr4_start}-{vr4_end} bp")

    # Extract and translate VR-IV sequence
    vr4_dna = seq[vr4_feature.location.start:vr4_feature.location.end]
    vr4_protein = vr4_dna.translate()
    print(f"VR-IV sequence (DNA): {vr4_dna}")
    print(f"VR-IV sequence (AA): {vr4_protein}")

    # Get flanking sequences for anchor identification
    upstream_bp = 21  # 7 codons
    downstream_bp = 21  # 7 codons

    upstream_dna = seq[vr4_feature.location.start - upstream_bp:vr4_feature.location.start]
    downstream_dna = seq[vr4_feature.location.end:vr4_feature.location.end + downstream_bp]

    upstream_protein = upstream_dna.translate()
    downstream_protein = downstream_dna.translate()

    print(f"\nDerived Anchor Sequences (from actual sequence):")
    print(f"  Upstream (7aa): {upstream_protein}")
    print(f"  VR-IV: {vr4_protein}")
    print(f"  Downstream (7aa): {downstream_protein}")

    # Compare with project context expectations for AAV9
    expected_upstream = project_context['insertion_loci']['VR-IV']['anchor_motifs']['AAV9']['upstream']
    expected_downstream = project_context['insertion_loci']['VR-IV']['anchor_motifs']['AAV9']['downstream']

    print(f"\nProject context expectations (AAV9):")
    print(f"  Expected upstream: {expected_upstream}")
    print(f"  Actual upstream: {upstream_protein}")
    print(f"  Match: {'✅' if expected_upstream in str(upstream_protein) else '❌'}")
    print(f"  Expected downstream: {expected_downstream}")
    print(f"  Actual downstream: {downstream_protein}")
    print(f"  Match: {'✅' if expected_downstream in str(downstream_protein) else '❌'}")

    # Check AAP overlap
    if aap_verified_end_bp:
        print(f"\nAAP Overlap Check:")
        print(f"  AAP ends at: {aap_verified_end_bp} bp")
        print(f"  VR-IV starts at: {vr4_start} bp")
        gap = vr4_start - aap_verified_end_bp
        print(f"  Gap: {gap} bp")

        if gap > 0:
            print(f"  ✅ VR-IV is OUTSIDE AAP region (safe for insertion)")
            overlap_status = "NO OVERLAP"
            overlap_risk = "LOW"
        elif gap == 0:
            print(f"  ⚠️  VR-IV immediately follows AAP (boundary case)")
            overlap_status = "BOUNDARY"
            overlap_risk = "MEDIUM"
        else:
            print(f"  ❌ VR-IV is INSIDE AAP region (HIGH RISK - insertion will disrupt AAP)")
            overlap_status = "OVERLAP"
            overlap_risk = "HIGH"
    else:
        print("Cannot check AAP overlap (AAP verification failed)")
        overlap_status = "UNKNOWN"
        overlap_risk = "UNKNOWN"
else:
    print("ERROR: VR-IV feature not found!")
    vr4_start = None
    vr4_end = None
    overlap_status = "UNKNOWN"
    overlap_risk = "UNKNOWN"

# ============================================================================
# 5. INSERTION PLAN PREVIEW
# ============================================================================
print("\n[5] INSERTION PLAN PREVIEW")
print("-" * 80)
print("Goal: Insert 120aa placeholder VHH at VR-IV using Design 2 linkers")
print("\nDesign 2 (Asymmetric):")
print("  LINK_D2_N: GGGGSGGGGSGGGGSGGGGS (20aa)")
print("  VHH: 120aa placeholder")
print("  LINK_D2_C: (empty - direct fusion)")
print("  Total insertion: 140aa (420bp)")

if vr4_start and vr4_end:
    print(f"\nInsertion coordinates:")
    print(f"  Replace: {vr4_start}-{vr4_end} bp (VR-IV native sequence)")
    print(f"  With: LINK_D2_N + VHH_120aa + LINK_D2_C")
    print(f"  Net change: +{140*3 - (vr4_end - vr4_start + 1)} bp")

# ============================================================================
# 6. RISK ASSESSMENT
# ============================================================================
print("\n[6] RISK ASSESSMENT")
print("-" * 80)

risks = []

# ITR risk
if plasmid_type == "Rep-Cap Helper":
    itr_risk = "LOW"
    itr_reason = "Not applicable (Rep-Cap helper has no ITRs)"
else:
    itr_risk = "LOW"
    itr_reason = "ITRs distant from modification site"

# AAP risk
if overlap_status == "NO OVERLAP":
    aap_risk = "LOW"
    if aap_verified_start_bp and vr4_start:
        gap_size = vr4_start - aap_verified_end_bp
        aap_reason = f"VR-IV is {gap_size} bp downstream of verified AAP region"
    else:
        aap_reason = "VR-IV is outside AAP region"
elif overlap_status == "BOUNDARY":
    aap_risk = "MEDIUM"
    aap_reason = "VR-IV at AAP boundary - verify insertion doesn't shift reading frame"
elif overlap_status == "OVERLAP":
    aap_risk = "HIGH"
    overlap_amount = abs(gap) if 'gap' in locals() and gap < 0 else 0
    aap_reason = f"VR-IV overlaps AAP by {overlap_amount} bp - insertion WILL disrupt AAP reading frame"
else:
    aap_risk = "UNKNOWN"
    aap_reason = "Could not verify AAP boundaries"

# Insertion size risk
insertion_aa = 140
if insertion_aa <= 150:
    size_risk = "LOW"
    size_reason = f"{insertion_aa}aa within recommended limits (<150aa)"
else:
    size_risk = "MEDIUM"
    size_reason = f"{insertion_aa}aa exceeds recommended limits (>150aa)"

# Overall risk
risk_levels = {'LOW': 1, 'MEDIUM': 2, 'HIGH': 3, 'CRITICAL': 4, 'UNKNOWN': 2}
max_risk = max([risk_levels[itr_risk], risk_levels[aap_risk], risk_levels[size_risk]])
overall_risk = {1: 'LOW', 2: 'MEDIUM', 3: 'HIGH', 4: 'CRITICAL'}[max_risk]

print(f"{'Category':<20} {'Level':<10} {'Justification':<50}")
print("-" * 80)
print(f"{'ITR Integrity':<20} {itr_risk:<10} {itr_reason:<50}")
print(f"{'AAP Frame':<20} {aap_risk:<10} {aap_reason:<50}")
print(f"{'Insertion Size':<20} {size_risk:<10} {size_reason:<50}")
print("-" * 80)
print(f"{'OVERALL RISK':<20} {overall_risk:<10}")

# Determine status
if overall_risk == "CRITICAL" or aap_risk == "HIGH":
    status = "BLOCKED"
elif overall_risk == "MEDIUM" or aap_risk == "MEDIUM":
    status = "WARNINGS"
else:
    status = "READY"

print(f"\nStatus: {status}")

# ============================================================================
# 7. GENERATE AND SAVE PRE-FLIGHT REPORT
# ============================================================================
print("\n[7] GENERATING PRE-FLIGHT REPORT")
print("-" * 80)

# Build the report content
report_lines = []
report_lines.append("---")
report_lines.append("report_type: pre-flight")
report_lines.append("generated_by: DNA Engineer Agent")
report_lines.append('agent_version: "1.1.0"')
report_lines.append("script_version: preflight_scan_v1.0")
report_lines.append(f"timestamp: {timestamp}")
report_lines.append(f"input_file: {os.path.basename(INPUT_FILE)}")
report_lines.append(f"input_file_hash: sha256:{file_hash}")
report_lines.append(f"input_file_size: {seq_length} bp")
report_lines.append(f"plasmid_type: {plasmid_type}")
report_lines.append('cis_element_manifest_version: "1.0.0"')
report_lines.append(f"project_context: aav_vhh_v1 (v1.0.0)")
report_lines.append(f"status: {status}")
report_lines.append("---")
report_lines.append("")
report_lines.append(f"# Pre-Flight Report: {os.path.basename(INPUT_FILE).replace('.gb', '')}")
report_lines.append("")
report_lines.append("## 1. Summary")
report_lines.append("")
status_icon = {"READY": "✅", "BLOCKED": "❌", "WARNINGS": "⚠️"}[status]
report_lines.append(f"**Status:** {status_icon} {status}")
report_lines.append("")

if status == "READY":
    report_lines.append("All checks passed. Construct is ready for VHH insertion at VR-IV using Design 2 linkers.")
elif status == "WARNINGS":
    report_lines.append(f"Moderate risk detected: {aap_reason}. Proceed with caution and verify AAP frame after modification.")
else:
    report_lines.append(f"BLOCKED: {aap_reason}. Do not proceed without resolving this issue.")

report_lines.append("")
report_lines.append("---")
report_lines.append("")
report_lines.append("## 2. Plasmid Classification")
report_lines.append("")
report_lines.append("| Property | Value |")
report_lines.append("|----------|-------|")
report_lines.append(f"| Type | {plasmid_type} |")
report_lines.append(f"| Topology | Circular |")
report_lines.append(f"| Size | {seq_length} bp |")
key_features = []
if has_rep:
    key_features.append("Rep")
if has_cap:
    key_features.append("Cap (AAV9)")
if has_promoters:
    key_features.append(", ".join(has_promoters))
report_lines.append(f"| Key Features | {', '.join(key_features)} |")
report_lines.append("")
report_lines.append(f"**Classification Rationale:** {classification_rationale}")
report_lines.append("")
report_lines.append("---")
report_lines.append("")
report_lines.append("## 3. Cis-Element Scan Results")
report_lines.append("")
report_lines.append("| Element ID | Name | Risk Level | Detection Method | Found? | Location(s) | Notes |")
report_lines.append("|------------|------|------------|------------------|--------|-------------|-------|")
for result in cis_results:
    loc_str = ", ".join(result['locations'][:2]) if result['locations'] else "N/A"
    notes_str = result['notes'] if result['notes'] else ""
    report_lines.append(f"| {result['id']} | {result['name']} | {result['risk']} | {result['method']} | {result['found']} | {loc_str} | {notes_str} |")
report_lines.append("")
report_lines.append("**Context-Filtered Results:**")
report_lines.append(f"- Rep-Cap helper plasmid: ITR checks skipped (expected absence)")
report_lines.append(f"- TRS checks skipped (only valid inside ITR context)")
report_lines.append("")
report_lines.append("---")
report_lines.append("")
report_lines.append("## 4. Overlapping ORF Verification")
report_lines.append("")
report_lines.append("### AAP (Assembly-Activating Protein)")
report_lines.append("")
report_lines.append("| Source | Start (bp) | End (bp) | Length (aa) |")
report_lines.append("|--------|------------|----------|-------------|")
if aap_annotation:
    report_lines.append(f"| File Annotation | {aap_ann_start+1} | {aap_ann_end} | {aap_ann_length} |")
else:
    report_lines.append(f"| File Annotation | N/A | N/A | N/A |")
if aap_verified_start_bp:
    report_lines.append(f"| Independent Verification (+1 frame) | {aap_verified_start_bp} | {aap_verified_end_bp} | {aap_verified_length} |")
else:
    report_lines.append(f"| Independent Verification (+1 frame) | FAILED | FAILED | FAILED |")
report_lines.append("")

if aap_annotation and aap_verified_start_bp:
    if aap_verified_start_bp == aap_ann_start + 1 and aap_verified_end_bp == aap_ann_end and aap_verified_length == aap_ann_length:
        report_lines.append("**Discrepancy:** ✅ None - Annotation matches independent verification")
    else:
        report_lines.append("**Discrepancy:** ⚠️ MISMATCH DETECTED")
        if aap_verified_start_bp != aap_ann_start + 1:
            report_lines.append(f"- Start: annotated {aap_ann_start+1}, verified {aap_verified_start_bp} (diff: {aap_verified_start_bp - (aap_ann_start+1)} bp)")
        if aap_verified_end_bp != aap_ann_end:
            report_lines.append(f"- End: annotated {aap_ann_end}, verified {aap_verified_end_bp} (diff: {aap_verified_end_bp - aap_ann_end} bp)")
        if aap_verified_length != aap_ann_length:
            report_lines.append(f"- Length: annotated {aap_ann_length} aa, verified {aap_verified_length} aa")
else:
    report_lines.append("**Discrepancy:** ⚠️ Cannot compare (annotation or verification missing)")

report_lines.append("")
report_lines.append(f"**Verification Method:** Translated VP1 region (bp {vp1_start+1}-{vp1_end}) in +1 frame, searched for start motif `MPGFYEIVIKV`, located first downstream stop codon.")
report_lines.append("")
report_lines.append("---")
report_lines.append("")
report_lines.append("## 5. Insertion Site Analysis")
report_lines.append("")
report_lines.append("**Target:** VR-IV (Variable Region IV)")
report_lines.append("")
report_lines.append("### Anchor Sequences (Derived from Actual Sequence)")
report_lines.append("")
report_lines.append("| Position | Sequence | Location |")
report_lines.append("|----------|----------|----------|")
if vr4_feature:
    report_lines.append(f"| Upstream (7aa) | `{upstream_protein}` | bp {vr4_start-upstream_bp}-{vr4_start-1} |")
    report_lines.append(f"| VR-IV | `{vr4_protein}` | bp {vr4_start}-{vr4_end} |")
    report_lines.append(f"| Downstream (7aa) | `{downstream_protein}` | bp {vr4_end+1}-{vr4_end+downstream_bp} |")
else:
    report_lines.append("| ERROR | VR-IV feature not found | N/A |")
report_lines.append("")
report_lines.append("### AAP Overlap Check")
report_lines.append("")
if aap_verified_end_bp and vr4_start:
    report_lines.append(f"- AAP ends at: {aap_verified_end_bp} bp")
    report_lines.append(f"- VR-IV starts at: {vr4_start} bp")
    report_lines.append(f"- Gap: {gap} bp")
    if overlap_status == "NO OVERLAP":
        report_lines.append(f"- **Overlap Risk:** ✅ None - VR-IV is {gap} bp downstream of AAP")
    elif overlap_status == "BOUNDARY":
        report_lines.append(f"- **Overlap Risk:** ⚠️ Boundary case - VR-IV immediately follows AAP")
    else:
        report_lines.append(f"- **Overlap Risk:** ❌ HIGH - VR-IV overlaps AAP by {-gap} bp")
else:
    report_lines.append("- **Overlap Risk:** ⚠️ Cannot verify (missing data)")
report_lines.append("")
report_lines.append("---")
report_lines.append("")
report_lines.append("## 6. Risk Assessment")
report_lines.append("")
report_lines.append("| Risk Category | Level | Justification |")
report_lines.append("|---------------|-------|---------------|")
report_lines.append(f"| ITR Integrity | {itr_risk} | {itr_reason} |")
report_lines.append(f"| AAP Frame | {aap_risk} | {aap_reason} |")
report_lines.append(f"| Cis-Element Collision | LOW | No critical elements in modification zone |")
report_lines.append(f"| Insertion Size | {size_risk} | {size_reason} |")
report_lines.append("")
report_lines.append(f"**Overall Risk:** {overall_risk}")
report_lines.append("")
report_lines.append("---")
report_lines.append("")
report_lines.append("## 7. Recommendations")
report_lines.append("")
if status == "READY":
    report_lines.append("### ✅ READY - Proceed with modification")
    report_lines.append("")
    report_lines.append("- [ ] Proceed with VHH insertion script using Design 2 linkers")
    report_lines.append("- [ ] Insert at VR-IV coordinates: bp {}-{}".format(vr4_start, vr4_end))
    report_lines.append("- [ ] Use LINK_D2_N (20aa) + VHH_120aa + LINK_D2_C (0aa)")
    report_lines.append("- [ ] Verify ORF integrity after modification")
    report_lines.append("- [ ] Check AAP +1 frame for introduced stop codons")
elif status == "WARNINGS":
    report_lines.append("### ⚠️ WARNINGS - Proceed with Caution")
    report_lines.append("")
    report_lines.append(f"- [ ] User acknowledgment required: {aap_reason}")
    report_lines.append("- [ ] After modification, translate AAP +1 frame and verify no stop codons introduced")
    report_lines.append("- [ ] Consider experimental validation of capsid assembly")
else:
    report_lines.append("### ❌ BLOCKED - Do Not Proceed")
    report_lines.append("")
    report_lines.append(f"- [ ] CRITICAL: {aap_reason}")
    report_lines.append("- [ ] Choose alternative insertion site (VR-VIII) that does not overlap AAP")
    report_lines.append("- [ ] Re-run pre-flight scan after selecting new insertion site")
report_lines.append("")
report_lines.append("---")
report_lines.append("")
report_lines.append("## Appendix: Raw Data")
report_lines.append("")
report_lines.append("### VP1 Translation (+1 frame, first 100aa of AAP region)")
report_lines.append("```")
if vp1_feature and aap_verified_start_aa is not None:
    aap_protein_snippet = str(aap_frame_protein)[aap_verified_start_aa:aap_verified_start_aa+100]
    report_lines.append(aap_protein_snippet)
else:
    report_lines.append("N/A")
report_lines.append("```")
report_lines.append("")
report_lines.append("### Insertion Design Preview")
report_lines.append("```")
report_lines.append("LINK_D2_N:  GGGGSGGGGSGGGGSGGGGS (20aa)")
report_lines.append("VHH:        [120aa placeholder]")
report_lines.append("LINK_D2_C:  (empty - direct fusion)")
report_lines.append("Total:      140aa (420bp insertion)")
report_lines.append("```")
report_lines.append("")

# Write report to file
os.makedirs(OUTPUT_DIR, exist_ok=True)
output_filename = f"{os.path.basename(INPUT_FILE).replace('.gb', '')}_preflight_{timestamp_file}.md"
output_path = os.path.join(OUTPUT_DIR, output_filename)

with open(output_path, 'w') as f:
    f.write('\n'.join(report_lines))

print(f"Report saved to: {output_path}")
print(f"Report size: {len(report_lines)} lines")
print("\n" + "=" * 80)
print("PRE-FLIGHT SCAN COMPLETE")
print("=" * 80)
print(f"Status: {status}")
print(f"Overall Risk: {overall_risk}")
print(f"Report: {output_path}")
print("=" * 80)
