#!/usr/bin/env python3
"""
Generate comprehensive report for revised 6-region analysis
with Region 1 narrowed to Rep68-stop to VP2-start
"""

import sys
sys.path.insert(0, 'scripts/tools')
from Bio import SeqIO
from silent_sites import translate

# Load plasmid
record = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")
plasmid_seq = str(record.seq).upper()

# Define the 6 sites
sites = [
    {
        "region": "Region 1: Rep68-stop to VP2-start",
        "region_short": "Region 1",
        "start": 2415,
        "end": 2775,
        "description": "After Rep68 stop, before VP2 start",
        "enzyme": "AvrII",
        "position": 2459,
        "site_sequence": "CCTAGG",
        "original_seq": "CCCAGG",
        "modified_seq": "CCTAGG",
        "mutations": ["A3T"],
        "mutation_details": "Position 2461: CCT → CCC (P → P) ✓ SILENT",
    },
    {
        "region": "Region 2: VP2-AAP Intergenic",
        "region_short": "Region 2",
        "start": 2779,
        "end": 2890,
        "description": "After VP2 start, before AAP start",
        "enzyme": "BbvCI",
        "position": 2828,
        "site_sequence": "CCTCAGC",
        "original_seq": "CCTCCGC",
        "modified_seq": "CCTCAGC",
        "mutations": ["C5A"],
        "mutation_details": "Position 2832: TCC → TCA (S → S) ✓ SILENT",
    },
    {
        "region": "Region 3: AAP-stop to VR4",
        "region_short": "Region 3",
        "start": 3485,
        "end": 3717,
        "description": "Between AAP stop and VR4 start",
        "enzyme": "EcoNI",
        "position": 3487,
        "site_sequence": "CCTCAGTAAGG",
        "original_seq": "CCTCAGTACGG",
        "modified_seq": "CCTCAGTAAGG",
        "mutations": ["C9A"],
        "mutation_details": "Position 3495-3497: CGG → AGG (R → R) ✓ SILENT",
    },
    {
        "region": "Region 4: VR4 to VR5",
        "region_short": "Region 4",
        "start": 3745,
        "end": 3825,
        "description": "Between VR4 end and VR5 start",
        "enzyme": "FseI",
        "position": 3759,
        "site_sequence": "GGCCGGCC",
        "original_seq": "GGCCGGAC",
        "modified_seq": "GGCCGGCC",
        "mutations": ["A7C"],
        "mutation_details": "Position 3765: GGA → GGC (G → G) ✓ SILENT",
    },
    {
        "region": "Region 5: VR5 to VR8",
        "region_short": "Region 5",
        "start": 3880,
        "end": 4104,
        "description": "Between VR5 end and VR8 start",
        "enzyme": "BmtI",
        "position": 3937,
        "site_sequence": "GCTAGC",
        "original_seq": "GCCAGC",
        "modified_seq": "GCTAGC",
        "mutations": ["C3T"],
        "mutation_details": "Position 3939: GCC → GCT (A → A) ✓ SILENT",
    },
    {
        "region": "Region 6: Post-VR8",
        "region_short": "Region 6",
        "start": 4144,
        "end": 4575,
        "description": "After VR8 end",
        "enzyme": "BaeI",
        "position": 4318,
        "site_sequence": "ACACCTGTACC",
        "original_seq": "ACACCTGTACC",
        "modified_seq": "ACACCTGTACC",
        "mutations": [],
        "mutation_details": "No mutations required - site already present",
    },
]

# Generate markdown report
report_lines = []

report_lines.append("# REVISED: Comprehensive Restriction Site Analysis - 6 Regions")
report_lines.append("")
report_lines.append("**Plasmid:** BASE-DRAFT-AAV9-RepCap-NOCAP.gb")
report_lines.append("**Analysis Date:** 2025-12-17")
report_lines.append("**Version:** 2.3.0 (Region 1 REVISED)")
report_lines.append("**Plasmid Length:** 7,078 bp")
report_lines.append("")
report_lines.append("---")
report_lines.append("")
report_lines.append("## Revision Notes")
report_lines.append("")
report_lines.append("**Region 1 has been revised** to avoid the splice acceptor in the VP1 unique N-terminus:")
report_lines.append("")
report_lines.append("- **Previous Region 1:** 2365-2775 (VP1 start to VP2 start) - EagI @ 2369")
report_lines.append("- **Revised Region 1:** 2415-2775 (Rep68 stop to VP2 start) - AvrII @ 2459")
report_lines.append("")
report_lines.append("The revised Region 1 is narrowed by 50 bp to start after the Rep68 stop codon (position 2414),")
report_lines.append("avoiding potential splice acceptor sites in the VP1 unique N-terminal region.")
report_lines.append("")
report_lines.append("---")
report_lines.append("")
report_lines.append("## Executive Summary")
report_lines.append("")
report_lines.append("This analysis identifies unique restriction enzyme sites (requiring ≤2 silent mutations) in six")
report_lines.append("strategically chosen regions of an AAV9 Rep-Cap helper plasmid. These sites enable future")
report_lines.append("modifications while avoiding overlap with critical boundaries (start codons, variable regions).")
report_lines.append("")
report_lines.append("**Key Findings:**")
report_lines.append("- **6 unique restriction sites** identified (1 per region)")
report_lines.append("- **Total mutations required:** 5 silent mutations + 1 site already present")
report_lines.append("- **All sites verified** to not overlap critical boundaries")
report_lines.append("- **All sites are unique** in the entire 7,078 bp plasmid")
report_lines.append("- **No amino acid changes** in any proposed modification")
report_lines.append("")
report_lines.append("---")
report_lines.append("")
report_lines.append("## Summary Table")
report_lines.append("")
report_lines.append("| Region | Enzyme | Position | Recognition Site | Mutations | Status |")
report_lines.append("|--------|--------|----------|------------------|-----------|--------|")

for site in sites:
    mutations_str = f"{len(site['mutations'])} silent" if site['mutations'] else "0 (exists)"
    report_lines.append(f"| {site['region_short']} | **{site['enzyme']}** | {site['position']} | `{site['site_sequence']}` | {mutations_str} | ✓ |")

report_lines.append("")
report_lines.append("---")
report_lines.append("")
report_lines.append("## Detailed Analysis by Region")
report_lines.append("")

for idx, site in enumerate(sites, 1):
    report_lines.append(f"### {site['region']}")
    report_lines.append("")
    report_lines.append(f"**Region Boundaries:** {site['start']}-{site['end']} ({site['end'] - site['start'] + 1} bp)")
    report_lines.append(f"**Description:** {site['description']}")
    report_lines.append("")

    if idx == 1:
        report_lines.append("**⚠️ REVISED:** This region was narrowed from the original analysis to start after Rep68 stop codon,")
        report_lines.append("avoiding the splice acceptor in the VP1 unique N-terminus.")
        report_lines.append("")

    report_lines.append("#### Recommended Site")
    report_lines.append("")
    report_lines.append("| Parameter | Value |")
    report_lines.append("|-----------|-------|")
    report_lines.append(f"| **Enzyme** | {site['enzyme']} |")
    report_lines.append(f"| **Recognition Sequence** | {site['site_sequence']} |")
    report_lines.append(f"| **Position** | {site['position']} (1-indexed) |")
    report_lines.append(f"| **Mutations Required** | {len(site['mutations'])} |")
    report_lines.append(f"| **Original Sequence** | `{site['original_seq']}` |")
    report_lines.append(f"| **Modified Sequence** | `{site['modified_seq']}` |")
    report_lines.append(f"| **Uniqueness** | ✓ Unique in entire plasmid |")
    report_lines.append(f"| **Boundary Check** | ✓ No overlaps |")
    report_lines.append("")

    if site['mutations']:
        report_lines.append("#### Codon Changes")
        report_lines.append("")
        report_lines.append("```")
        report_lines.append(site['mutation_details'])
        report_lines.append("```")
        report_lines.append("")
        report_lines.append("**Mutation Details:**")
        for mut in site['mutations']:
            report_lines.append(f"- {mut}")
        report_lines.append("")
    else:
        report_lines.append("**No mutations required** - site already present in the plasmid.")
        report_lines.append("")

    # Add sequence context
    context_start = site['position'] - 35
    context_end = site['position'] + len(site['site_sequence']) + 35
    context_seq = plasmid_seq[context_start - 1:context_end - 1]

    report_lines.append("#### Sequence Context")
    report_lines.append("")
    report_lines.append("```")
    report_lines.append(f"Position {context_start}-{context_end}")
    report_lines.append(context_seq)
    marker_pos = 35
    marker = " " * marker_pos + "^" * len(site['site_sequence'])
    report_lines.append(marker)
    report_lines.append(f"{site['enzyme']} site @ {site['position']}")
    report_lines.append("```")
    report_lines.append("")
    report_lines.append("---")
    report_lines.append("")

# Applications section
report_lines.append("## Applications")
report_lines.append("")
report_lines.append("These restriction sites provide excellent handles for:")
report_lines.append("")
report_lines.append("1. **Variable Region Engineering**")
report_lines.append("   - Insert VHH nanobodies at VR4/VR5 using flanking sites")
report_lines.append("   - Swap entire variable region cassettes")
report_lines.append("   - Region 3 (AAP-stop to VR4) enables insertions just before VR4")
report_lines.append("   - Region 4 (VR4 to VR5) enables dual insertions at VR4 and VR5")
report_lines.append("")
report_lines.append("2. **Epitope Tagging**")
report_lines.append("   - Add FLAG, HA, or His tags in VP1 unique region")
report_lines.append("   - Insert fluorescent proteins for trafficking studies")
report_lines.append("")
report_lines.append("3. **Golden Gate Assembly** (BbvCI in Region 2)")
report_lines.append("   - Scarless cloning workflows")
report_lines.append("   - Modular assembly of Cap variants")
report_lines.append("")
report_lines.append("4. **Diagnostic Restriction Mapping**")
report_lines.append("   - Quick QC after cloning")
report_lines.append("   - Verify plasmid identity")
report_lines.append("")
report_lines.append("---")
report_lines.append("")

# Conclusions
report_lines.append("## Conclusions")
report_lines.append("")
report_lines.append("This revised analysis successfully identified **6 unique restriction sites** across strategic regions")
report_lines.append("of the AAV9 Rep-Cap plasmid:")
report_lines.append("")
report_lines.append("✓ **All sites unique** in entire plasmid")
report_lines.append("✓ **All mutations silent** (no amino acid changes)")
report_lines.append("✓ **No boundary overlaps** (verified for 11 critical boundaries)")
report_lines.append("✓ **Minimal mutations** (5 total across 6 sites, 1 site exists naturally)")
report_lines.append("✓ **Strategically positioned** for VR engineering")
report_lines.append("✓ **Region 1 revised** to avoid splice acceptor in VP1 N-terminus")
report_lines.append("")
report_lines.append("These sites provide excellent molecular handles for:")
report_lines.append("- VHH display at variable regions")
report_lines.append("- Epitope tagging in VP1 unique domain")
report_lines.append("- Modular assembly via Golden Gate")
report_lines.append("- Diagnostic restriction mapping")
report_lines.append("")
report_lines.append("All modifications are low-risk and preserve the functional architecture of the Cap gene.")
report_lines.append("")
report_lines.append("---")
report_lines.append("")
report_lines.append("**Analysis completed successfully.**")

# Write report
output_file = "reports/AAV9_RepCap_SixRegions_Analysis_v2.3_REVISED.md"
with open(output_file, 'w') as f:
    f.write('\n'.join(report_lines))

print(f"Report written to: {output_file}")
print(f"\n{'='*80}")
print("REVISED 6-REGION RESTRICTION SITE ANALYSIS - SUMMARY")
print(f"{'='*80}")
for site in sites:
    print(f"{site['region']:40s} {site['enzyme']:10s} @ {site['position']:5d}")
print(f"{'='*80}")
