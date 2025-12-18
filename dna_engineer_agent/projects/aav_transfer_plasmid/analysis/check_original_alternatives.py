#!/usr/bin/env python3
"""
Check Original Alternative Enzyme Candidates in Both Contexts

This script tests the ORIGINAL alternative enzyme candidates that were identified
during the six-regions analysis for the RepCap plasmid. It checks if these
alternatives are unique in BOTH:
1. RepCap source plasmid (AAV9-RepCap-NOCAP-WITH-RESTRICTION-SITES.gb)
2. Transfer plasmids (ssAAV and scAAV v02)

Based on alternatives documented in:
- projects/nocap_restriction_sites/verification/verify_region3_alternatives.py
- projects/nocap_restriction_sites/verification/verify_region4_alternatives.py
- reports/AAV9_RepCap_SixRegions_Analysis_v2.3_REVISED.md
"""

from Bio import SeqIO
from Bio.Restriction import Analysis, RestrictionBatch
from Bio.Restriction import Restriction
from pathlib import Path
import sys

# Original alternative candidates identified during RepCap analysis
# Format: (region_num, region_name, enzyme, position_in_repcap, notes)
ORIGINAL_ALTERNATIVES = {
    1: [
        # Region 1: Rep68-stop to VP2-start (2415-2775 in RepCap)
        # From FINAL summary: Recommended SmaI @ 2505
        # From revised report: AvrII @ 2459
        ('SmaI', 2505, 'Current recommended'),
        ('AvrII', 2459, 'Alternative in v2.3 report'),
        ('EagI', 2369, 'Original v2.0, too close to VP1 start'),
    ],

    2: [
        # Region 2: VP2-AAP Intergenic (2779-2890 in RepCap)
        # From FINAL summary: BbvCI @ 2828
        ('BbvCI', 2828, 'Current recommended - Type IIS'),
    ],

    3: [
        # Region 3: AAP-stop to VR4 (3485-3717 in RepCap)
        # From FINAL summary: AgeI @ 3583
        # From verification script: PshAI @ 3533, BsmBI @ 3566
        # From v2.3 report: EcoNI @ 3487
        ('AgeI', 3583, 'Current recommended in FINAL'),
        ('EcoNI', 3487, 'Alternative in v2.3 report'),
        ('PshAI', 3533, 'Alternative - Type II'),
        ('BsmBI', 3566, 'Alternative - Type IIS'),
        ('EcoRV', None, 'Excluded - used by Genscript'),
    ],

    4: [
        # Region 4: VR4 to VR5 (3745-3825 in RepCap)
        # From FINAL summary: BsrGI @ 3780 (WORKS in v02!)
        # From verification: NaeI @ 3760, RsrII @ 3762
        # From v2.3 report: FseI @ 3759
        ('BsrGI', 3780, 'Current recommended - WORKS in v02 ✓'),
        ('FseI', 3759, 'Alternative in v2.3 report - fussy enzyme'),
        ('NaeI', 3760, 'Alternative - common enzyme'),
        ('NgoMIV', 3760, 'Alternative - same site as NaeI'),
        ('RsrII', 3762, 'Alternative - IUPAC ambiguity'),
    ],

    5: [
        # Region 5: VR5 to VR8 (3880-4104 in RepCap)
        # From FINAL summary: BmtI @ 3937 (WORKS in v02!)
        ('BmtI', 3937, 'Current recommended - WORKS in v02 ✓'),
        ('NheI', 3937, 'Same recognition sequence as BmtI'),
    ],

    6: [
        # Region 6: Post-VR8 (4144-4575 in RepCap, narrowed to 4144-4343)
        # From FINAL summary: BstZ17I @ 4163 (WORKS in v02!)
        # From v2.3 report: BaeI @ 4318
        ('BstZ17I', 4163, 'Current recommended - WORKS in v02 ✓'),
        ('BaeI', 4318, 'Alternative in v2.3 report - double cutter'),
        ('AfeI', 4442, 'Alternative - beyond 200bp limit'),
    ],
}


def load_plasmids():
    """Load all three plasmids for comparison"""
    plasmids = {}

    base_path = Path('test_data')

    files = {
        'repcap': base_path / 'AAV9-RepCap-NOCAP-WITH-RESTRICTION-SITES.gb',
        'ssaav': base_path / 'pGS-ssAAV-EF1A-VP1-rBG_v02.gb',
        'scaav': base_path / 'pGS-scAAV-EF1A-VP1-rBG_v02.gb'
    }

    for name, filepath in files.items():
        if filepath.exists():
            plasmids[name] = SeqIO.read(str(filepath), 'genbank')
            print(f"✓ Loaded {name}: {len(plasmids[name])} bp")
        else:
            print(f"❌ Missing: {filepath}")
            return None

    return plasmids


def check_enzyme_uniqueness(enzyme_name, seq):
    """
    Check if enzyme cuts exactly once in the sequence.
    Returns (is_unique, count, positions)
    """
    try:
        enzyme_obj = getattr(Restriction, enzyme_name)
        analysis = Analysis(RestrictionBatch([enzyme_obj]), seq)
        sites = analysis.full()

        if enzyme_obj in sites:
            positions = sites[enzyme_obj]
            count = len(positions)
            return (count == 1, count, positions)
        else:
            return (False, 0, [])

    except AttributeError:
        # Enzyme not found in Biopython
        return (False, -1, [])
    except Exception as e:
        return (False, -2, [])


def main():
    print("="*80)
    print("CHECKING ORIGINAL ALTERNATIVE ENZYME CANDIDATES")
    print("="*80)
    print("\nTesting alternatives from original six-regions RepCap analysis")
    print("in BOTH RepCap source and assembled transfer plasmid contexts\n")

    # Load plasmids
    plasmids = load_plasmids()
    if not plasmids:
        print("\n❌ Failed to load required plasmids")
        sys.exit(1)

    print()

    # Summary data structure
    viable_by_region = {}

    # Test each region
    for region_num in sorted(ORIGINAL_ALTERNATIVES.keys()):
        alternatives = ORIGINAL_ALTERNATIVES[region_num]

        print("\n" + "="*80)
        print(f"REGION {region_num}")
        print("="*80)

        viable_alternatives = []

        for enzyme_name, position, notes in alternatives:
            if position is None:
                print(f"\n  {enzyme_name}: {notes} - SKIPPED")
                continue

            print(f"\n  {enzyme_name} @ position {position} ({notes})")

            # Check uniqueness in all three contexts
            repcap_unique, repcap_count, repcap_pos = check_enzyme_uniqueness(
                enzyme_name, plasmids['repcap'].seq
            )

            ssaav_unique, ssaav_count, ssaav_pos = check_enzyme_uniqueness(
                enzyme_name, plasmids['ssaav'].seq
            )

            scaav_unique, scaav_count, scaav_pos = check_enzyme_uniqueness(
                enzyme_name, plasmids['scaav'].seq
            )

            # Handle enzyme not found
            if repcap_count == -1:
                print(f"    ⚠️  Enzyme not found in Biopython Restriction library")
                continue

            # Report counts
            print(f"    RepCap:  {repcap_count} site(s)", end="")
            if repcap_unique:
                print(" ✓")
            else:
                print(f" ❌ (positions: {[p+1 for p in repcap_pos]})")

            print(f"    ssAAV:   {ssaav_count} site(s)", end="")
            if ssaav_unique:
                print(" ✓")
            else:
                if ssaav_count == 0:
                    print(" (none)")
                else:
                    print(f" ❌ (positions: {[p+1 for p in ssaav_pos]})")

            print(f"    scAAV:   {scaav_count} site(s)", end="")
            if scaav_unique:
                print(" ✓")
            else:
                if scaav_count == 0:
                    print(" (none)")
                else:
                    print(f" ❌ (positions: {[p+1 for p in scaav_pos]})")

            # Check if viable in ALL contexts
            all_unique = repcap_unique and ssaav_unique and scaav_unique

            if all_unique:
                print(f"    >>> ✅ VIABLE IN ALL CONTEXTS")
                viable_alternatives.append({
                    'enzyme': enzyme_name,
                    'repcap_pos': position,
                    'ssaav_pos': ssaav_pos[0] + 1 if ssaav_pos else None,
                    'scaav_pos': scaav_pos[0] + 1 if scaav_pos else None,
                    'notes': notes
                })
            else:
                reasons = []
                if not repcap_unique:
                    reasons.append("RepCap")
                if not ssaav_unique:
                    reasons.append("ssAAV")
                if not scaav_unique:
                    reasons.append("scAAV")
                print(f"    >>> ❌ NOT VIABLE - Non-unique in: {', '.join(reasons)}")

        viable_by_region[region_num] = viable_alternatives

    # Final summary
    print("\n\n" + "="*80)
    print("FINAL SUMMARY: VIABLE ALTERNATIVES BY REGION")
    print("="*80)

    total_viable = 0

    for region_num in sorted(viable_by_region.keys()):
        viable = viable_by_region[region_num]
        total_viable += len(viable)

        if viable:
            print(f"\n✅ Region {region_num}: {len(viable)} viable alternative(s)")
            for alt in viable:
                print(f"   - {alt['enzyme']:12} @ RepCap:{alt['repcap_pos']:4} | "
                      f"ssAAV:{alt['ssaav_pos']:4} | scAAV:{alt['scaav_pos']:4}")
                print(f"     Notes: {alt['notes']}")
        else:
            print(f"\n❌ Region {region_num}: NO viable alternatives")
            print(f"   All original candidates fail uniqueness in transfer plasmid context")

    print(f"\n{'='*80}")
    print(f"TOTAL VIABLE SITES ACROSS ALL REGIONS: {total_viable}")
    print(f"{'='*80}")

    # Recommendations
    print("\n\n" + "="*80)
    print("RECOMMENDATIONS")
    print("="*80)

    if total_viable >= 3:
        print(f"\n✅ SUFFICIENT: {total_viable} viable sites found across all regions")
        print("   This is adequate for basic modular cloning strategies.")
    else:
        print(f"\n⚠️  LIMITED: Only {total_viable} viable sites found")
        print("   Consider:")
        print("   1. Using a different promoter (without problematic sites)")
        print("   2. Redesigning problematic regions in RepCap source")
        print("   3. Using partial digest strategies for non-unique sites")

    # Identify problematic regions
    no_viable = [r for r, v in viable_by_region.items() if not v]
    if no_viable:
        print(f"\n⚠️  Regions {no_viable} have NO viable alternatives")
        print("   These regions require manual review and redesign")

    print("\n" + "="*80)


if __name__ == '__main__':
    main()
