#!/usr/bin/env python3
"""
Find Restriction Sites Unique in BOTH RepCap and Transfer Plasmid Contexts

This script searches for alternative restriction sites in the 6 regions of VP1
that are unique in BOTH:
1. The RepCap source plasmid (AAV9-RepCap-NOCAP-WITH-RESTRICTION-SITES.gb)
2. The assembled transfer plasmids (ssAAV and scAAV)

Critical for ensuring sites work for cloning in both contexts.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Restriction import Analysis, RestrictionBatch, AllEnzymes
import sys
from pathlib import Path

# Curated list of 61 reliable restriction enzymes
CURATED_ENZYMES = [
    # Standard 6-bp cutters
    'AfeI', 'AflII', 'AgeI', 'ApaI', 'ApaLI', 'AvrII', 'BamHI', 'BglII',
    'BmtI', 'BbvCI', 'BsaI', 'BsmBI', 'BspEI', 'BsrGI', 'BstBI', 'BstZ17I',
    'ClaI', 'EagI', 'EcoRI', 'EcoRV', 'Esp3I', 'HindIII', 'HpaI', 'KpnI',
    'MfeI', 'MluI', 'MscI', 'NcoI', 'NdeI', 'NheI', 'NruI', 'NsiI',
    'PciI', 'PmlI', 'PstI', 'PvuI', 'PvuII', 'SacI', 'SacII', 'SalI',
    'ScaI', 'SmaI', 'SnaBI', 'SpeI', 'SphI', 'StuI', 'XbaI', 'XhoI', 'XmaI',
    # 8-bp cutters
    'AscI', 'AsiSI', 'FseI', 'NotI', 'PacI', 'PmeI', 'SbfI', 'SrfI', 'SwaI',
    # Minimal ambiguity
    'AflIII', 'BsaWI', 'HincII'
]


def load_plasmids():
    """Load all relevant plasmids"""
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


def get_vp1_region_in_context(plasmid, vp1_label='VP1'):
    """
    Find VP1 CDS coordinates in a plasmid.
    Returns (start, end, sequence) or None if not found.
    """
    for feature in plasmid.features:
        if feature.qualifiers.get('label', [''])[0] == vp1_label:
            vp1_start = feature.location.start
            vp1_end = feature.location.end
            vp1_seq = feature.extract(plasmid.seq)
            return (vp1_start, vp1_end, vp1_seq)
    return None


def check_enzyme_uniqueness(enzyme, seq):
    """
    Check if enzyme cuts exactly once in the sequence.
    Returns (is_unique, count, positions)
    """
    try:
        # Import enzyme dynamically
        from Bio.Restriction import Restriction
        enzyme_obj = getattr(Restriction, enzyme)

        analysis = Analysis(RestrictionBatch([enzyme_obj]), seq)
        sites = analysis.full()

        if enzyme_obj in sites:
            positions = sites[enzyme_obj]
            count = len(positions)
            return (count == 1, count, positions)
        else:
            return (False, 0, [])

    except Exception as e:
        return (False, 0, [])


def find_sites_in_region(region_seq, region_abs_start, enzymes):
    """
    Find all restriction sites in a region.
    Returns dict of {enzyme: [(abs_position, rel_position, site_seq), ...]}
    """
    sites_found = {}

    for enzyme_name in enzymes:
        try:
            from Bio.Restriction import Restriction
            enzyme = getattr(Restriction, enzyme_name)

            # Search for sites
            site_pattern = str(enzyme.site)
            region_str = str(region_seq)

            # Find all occurrences
            pos = 0
            while True:
                pos = region_str.find(site_pattern, pos)
                if pos == -1:
                    break

                abs_pos = region_abs_start + pos
                rel_pos = pos

                if enzyme_name not in sites_found:
                    sites_found[enzyme_name] = []

                sites_found[enzyme_name].append((abs_pos, rel_pos, site_pattern))
                pos += 1

        except:
            continue

    return sites_found


def analyze_region(region_num, region_name, vp1_rel_start, vp1_rel_end, plasmids):
    """
    Analyze a specific region for alternative restriction sites.

    Returns list of viable alternatives: [(enzyme, details), ...]
    """
    print(f"\n{'='*80}")
    print(f"REGION {region_num}: {region_name}")
    print(f"VP1 Relative Position: {vp1_rel_start}-{vp1_rel_end}")
    print(f"{'='*80}")

    # Get VP1 coordinates in each context
    repcap_vp1 = get_vp1_region_in_context(plasmids['repcap'])
    ssaav_vp1 = get_vp1_region_in_context(plasmids['ssaav'])
    scaav_vp1 = get_vp1_region_in_context(plasmids['scaav'])

    if not all([repcap_vp1, ssaav_vp1, scaav_vp1]):
        print("❌ Could not find VP1 in one or more plasmids!")
        return []

    # Extract region sequence (same in all contexts - it's from VP1)
    region_seq = repcap_vp1[2][vp1_rel_start:vp1_rel_end]
    print(f"\nRegion length: {len(region_seq)} bp")
    print(f"Region sequence: {region_seq[:50]}...{region_seq[-50:]}")

    # Find all potential sites in this region
    region_sites = find_sites_in_region(region_seq, 0, CURATED_ENZYMES)

    if not region_sites:
        print("\n❌ No restriction sites found in this region!")
        return []

    print(f"\n✓ Found {len(region_sites)} enzymes with sites in this region")

    # Check each site for uniqueness in all three contexts
    viable_alternatives = []

    for enzyme in sorted(region_sites.keys()):
        # Check uniqueness in RepCap
        repcap_unique, repcap_count, repcap_pos = check_enzyme_uniqueness(
            enzyme, plasmids['repcap'].seq
        )

        # Check uniqueness in ssAAV
        ssaav_unique, ssaav_count, ssaav_pos = check_enzyme_uniqueness(
            enzyme, plasmids['ssaav'].seq
        )

        # Check uniqueness in scAAV
        scaav_unique, scaav_count, scaav_pos = check_enzyme_uniqueness(
            enzyme, plasmids['scaav'].seq
        )

        # Only viable if unique in ALL contexts
        all_unique = repcap_unique and ssaav_unique and scaav_unique

        if all_unique:
            site_info = region_sites[enzyme][0]  # First occurrence
            viable_alternatives.append({
                'enzyme': enzyme,
                'rel_pos': vp1_rel_start + site_info[1],
                'site_seq': site_info[2],
                'repcap_pos': repcap_pos[0] if repcap_pos else None,
                'ssaav_pos': ssaav_pos[0] if ssaav_pos else None,
                'scaav_pos': scaav_pos[0] if scaav_pos else None
            })

    # Report results
    print(f"\n{'─'*80}")
    print(f"RESULTS FOR REGION {region_num}")
    print(f"{'─'*80}")

    if viable_alternatives:
        print(f"\n✅ Found {len(viable_alternatives)} VIABLE alternative(s):\n")

        for alt in viable_alternatives:
            print(f"  {alt['enzyme']:12} ({alt['site_seq']:10})")
            print(f"    VP1 relative position: {alt['rel_pos']}")
            print(f"    RepCap position: {alt['repcap_pos'] + 1}")
            print(f"    ssAAV position:  {alt['ssaav_pos'] + 1}")
            print(f"    scAAV position:  {alt['scaav_pos'] + 1}")
            print(f"    Status: ✅ UNIQUE in all contexts\n")
    else:
        print("\n❌ NO VIABLE ALTERNATIVES FOUND")
        print("    No enzymes are unique in all three contexts for this region.")
        print("    ⚠️  FLAG FOR MANUAL REVIEW\n")

        # Show why common enzymes failed
        print("  Common enzymes in region but NOT unique:")
        for enzyme in sorted(region_sites.keys())[:5]:  # Show first 5
            repcap_unique, repcap_count, _ = check_enzyme_uniqueness(
                enzyme, plasmids['repcap'].seq
            )
            ssaav_unique, ssaav_count, _ = check_enzyme_uniqueness(
                enzyme, plasmids['ssaav'].seq
            )
            scaav_unique, scaav_count, _ = check_enzyme_uniqueness(
                enzyme, plasmids['scaav'].seq
            )

            print(f"    {enzyme:12} RepCap:{repcap_count:2} ssAAV:{ssaav_count:2} scAAV:{scaav_count:2}")

    return viable_alternatives


def main():
    """Main analysis workflow"""
    print("="*80)
    print("DUAL-CONTEXT RESTRICTION SITE ANALYSIS")
    print("="*80)
    print("\nSearching for sites unique in BOTH RepCap and Transfer plasmids\n")

    # Load plasmids
    plasmids = load_plasmids()
    if not plasmids:
        print("\n❌ Failed to load required plasmids")
        sys.exit(1)

    # Define the 6 regions (relative to VP1 start, 0-based)
    # These match the original six-regions analysis
    regions = [
        (1, "Rep68-stop to VP2-start", 0, 280),  # Around position 141 (SmaI)
        (2, "VP2-AAP Intergenic", 400, 550),     # Around position 464 (BbvCI)
        (3, "AAP-stop to VR4", 1150, 1300),      # Around position 1219 (AgeI)
        (4, "VR4 to VR5", 1350, 1480),           # Around position 1416 (BsrGI)
        (5, "VR5 to VR8", 1500, 1650),           # Around position 1573 (BmtI) - Already OK
        (6, "Post-VR8", 1750, 1850),             # Around position 1799 (BstZ17I) - Already OK
    ]

    # Analyze each region
    all_results = {}

    for region_num, region_name, start, end in regions:
        results = analyze_region(region_num, region_name, start, end, plasmids)
        all_results[region_num] = results

    # Summary
    print("\n" + "="*80)
    print("FINAL SUMMARY")
    print("="*80)

    for region_num in sorted(all_results.keys()):
        results = all_results[region_num]
        status = "✅ OK" if results else "❌ FLAG"
        count = len(results)
        print(f"\nRegion {region_num}: {count} viable alternative(s) {status}")

        if results:
            for alt in results[:3]:  # Show up to 3 alternatives
                print(f"  - {alt['enzyme']} ({alt['site_seq']})")

    # Regions needing attention
    flagged = [r for r, results in all_results.items() if not results]

    if flagged:
        print(f"\n⚠️  FLAGGED REGIONS NEEDING MANUAL REVIEW: {flagged}")
        print("    These regions have no enzymes unique in all three contexts.")
        print("    Consider:")
        print("      1. Using sites from other regions")
        print("      2. Changing the promoter to avoid conflicts")
        print("      3. Accepting non-unique sites with careful digest planning")
    else:
        print("\n✅ ALL REGIONS HAVE VIABLE ALTERNATIVES!")

    print("\n" + "="*80)


if __name__ == '__main__':
    main()
