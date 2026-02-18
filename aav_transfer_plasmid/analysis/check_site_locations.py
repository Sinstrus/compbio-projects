#!/usr/bin/env python3
"""
Determine if "1 site" occurrences are in VP1 or backbone

For sites showing 1 occurrence in transfer plasmid:
- If position is within VP1 coordinates → that's our engineered site → VIABLE
- If position is outside VP1 coordinates → backbone conflict → NOT VIABLE
"""

from Bio import SeqIO
from Bio.Restriction import Analysis, RestrictionBatch, Restriction
from pathlib import Path

# Sites that showed "1 occurrence" and need checking
SITES_TO_CHECK = {
    'EagI': 2369,
    'BsrGI': 3780,
    'FseI': 3759,
    'NaeI': 3760,
    'NgoMIV': 3760,
    'BmtI': 3937,
    'NheI': 3937,
    'BstZ17I': 4163,
    'AfeI': 4442,
}


def get_vp1_coords(plasmid):
    """Find VP1 CDS coordinates in plasmid"""
    for feature in plasmid.features:
        if feature.qualifiers.get('label', [''])[0] == 'VP1':
            return (feature.location.start, feature.location.end)
    return None


def check_enzyme_position(enzyme_name, seq):
    """Get position of enzyme site"""
    try:
        enzyme_obj = getattr(Restriction, enzyme_name)
        analysis = Analysis(RestrictionBatch([enzyme_obj]), seq)
        sites = analysis.full()

        if enzyme_obj in sites:
            positions = sites[enzyme_obj]
            return positions
        else:
            return []

    except:
        return []


def main():
    # Load transfer plasmids
    base_path = Path('test_data')
    ssaav = SeqIO.read(str(base_path / 'pGS-ssAAV-EF1A-VP1-rBG_v02.gb'), 'genbank')
    scaav = SeqIO.read(str(base_path / 'pGS-scAAV-EF1A-VP1-rBG_v02.gb'), 'genbank')

    # Get VP1 coordinates
    ssaav_vp1 = get_vp1_coords(ssaav)
    scaav_vp1 = get_vp1_coords(scaav)

    print("="*80)
    print("CHECKING SITE LOCATIONS: VP1 vs BACKBONE")
    print("="*80)
    print(f"\nssAAV VP1 coordinates: {ssaav_vp1[0]}-{ssaav_vp1[1]} ({ssaav_vp1[1] - ssaav_vp1[0]} bp)")
    print(f"scAAV VP1 coordinates: {scaav_vp1[0]}-{scaav_vp1[1]} ({scaav_vp1[1] - scaav_vp1[0]} bp)")

    print("\n" + "="*80)

    for enzyme_name, repcap_pos in sorted(SITES_TO_CHECK.items()):
        print(f"\n{enzyme_name} (would be at RepCap position {repcap_pos}):")

        # Check ssAAV
        ssaav_positions = check_enzyme_position(enzyme_name, ssaav.seq)
        if ssaav_positions:
            for pos in ssaav_positions:
                pos_1based = pos + 1
                in_vp1 = ssaav_vp1[0] <= pos < ssaav_vp1[1]
                location = "VP1" if in_vp1 else "BACKBONE"
                status = "✅" if in_vp1 else "❌"
                print(f"  ssAAV: position {pos_1based:4} → {location:8} {status}")
        else:
            print(f"  ssAAV: 0 sites → Can engineer ✅")

        # Check scAAV
        scaav_positions = check_enzyme_position(enzyme_name, scaav.seq)
        if scaav_positions:
            for pos in scaav_positions:
                pos_1based = pos + 1
                in_vp1 = scaav_vp1[0] <= pos < scaav_vp1[1]
                location = "VP1" if in_vp1 else "BACKBONE"
                status = "✅" if in_vp1 else "❌"
                print(f"  scAAV: position {pos_1based:4} → {location:8} {status}")
        else:
            print(f"  scAAV: 0 sites → Can engineer ✅")

        # Determine viability
        ssaav_ok = not ssaav_positions or all(ssaav_vp1[0] <= pos < ssaav_vp1[1] for pos in ssaav_positions)
        scaav_ok = not scaav_positions or all(scaav_vp1[0] <= pos < scaav_vp1[1] for pos in scaav_positions)

        if ssaav_ok and scaav_ok:
            print(f"  >>> ✅ VIABLE - All sites in VP1 or no conflicts")
        else:
            conflicts = []
            if not ssaav_ok:
                conflicts.append("ssAAV backbone")
            if not scaav_ok:
                conflicts.append("scAAV backbone")
            print(f"  >>> ❌ NOT VIABLE - Conflicts in {', '.join(conflicts)}")

    print("\n" + "="*80)


if __name__ == '__main__':
    main()
