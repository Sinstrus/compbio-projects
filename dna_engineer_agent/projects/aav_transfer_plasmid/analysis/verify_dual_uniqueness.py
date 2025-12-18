#!/usr/bin/env python3
"""
Verify that "transfer plasmid unique" sites would ALSO be unique in RepCap

Logic:
- If enzyme recognition sequence appears 0 times in current RepCap
- AND 0 times in transfer plasmid backbones
- THEN engineering it would make it unique in BOTH contexts (1 site each)

This verifies our top candidates work everywhere.
"""

from Bio import SeqIO
from Bio.Restriction import Analysis, RestrictionBatch, Restriction
from pathlib import Path

# Top candidates from transfer plasmid analysis (1-mutation, perfectly unique)
CANDIDATES = {
    'Region 1': ['AvrII', 'BsaI', 'BsmBI', 'BspEI', 'Esp3I', 'HindIII', 'HpaI', 'KpnI',
                 'MluI', 'NcoI', 'SalI', 'SnaBI', 'XbaI', 'BstBI'],
    'Region 2': ['BsaI', 'BsmBI', 'BspEI', 'Esp3I', 'KpnI', 'SalI', 'XbaI'],
    'Region 3': ['BsaI', 'BsmBI', 'ClaI', 'EcoRV', 'Esp3I', 'HindIII', 'HpaI',
                 'KpnI', 'MluI', 'NcoI', 'NsiI', 'SalI', 'SnaBI']
}


def count_enzyme_sites(enzyme_name, seq):
    """Count enzyme sites in sequence"""
    try:
        enzyme = getattr(Restriction, enzyme_name)
        analysis = Analysis(RestrictionBatch([enzyme]), seq)
        sites = analysis.full()
        if enzyme in sites:
            return len(sites[enzyme])
        return 0
    except:
        return -1


def main():
    # Load plasmids
    base_path = Path('test_data')
    repcap = SeqIO.read(str(base_path / 'AAV9-RepCap-NOCAP-WITH-RESTRICTION-SITES.gb'), 'genbank')
    ssaav = SeqIO.read(str(base_path / 'pGS-ssAAV-EF1A-VP1-rBG_v02.gb'), 'genbank')
    scaav = SeqIO.read(str(base_path / 'pGS-scAAV-EF1A-VP1-rBG_v02.gb'), 'genbank')

    print("="*80)
    print("VERIFICATION: Do Transfer-Unique Sites Also Work in RepCap?")
    print("="*80)
    print("\nChecking if our top candidates would be unique in RepCap too...\n")

    for region_name, enzymes in CANDIDATES.items():
        print(f"\n{region_name}:")
        print("-" * 60)
        print(f"{'Enzyme':<12} {'RepCap':<10} {'ssAAV':<10} {'scAAV':<10} {'Both Unique?'}")
        print("-" * 60)

        for enzyme in enzymes:
            repcap_count = count_enzyme_sites(enzyme, repcap.seq)
            ssaav_count = count_enzyme_sites(enzyme, ssaav.seq)
            scaav_count = count_enzyme_sites(enzyme, scaav.seq)

            if repcap_count == -1:
                continue

            # Site is unique in BOTH contexts if:
            # - 0 in RepCap (we can engineer it → 1 site)
            # - 0 in both transfer plasmids (no backbone conflicts)
            both_unique = (repcap_count == 0 and ssaav_count == 0 and scaav_count == 0)

            status = "✅ YES" if both_unique else "❌ NO"
            if repcap_count > 0:
                status += f" (RepCap conflict)"
            elif ssaav_count > 0 or scaav_count > 0:
                status += f" (Transfer conflict)"

            print(f"{enzyme:<12} {repcap_count:<10} {ssaav_count:<10} {scaav_count:<10} {status}")

    # Summary
    print("\n\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    all_both_unique = []

    for region_name, enzymes in CANDIDATES.items():
        both_unique_in_region = []

        for enzyme in enzymes:
            repcap_count = count_enzyme_sites(enzyme, repcap.seq)
            ssaav_count = count_enzyme_sites(enzyme, ssaav.seq)
            scaav_count = count_enzyme_sites(enzyme, scaav.seq)

            if repcap_count == 0 and ssaav_count == 0 and scaav_count == 0:
                both_unique_in_region.append(enzyme)
                all_both_unique.append(enzyme)

        print(f"\n{region_name}:")
        print(f"  Unique in BOTH contexts: {len(both_unique_in_region)}/{len(enzymes)}")
        if both_unique_in_region:
            print(f"  Enzymes: {', '.join(both_unique_in_region)}")

    print(f"\n{'='*80}")
    print(f"TOTAL ENZYMES UNIQUE IN BOTH CONTEXTS: {len(set(all_both_unique))}")
    print(f"{'='*80}")

    # List unique enzymes across all regions
    unique_all = set(all_both_unique)
    print(f"\nEnzymes available in at least one region (unique in both contexts):")
    for enzyme in sorted(unique_all):
        regions_available = [r for r, enzs in CANDIDATES.items()
                           if enzyme in enzs and
                           count_enzyme_sites(enzyme, repcap.seq) == 0 and
                           count_enzyme_sites(enzyme, ssaav.seq) == 0 and
                           count_enzyme_sites(enzyme, scaav.seq) == 0]
        print(f"  {enzyme:<12} → {', '.join(regions_available)}")

    print("\n" + "="*80)


if __name__ == '__main__':
    main()
