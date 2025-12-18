#!/usr/bin/env python3
"""
Correct Interpretation: Which Alternatives Are Viable Design Choices?

For alternatives that show "0 in RepCap":
- This means we haven't engineered them YET
- We CAN still engineer them (add silent mutation)
- The question is: would they be unique in the TRANSFER PLASMID?

Key insight:
- If alternative shows 0 in ssAAV/scAAV → sequence NOT in backbone → VIABLE choice
- If alternative shows 1+ in ssAAV/scAAV → sequence IN backbone → NOT VIABLE

We only care about transfer plasmid backbone conflicts!
"""

from Bio import SeqIO
from Bio.Restriction import Analysis, RestrictionBatch, Restriction
from pathlib import Path

ORIGINAL_ALTERNATIVES = {
    1: [
        ('SmaI', 2505, 'Current choice in RepCap'),
        ('AvrII', 2459, 'Alternative candidate'),
        ('EagI', 2369, 'Alternative candidate'),
    ],
    2: [
        ('BbvCI', 2828, 'Current choice in RepCap'),
    ],
    3: [
        ('AgeI', 3583, 'Current choice in RepCap'),
        ('EcoNI', 3487, 'Alternative candidate'),
        ('PshAI', 3533, 'Alternative candidate'),
        ('BsmBI', 3566, 'Alternative candidate'),
    ],
    4: [
        ('BsrGI', 3780, 'Current choice in RepCap'),
        ('FseI', 3759, 'Alternative candidate'),
        ('NaeI', 3760, 'Alternative candidate'),
        ('NgoMIV', 3760, 'Alternative candidate'),
    ],
    5: [
        ('BmtI', 3937, 'Current choice in RepCap'),
        ('NheI', 3937, 'Same site as BmtI'),
    ],
    6: [
        ('BstZ17I', 4163, 'Current choice in RepCap'),
        ('BaeI', 4318, 'Alternative - naturally occurs 2x'),
        ('AfeI', 4442, 'Alternative candidate'),
    ],
}


def check_enzyme_uniqueness(enzyme_name, seq):
    """Check if enzyme cuts exactly once in the sequence."""
    try:
        enzyme_obj = getattr(Restriction, enzyme_name)
        analysis = Analysis(RestrictionBatch([enzyme_obj]), seq)
        sites = analysis.full()

        if enzyme_obj in sites:
            positions = sites[enzyme_obj]
            count = len(positions)
            return (count == 1, count, positions)
        else:
            return (True, 0, [])  # 0 sites = unique if we add it

    except AttributeError:
        return (False, -1, [])
    except Exception as e:
        return (False, -2, [])


def main():
    # Load transfer plasmids (we check backbone conflicts)
    base_path = Path('test_data')

    ssaav = SeqIO.read(str(base_path / 'pGS-ssAAV-EF1A-VP1-rBG_v02.gb'), 'genbank')
    scaav = SeqIO.read(str(base_path / 'pGS-scAAV-EF1A-VP1-rBG_v02.gb'), 'genbank')

    print("="*80)
    print("CORRECTED ANALYSIS: VIABLE DESIGN CHOICES")
    print("="*80)
    print("\nWhich original alternatives would work IF we engineer them?\n")
    print("Key: If site count is 0 in transfer plasmid, we CAN engineer it.")
    print("     If site count is 1+, sequence exists in backbone → conflict.\n")

    viable_by_region = {}

    for region_num in sorted(ORIGINAL_ALTERNATIVES.keys()):
        alternatives = ORIGINAL_ALTERNATIVES[region_num]

        print("\n" + "="*80)
        print(f"REGION {region_num}")
        print("="*80)

        viable = []

        for enzyme_name, position, notes in alternatives:
            # Check if recognition sequence appears in transfer plasmid backbones
            ssaav_unique, ssaav_count, ssaav_pos = check_enzyme_uniqueness(
                enzyme_name, ssaav.seq
            )
            scaav_unique, scaav_count, scaav_pos = check_enzyme_uniqueness(
                enzyme_name, scaav.seq
            )

            if ssaav_count == -1:
                print(f"\n  {enzyme_name}: Enzyme not in Biopython - SKIP")
                continue

            print(f"\n  {enzyme_name} @ position {position}")
            print(f"    {notes}")
            print(f"    ssAAV backbone: {ssaav_count} site(s)", end="")
            if ssaav_count == 0:
                print(" ✓ (can engineer)")
            else:
                print(f" ❌ (conflicts at {[p+1 for p in ssaav_pos]})")

            print(f"    scAAV backbone: {scaav_count} site(s)", end="")
            if scaav_count == 0:
                print(" ✓ (can engineer)")
            else:
                print(f" ❌ (conflicts at {[p+1 for p in scaav_pos]})")

            # Viable if 0 occurrences in BOTH transfer plasmid backbones
            # (or 1 occurrence which would be the engineered site in VP1)
            if ssaav_count <= 1 and scaav_count <= 1:
                # Need to check if the 1 occurrence is in VP1 region or backbone
                # For now, assume 0 = definitely viable, 1 = need to check
                if ssaav_count == 0 and scaav_count == 0:
                    print(f"    >>> ✅ VIABLE - Can engineer without conflicts")
                    viable.append(enzyme_name)
                else:
                    # 1 site exists - need to determine if it's in VP1 or backbone
                    print(f"    >>> ⚠️  MAYBE - Need to check if site is in backbone or VP1")
            else:
                print(f"    >>> ❌ NOT VIABLE - Multiple sites in backbone")

        viable_by_region[region_num] = viable

    # Summary
    print("\n\n" + "="*80)
    print("SUMMARY: VIABLE ALTERNATIVES TO ENGINEER")
    print("="*80)

    for region_num in sorted(viable_by_region.keys()):
        viable = viable_by_region[region_num]
        if viable:
            print(f"\n✅ Region {region_num}: {len(viable)} viable option(s)")
            for enzyme in viable:
                print(f"   - {enzyme}")
        else:
            print(f"\n❌ Region {region_num}: No alternatives with 0 backbone conflicts")

    print("\n" + "="*80)


if __name__ == '__main__':
    main()
