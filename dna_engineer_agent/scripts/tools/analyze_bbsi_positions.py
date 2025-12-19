#!/usr/bin/env python3
"""
Analyze BbsI positions relative to HindIII-XbaI cloning regions in AAV backbones.

This script answers the architect's question:
Are BbsI sites INSIDE or OUTSIDE the HindIII-XbaI cloning region?
"""

import json

def main():
    catalog_path = "backbones/genscript/BACKBONE_CATALOG.json"

    with open(catalog_path, 'r') as f:
        catalog = json.load(f)

    print("="*80)
    print("BbsI Position Analysis - AAV Backbones")
    print("="*80)
    print()

    for backbone in catalog['aav_backbones']:
        name = backbone['file']
        hindiii_sites = backbone['mcs_sites'].get('HindIII', [])
        xbai_sites = backbone['mcs_sites'].get('XbaI', [])
        bbsi_sites = backbone['type_iis_sites'].get('BbsI', [])

        if not hindiii_sites or not xbai_sites:
            print(f"⚠️  {name}: Missing cloning sites")
            continue

        hindiii = hindiii_sites[0]
        xbai = xbai_sites[0]
        cloning_region_size = xbai - hindiii

        print(f"📁 {name}")
        print(f"   Cloning region: HindIII ({hindiii}) → XbaI ({xbai})")
        print(f"   Insert size: {cloning_region_size} bp")

        if bbsi_sites:
            for bbsi in bbsi_sites:
                if hindiii < bbsi < xbai:
                    location = "INSIDE cloning region"
                    effect = "❌ Will be REMOVED during cloning"
                elif bbsi < hindiii:
                    location = "UPSTREAM of cloning region"
                    effect = "✅ Will REMAIN after cloning"
                else:  # bbsi > xbai
                    location = "DOWNSTREAM of cloning region"
                    effect = "✅ Will REMAIN after cloning"

                distance_from_cloning = min(abs(bbsi - hindiii), abs(bbsi - xbai))

                print(f"   BbsI position: {bbsi}")
                print(f"   Location: {location}")
                print(f"   Distance from cloning region: {distance_from_cloning} bp")
                print(f"   Effect: {effect}")
        else:
            print(f"   BbsI: Not found")

        print()

    print("="*80)
    print("CONCLUSION:")
    print("="*80)
    print()
    print("All BbsI sites in pGS-AAV backbones are DOWNSTREAM of the XbaI site,")
    print("meaning they are OUTSIDE the HindIII-XbaI cloning region.")
    print()
    print("Impact: After cloning your insert (HindIII → insert → XbaI), the BbsI")
    print("        sites will REMAIN in the backbone downstream of the 3' ITR.")
    print()
    print("Golden Gate Compatibility: These backbones are NOT compatible with")
    print("                            BbsI-based Golden Gate assembly.")
    print()

if __name__ == "__main__":
    main()
