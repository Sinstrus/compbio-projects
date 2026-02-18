#!/usr/bin/env python3
"""
Assemble v04 ssAAV transfer plasmid using cloning strategy.

Instead of replacing stuffer region end-to-end, this script:
1. Extracts VP1 from v04 RepCap
2. Builds EF1a-VP1-polyA cassette
3. Adds HindIII/XbaI flanking sites for cloning
4. Inserts cassette into backbone between HindIII (187) and XbaI (2278)
5. Preserves ITR annotations from v02 backbone

Only generates ssAAV version.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Restriction import Analysis, RestrictionBatch, AvrII, BspEI, BsmBI, BsrGI, BmtI, BstZ17I
from pathlib import Path
import sys

def main():
    print("="*80)
    print("ASSEMBLING v04 ssAAV TRANSFER PLASMID - CLONING STRATEGY")
    print("="*80)

    # Load files
    repcap_v04 = SeqIO.read('test_data/AAV9-RepCap-NOCAP-v04.gb', 'genbank')
    transfer_v03 = SeqIO.read('test_data/pGS-ssAAV-EF1A-VP1-rBG_v03.gb', 'genbank')
    backbone_v02 = SeqIO.read('test_data/pGS-ssAAV-ITR128-Amp-empty_v02.gb', 'genbank')

    print(f"\n✓ Loaded RepCap v04: {len(repcap_v04.seq)} bp")
    print(f"✓ Loaded transfer v03 (for EF1a/polyA): {len(transfer_v03.seq)} bp")
    print(f"✓ Loaded backbone v02: {len(backbone_v02.seq)} bp")

    # Extract VP1 from RepCap v04
    print("\n" + "="*80)
    print("STEP 1: EXTRACTING CASSETTE COMPONENTS")
    print("="*80)

    # Get VP1 from v04 RepCap
    vp1_v04_feature = None
    for feature in repcap_v04.features:
        if feature.qualifiers.get('label', [''])[0] == 'VP1':
            vp1_v04_feature = feature
            print(f"✓ Found VP1 in v04 RepCap: {feature.location.start+1}-{feature.location.end}")
            break

    if not vp1_v04_feature:
        print("❌ ERROR: Could not find VP1 in v04 RepCap")
        sys.exit(1)

    # Get EF1a (1194 bp version from v01) and polyA from transfer plasmids
    transfer_v01 = SeqIO.read('test_data/pGS-ssAAV-EF1A-VP1-rBG_v01.gb', 'genbank')

    ef1a_feature = None
    polya_feature = None

    # Get correct 1194 bp EF1α from v01
    for feature in transfer_v01.features:
        label = feature.qualifiers.get('label', [''])[0]
        if 'EF1' in label or 'EF-1' in label or 'EF1a' in label:
            ef1a_feature = feature
            print(f"✓ Found EF1α (1194 bp) in v01 transfer: {feature.location.start+1}-{feature.location.end}")
            break

    # Get polyA from v03
    for feature in transfer_v03.features:
        label = feature.qualifiers.get('label', [''])[0]
        if 'polyA' in label or 'rBG' in label:
            polya_feature = feature
            print(f"✓ Found polyA in v03 transfer: {feature.location.start+1}-{feature.location.end}")
            break

    if not all([ef1a_feature, vp1_v04_feature, polya_feature]):
        print("❌ ERROR: Could not find all required features")
        sys.exit(1)

    # Extract sequences
    ef1a_seq = ef1a_feature.extract(transfer_v01.seq)  # 1194 bp version from v01
    vp1_seq = vp1_v04_feature.extract(repcap_v04.seq)
    polya_seq = polya_feature.extract(transfer_v03.seq)

    print(f"\n✓ EF1α: {len(ef1a_seq)} bp")
    print(f"✓ VP1: {len(vp1_seq)} bp")
    print(f"✓ polyA: {len(polya_seq)} bp")

    # Build cassette with HindIII/XbaI flanking sites
    print("\n" + "="*80)
    print("STEP 2: BUILDING CASSETTE WITH CLONING SITES")
    print("="*80)

    # HindIII: AAGCTT, XbaI: TCTAGA
    hindiii_site = "AAGCTT"
    xbai_site = "TCTAGA"

    # Build: HindIII - EF1a - VP1 - polyA - XbaI
    cassette_seq = Seq(hindiii_site) + ef1a_seq + vp1_seq + polya_seq + Seq(xbai_site)

    print(f"✓ Cassette built: {len(cassette_seq)} bp")
    print(f"  Structure: HindIII-EF1α-VP1-polyA-XbaI")
    print(f"  HindIII at position 1")
    print(f"  XbaI at position {len(cassette_seq)-5}")

    # Clone into backbone
    print("\n" + "="*80)
    print("STEP 3: CLONING CASSETTE INTO BACKBONE")
    print("="*80)

    # Backbone HindIII at 187 (1-indexed), XbaI at 2278 (1-indexed)
    # 0-indexed: HindIII at 186, XbaI at 2277
    hindiii_pos_0idx = 186
    xbai_pos_0idx = 2277

    print(f"Backbone HindIII: position {hindiii_pos_0idx+1} (1-indexed)")
    print(f"Backbone XbaI: position {xbai_pos_0idx+1} (1-indexed)")

    # Extract backbone fragments
    # In proper cloning:
    # - Digest backbone with HindIII + XbaI, removing region between them
    # - Insert cassette with HindIII and XbaI sites
    # - Both restriction sites are regenerated in the final plasmid

    # Left fragment: everything BEFORE HindIII site
    left_fragment = backbone_v02.seq[:hindiii_pos_0idx]

    # Right fragment: everything AFTER XbaI site
    right_fragment = backbone_v02.seq[xbai_pos_0idx + 6:]

    # Cassette: full insert WITH HindIII and XbaI sites
    # (This simulates the ligation where both sites are regenerated)
    cassette_with_sites = cassette_seq  # Already has HindIII-insert-XbaI

    print(f"\n✓ Left fragment: 1-{len(left_fragment)} ({len(left_fragment)} bp, before HindIII)")
    print(f"✓ Cassette: {len(cassette_with_sites)} bp (HindIII-EF1α-VP1-polyA-XbaI)")
    print(f"✓ Right fragment: {len(right_fragment)} bp (after XbaI)")

    # Assemble: left + cassette + right
    # This recreates both HindIII and XbaI sites in the final plasmid
    final_seq = left_fragment + cassette_with_sites + right_fragment

    print(f"\n✓ Final plasmid: {len(final_seq)} bp")

    # Build features for the transfer plasmid
    print("\n" + "="*80)
    print("STEP 4: ANNOTATING FEATURES")
    print("="*80)

    features = []

    # Copy ITRs and backbone features from v02 backbone
    # We need to adjust positions for features after the insertion
    insertion_point = hindiii_pos_0idx  # Start of HindIII site (where we cut)
    insertion_length = len(cassette_with_sites)  # Full cassette with sites
    removed_length = (xbai_pos_0idx + 6) - insertion_point  # Everything from HindIII to after XbaI
    shift = insertion_length - removed_length

    print(f"Insertion point: {insertion_point} (HindIII position)")
    print(f"Insertion length: {insertion_length} bp (cassette with sites)")
    print(f"Removed length: {removed_length} bp (HindIII to XbaI region)")
    print(f"Net shift: {shift:+d} bp")

    for feature in backbone_v02.features:
        # If feature is entirely before insertion point, keep as-is
        if feature.location.end <= insertion_point:
            features.append(feature)
        # If feature is entirely after the removed region, shift it
        elif feature.location.start >= xbai_pos_0idx + 6:
            new_start = feature.location.start + shift
            new_end = feature.location.end + shift
            new_feature = SeqFeature(
                FeatureLocation(new_start, new_end, strand=feature.location.strand),
                type=feature.type,
                qualifiers=feature.qualifiers.copy()
            )
            features.append(new_feature)
        # If feature spans insertion, skip it (it's the removed stuffer region)

    # Add cassette features
    # Cassette starts right where HindIII begins (at insertion_point)
    cassette_start = len(left_fragment)

    # HindIII site (first 6 bp of cassette)
    hindiii_site_start = cassette_start
    hindiii_site_end = hindiii_site_start + 6
    features.append(SeqFeature(
        FeatureLocation(hindiii_site_start, hindiii_site_end, strand=1),
        type='misc_feature',
        qualifiers={'label': ['HindIII'], 'note': ['Cloning site (5\' junction)']}
    ))
    print(f"✓ HindIII site: {hindiii_site_start+1}-{hindiii_site_end}")

    # EF1α promoter (after HindIII)
    ef1a_start = hindiii_site_end
    ef1a_end = ef1a_start + len(ef1a_seq)
    features.append(SeqFeature(
        FeatureLocation(ef1a_start, ef1a_end, strand=1),
        type='promoter',
        qualifiers={
            'label': ['EF1α promoter'],
            'note': ['Human EEF1A1 promoter with first intron (1194 bp)']
        }
    ))
    print(f"✓ EF1α promoter: {ef1a_start+1}-{ef1a_end} (1194 bp)")

    # VP1
    vp1_start = ef1a_end
    vp1_end = vp1_start + len(vp1_seq)
    features.append(SeqFeature(
        FeatureLocation(vp1_start, vp1_end, strand=1),
        type='CDS',
        qualifiers={
            'label': ['VP1'],
            'note': ['AAV9 VP1 capsid protein with 6 unique restriction sites']
        }
    ))
    print(f"✓ VP1 CDS: {vp1_start+1}-{vp1_end}")

    # polyA
    polya_start = vp1_end
    polya_end = polya_start + len(polya_seq)
    features.append(SeqFeature(
        FeatureLocation(polya_start, polya_end, strand=1),
        type='polyA_signal',
        qualifiers={'label': ['bGH polyA']}
    ))
    print(f"✓ polyA: {polya_start+1}-{polya_end}")

    # XbaI site (last 6 bp of cassette)
    xbai_site_start = polya_end
    xbai_site_end = xbai_site_start + 6
    features.append(SeqFeature(
        FeatureLocation(xbai_site_start, xbai_site_end, strand=1),
        type='misc_feature',
        qualifiers={'label': ['XbaI'], 'note': ['Cloning site (3\' junction)']}
    ))
    print(f"✓ XbaI site: {xbai_site_start+1}-{xbai_site_end}")

    # Create final record
    transfer_plasmid = SeqRecord(
        final_seq,
        id="pGS-ssAAV-v04",
        name="pGS-ssAAV-v04",
        description="ssAAV transfer plasmid with EF1a-VP1-polyA (6 unique restriction sites)",
        features=features,
        annotations={'molecule_type': 'DNA', 'topology': 'circular'}
    )

    # Verify restriction sites
    print("\n" + "="*80)
    print("STEP 5: VERIFYING RESTRICTION SITES")
    print("="*80)

    # First verify cloning sites
    from Bio.Restriction import HindIII, XbaI
    cloning_enzymes = [HindIII, XbaI]
    cloning_batch = RestrictionBatch(cloning_enzymes)
    cloning_analysis = Analysis(cloning_batch, final_seq)
    cloning_sites = cloning_analysis.full()

    print("Cloning sites:")
    for enzyme in cloning_enzymes:
        enzyme_sites = cloning_sites.get(enzyme, [])
        count = len(enzyme_sites)
        positions = [p+1 for p in enzyme_sites]

        if count == 1:
            print(f"✅ {enzyme.__name__:<10s} {count} site @ position {positions[0]}")
        else:
            print(f"❌ {enzyme.__name__:<10s} {count} sites @ {positions} - SHOULD BE 1!")
            sys.exit(1)

    # Now verify the 6 unique sites in VP1
    print("\nVP1 unique restriction sites:")
    enzymes = [AvrII, BspEI, BsmBI, BsrGI, BmtI, BstZ17I]
    batch = RestrictionBatch(enzymes)
    analysis = Analysis(batch, final_seq)
    sites = analysis.full()

    all_unique = True
    for enzyme in enzymes:
        enzyme_sites = sites.get(enzyme, [])
        count = len(enzyme_sites)
        positions = [p+1 for p in enzyme_sites]

        if count == 1:
            print(f"✅ {enzyme.__name__:<10s} {count} site @ position {positions[0]}")
        else:
            print(f"❌ {enzyme.__name__:<10s} {count} sites @ {positions}")
            all_unique = False

    if not all_unique:
        print("\n❌ ERROR: Not all sites are unique!")
        sys.exit(1)

    print("\n✅ All cloning sites and VP1 sites verified!")

    # Save
    output_path = 'test_data/pGS-ssAAV-EF1A-VP1-rBG_v04.gb'
    SeqIO.write(transfer_plasmid, output_path, 'genbank')

    print("\n" + "="*80)
    print("✅ ASSEMBLY COMPLETE")
    print("="*80)
    print(f"✓ Saved: {output_path}")
    print(f"✓ Size: {len(final_seq)} bp")
    print(f"✓ All 6 restriction sites unique")
    print(f"✓ ITR annotations preserved from backbone")

if __name__ == '__main__':
    main()
