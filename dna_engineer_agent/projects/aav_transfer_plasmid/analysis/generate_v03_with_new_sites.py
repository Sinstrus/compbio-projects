#!/usr/bin/env python3
"""
Generate v03 RepCap and Transfer Plasmids with New Site Choices

NEW SITES (Regions 1-3):
- Region 1: SnaBI
- Region 2: XbaI
- Region 3: MluI

KEEP EXISTING (Regions 4-6):
- Region 4: BsrGI @ 3780
- Region 5: BmtI @ 3937
- Region 6: BstZ17I @ 4163

Process:
1. Load ORIGINAL RepCap (no mutations)
2. Apply mutations for all 6 sites (3 new + 3 kept from first round)
3. Generate transfer plasmids
4. Verify uniqueness
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Restriction import Analysis, RestrictionBatch, Restriction
from pathlib import Path
import sys

# Genetic code
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def find_enzyme_site_location(enzyme_name, region_seq, vp1_start_in_plasmid, region_start_in_vp1):
    """
    Find where an enzyme recognition sequence should be placed in a region.
    Returns (position_in_plasmid, mutation_info) or None
    """
    try:
        enzyme = getattr(Restriction, enzyme_name)
        site_seq = str(enzyme.site)

        print(f"\n  Searching for {enzyme_name} site: {site_seq}")

        # Search for best match (minimum Hamming distance)
        best_pos = None
        best_dist = float('inf')
        best_window = None

        for i in range(len(region_seq) - len(site_seq) + 1):
            window = region_seq[i:i+len(site_seq)]
            dist = sum(c1 != c2 for c1, c2 in zip(window.upper(), site_seq.upper()))

            if dist < best_dist:
                best_dist = dist
                best_pos = i
                best_window = window

        if best_dist > 2:
            print(f"    ❌ No good match found (min distance: {best_dist})")
            return None

        # Calculate absolute position in plasmid
        abs_pos_in_vp1 = region_start_in_vp1 + best_pos
        abs_pos_in_plasmid = vp1_start_in_plasmid + abs_pos_in_vp1

        print(f"    Found at VP1 position {abs_pos_in_vp1} (plasmid position {abs_pos_in_plasmid})")
        print(f"    Current: {best_window}")
        print(f"    Target:  {site_seq}")
        print(f"    Mutations needed: {best_dist}")

        # Identify specific mutations
        mutations = []
        for i, (current, target) in enumerate(zip(best_window.upper(), site_seq.upper())):
            if current != target:
                mut_pos = abs_pos_in_plasmid + i
                mutations.append({
                    'position': mut_pos,
                    'from': current,
                    'to': target,
                    'offset_in_site': i
                })
                print(f"      Mutation {len(mutations)}: position {mut_pos} ({current}→{target})")

        return {
            'enzyme': enzyme_name,
            'site': site_seq,
            'position': abs_pos_in_plasmid,
            'vp1_position': abs_pos_in_vp1,
            'mutations': mutations,
            'window': best_window
        }

    except Exception as e:
        print(f"    ❌ Error: {e}")
        return None


def verify_silent_mutations(plasmid_seq, vp1_start, mutations):
    """Verify that mutations are silent in VP1 frame"""
    print(f"\n  Verifying mutations are silent:")

    all_silent = True

    for mut in mutations:
        pos = mut['position']
        vp1_offset = pos - vp1_start

        # Get codon
        codon_start = vp1_start + (vp1_offset // 3) * 3
        codon_end = codon_start + 3

        original_codon = str(plasmid_seq[codon_start:codon_end]).upper()

        # Apply mutation
        mutated_seq = list(str(plasmid_seq).upper())
        mutated_seq[pos] = mut['to']
        mutated_codon = ''.join(mutated_seq[codon_start:codon_end])

        original_aa = CODON_TABLE.get(original_codon, '?')
        mutated_aa = CODON_TABLE.get(mutated_codon, '?')

        is_silent = (original_aa == mutated_aa)
        status = "✓ Silent" if is_silent else "✗ NOT SILENT"

        print(f"    Position {pos}: {original_codon} → {mutated_codon} ({original_aa} → {mutated_aa}) {status}")

        if not is_silent:
            all_silent = False

    return all_silent


def apply_mutations(plasmid_seq, all_mutations):
    """Apply all mutations to plasmid sequence"""
    seq_list = list(str(plasmid_seq).upper())

    for mut in all_mutations:
        pos = mut['position']
        seq_list[pos] = mut['to']

    return Seq(''.join(seq_list))


def main():
    print("="*80)
    print("GENERATING v03 REPCAP AND TRANSFER PLASMIDS")
    print("="*80)
    print("\nNew sites for Regions 1-3:")
    print("  Region 1: SnaBI")
    print("  Region 2: XbaI")
    print("  Region 3: MluI")
    print("\nKeeping sites from Regions 4-6:")
    print("  Region 4: BsrGI @ 3780")
    print("  Region 5: BmtI @ 3937")
    print("  Region 6: BstZ17I @ 4163")

    # Load ORIGINAL RepCap plasmid (no mutations)
    base_path = Path('test_data')
    original_file = base_path / 'BASE-DRAFT-AAV9-RepCap-NOCAP.gb'

    if not original_file.exists():
        print(f"\n❌ Original RepCap not found: {original_file}")
        sys.exit(1)

    print(f"\n✓ Loading original RepCap: {original_file}")
    original_repcap = SeqIO.read(str(original_file), 'genbank')

    # Find VP1
    vp1_feature = None
    for feature in original_repcap.features:
        if feature.qualifiers.get('label', [''])[0] == 'VP1':
            vp1_feature = feature
            break

    if not vp1_feature:
        print("❌ Could not find VP1 in RepCap")
        sys.exit(1)

    vp1_start = vp1_feature.location.start
    vp1_seq = str(vp1_feature.extract(original_repcap.seq))

    print(f"✓ Found VP1: positions {vp1_start}-{vp1_feature.location.end} ({len(vp1_seq)} bp)")

    # Define regions (0-indexed relative to VP1 start)
    regions_new = [
        (1, "Region 1: Rep68-stop to VP2-start", 50, 410, "SnaBI"),
        (2, "Region 2: VP2-AAP Intergenic", 400, 550, "XbaI"),
        (3, "Region 3: AAP-stop to VR4", 1120, 1350, "MluI"),
    ]

    # Regions to keep (with known positions from first round)
    regions_keep = [
        (4, "Region 4: VR4 to VR5", "BsrGI", 3780, [(3783, 'C', 'A')]),  # GTC → GTA
        (5, "Region 5: VR5 to VR8", "BmtI", 3937, [(3939, 'C', 'T')]),   # GCC → GCT
        (6, "Region 6: Post-VR8", "BstZ17I", 4163, [(4164, 'A', 'T')]),  # GGA → GGT
    ]

    # Find sites for new regions
    print("\n" + "="*80)
    print("STEP 1: FINDING NEW RESTRICTION SITES")
    print("="*80)

    all_site_info = []
    all_mutations = []

    for region_num, region_name, start, end, enzyme_name in regions_new:
        print(f"\n{region_name} → {enzyme_name}")

        region_seq = vp1_seq[start:end]
        site_info = find_enzyme_site_location(enzyme_name, region_seq, vp1_start, start)

        if not site_info:
            print(f"  ❌ Failed to find {enzyme_name} site")
            continue

        # Verify silent
        if not verify_silent_mutations(original_repcap.seq, vp1_start, site_info['mutations']):
            print(f"  ❌ {enzyme_name} mutations are NOT silent!")
            continue

        all_site_info.append(site_info)
        all_mutations.extend(site_info['mutations'])
        print(f"  ✓ {enzyme_name} site verified")

    # Add kept regions
    print("\n" + "="*80)
    print("STEP 2: ADDING KEPT RESTRICTION SITES")
    print("="*80)

    for region_num, region_name, enzyme_name, position, mutations_list in regions_keep:
        print(f"\n{region_name} → {enzyme_name} @ {position}")

        # Add mutations
        for pos, from_base, to_base in mutations_list:
            all_mutations.append({
                'position': pos,
                'from': from_base,
                'to': to_base
            })
            print(f"  Mutation: position {pos} ({from_base}→{to_base})")

        # Verify silent
        verify_silent_mutations(original_repcap.seq, vp1_start,
                              [{'position': p, 'from': f, 'to': t} for p, f, t in mutations_list])

    # Apply all mutations
    print("\n" + "="*80)
    print("STEP 3: APPLYING MUTATIONS TO REPCAP")
    print("="*80)

    print(f"\nApplying {len(all_mutations)} total mutations...")
    new_repcap_seq = apply_mutations(original_repcap.seq, all_mutations)

    # Create new RepCap record
    new_repcap = SeqRecord(
        new_repcap_seq,
        id=original_repcap.id,
        name="AAV9-RepCap-v03",
        description="AAV9 RepCap with 6 restriction sites (SnaBI, XbaI, MluI, BsrGI, BmtI, BstZ17I)",
        annotations=original_repcap.annotations
    )
    new_repcap.features = original_repcap.features[:]

    # Save new RepCap
    output_file = base_path / 'AAV9-RepCap-NOCAP-v03.gb'
    SeqIO.write(new_repcap, str(output_file), 'genbank')
    print(f"✓ Saved new RepCap: {output_file}")

    # Verify all 6 sites are unique in new RepCap
    print("\n" + "="*80)
    print("STEP 4: VERIFYING UNIQUENESS IN REPCAP")
    print("="*80)

    all_enzymes = ['SnaBI', 'XbaI', 'MluI', 'BsrGI', 'BmtI', 'BstZ17I']

    for enzyme_name in all_enzymes:
        try:
            enzyme = getattr(Restriction, enzyme_name)
            analysis = Analysis(RestrictionBatch([enzyme]), new_repcap_seq)
            sites = analysis.full()

            if enzyme in sites:
                count = len(sites[enzyme])
                status = "✓ Unique" if count == 1 else f"✗ {count} sites"
                print(f"  {enzyme_name:<12} {count} site(s)  {status}")
            else:
                print(f"  {enzyme_name:<12} 0 sites  ✗ NOT FOUND")
        except:
            print(f"  {enzyme_name:<12} ERROR")

    print("\n✓ RepCap generation complete!")
    print("\n" + "="*80)
    print("Next: Generate transfer plasmids using this RepCap")
    print("="*80)


if __name__ == '__main__':
    main()
