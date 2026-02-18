#!/usr/bin/env python3
"""
Generate v03 RepCap with Silent Mutations Only

Search for positions where the 1-mutation change is actually silent.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import Restriction
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


def hamming_distance(s1, s2):
    """Calculate Hamming distance"""
    if len(s1) != len(s2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(s1.upper(), s2.upper()))


def is_silent_at_position(vp1_seq, position, new_base):
    """Check if mutation at VP1 position is silent"""
    codon_idx = position // 3
    pos_in_codon = position % 3

    if codon_idx * 3 + 2 >= len(vp1_seq):
        return False

    original_codon = vp1_seq[codon_idx * 3 : codon_idx * 3 + 3].upper()
    if len(original_codon) != 3:
        return False

    mutated_codon = list(original_codon)
    mutated_codon[pos_in_codon] = new_base.upper()
    mutated_codon = ''.join(mutated_codon)

    original_aa = CODON_TABLE.get(original_codon, '?')
    mutated_aa = CODON_TABLE.get(mutated_codon, '?')

    return original_aa == mutated_aa and original_aa != '?'


def find_silent_site(enzyme_name, region_seq, region_start_in_vp1, vp1_seq):
    """
    Find a position in the region where the enzyme site can be created
    with a silent mutation.
    """
    try:
        enzyme = getattr(Restriction, enzyme_name)
        site_seq = str(enzyme.site).upper()

        print(f"\n  Searching for {enzyme_name} (recognition: {site_seq})...")

        # Find all potential positions
        candidates = []

        for i in range(len(region_seq) - len(site_seq) + 1):
            window = region_seq[i:i+len(site_seq)].upper()
            dist = hamming_distance(window, site_seq)

            if dist <= 2:  # Max 2 mutations
                vp1_pos = region_start_in_vp1 + i

                # Check each mutation
                mutations = []
                all_silent = True

                for j, (curr, target) in enumerate(zip(window, site_seq)):
                    if curr != target:
                        is_silent = is_silent_at_position(vp1_seq, vp1_pos + j, target)

                        mutations.append({
                            'vp1_pos': vp1_pos + j,
                            'offset': j,
                            'from': curr,
                            'to': target,
                            'silent': is_silent
                        })

                        if not is_silent:
                            all_silent = False

                if all_silent and len(mutations) > 0:
                    candidates.append({
                        'vp1_pos': vp1_pos,
                        'window': window,
                        'mutations': mutations,
                        'num_mutations': len(mutations)
                    })

        if not candidates:
            print(f"    ❌ No silent positions found")
            return None

        # Sort by fewest mutations
        candidates.sort(key=lambda x: x['num_mutations'])
        best = candidates[0]

        print(f"    ✓ Found at VP1 position {best['vp1_pos']}")
        print(f"      Current: {best['window']}")
        print(f"      Target:  {site_seq}")
        print(f"      Mutations: {best['num_mutations']} (all silent)")

        for mut in best['mutations']:
            print(f"        Position {mut['vp1_pos']}: {mut['from']}→{mut['to']}")

        return best

    except Exception as e:
        print(f"    ❌ Error: {e}")
        return None


def main():
    print("="*80)
    print("GENERATING v03 REPCAP WITH SILENT MUTATIONS")
    print("="*80)

    # Load original RepCap
    base_path = Path('test_data')
    original_file = base_path / 'BASE-DRAFT-AAV9-RepCap-NOCAP.gb'

    print(f"\n✓ Loading: {original_file}")
    original_repcap = SeqIO.read(str(original_file), 'genbank')

    # Find VP1
    vp1_feature = None
    for feature in original_repcap.features:
        if feature.qualifiers.get('label', [''])[0] == 'VP1':
            vp1_feature = feature
            break

    vp1_start = vp1_feature.location.start  # 0-indexed
    vp1_seq = str(vp1_feature.extract(original_repcap.seq))

    print(f"✓ VP1: positions {vp1_start}-{vp1_feature.location.end} ({len(vp1_seq)} bp)")

    # Define regions
    new_sites = [
        ("Region 1: Rep68-stop to VP2-start", 50, 410, "SnaBI"),
        ("Region 2: VP2-AAP Intergenic", 400, 550, "XbaI"),
        ("Region 3: AAP-stop to VR4", 1120, 1350, "MluI"),
    ]

    # Known mutations from first round (Regions 4-6)
    # These are from the original FINAL summary - positions are 1-indexed
    kept_sites = [
        ("BsrGI", 3780, [(3783, 'C', 'A')]),   # GTC→GTA (V→V)
        ("BmtI", 3937, [(3939, 'C', 'T')]),    # GCC→GCT (A→A)
        ("BstZ17I", 4163, [(4164, 'A', 'T')]), # GGA→GGT (G→G)
    ]

    # Find new sites
    print("\n" + "="*80)
    print("FINDING NEW RESTRICTION SITES (WITH SILENT MUTATIONS)")
    print("="*80)

    all_mutations = []
    found_sites = {}

    for region_name, start, end, enzyme_name in new_sites:
        print(f"\n{region_name} → {enzyme_name}")

        region_seq = vp1_seq[start:end]
        site_info = find_silent_site(enzyme_name, region_seq, start, vp1_seq)

        if site_info:
            found_sites[enzyme_name] = site_info
            # Convert VP1 positions to plasmid positions
            for mut in site_info['mutations']:
                all_mutations.append({
                    'position': vp1_start + mut['vp1_pos'],  # Convert to plasmid coords
                    'from': mut['from'],
                    'to': mut['to']
                })
        else:
            print(f"  ⚠️  Could not find silent mutation for {enzyme_name}")

    # Add kept sites
    print("\n" + "="*80)
    print("ADDING KEPT SITES FROM REGIONS 4-6")
    print("="*80)

    for enzyme_name, site_pos, mutations in kept_sites:
        print(f"\n{enzyme_name} @ position {site_pos}")
        for pos_1idx, from_base, to_base in mutations:
            # Positions in kept_sites are 1-indexed, convert to 0-indexed
            pos_0idx = pos_1idx - 1
            all_mutations.append({
                'position': pos_0idx,
                'from': from_base,
                'to': to_base
            })
            print(f"  Position {pos_1idx} (0-idx: {pos_0idx}): {from_base}→{to_base}")

    # Apply mutations
    print("\n" + "="*80)
    print(f"APPLYING {len(all_mutations)} MUTATIONS")
    print("="*80)

    seq_list = list(str(original_repcap.seq).upper())

    for mut in all_mutations:
        pos = mut['position']
        old_val = seq_list[pos]

        if old_val != mut['from'].upper():
            print(f"⚠️  Warning: Position {pos} has {old_val}, expected {mut['from']}")

        seq_list[pos] = mut['to'].upper()
        print(f"  Position {pos}: {mut['from']}→{mut['to']}")

    new_seq = Seq(''.join(seq_list))

    # Create new record
    new_repcap = SeqRecord(
        new_seq,
        id="AAV9-RepCap-v03",
        name="AAV9-RepCap-v03",
        description="AAV9 RepCap with 6 restriction sites (SnaBI, XbaI, MluI, BsrGI, BmtI, BstZ17I)",
        annotations=original_repcap.annotations
    )
    new_repcap.features = original_repcap.features[:]

    # Save
    output_file = base_path / 'AAV9-RepCap-NOCAP-v03.gb'
    SeqIO.write(new_repcap, str(output_file), 'genbank')
    print(f"\n✓ Saved: {output_file}")

    print("\n" + "="*80)
    print("VERIFICATION")
    print("="*80)

    # Check each site
    from Bio.Restriction import Analysis, RestrictionBatch

    all_enzymes = ['SnaBI', 'XbaI', 'MluI', 'BsrGI', 'BmtI', 'BstZ17I']

    for enzyme_name in all_enzymes:
        try:
            enzyme = getattr(Restriction, enzyme_name)
            analysis = Analysis(RestrictionBatch([enzyme]), new_seq)
            sites = analysis.full()

            count = len(sites.get(enzyme, []))
            status = "✓" if count == 1 else "✗"
            print(f"  {enzyme_name:<12} {count} site(s)  {status}")

        except Exception as e:
            print(f"  {enzyme_name:<12} ERROR: {e}")

    print("\n✓ Done!")


if __name__ == '__main__':
    main()
