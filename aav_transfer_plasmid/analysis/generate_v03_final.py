#!/usr/bin/env python3
"""
Generate v03 RepCap and Transfer Plasmids - Final Version

Sites chosen:
- Region 1: XbaI (2 mutations)
- Region 2: BspEI (2 mutations)
- Region 3: BsmBI (1 mutation)
- Region 4: BsrGI (keep from v02)
- Region 5: BmtI (keep from v02)
- Region 6: BstZ17I (keep from v02)
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
    """Find where enzyme can be created with silent mutations"""
    try:
        enzyme = getattr(Restriction, enzyme_name)
        site_seq = str(enzyme.site).upper()

        for i in range(len(region_seq) - len(site_seq) + 1):
            window = region_seq[i:i+len(site_seq)].upper()

            # Check differences
            mutations = []
            for j, (curr, target) in enumerate(zip(window, site_seq)):
                if curr != target:
                    mutations.append((j, curr, target))

            if len(mutations) == 0 or len(mutations) > 2:
                continue

            # Check if all mutations are silent
            vp1_pos = region_start_in_vp1 + i
            all_silent = True
            silent_mutations = []

            for offset, curr, target in mutations:
                pos = vp1_pos + offset
                is_sil = is_silent_at_position(vp1_seq, pos, target)

                if not is_sil:
                    all_silent = False
                    break
                else:
                    # Get codon info for reporting
                    codon_idx = pos // 3
                    orig_codon = vp1_seq[codon_idx*3:codon_idx*3+3].upper()
                    mut_codon = list(orig_codon)
                    mut_codon[pos % 3] = target.upper()
                    mut_codon = ''.join(mut_codon)

                    silent_mutations.append({
                        'vp1_pos': pos,
                        'from': curr,
                        'to': target,
                        'orig_codon': orig_codon,
                        'mut_codon': mut_codon,
                        'aa': CODON_TABLE.get(orig_codon, '?')
                    })

            if all_silent:
                return {
                    'vp1_pos': vp1_pos,
                    'window': window,
                    'site': site_seq,
                    'mutations': silent_mutations,
                    'num_muts': len(mutations)
                }

        return None
    except Exception as e:
        print(f"Error finding {enzyme_name}: {e}")
        return None


def main():
    print("="*80)
    print("GENERATING v03 REPCAP AND TRANSFER PLASMIDS")
    print("="*80)
    print("\nSites to engineer:")
    print("  Region 1: XbaI")
    print("  Region 2: BspEI")
    print("  Region 3: BsmBI")
    print("  Region 4: BsrGI (keep)")
    print("  Region 5: BmtI (keep)")
    print("  Region 6: BstZ17I (keep)")

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

    print(f"✓ VP1: positions {vp1_start}-{vp1_feature.location.end} (0-indexed, {len(vp1_seq)} bp)")

    # Define new sites to find
    new_sites = [
        ("Region 1", 50, 410, "XbaI"),
        ("Region 2", 400, 550, "BspEI"),
        ("Region 3", 1120, 1350, "BsmBI"),
    ]

    # Known sites to keep (from original analysis, 1-indexed positions)
    kept_sites = [
        ("BsrGI", 3780, [(3783, 'C', 'A')]),   # GTC→GTA (V→V)
        ("BmtI", 3937, [(3939, 'C', 'T')]),    # GCC→GCT (A→A)
        ("BstZ17I", 4163, [(4164, 'A', 'T')]), # GGA→GGT (G→G)
    ]

    # Find new sites
    print("\n" + "="*80)
    print("STEP 1: FINDING NEW SITES WITH SILENT MUTATIONS")
    print("="*80)

    all_mutations = []
    found_sites = {}

    for region_name, start, end, enzyme_name in new_sites:
        print(f"\n{region_name}: {enzyme_name}")

        region_seq = vp1_seq[start:end]
        site_info = find_silent_site(enzyme_name, region_seq, start, vp1_seq)

        if not site_info:
            print(f"  ❌ Could not find silent mutations for {enzyme_name}")
            sys.exit(1)

        found_sites[enzyme_name] = site_info

        print(f"  ✓ Found at VP1 position {site_info['vp1_pos']}")
        print(f"    Current: {site_info['window']}")
        print(f"    Target:  {site_info['site']}")
        print(f"    Mutations: {site_info['num_muts']}")

        for mut in site_info['mutations']:
            plasmid_pos = vp1_start + mut['vp1_pos']
            print(f"      Position {mut['vp1_pos']} (plasmid: {plasmid_pos}): "
                  f"{mut['from']}→{mut['to']} | "
                  f"{mut['orig_codon']}→{mut['mut_codon']} ({mut['aa']}→{mut['aa']})")

            all_mutations.append({
                'position': plasmid_pos,
                'from': mut['from'],
                'to': mut['to']
            })

    # Add kept sites
    print("\n" + "="*80)
    print("STEP 2: ADDING KEPT SITES FROM REGIONS 4-6")
    print("="*80)

    for enzyme_name, site_pos, mutations in kept_sites:
        print(f"\n{enzyme_name} @ position {site_pos}")
        for pos_1idx, from_base, to_base in mutations:
            pos_0idx = pos_1idx - 1
            all_mutations.append({
                'position': pos_0idx,
                'from': from_base,
                'to': to_base
            })
            print(f"  Position {pos_1idx} (0-idx: {pos_0idx}): {from_base}→{to_base}")

    # Apply mutations
    print("\n" + "="*80)
    print(f"STEP 3: APPLYING {len(all_mutations)} MUTATIONS TO REPCAP")
    print("="*80)

    seq_list = list(str(original_repcap.seq).upper())

    for mut in all_mutations:
        pos = mut['position']
        old_val = seq_list[pos]

        if old_val != mut['from'].upper():
            print(f"⚠️  Warning: Position {pos} has {old_val}, expected {mut['from']}")

        seq_list[pos] = mut['to'].upper()

    new_seq = Seq(''.join(seq_list))

    # Create new RepCap record
    new_repcap = SeqRecord(
        new_seq,
        id="AAV9-RepCap-v03",
        name="AAV9-RepCap-v03",
        description="AAV9 RepCap with 6 restriction sites (XbaI, BspEI, BsmBI, BsrGI, BmtI, BstZ17I)",
        annotations=original_repcap.annotations
    )
    new_repcap.features = original_repcap.features[:]

    # Save RepCap
    output_file = base_path / 'AAV9-RepCap-NOCAP-v03.gb'
    SeqIO.write(new_repcap, str(output_file), 'genbank')
    print(f"\n✓ Saved RepCap: {output_file}")

    # Verify sites
    print("\n" + "="*80)
    print("STEP 4: VERIFYING SITES IN REPCAP")
    print("="*80)

    all_enzymes = ['XbaI', 'BspEI', 'BsmBI', 'BsrGI', 'BmtI', 'BstZ17I']

    for enzyme_name in all_enzymes:
        try:
            enzyme = getattr(Restriction, enzyme_name)
            analysis = Analysis(RestrictionBatch([enzyme]), new_seq)
            sites = analysis.full()

            count = len(sites.get(enzyme, []))
            positions = sites.get(enzyme, [])

            if count == 1:
                print(f"  ✅ {enzyme_name:<12} 1 site @ position {positions[0] + 1}")
            else:
                print(f"  ❌ {enzyme_name:<12} {count} site(s)")

        except Exception as e:
            print(f"  ❌ {enzyme_name:<12} ERROR: {e}")

    print("\n✓ RepCap v03 complete!")
    print(f"✓ Saved: {output_file}")


if __name__ == '__main__':
    main()
