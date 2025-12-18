#!/usr/bin/env python3
"""
Find enzymes that are:
1. Unique in both RepCap and transfer plasmids (0 sites currently)
2. Can be created with SILENT mutations only
"""

from Bio import SeqIO
from Bio.Restriction import Restriction
from pathlib import Path

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

# Candidates that are unique in both contexts
DUAL_UNIQUE = {
    'Region 1': ['AvrII', 'BsmBI', 'BspEI', 'Esp3I', 'HpaI', 'MluI', 'SnaBI', 'XbaI'],
    'Region 2': ['BsmBI', 'BspEI', 'Esp3I', 'XbaI'],
    'Region 3': ['BsmBI', 'EcoRV', 'Esp3I', 'HpaI', 'MluI', 'NsiI', 'SnaBI']
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
    """Find if enzyme can be created with silent mutations"""
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

            for offset, curr, target in mutations:
                if not is_silent_at_position(vp1_seq, vp1_pos + offset, target):
                    all_silent = False
                    break

            if all_silent:
                return {
                    'vp1_pos': vp1_pos,
                    'window': window,
                    'site': site_seq,
                    'mutations': mutations,
                    'num_muts': len(mutations)
                }

        return None
    except:
        return None


def main():
    # Load RepCap
    base_path = Path('test_data')
    original = SeqIO.read(str(base_path / 'BASE-DRAFT-AAV9-RepCap-NOCAP.gb'), 'genbank')

    # Find VP1
    vp1_feature = None
    for feature in original.features:
        if feature.qualifiers.get('label', [''])[0] == 'VP1':
            vp1_feature = feature
            break

    vp1_seq = str(vp1_feature.extract(original.seq))

    # Define regions
    regions = {
        'Region 1': (50, 410),
        'Region 2': (400, 550),
        'Region 3': (1120, 1350)
    }

    print("="*80)
    print("FINDING DUAL-UNIQUE ENZYMES WITH SILENT MUTATIONS")
    print("="*80)

    for region_name in ['Region 1', 'Region 2', 'Region 3']:
        start, end = regions[region_name]
        enzymes = DUAL_UNIQUE[region_name]
        region_seq = vp1_seq[start:end]

        print(f"\n{region_name} ({len(enzymes)} candidates):")
        print("-" * 60)

        viable = []

        for enzyme_name in enzymes:
            result = find_silent_site(enzyme_name, region_seq, start, vp1_seq)

            if result:
                viable.append((enzyme_name, result))
                print(f"  ✅ {enzyme_name:<12} {result['num_muts']} mutation(s) @ VP1 pos {result['vp1_pos']}")
            else:
                print(f"  ❌ {enzyme_name:<12} No silent option found")

        if viable:
            print(f"\n  VIABLE OPTIONS FOR {region_name}:")
            for enzyme_name, result in viable:
                print(f"    → {enzyme_name}")

    print("\n" + "="*80)
    print("RECOMMENDED COMBINATION")
    print("="*80)

    # Find best combination
    for region_name in ['Region 1', 'Region 2', 'Region 3']:
        start, end = regions[region_name]
        enzymes = DUAL_UNIQUE[region_name]
        region_seq = vp1_seq[start:end]

        for enzyme_name in enzymes:
            result = find_silent_site(enzyme_name, region_seq, start, vp1_seq)
            if result:
                print(f"\n{region_name}: {enzyme_name}")
                break


if __name__ == '__main__':
    main()
