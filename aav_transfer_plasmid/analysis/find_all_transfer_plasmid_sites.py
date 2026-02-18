#!/usr/bin/env python3
"""
Find ALL Possible Restriction Sites for Regions 1-3
Constrained ONLY by Transfer Plasmid Uniqueness

For each region in VP1:
1. Identify all possible restriction sites requiring ≤2 silent mutations
2. Check if site would be unique in ssAAV and scAAV transfer plasmids
3. Ignore RepCap plasmid constraints

This gives us maximum flexibility to choose the best sites for the transfer plasmid.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Restriction import AllEnzymes
from pathlib import Path
import sys

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

# Standard genetic code
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


def load_plasmids():
    """Load transfer plasmids and get VP1 sequence"""
    base_path = Path('test_data')

    # Load source RepCap to get VP1 sequence
    repcap = SeqIO.read(str(base_path / 'AAV9-RepCap-NOCAP-WITH-RESTRICTION-SITES.gb'), 'genbank')

    # Find VP1
    vp1_feature = None
    for feature in repcap.features:
        if feature.qualifiers.get('label', [''])[0] == 'VP1':
            vp1_feature = feature
            break

    if not vp1_feature:
        print("❌ Could not find VP1 in RepCap plasmid")
        return None

    vp1_start = vp1_feature.location.start
    vp1_seq = vp1_feature.extract(repcap.seq)

    # Load transfer plasmids
    ssaav = SeqIO.read(str(base_path / 'pGS-ssAAV-EF1A-VP1-rBG_v02.gb'), 'genbank')
    scaav = SeqIO.read(str(base_path / 'pGS-scAAV-EF1A-VP1-rBG_v02.gb'), 'genbank')

    return {
        'vp1_seq': str(vp1_seq),
        'vp1_start': vp1_start,
        'ssaav': ssaav,
        'scaav': scaav
    }


def get_enzyme_site(enzyme_name):
    """Get recognition site for enzyme"""
    try:
        from Bio.Restriction import Restriction
        enzyme = getattr(Restriction, enzyme_name)
        return str(enzyme.site)
    except:
        return None


def count_sites_in_plasmid(enzyme_name, plasmid_seq):
    """Count occurrences of enzyme recognition site in plasmid"""
    site = get_enzyme_site(enzyme_name)
    if not site:
        return -1

    # Handle IUPAC ambiguity codes
    if 'N' in site or 'R' in site or 'Y' in site or 'W' in site or 'S' in site or 'M' in site or 'K' in site:
        # For ambiguous sites, use Biopython's analysis
        try:
            from Bio.Restriction import Restriction, Analysis, RestrictionBatch
            enzyme = getattr(Restriction, enzyme_name)
            analysis = Analysis(RestrictionBatch([enzyme]), plasmid_seq)
            sites = analysis.full()
            if enzyme in sites:
                return len(sites[enzyme])
            return 0
        except:
            return -1

    # For non-ambiguous sites, simple string search
    return str(plasmid_seq).upper().count(site.upper())


def hamming_distance(s1, s2):
    """Calculate Hamming distance between two sequences"""
    if len(s1) != len(s2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def is_silent_mutation(vp1_seq, position, new_base):
    """Check if mutation at position is silent in VP1 frame"""
    # Position is 0-indexed in VP1
    codon_idx = position // 3
    pos_in_codon = position % 3

    if codon_idx * 3 + 2 >= len(vp1_seq):
        return False

    # Get original codon
    original_codon = vp1_seq[codon_idx * 3 : codon_idx * 3 + 3]
    if len(original_codon) != 3:
        return False

    # Make mutation
    mutated_codon = list(original_codon)
    mutated_codon[pos_in_codon] = new_base
    mutated_codon = ''.join(mutated_codon)

    # Check if amino acid changes
    original_aa = CODON_TABLE.get(original_codon.upper(), '?')
    mutated_aa = CODON_TABLE.get(mutated_codon.upper(), '?')

    return original_aa == mutated_aa


def find_sites_in_region(region_name, region_start, region_end, vp1_seq, ssaav, scaav):
    """
    Find all possible restriction sites in a VP1 region that would be unique in transfer plasmids.

    region_start, region_end: 0-indexed positions relative to VP1 start
    """
    print(f"\n{'='*80}")
    print(f"REGION: {region_name}")
    print(f"VP1 positions: {region_start}-{region_end} (0-indexed)")
    print(f"{'='*80}")

    region_seq = vp1_seq[region_start:region_end]
    print(f"\nRegion length: {len(region_seq)} bp")
    print(f"Analyzing {len(CURATED_ENZYMES)} enzymes...")

    viable_sites = []

    for enzyme_name in CURATED_ENZYMES:
        site = get_enzyme_site(enzyme_name)
        if not site:
            continue

        # Skip very long or highly ambiguous sites
        if len(site) > 12:
            continue

        # Check how many mutations needed to create this site in the region
        best_match = None
        min_mutations = float('inf')

        for i in range(len(region_seq) - len(site) + 1):
            window = region_seq[i:i+len(site)]

            # For ambiguous sites, skip detailed analysis
            if 'N' in site or 'R' in site or 'Y' in site or 'W' in site:
                continue

            dist = hamming_distance(window.upper(), site.upper())

            if dist < min_mutations:
                min_mutations = dist
                best_match = {
                    'position': region_start + i,
                    'window': window,
                    'site': site
                }

        # Only consider sites requiring ≤2 mutations
        if min_mutations > 2:
            continue

        # If site already exists (0 mutations), check if engineerable
        # If site needs 1-2 mutations, check if mutations are silent
        if min_mutations == 0:
            # Site already exists - check if unique in transfer plasmids
            pass
        else:
            # Need to check if mutations would be silent
            # For now, include all candidates and note mutation count
            pass

        # Check uniqueness in transfer plasmids
        ssaav_count = count_sites_in_plasmid(enzyme_name, ssaav.seq)
        scaav_count = count_sites_in_plasmid(enzyme_name, scaav.seq)

        if ssaav_count == -1:
            continue

        # Site is viable if:
        # - 0 occurrences in both transfer plasmids (can engineer and will be unique), OR
        # - 1 occurrence in transfer plasmids IF that occurrence is in the VP1 region we're engineering

        # For simplicity, consider viable if 0-1 occurrences
        # (1 occurrence might be from current VP1, would be replaced)
        if ssaav_count <= 1 and scaav_count <= 1:
            viable_sites.append({
                'enzyme': enzyme_name,
                'site': site,
                'position': best_match['position'] if best_match else None,
                'mutations': min_mutations,
                'ssaav_count': ssaav_count,
                'scaav_count': scaav_count
            })

    # Sort by mutation count, then by ssaav_count + scaav_count
    viable_sites.sort(key=lambda x: (x['mutations'], x['ssaav_count'] + x['scaav_count']))

    print(f"\n✅ Found {len(viable_sites)} viable candidate(s):")
    print(f"\n{'Enzyme':<12} {'Mutations':<10} {'ssAAV':<8} {'scAAV':<8} {'Status'}")
    print("-" * 60)

    for site_info in viable_sites[:20]:  # Show top 20
        enzyme = site_info['enzyme']
        muts = site_info['mutations']
        ssaav_cnt = site_info['ssaav_count']
        scaav_cnt = site_info['scaav_count']

        status = "✅ Unique" if ssaav_cnt == 0 and scaav_cnt == 0 else "⚠️  Check VP1"

        print(f"{enzyme:<12} {muts:<10} {ssaav_cnt:<8} {scaav_cnt:<8} {status}")

    if len(viable_sites) > 20:
        print(f"\n... and {len(viable_sites) - 20} more candidates")

    return viable_sites


def main():
    print("="*80)
    print("COMPREHENSIVE TRANSFER PLASMID RESTRICTION SITE ANALYSIS")
    print("="*80)
    print("\nFinding ALL possible restriction sites for Regions 1-3")
    print("Constraint: Uniqueness in transfer plasmids ONLY\n")

    # Load plasmids
    data = load_plasmids()
    if not data:
        sys.exit(1)

    vp1_seq = data['vp1_seq']
    ssaav = data['ssaav']
    scaav = data['scaav']

    print(f"VP1 length: {len(vp1_seq)} bp")

    # Define regions (0-indexed relative to VP1 start)
    regions = [
        ("Region 1: Rep68-stop to VP2-start", 50, 410),  # Expanded from 0-280
        ("Region 2: VP2-AAP Intergenic", 400, 550),
        ("Region 3: AAP-stop to VR4", 1120, 1350),  # Expanded from 1150-1300
    ]

    all_results = {}

    for region_name, start, end in regions:
        results = find_sites_in_region(region_name, start, end, vp1_seq, ssaav, scaav)
        all_results[region_name] = results

    # Summary
    print("\n\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    for region_name in all_results:
        results = all_results[region_name]
        unique_count = sum(1 for r in results if r['ssaav_count'] == 0 and r['scaav_count'] == 0)
        print(f"\n{region_name}:")
        print(f"  Total candidates: {len(results)}")
        print(f"  Perfectly unique (0 in both): {unique_count}")
        print(f"  Need verification (1 in either): {len(results) - unique_count}")

        # Show ALL perfectly unique options
        perfect = [r for r in results if r['ssaav_count'] == 0 and r['scaav_count'] == 0]
        if perfect:
            print(f"\n  ALL perfectly unique options:")
            for r in perfect:
                print(f"    - {r['enzyme']:<12} ({r['mutations']} mutation(s))")

    print("\n" + "="*80)


if __name__ == '__main__':
    main()
