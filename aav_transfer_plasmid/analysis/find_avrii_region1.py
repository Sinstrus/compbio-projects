#!/usr/bin/env python3
"""
Find the exact position and silent mutation to create AvrII site in VP1 Region 1.
AvrII recognition: CCTAGG
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Restriction import AvrII
import sys

def get_codon_table():
    """Standard genetic code."""
    return {
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

def find_vp1_start(record):
    """Find VP1 CDS start position."""
    for feature in record.features:
        if feature.type == 'CDS':
            if 'VP1' in feature.qualifiers.get('note', [''])[0] or \
               'VP1' in feature.qualifiers.get('label', [''])[0]:
                return feature.location.start
    return None

# Load RepCap
repcap = SeqIO.read('test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb', 'genbank')
vp1_start = find_vp1_start(repcap)

if not vp1_start:
    print("ERROR: Could not find VP1 start")
    sys.exit(1)

print(f"VP1 starts at position {vp1_start + 1}")
print()

# Region 1 is around VP1 aa 59 (based on v03 XbaI at VP1 aa 59)
# Search VP1 positions 1-150 to find all possible AvrII sites
region1_start_aa = 1
region1_end_aa = 150

target_seq = "CCTAGG"  # AvrII recognition
codon_table = get_codon_table()

print(f"Searching for AvrII sites in Region 1 (VP1 aa {region1_start_aa}-{region1_end_aa})...")
print()

candidates = []

for vp1_aa_pos in range(region1_start_aa, region1_end_aa):
    # Convert to nucleotide position (0-indexed)
    nt_pos = vp1_start + (vp1_aa_pos * 3)

    # Check 6bp window
    for offset in [0, 1, 2]:  # Check all 3 reading frame positions
        start = nt_pos + offset
        if start + 6 > len(repcap.seq):
            continue

        current = str(repcap.seq[start:start+6])

        # Calculate Hamming distance
        hamming = sum(1 for a, b in zip(current, target_seq) if a != b)

        if hamming <= 2:  # Only consider if 1-2 mutations needed
            # Figure out which codons are affected
            # Need to check if mutations are silent

            # Determine codon boundaries relative to this window
            reading_frame_offset = (start - vp1_start) % 3

            # Build mutation info
            mutations = []
            is_silent = True

            for i, (orig, new) in enumerate(zip(current, target_seq)):
                if orig != new:
                    abs_pos = start + i

                    # Find which codon this position is in
                    codon_start = vp1_start + ((abs_pos - vp1_start) // 3) * 3
                    codon_pos_in_codon = (abs_pos - codon_start)

                    # Get original and new codons
                    orig_codon = str(repcap.seq[codon_start:codon_start+3])
                    new_codon_list = list(orig_codon)
                    new_codon_list[codon_pos_in_codon] = new
                    new_codon = ''.join(new_codon_list)

                    # Get amino acids
                    orig_aa = codon_table.get(orig_codon, '?')
                    new_aa = codon_table.get(new_codon, '?')

                    vp1_aa = (codon_start - vp1_start) // 3 + 1

                    mutations.append({
                        'abs_pos': abs_pos + 1,  # 1-indexed
                        'vp1_aa': vp1_aa,
                        'orig': orig,
                        'new': new,
                        'orig_codon': orig_codon,
                        'new_codon': new_codon,
                        'orig_aa': orig_aa,
                        'new_aa': new_aa,
                        'silent': orig_aa == new_aa
                    })

                    if orig_aa != new_aa:
                        is_silent = False

            if is_silent and hamming > 0:
                candidates.append({
                    'vp1_aa_region': vp1_aa_pos,
                    'start_pos': start + 1,  # 1-indexed
                    'current': current,
                    'target': target_seq,
                    'hamming': hamming,
                    'mutations': mutations
                })

# Sort by position
candidates.sort(key=lambda x: x['start_pos'])

print(f"Found {len(candidates)} candidate positions for AvrII with silent mutations:")
print()

for cand in candidates:
    print(f"Position {cand['start_pos']} (VP1 aa region ~{cand['vp1_aa_region']})")
    print(f"  Current: {cand['current']} → Target: {cand['target']}")
    print(f"  Mutations needed: {cand['hamming']}")
    for mut in cand['mutations']:
        print(f"    - Pos {mut['abs_pos']} (VP1 aa {mut['vp1_aa']}): "
              f"{mut['orig']}→{mut['new']} | "
              f"{mut['orig_codon']}→{mut['new_codon']} ({mut['orig_aa']}→{mut['new_aa']}) "
              f"{'✅ SILENT' if mut['silent'] else '❌ NOT SILENT'}")
    print()

# Pick the best one (prefer 1 mutation, earlier position)
if candidates:
    best = min(candidates, key=lambda x: (x['hamming'], x['start_pos']))
    print("=" * 70)
    print("RECOMMENDED CHOICE:")
    print(f"  Position: {best['start_pos']} (VP1 aa region ~{best['vp1_aa_region']})")
    print(f"  Mutations: {best['hamming']}")
    for mut in best['mutations']:
        print(f"    - Pos {mut['abs_pos']} (VP1 aa {mut['vp1_aa']}): "
              f"{mut['orig']}→{mut['new']} | "
              f"{mut['orig_codon']}→{mut['new_codon']} ({mut['orig_aa']}→{mut['new_aa']})")
else:
    print("No candidates found!")
