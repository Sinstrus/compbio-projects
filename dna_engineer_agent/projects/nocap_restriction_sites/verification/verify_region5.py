#!/usr/bin/env python3
"""
Verify Region 5 BmtI site is actually silent
"""

from Bio import SeqIO

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

# Load plasmid
record = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")
plasmid_seq = str(record.seq).upper()

# BmtI site details
site_position = 3937  # 1-indexed
recognition = "GCTAGC"
mutation = "C3T"  # Position 3 in recognition sequence

# Calculate mutation position
mut_offset = int(mutation[1]) - 1  # 0-indexed
abs_position = site_position + mut_offset  # 1-indexed

print("="*80)
print("VERIFYING REGION 5: BmtI @ 3937")
print("="*80)
print(f"\nRecognition sequence: {recognition}")
print(f"Mutation: {mutation} (position {mut_offset + 1} in site)")
print(f"Absolute mutation position: {abs_position}")

# Get sequence context
original_site = plasmid_seq[site_position-1:site_position-1+len(recognition)]
print(f"\nOriginal site sequence: {original_site}")

# Apply mutation
modified_site = list(original_site)
modified_site[mut_offset] = 'T'
modified_site = ''.join(modified_site)
print(f"Modified site sequence: {modified_site}")

# Region 5 protein starts at 3880 according to analysis
protein_start = 3880

# Calculate which codon contains position abs_position
offset_from_protein_start = abs_position - protein_start  # In 1-indexed coords
print(f"\nOffset from protein start ({protein_start}): {offset_from_protein_start}")

# Find codon
codon_number = offset_from_protein_start // 3
position_in_codon = (offset_from_protein_start - 1) % 3  # 0-indexed
codon_start = protein_start + (codon_number * 3) - 1  # Adjust for 0-indexed

print(f"Codon number: {codon_number + 1}")
print(f"Position in codon: {position_in_codon + 1} (of 3)")
print(f"Codon starts at position: {codon_start + 1} (1-indexed)")

# Get original codon
original_codon = plasmid_seq[codon_start:codon_start+3]
print(f"\nOriginal codon: {original_codon}")
print(f"Original amino acid: {CODON_TABLE.get(original_codon, '?')}")

# Get mutated codon
mutated_codon = list(original_codon)
pos_in_codon_for_mutation = (abs_position - 1) - codon_start  # 0-indexed
mutated_codon[pos_in_codon_for_mutation] = 'T'
mutated_codon = ''.join(mutated_codon)
print(f"\nMutated codon: {mutated_codon}")
print(f"Mutated amino acid: {CODON_TABLE.get(mutated_codon, '?')}")

# Check if silent
original_aa = CODON_TABLE.get(original_codon, '?')
mutated_aa = CODON_TABLE.get(mutated_codon, '?')

print(f"\n{'='*80}")
if original_aa == mutated_aa:
    print(f"✓ SILENT MUTATION: {original_codon} → {mutated_codon} ({original_aa} → {mutated_aa})")
    print("This site is safe to use.")
else:
    print(f"✗ NOT SILENT: {original_codon} → {mutated_codon} ({original_aa} → {mutated_aa})")
    print("This mutation changes the amino acid!")
    print("\nSuggested action: Select an alternative site for Region 5")

print(f"{'='*80}")
