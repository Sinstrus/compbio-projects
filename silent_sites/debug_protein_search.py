#!/usr/bin/env python3
"""
Debug script to find protein sequence in DNA
"""

import sys

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

def translate(dna_seq, frame=0, include_stops=False):
    """Translate DNA to protein"""
    dna_seq = dna_seq.upper()[frame:]
    protein = []
    
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        if len(codon) == 3:
            aa = CODON_TABLE.get(codon, 'X')
            if aa == '*':
                if include_stops:
                    protein.append('*')
                else:
                    break
            else:
                protein.append(aa)
    
    return ''.join(protein)

def find_protein_matches(dna, protein, min_match_length=10):
    """Find partial matches of protein in all frames"""
    protein = protein.upper()
    dna = dna.upper()
    
    print(f"\n{'='*80}")
    print(f"Searching for protein: {protein}")
    print(f"Protein length: {len(protein)} aa")
    print(f"DNA length: {len(dna)} bp")
    print(f"{'='*80}\n")
    
    # Try all 3 frames
    for frame in [0, 1, 2]:
        print(f"\n--- Frame {frame} ---")
        translated = translate(dna, frame, include_stops=True)
        print(f"Full translation length: {len(translated)} aa")
        
        # Check for exact match
        if protein in translated:
            pos = translated.index(protein)
            dna_start = frame + (pos * 3)
            dna_end = dna_start + (len(protein) * 3)
            print(f"✓ EXACT MATCH FOUND!")
            print(f"  Protein position: {pos+1}-{pos+len(protein)}")
            print(f"  DNA position: {dna_start+1}-{dna_end}")
            return
        
        # Look for partial matches
        best_match_len = 0
        best_match_pos = -1
        
        # Check all positions
        for i in range(len(translated) - min_match_length + 1):
            for length in range(min_match_length, len(protein) + 1):
                if i + length > len(translated):
                    break
                
                window = translated[i:i+length]
                if window in protein:
                    if length > best_match_len:
                        best_match_len = length
                        best_match_pos = i
        
        if best_match_len > 0:
            print(f"✓ Partial match found: {best_match_len} aa")
            print(f"  Position in translation: {best_match_pos+1}")
            dna_pos = frame + (best_match_pos * 3)
            print(f"  DNA position: {dna_pos+1}")
            match_seq = translated[best_match_pos:best_match_pos+best_match_len]
            print(f"  Matched sequence: {match_seq}")
            protein_pos = protein.index(match_seq)
            print(f"  Position in target protein: {protein_pos+1}-{protein_pos+best_match_len}")
        else:
            print(f"✗ No partial matches found (min length: {min_match_length})")
        
        # Show first 100 aa of translation
        preview = translated[:100]
        print(f"\n  Translation preview (first 100 aa):")
        for i in range(0, len(preview), 60):
            print(f"  {preview[i:i+60]}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python debug_protein_search.py <DNA> <PROTEIN>")
        sys.exit(1)
    
    dna = sys.argv[1]
    protein = sys.argv[2]
    
    find_protein_matches(dna, protein)
