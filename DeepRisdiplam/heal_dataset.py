import pandas as pd
import gzip
import csv
import sys
import os

# --- Configuration ---
INPUT_FILE = "output_sequences_full.tsv.gz"
OUTPUT_FILE = "output_sequences_healed.tsv.gz"

csv.field_size_limit(sys.maxsize)

def heal_junction(seq_pseudo, seq_intron2):
    """
    Attempts to fix the boundary between Pseudoexon and Intron 2.
    Goal: Make Intron 2 start with 'GT'.
    """
    # Combine them to see the full context
    # We take the last 5bp of exon and first 5bp of intron
    if len(seq_pseudo) < 5 or len(seq_intron2) < 5:
        return seq_pseudo, seq_intron2, "TOO_SHORT"
        
    junction = (seq_pseudo[-5:] + seq_intron2[:5]).upper()
    # Canonical: ...AG|GT... (Index 5 is start of intron)
    
    # Case 0: Already Good (Starts with GT)
    if seq_intron2[:2].upper() == "GT":
        return seq_pseudo, seq_intron2, "OK"
        
    # Case 1: Shifted Left (The 'G' is in the Exon)
    # Pattern: ...G|T... (Intron starts with T)
    if seq_pseudo[-1].upper() == 'G' and seq_intron2[0].upper() == 'T':
        # Move 'G' from Exon to Intron
        new_pseudo = seq_pseudo[:-1]
        new_intron2 = seq_pseudo[-1] + seq_intron2
        return new_pseudo, new_intron2, "HEALED_LEFT"
        
    # Case 2: Shifted Right (The 'A' or 'G' of exon is in Intron)
    # Pattern: ...|GGT... or ...|AGT...
    # This is rare but possible if 1-based coord was exclusive
    if seq_intron2[1:3].upper() == "GT":
        # Move first char of Intron to Exon
        new_pseudo = seq_pseudo + seq_intron2[0]
        new_intron2 = seq_intron2[1:]
        return new_pseudo, new_intron2, "HEALED_RIGHT"

    return seq_pseudo, seq_intron2, "UNFIXABLE"

def main():
    print("--- Auto-Healing Dataset ---")
    
    stats = {"OK": 0, "HEALED_LEFT": 0, "HEALED_RIGHT": 0, "UNFIXABLE": 0, "TOO_SHORT": 0}
    
    with gzip.open(INPUT_FILE, 'rt') as f_in, gzip.open(OUTPUT_FILE, 'wt') as f_out:
        reader = csv.DictReader(f_in, delimiter='\t')
        
        # Write Header
        writer = csv.DictWriter(f_out, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()
        
        count = 0
        for row in reader:
            # fix Pseudo -> Intron 2
            p, i2, status = heal_junction(row['Seq3_Pseudoexon'], row['Seq4_Intron2'])
            
            row['Seq3_Pseudoexon'] = p
            row['Seq4_Intron2'] = i2
            
            stats[status] += 1
            writer.writerow(row)
            
            count += 1
            if count % 50000 == 0:
                print(f"Processed {count}...", end='\r')
                
    print(f"\nDone! Processed {count} sequences.")
    print("-" * 30)
    print("Healing Statistics:")
    for k, v in stats.items():
        print(f"  {k}: {v} ({v/count*100:.1f}%)")
    print("-" * 30)
    print(f"Saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()