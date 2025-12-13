import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
import os
import gzip
import csv
import sys

# --- Configuration ---
INPUT_FILE = "output_sequences_full.tsv.gz"
INPUT_LABELS = "dataset_y_labels.tsv" 

csv.field_size_limit(sys.maxsize)

def get_donor_junction_corrected(seq_pseudo, seq_intron2):
    """
    Extracts the junction with a 1bp correction for 0-based coordinates.
    The 'End' of the pseudoexon column actually contains the first base ('G') of the intron.
    """
    # 1. Reconstruct the raw sequence around the cut
    # We take the last 11 bases of the 'Exon' column and first 9 of 'Intron'
    # Why? Because the Exon column stole the first base of the intron.
    
    if len(seq_pseudo) < 11 or len(seq_intron2) < 9:
        return None
        
    # Exon Part: Take -11 to -1 (Drop the last char, which is the Intron's 'G')
    real_exon_end = seq_pseudo[-11:-1]
    
    # Intron Part: Take the last char of 'Exon' + start of 'Intron'
    stolen_g = seq_pseudo[-1]
    real_intron_start = stolen_g + seq_intron2[:9]
    
    # Concatenate: 10bp Exon + 10bp Intron
    full_seq = real_exon_end + real_intron_start
    return full_seq.upper()

def main():
    print("--- Generating Corrected Risdiplam Motif ---")
    
    try:
        labels_df = pd.read_csv(INPUT_LABELS, sep='\t')
        positive_ids = set(labels_df[labels_df['label'] == 1]['event_id'])
    except FileNotFoundError:
        print(f"Error: {INPUT_LABELS} not found.")
        return

    sequences = []
    
    try:
        with gzip.open(INPUT_FILE, 'rt') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row['event_id'] in positive_ids:
                    # Use the CORRECTED extractor
                    seq = get_donor_junction_corrected(
                        row['Seq3_Pseudoexon'], 
                        row['Seq4_Intron2']
                    )
                    if seq and len(seq) == 20:
                        sequences.append(seq)
    except FileNotFoundError:
        print(f"Error: {INPUT_FILE} not found.")
        return

    print(f"Aligned {len(sequences)} sequences (Corrected for 0-based shift).")

    # Plot
    counts_df = logomaker.alignment_to_matrix(sequences)
    
    plt.figure(figsize=(12, 5))
    logo = logomaker.Logo(counts_df, shade_below=.5, fade_below=.5, color_scheme='classic')
    logo.style_spines(visible=False)
    logo.style_spines(spines=['bottom', 'left'], visible=True)
    logo.ax.set_xticks(range(20))
    logo.ax.set_xticklabels(np.arange(-10, 10))
    logo.ax.set_title(f"Corrected 5' Splice Site Motif (n={len(sequences)})\n(Shifted 1bp to fix coordinate error)")
    plt.axvline(9.5, color='black', linestyle='--', linewidth=1.5)
    
    plt.savefig("risdiplam_motif_logo_corrected.png")
    print("Saved logo to risdiplam_motif_logo_corrected.png")

if __name__ == "__main__":
    main()