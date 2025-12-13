import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
import os
import gzip
import csv
import sys

# --- Configuration ---
# We use the raw source file to get exact boundaries
INPUT_FILE = "output_sequences_full.tsv.gz"
INPUT_LABELS = "dataset_y_labels.tsv" # To filter for positives

# Increase CSV limit for large fields
csv.field_size_limit(sys.maxsize)

def get_donor_junction_from_source(seq_pseudo, seq_intron2):
    """
    Extracts the junction directly from the separate components.
    Target: 10bp Exon (End of Pseudo) + 10bp Intron (Start of Intron 2)
    """
    # Exon side (Left of junction)
    if len(seq_pseudo) < 10:
        return None
    exon_part = seq_pseudo[-10:]
    
    # Intron side (Right of junction)
    if len(seq_intron2) < 10:
        return None
    intron_part = seq_intron2[:10]
    
    # Concatenate
    full_seq = exon_part + intron_part
    return full_seq.upper() # Force uppercase to avoid logo artifacts

def main():
    print("--- Generating Ground Truth Risdiplam Motif ---")
    
    # 1. Load Labels to identify Positives
    print("Loading labels...")
    try:
        labels_df = pd.read_csv(INPUT_LABELS, sep='\t')
        # Get set of positive IDs
        positive_ids = set(labels_df[labels_df['label'] == 1]['event_id'])
        print(f"Found {len(positive_ids)} positive events.")
    except FileNotFoundError:
        print(f"Error: {INPUT_LABELS} not found.")
        return

    # 2. Stream Source File and Extract
    print(f"Streaming {INPUT_FILE}...")
    sequences = []
    
    try:
        with gzip.open(INPUT_FILE, 'rt') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            for row in reader:
                eid = row['event_id']
                
                if eid in positive_ids:
                    # Extract directly from columns
                    seq = get_donor_junction_from_source(
                        row['Seq3_Pseudoexon'], 
                        row['Seq4_Intron2']
                    )
                    
                    if seq and len(seq) == 20:
                        sequences.append(seq)
                        
    except FileNotFoundError:
        print(f"Error: {INPUT_FILE} not found.")
        return

    print(f"Extracted {len(sequences)} valid sequences aligned at the donor.")

    # 3. Create Matrix
    print("Calculating frequencies...")
    counts_df = logomaker.alignment_to_matrix(sequences)
    
    # 4. Plot
    print("Plotting...")
    plt.figure(figsize=(12, 5))
    
    logo = logomaker.Logo(counts_df, shade_below=.5, fade_below=.5, color_scheme='classic')
    
    logo.style_spines(visible=False)
    logo.style_spines(spines=['bottom', 'left'], visible=True)
    
    # Set X-axis labels (-10 to +10)
    logo.ax.set_xticks(range(20))
    logo.ax.set_xticklabels(np.arange(-10, 10))
    logo.ax.set_xlabel("Position relative to Splice Donor (bp)")
    logo.ax.set_ylabel("Information Content (Bits)")
    logo.ax.set_title(f"True 5' Splice Site Motif of Risdiplam Targets (n={len(sequences)})")
    
    # Highlight the Splice Line (Between -1 and 0, which is Index 9.5)
    plt.axvline(9.5, color='black', linestyle='--', linewidth=1.5, label='Exon/Intron Boundary')
    
    save_path = "risdiplam_motif_logo_truth.png"
    plt.savefig(save_path)
    print(f"Saved logo to {save_path}")

if __name__ == "__main__":
    main()