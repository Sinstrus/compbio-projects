import torch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
import os

# --- Configuration ---
VERSION = "v4"
INPUT_X_FILE = "dataset_X_sequences.tsv"
INPUT_Y_FILE = "dataset_y_labels.tsv"

# Constants
INTRON_1_PART_LEN = 250
INTRON_2_PART_LEN = 250
EXON_B_LEN = 50
OFFSET_END_FROM_RIGHT = INTRON_2_PART_LEN * 2 + EXON_B_LEN 

def get_donor_window(seq):
    """Extracts -10 to +10 around the Pseudoexon-Intron2 Junction"""
    idx_donor = len(seq) - OFFSET_END_FROM_RIGHT
    start = idx_donor - 10
    end = idx_donor + 10
    
    # Extract and FORCE UPPERCASE to fix the visualization artifact
    return seq[start:end].upper()

def main():
    print("--- Generating Risdiplam Motif Logo (Fixed Case) ---")
    
    # 1. Load Data
    try:
        df_x = pd.read_csv(INPUT_X_FILE, sep='\t')
        df_y = pd.read_csv(INPUT_Y_FILE, sep='\t')
        df = pd.merge(df_x, df_y, on="event_id")
    except FileNotFoundError:
        print("Data files not found.")
        return

    positives = df[df['label'] == 1]
    print(f"Aligning {len(positives)} positive sequences...")

    # 2. Extract Windows
    sequences = []
    for seq in positives['minigene_sequence']:
        window = get_donor_window(seq)
        if len(window) == 20:
            sequences.append(window)
            
    if not sequences:
        print("Error extracting sequences.")
        return

    # 3. Create Matrix
    print("Calculating frequencies...")
    counts_df = logomaker.alignment_to_matrix(sequences)
    
    # 4. Plot
    print("Plotting...")
    plt.figure(figsize=(10, 4))
    
    logo = logomaker.Logo(counts_df, shade_below=.5, fade_below=.5)
    
    logo.style_spines(visible=False)
    logo.style_spines(spines=['bottom'], visible=True)
    logo.ax.set_xticks(range(20))
    logo.ax.set_xticklabels(np.arange(-10, 10))
    logo.ax.set_ylabel("Information Content")
    logo.ax.set_title(f"Learned 5' Splice Site Motif (Risdiplam Responsive)\n(n={len(sequences)})")
    
    # Highlight the GT position (Index 10/11)
    plt.axvline(9.5, color='black', linestyle='--', linewidth=1)
    
    save_path = f"risdiplam_motif_logo_{VERSION}_fixed.png"
    plt.savefig(save_path)
    print(f"Saved logo to {save_path}")

if __name__ == "__main__":
    main()