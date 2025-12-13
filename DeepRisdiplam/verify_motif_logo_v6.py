import torch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
import os
import torch.nn as nn

# --- Configuration ---
VERSION = "v6"
INPUT_X_FILE = "dataset_X_sequences_strict.tsv"
INPUT_Y_FILE = "dataset_y_labels_strict.tsv"
OFFSET_END_FROM_RIGHT = 550 

def get_donor_window(seq):
    idx_donor = len(seq) - OFFSET_END_FROM_RIGHT
    start = idx_donor - 10
    end = idx_donor + 10
    return seq[start:end].upper()

def main():
    print("--- Verifying V6 Motif (Strict GT Only) ---")
    
    try:
        df_x = pd.read_csv(INPUT_X_FILE, sep='\t')
        df_y = pd.read_csv(INPUT_Y_FILE, sep='\t')
        df = pd.merge(df_x, df_y, on="event_id")
    except FileNotFoundError:
        print("Data files not found.")
        return

    positives = df[df['label'] == 1]
    sequences = []
    
    for seq in positives['minigene_sequence']:
        window = get_donor_window(seq)
        if len(window) == 20:
            sequences.append(window)
            
    print(f"Aligning {len(sequences)} sequences...")
    counts_df = logomaker.alignment_to_matrix(sequences)
    
    plt.figure(figsize=(10, 4))
    logo = logomaker.Logo(counts_df, shade_below=.5, fade_below=.5, color_scheme='classic')
    logo.style_spines(visible=False)
    logo.style_spines(spines=['bottom'], visible=True)
    logo.ax.set_xticks(range(20))
    logo.ax.set_xticklabels(np.arange(-10, 10))
    logo.ax.set_title(f"Strict V6 Motif (n={len(sequences)})\n(Only confirmed GT start sites)")
    
    plt.axvline(9.5, color='black', linestyle='--', linewidth=1)
    plt.savefig(f"risdiplam_motif_logo_{VERSION}.png")
    print(f"Saved logo to risdiplam_motif_logo_{VERSION}.png")

if __name__ == "__main__":
    main()