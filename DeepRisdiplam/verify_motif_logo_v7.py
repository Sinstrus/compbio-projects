import torch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
import os
import torch.nn as nn

# --- Configuration ---
VERSION = "v7_components"
INPUT_FILE = "dataset_components.tsv"

def main():
    print(f"--- Verifying V7 Motif (Component-Based) ---")
    
    try:
        df = pd.read_csv(INPUT_FILE, sep='\t')
    except FileNotFoundError:
        print("Data file not found.")
        return

    positives = df[df['label'] == 1]
    sequences = []
    
    print(f"Aligning {len(positives)} sequences...")
    
    for idx, row in positives.iterrows():
        # Get the junction directly from components
        # Pseudo End (last 10) + Intron 2 Start (first 10)
        
        ps = row['Pseudoexon']
        i2 = row['Intron2']
        
        if len(ps) < 10 or len(i2) < 10: continue
            
        junction = ps[-10:] + i2[:10]
        sequences.append(junction.upper())
            
    # Plot
    counts_df = logomaker.alignment_to_matrix(sequences)
    
    plt.figure(figsize=(10, 4))
    logo = logomaker.Logo(counts_df, shade_below=.5, fade_below=.5, color_scheme='classic')
    logo.style_spines(visible=False)
    logo.style_spines(spines=['bottom'], visible=True)
    logo.ax.set_xticks(range(20))
    logo.ax.set_xticklabels(np.arange(-10, 10))
    logo.ax.set_title(f"V7 Component Motif (n={len(sequences)})\n(Guaranteed Alignment)")
    
    plt.axvline(9.5, color='black', linestyle='--', linewidth=1)
    
    save_path = f"risdiplam_motif_logo_{VERSION}.png"
    plt.savefig(save_path)
    print(f"Saved logo to {save_path}")

if __name__ == "__main__":
    main()