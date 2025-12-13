import torch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
import os
import torch.nn as nn

# --- Configuration ---
VERSION = "v5"
MODEL_PATH = f"deep_risdiplam_model_{VERSION}.pth"
INPUT_X_FILE = "dataset_X_sequences_healed.tsv" # Use Healed
INPUT_Y_FILE = "dataset_y_labels_healed.tsv"

OFFSET_END_FROM_RIGHT = 550 # Intron 2 (500) + Exon B (50)
MAX_SEQ_LEN = 1500
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

class RisdiplamModel(nn.Module):
    def __init__(self):
        super(RisdiplamModel, self).__init__()
        self.conv1 = nn.Conv1d(4, 64, 11, padding='same')
        self.bn1 = nn.BatchNorm1d(64)
        self.relu = nn.ReLU()
        self.pool = nn.MaxPool1d(2)
        self.conv2 = nn.Conv1d(64, 128, 11, padding='same', dilation=2)
        self.bn2 = nn.BatchNorm1d(128)
        self.conv3 = nn.Conv1d(128, 128, 11, padding='same', dilation=4)
        self.bn3 = nn.BatchNorm1d(128)
        self.global_pool = nn.AdaptiveAvgPool1d(1)
        self.flatten = nn.Flatten()
        self.fc1 = nn.Linear(128, 64)
        self.dropout = nn.Dropout(0.6)
        self.fc2 = nn.Linear(64, 1) 
    def forward(self, x):
        x = self.pool(self.relu(self.bn1(self.conv1(x))))
        x = self.pool(self.relu(self.bn2(self.conv2(x))))
        x = self.relu(self.bn3(self.conv3(x)))
        x = self.fc2(self.dropout(self.relu(self.fc1(self.flatten(self.global_pool(x))))))
        return x

def get_donor_window(seq):
    # Now that data is healed, the donor should be exactly here
    idx_donor = len(seq) - OFFSET_END_FROM_RIGHT
    start = idx_donor - 10
    end = idx_donor + 10
    return seq[start:end].upper()

def main():
    print("--- Verifying V5 Motif (Healed Data) ---")
    
    # 1. Load Data
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
            
    # 2. Plot
    print(f"Aligning {len(sequences)} sequences...")
    counts_df = logomaker.alignment_to_matrix(sequences)
    
    plt.figure(figsize=(10, 4))
    logo = logomaker.Logo(counts_df, shade_below=.5, fade_below=.5, color_scheme='classic')
    logo.style_spines(visible=False)
    logo.style_spines(spines=['bottom'], visible=True)
    logo.ax.set_xticks(range(20))
    logo.ax.set_xticklabels(np.arange(-10, 10))
    logo.ax.set_title(f"Healed V5 Motif (n={len(sequences)})")
    
    # Expect GT at index 10 and 11
    plt.axvline(9.5, color='black', linestyle='--', linewidth=1)
    
    plt.savefig(f"risdiplam_motif_logo_{VERSION}_healed.png")
    print(f"Saved logo to risdiplam_motif_logo_{VERSION}_healed.png")

if __name__ == "__main__":
    main()