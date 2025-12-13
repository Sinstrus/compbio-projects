import torch
import torch.nn as nn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
from sklearn.model_selection import train_test_split
import sys
import os

# --- CONFIGURATION ---
VERSION = "v11_hybrid_franken"
MODEL_PATH = f"deep_risdiplam_model_{VERSION}.pth"
INPUT_FILE = "dataset_components_v8.tsv"
MAX_SEQ_LEN = 1500
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# --- MODEL (V11) ---
class RisdiplamHybrid(nn.Module):
    def __init__(self):
        super(RisdiplamHybrid, self).__init__()
        self.conv1 = nn.Conv1d(4, 64, 11, padding='same')
        self.bn1 = nn.BatchNorm1d(64)
        self.relu = nn.ReLU()
        self.pool = nn.MaxPool1d(2) 
        self.conv2 = nn.Conv1d(64, 128, 11, padding='same', dilation=2)
        self.bn2 = nn.BatchNorm1d(128)
        self.d_model = 128
        self.nhead = 4
        self.attention = nn.MultiheadAttention(embed_dim=self.d_model, num_heads=self.nhead, batch_first=True)
        self.global_pool = nn.AdaptiveAvgPool1d(1)
        self.flatten = nn.Flatten()
        self.fc1 = nn.Linear(128, 64)
        self.dropout = nn.Dropout(0.65) 
        self.fc2 = nn.Linear(64, 1)

    def _attn_block(self, x):
        attn_out, _ = self.attention(x, x, x)
        return x + attn_out

    def forward(self, x):
        x = self.relu(self.bn1(self.conv1(x)))
        x = self.pool(x)
        x = self.relu(self.bn2(self.conv2(x)))
        x = x.permute(0, 2, 1) 
        x = self._attn_block(x)
        x = x.permute(0, 2, 1)
        x = self.global_pool(x)
        x = self.flatten(x)
        x = self.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.fc2(x)
        return x

# --- HELPER FUNCTIONS ---
def vectorized_one_hot_encode(seqs, max_len=MAX_SEQ_LEN):
    batch_size = len(seqs)
    mapping = np.zeros(128, dtype=np.int8) + 4 
    mapping[ord('A')] = 0; mapping[ord('a')] = 0
    mapping[ord('C')] = 1; mapping[ord('c')] = 1
    mapping[ord('G')] = 2; mapping[ord('g')] = 2
    mapping[ord('T')] = 3; mapping[ord('t')] = 3
    arr = np.zeros((batch_size, max_len), dtype=np.int8) + 4
    for i, seq in enumerate(seqs):
        l = min(len(seq), max_len)
        arr[i, :l] = mapping[np.frombuffer(seq[:l].encode('latin1'), dtype=np.uint8)]
    one_hot = np.zeros((batch_size, 4, max_len), dtype=np.float32)
    for base_val in range(4):
        one_hot[:, base_val, :] = (arr == base_val)
    return torch.tensor(one_hot)

def stitch(row):
    return row['ExonA'] + row['Intron1'] + row['Pseudoexon'] + row['Intron2'] + row['ExonB']

def extract_donor_window(row, window_size=6):
    """
    Extracts the window around the splice donor.
    The Donor is at the start of Intron 2.
    We want the last 3 bases of Pseudoexon and first 6 of Intron 2.
    """
    pseudo = row['Pseudoexon']
    intron2 = row['Intron2']
    
    # -3 from Exon, +6 from Intron (GT is at +0, +1)
    context_up = pseudo[-3:] if len(pseudo) >= 3 else "N" * 3
    context_down = intron2[:6] if len(intron2) >= 6 else "N" * 6
    
    return context_up + context_down

def create_logo_matrix(sequences):
    """
    Converts a list of strings into a probability matrix for Logomaker.
    """
    if not sequences: return None
    
    # Filter out length mismatches
    length = len(sequences[0])
    valid_seqs = [s for s in sequences if len(s) == length]
    
    counts = np.zeros((length, 4)) # A, C, G, T
    map_base = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3}
    
    for seq in valid_seqs:
        for i, char in enumerate(seq):
            if char in map_base:
                counts[i, map_base[char]] += 1
                
    # Normalize
    sums = counts.sum(axis=1, keepdims=True)
    sums[sums==0] = 1
    probs = counts / sums
    return pd.DataFrame(probs, columns=['A', 'C', 'G', 'T'])

# --- MAIN ---
def main():
    print("--- Analyzing Motifs of Hits vs Misses ---")
    
    df = pd.read_csv(INPUT_FILE, sep='\t')
    df = df.dropna(subset=['ExonA', 'Intron1', 'Pseudoexon', 'Intron2', 'ExonB'])
    _, val_df = train_test_split(df, test_size=0.1, random_state=42, stratify=df['label'])
    
    positives = val_df[val_df['label'] == 1]
    
    model = RisdiplamHybrid().to(device)
    model.load_state_dict(torch.load(MODEL_PATH, map_location=device))
    model.eval()
    
    tp_seqs = []
    fn_seqs = []
    
    print("Classifying Validation Positives...")
    with torch.no_grad():
        for _, row in positives.iterrows():
            full_seq = stitch(row)
            tensor = vectorized_one_hot_encode([full_seq]).to(device)
            prob = torch.sigmoid(model(tensor)).item()
            
            donor_seq = extract_donor_window(row)
            
            if prob > 0.5:
                tp_seqs.append(donor_seq)
            else:
                fn_seqs.append(donor_seq)
                
    print(f"Hits (TP): {len(tp_seqs)}")
    print(f"Misses (FN): {len(fn_seqs)}")
    
    # Plotting
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    
    # Logo 1: The Easy Ones
    df_tp = create_logo_matrix(tp_seqs)
    if df_tp is not None:
        logo_tp = logomaker.Logo(df_tp, ax=axes[0], shade_below=.5, fade_below=.5)
        axes[0].set_title(f"Motif of CAUGHT Positives (n={len(tp_seqs)})")
        axes[0].set_ylabel("Probability")
        # Mark the splice site (between index 2 and 3)
        axes[0].axvline(2.5, color='red', linestyle='--', linewidth=2)
    
    # Logo 2: The Hard Ones
    df_fn = create_logo_matrix(fn_seqs)
    if df_fn is not None:
        logo_fn = logomaker.Logo(df_fn, ax=axes[1], shade_below=.5, fade_below=.5)
        axes[1].set_title(f"Motif of MISSED Positives (n={len(fn_seqs)})")
        axes[1].set_ylabel("Probability")
        axes[1].axvline(2.5, color='red', linestyle='--', linewidth=2)
        
    plt.tight_layout()
    plt.savefig("motif_comparison.png")
    print("\nSaved 'motif_comparison.png'")
    print("Interpretation:")
    print("Indices 0-2: End of Exon (-3, -2, -1)")
    print("Indices 3-8: Start of Intron (+1, +2 ...)")
    print("Look at the +1 (Index 3) and +2 (Index 4). Are they both GT?")
    print("Look at -1 (Index 2). Risdiplam likes a 'G' here.")

if __name__ == "__main__":
    main()