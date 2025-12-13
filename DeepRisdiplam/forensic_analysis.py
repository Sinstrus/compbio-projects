import torch
import torch.nn as nn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import train_test_split
import sys
import os

# --- CONFIGURATION ---
VERSION = "v11_hybrid_franken"
MODEL_PATH = f"deep_risdiplam_model_{VERSION}.pth"
INPUT_FILE = "dataset_components_v8.tsv"
MAX_SEQ_LEN = 1500
BATCH_SIZE = 32
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# --- RECONSTRUCT MODEL (V11 Architecture) ---
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

def get_gc_content(seq):
    if len(seq) == 0: return 0
    g = seq.count('G') + seq.count('g')
    c = seq.count('C') + seq.count('c')
    return (g + c) / len(seq)

def stitch(row):
    return row['ExonA'] + row['Intron1'] + row['Pseudoexon'] + row['Intron2'] + row['ExonB']

# --- MAIN ANALYSIS ---
def main():
    print(f"--- Forensics on {VERSION} ---")
    
    # 1. Load Data & Model
    df = pd.read_csv(INPUT_FILE, sep='\t')
    df = df.dropna(subset=['ExonA', 'Intron1', 'Pseudoexon', 'Intron2', 'ExonB'])
    
    # Re-create the validation split EXACTLY as used in training
    _, val_df = train_test_split(df, test_size=0.1, random_state=42, stratify=df['label'])
    
    positives = val_df[val_df['label'] == 1]
    print(f"Total Positives in Validation Set: {len(positives)}")
    
    model = RisdiplamHybrid().to(device)
    model.load_state_dict(torch.load(MODEL_PATH, map_location=device))
    model.eval()
    
    # 2. Run Inference on Positives Only
    results = []
    
    print("Scanning Validation Positives...")
    with torch.no_grad():
        for _, row in positives.iterrows():
            full_seq = stitch(row)
            tensor = vectorized_one_hot_encode([full_seq]).to(device)
            logit = model(tensor).item()
            prob = torch.sigmoid(torch.tensor(logit)).item()
            
            # Classification
            pred = 1 if prob > 0.5 else 0
            status = "TP" if pred == 1 else "FN"
            
            results.append({
                "event_id": row['event_id'],
                "status": status,
                "prob": prob,
                "len_total": len(full_seq),
                "len_pseudo": len(row['Pseudoexon']),
                "len_i1": len(row['Intron1']),
                "len_i2": len(row['Intron2']),
                "gc_total": get_gc_content(full_seq),
                "gc_pseudo": get_gc_content(row['Pseudoexon']),
                # Check for Canonical Donor (GT) at start of Intron 2
                "has_gt_donor": row['Intron2'].upper().startswith("GT")
            })
            
    res_df = pd.DataFrame(results)
    
    # 3. Report
    print("\n--- CONFUSION REPORT ---")
    print(res_df['status'].value_counts())
    
    tps = res_df[res_df['status'] == "TP"]
    fns = res_df[res_df['status'] == "FN"]
    
    print(f"\nCaught (True Positives): {len(tps)}")
    print(f"Missed (False Negatives): {len(fns)}")
    
    # 4. Statistical Comparisons
    print("\n--- FEATURE COMPARISON (TP vs FN) ---")
    features = ['len_total', 'len_pseudo', 'gc_total', 'gc_pseudo', 'prob']
    
    print(f"{'Feature':<15} | {'TP Mean':<10} | {'FN Mean':<10} | {'Diff'}")
    print("-" * 50)
    for f in features:
        tp_m = tps[f].mean()
        fn_m = fns[f].mean()
        print(f"{f:<15} | {tp_m:<10.4f} | {fn_m:<10.4f} | {tp_m - fn_m:.4f}")

    # Check Donor Status
    print("\n--- DONOR SITE CHECK ---")
    print(f"TPs with GT Donor: {tps['has_gt_donor'].sum()}/{len(tps)}")
    print(f"FNs with GT Donor: {fns['has_gt_donor'].sum()}/{len(fns)}")

    # 5. Save the List of Hard Positives
    fns.to_csv("analysis_hard_positives_failures.csv", index=False)
    print("\nSaved list of failed events to 'analysis_hard_positives_failures.csv'.")
    print("ACTION: Open this CSV. Look at the 'event_id'. Check if they share a gene name or chromosomal region.")

    # 6. Visualization
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    sns.boxplot(x='status', y='len_pseudo', data=res_df)
    plt.title("Pseudoexon Length: TP vs FN")
    
    plt.subplot(1, 2, 2)
    sns.boxplot(x='status', y='gc_pseudo', data=res_df)
    plt.title("Pseudoexon GC%: TP vs FN")
    
    plt.savefig("forensic_analysis.png")
    print("Saved plots to 'forensic_analysis.png'")

if __name__ == "__main__":
    main()