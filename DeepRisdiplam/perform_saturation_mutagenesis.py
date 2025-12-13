import torch
import torch.nn as nn
import pandas as pd
import numpy as np
import csv
import matplotlib.pyplot as plt
from tqdm import tqdm
import copy

# --- Configuration ---
VERSION = "v8_components"
MODEL_PATH = f"deep_risdiplam_model_{VERSION}.pth"
INPUT_FILE = "dataset_components_v8.tsv"
MAX_SEQ_LEN = 3000
BATCH_SIZE = 512
TOP_N = 5
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# --- Model Architecture ---
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
        self.dropout = nn.Dropout(0.5) 
        self.fc2 = nn.Linear(64, 1)

    def forward(self, x):
        x = self.relu(self.bn1(self.conv1(x)))
        x = self.pool(x)
        x = self.relu(self.bn2(self.conv2(x)))
        x = self.relu(self.bn3(self.conv3(x)))
        x = self.global_pool(x)
        x = self.flatten(x)
        x = self.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.fc2(x)
        return x

def one_hot_encode(seq, max_len=MAX_SEQ_LEN):
    mapping = {'A':0, 'a':0, 'C':1, 'c':1, 'G':2, 'g':2, 'T':3, 't':3}
    seq_len = min(len(seq), max_len)
    tensor = np.zeros((4, max_len), dtype=np.float32)
    for i in range(seq_len):
        if seq[i] in mapping:
            tensor[mapping[seq[i]], i] = 1.0
    return torch.tensor(tensor)

def stitch(row):
    parts = [row['ExonA'], row['Intron1'], row['Pseudoexon'], row['Intron2'], row['ExonB']]
    full_seq = "".join(parts)
    
    boundaries = []
    current_pos = 0
    for p in parts:
        current_pos += len(p)
        boundaries.append(current_pos)
        
    return full_seq, boundaries

def get_top_responders(model, filename, n=10):
    print(f"Scanning for top {n} responders...")
    candidates = []
    scan_limit = 2000 
    
    with open(filename, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['label'] == '1':
                candidates.append(row)
            if len(candidates) >= scan_limit: break
    
    scored_candidates = []
    model.eval()
    with torch.no_grad():
        for i in range(0, len(candidates), BATCH_SIZE):
            batch = candidates[i:i+BATCH_SIZE]
            seqs = [stitch(r)[0] for r in batch]
            tensors = torch.stack([one_hot_encode(s) for s in seqs]).to(device)
            preds = torch.sigmoid(model(tensors)).cpu().numpy().flatten()
            
            for j, score in enumerate(preds):
                scored_candidates.append((score, batch[j]))
    
    scored_candidates.sort(key=lambda x: x[0], reverse=True)
    return scored_candidates[:n]

def mutate_sequence(seq):
    mutants = []
    bases = ['A', 'C', 'G', 'T']
    for i in range(len(seq)):
        original = seq[i]
        for b in bases:
            if b.upper() == original.upper(): continue
            mut_seq = seq[:i] + b + seq[i+1:]
            mutants.append((mut_seq, f"Sub {original}{i}{b}", i))
    return mutants

def analyze_mutations(model, top_candidates):
    print(f"\nPerforming End-to-End Mutagenesis on {len(top_candidates)} sequences...")
    
    results = [] 
    boundaries_map = {}
    
    for score, row in top_candidates:
        full_seq, boundaries = stitch(row)
        original_score = score
        event_id = row['event_id']
        boundaries_map[event_id] = boundaries
        
        # --- DEBUG PRINT FOR USER ---
        print(f"\nAnalysis for: {event_id[:30]}...")
        print(f"  ExonA:   {len(row['ExonA'])} bp")
        print(f"  Intron1: {len(row['Intron1'])} bp")
        print(f"  Pseudo:  {len(row['Pseudoexon'])} bp  <-- (Starts at {len(row['ExonA'])+len(row['Intron1'])})")
        print(f"  Intron2: {len(row['Intron2'])} bp")
        print(f"  ExonB:   {len(row['ExonB'])} bp")
        print(f"  Total:   {len(full_seq)} bp")
        print(f"  Base Score: {original_score:.4f}")
        # ----------------------------

        mutants = mutate_sequence(full_seq)
        batch_seqs = []
        batch_meta = []
        
        for m_seq, m_desc, m_pos in mutants:
            batch_seqs.append(one_hot_encode(m_seq))
            batch_meta.append((m_desc, m_pos))
            
            if len(batch_seqs) == BATCH_SIZE:
                t = torch.stack(batch_seqs).to(device)
                with torch.no_grad():
                    new_scores = torch.sigmoid(model(t)).cpu().numpy().flatten()
                
                for k, new_s in enumerate(new_scores):
                    delta = original_score - new_s
                    desc, pos = batch_meta[k]
                    if delta > 0.01:
                        results.append({
                            'id': event_id,
                            'pos': pos,
                            'desc': desc,
                            'delta': delta,
                            'score_before': original_score,
                            'score_after': new_s
                        })
                batch_seqs = []
                batch_meta = []
        
        # Cleanup final batch
        if batch_seqs:
            t = torch.stack(batch_seqs).to(device)
            with torch.no_grad():
                new_scores = torch.sigmoid(model(t)).cpu().numpy().flatten()
            for k, new_s in enumerate(new_scores):
                delta = original_score - new_s
                desc, pos = batch_meta[k]
                if delta > 0.01:
                    results.append({
                        'id': event_id,
                        'pos': pos,
                        'desc': desc,
                        'delta': delta,
                        'score_before': original_score,
                        'score_after': new_s
                    })

    return pd.DataFrame(results), boundaries_map

def plot_heatmap(df, boundaries_map):
    print("\nGenerating End-to-End Sensitivity Plot...")
    
    if df.empty:
        print("No significant mutations found.")
        return

    top_id = df['id'].unique()[0]
    df_single = df[df['id'] == top_id]
    bounds = boundaries_map[top_id]
    
    # Get component names for x-axis
    labels = ["Exon A", "Intron 1", "Pseudo", "Intron 2", "Exon B"]
    
    plt.figure(figsize=(15, 6))
    
    # 1. Plot the Deltas
    plt.scatter(df_single['pos'], df_single['delta'], alpha=0.6, c='crimson', s=15, label='Mutation Impact')
    
    # 2. Draw Regions
    prev_b = 0
    colors = ['#f0f0f0', '#e0e0e0', '#ffeb99', '#e0e0e0', '#f0f0f0'] 
    
    for i, b in enumerate(bounds):
        # Shaded background
        plt.axvspan(prev_b, b, color=colors[i], alpha=0.3, zorder=0)
        # Vertical Line
        plt.axvline(b, color='black', linestyle=':', alpha=0.5)
        # Label
        mid_point = (prev_b + b) / 2
        plt.text(mid_point, -0.02, labels[i], ha='center', fontsize=9, fontweight='bold', transform=plt.gca().get_xaxis_transform())
        prev_b = b

    plt.title(f"Mutational Impact Map: {top_id[:30]}...\n(Higher Point = Mutation Damages Prediction)")
    plt.xlabel("Position (bp)")
    plt.ylabel("Drop in Confidence (Delta)")
    plt.ylim(bottom=0)
    plt.legend()
    plt.tight_layout()
    plt.savefig("end_to_end_mutagenesis.png")
    print("Saved plot to 'end_to_end_mutagenesis.png'")

def main():
    print(f"--- Loading Model ---")
    model = RisdiplamModel().to(device)
    try:
        model.load_state_dict(torch.load(MODEL_PATH, map_location=device))
    except Exception as e:
        print(f"Error: {e}")
        return
    
    top_candidates = get_top_responders(model, INPUT_FILE, n=TOP_N)
    
    df_results, boundaries_map = analyze_mutations(model, top_candidates)
    
    print("\n" + "="*80)
    print("TOP 10 MOST DESTRUCTIVE MUTATIONS (Global)")
    print("="*80)
    if not df_results.empty:
        top10 = df_results.sort_values(by='delta', ascending=False).head(10)
        print(f"{'Position':<10} {'Mutation':<15} {'Before':<10} {'After':<10} {'Delta':<10}")
        print("-" * 60)
        for _, row in top10.iterrows():
            print(f"{row['pos']:<10} {row['desc']:<15} {row['score_before']:.4f}     {row['score_after']:.4f}     {row['delta']:.4f}")
    
    plot_heatmap(df_results, boundaries_map)

if __name__ == "__main__":
    main()