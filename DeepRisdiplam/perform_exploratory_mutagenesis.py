import torch
import torch.nn as nn
import pandas as pd
import numpy as np
import csv
import matplotlib.pyplot as plt
import itertools
from tqdm import tqdm
# [Optimization] Import AMP for Mixed Precision
from torch.cuda.amp import autocast

# --- Configuration ---
VERSION = "v8_components"
MODEL_PATH = f"deep_risdiplam_model_{VERSION}.pth"
INPUT_FILE = "dataset_components_v8.tsv"
MAX_SEQ_LEN = 3000

# [Optimization] Increased Batch Size due to FP16 memory savings
# RTX 3070 Ti (8GB) can likely handle 2048 or even 4096 in mixed precision for this model size.
BATCH_SIZE = 2048  
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

def get_top_responders(model, filename, n=TOP_N):
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
    
    # [Optimization] Use inference_mode instead of no_grad (faster, less overhead)
    with torch.inference_mode():
        for i in range(0, len(candidates), BATCH_SIZE):
            batch = candidates[i:i+BATCH_SIZE]
            seqs = [stitch(r)[0] for r in batch]
            tensors = torch.stack([one_hot_encode(s) for s in seqs]).to(device)
            
            # [Optimization] Autocast for Mixed Precision Inference
            with autocast():
                preds = torch.sigmoid(model(tensors)).cpu().numpy().flatten()
                
            for j, score in enumerate(preds):
                scored_candidates.append((score, batch[j]))
    
    scored_candidates.sort(key=lambda x: x[0], reverse=True)
    return scored_candidates[:n]

def generate_2bp_mutations(seq):
    """
    Generates:
    1. 2-bp Substitutions (at i, i+1)
    2. 2-bp Deletions (removing i, i+1)
    3. 2-bp Insertions (at i)
    """
    mutants = []
    bases = ['A', 'C', 'G', 'T']
    dimers = [''.join(p) for p in itertools.product(bases, repeat=2)]
    
    L = len(seq)
    
    # 1. Substitutions (i, i+1)
    for i in range(L - 1):
        original = seq[i:i+2]
        for dimer in dimers:
            if dimer == original: continue
            mut_seq = seq[:i] + dimer + seq[i+2:]
            mutants.append((mut_seq, f"Sub {original} at {i} -> {dimer}", i))
            
    # 2. Deletions (i, i+1)
    for i in range(L - 1):
        mut_seq = seq[:i] + seq[i+2:]
        mutants.append((mut_seq, f"Del 2bp at {i}", i))
        
    # 3. Insertions (at i)
    for i in range(L + 1):
        for dimer in dimers:
            mut_seq = seq[:i] + dimer + seq[i:]
            mutants.append((mut_seq, f"Ins {dimer} at {i}", i))
            
    return mutants

def analyze_mutations(model, top_candidates):
    print(f"\nPerforming Deep 2-bp Exploratory Mutagenesis on {len(top_candidates)} sequences...")
    
    results_map = {} 
    boundaries_map = {}
    
    for score, row in top_candidates:
        full_seq, boundaries = stitch(row)
        original_score = score
        event_id = row['event_id']
        boundaries_map[event_id] = boundaries
        
        print(f"Processing {event_id[:20]}... (Len: {len(full_seq)})")
        
        mutants = generate_2bp_mutations(full_seq)
        
        pos_impact = {} 
        
        batch_seqs = []
        batch_meta = []
        
        def flush_batch(b_seqs, b_meta):
            t = torch.stack(b_seqs).to(device)
            
            # [Optimization] Use inference_mode + autocast
            with torch.inference_mode():
                with autocast():
                    # AMP runs the convolution ops in FP16
                    new_scores = torch.sigmoid(model(t)).float().cpu().numpy().flatten()
            
            for k, new_s in enumerate(new_scores):
                change = new_s - original_score 
                desc, pos = b_meta[k]
                
                if pos not in pos_impact:
                    pos_impact[pos] = {'min_change': 0.0, 'max_change': 0.0, 'worst_mut': '', 'best_mut': ''}
                
                if change < pos_impact[pos]['min_change']:
                    pos_impact[pos]['min_change'] = change
                    pos_impact[pos]['worst_mut'] = desc
                
                if change > pos_impact[pos]['max_change']:
                    pos_impact[pos]['max_change'] = change
                    pos_impact[pos]['best_mut'] = desc

        for m_seq, m_desc, m_pos in tqdm(mutants, leave=False):
            batch_seqs.append(one_hot_encode(m_seq))
            batch_meta.append((m_desc, m_pos))
            
            if len(batch_seqs) == BATCH_SIZE:
                flush_batch(batch_seqs, batch_meta)
                batch_seqs = []
                batch_meta = []
                
        if batch_seqs:
            flush_batch(batch_seqs, batch_meta)
            
        results_map[event_id] = pos_impact

    return results_map, boundaries_map

def plot_exploratory(results_map, boundaries_map):
    print("\nGenerating Exploratory Plots...")
    
    top_id = list(results_map.keys())[0]
    data = results_map[top_id]
    bounds = boundaries_map[top_id]
    
    positions = sorted(data.keys())
    min_changes = [data[p]['min_change'] for p in positions]
    max_changes = [data[p]['max_change'] for p in positions]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10), sharex=True)
    
    labels = ["Exon A", "Intron 1", "Pseudo", "Intron 2", "Exon B"]
    colors = ['#f0f0f0', '#e0e0e0', '#ffeb99', '#e0e0e0', '#f0f0f0']
    
    # --- Plot 1: Negative Impact ---
    ax1.plot(positions, min_changes, color='crimson', linewidth=1, alpha=0.8)
    ax1.fill_between(positions, min_changes, 0, color='crimson', alpha=0.3)
    ax1.set_ylabel("Change (Decrease)")
    ax1.set_title(f"Destructive Potential (Max Decrease)\\n{top_id[:40]}")
    ax1.grid(True, alpha=0.2)
    
    prev_b = 0
    for i, b in enumerate(bounds):
        ax1.axvspan(prev_b, b, color=colors[i], alpha=0.3, zorder=0)
        ax1.axvline(b, color='black', linestyle=':', alpha=0.5)
        prev_b = b

    # --- Plot 2: Positive Impact ---
    ax2.plot(positions, max_changes, color='forestgreen', linewidth=1, alpha=0.8)
    ax2.fill_between(positions, max_changes, 0, color='forestgreen', alpha=0.3)
    ax2.set_ylabel("Change (Increase)")
    ax2.set_title("Enhancement Potential (Max Increase)")
    ax2.set_xlabel("Position (bp)")
    ax2.grid(True, alpha=0.2)
    
    prev_b = 0
    for i, b in enumerate(bounds):
        ax2.axvspan(prev_b, b, color=colors[i], alpha=0.3, zorder=0)
        ax2.axvline(b, color='black', linestyle=':', alpha=0.5)
        mid_point = (prev_b + b) / 2
        ax2.text(mid_point, -0.05, labels[i], ha='center', fontsize=9, fontweight='bold', transform=ax2.get_xaxis_transform())
        prev_b = b

    plt.tight_layout()
    plt.savefig("exploratory_mutagenesis_2bp.png")
    print("Saved plot to 'exploratory_mutagenesis_2bp.png'")
    
    print("\n" + "="*80)
    print(f"TOP 5 DESTRUCTIVE MUTATIONS for {top_id[:20]}")
    print("="*80)
    sorted_dest = sorted(data.items(), key=lambda x: x[1]['min_change'])
    for p, info in sorted_dest[:5]:
        print(f"Pos {p:<5} | {info['worst_mut']:<30} | Change: {info['min_change']:.4f}")
        
    print("\n" + "="*80)
    print(f"TOP 5 ENHANCING MUTATIONS for {top_id[:20]}")
    print("="*80)
    sorted_enh = sorted(data.items(), key=lambda x: x[1]['max_change'], reverse=True)
    real_enh = [x for x in sorted_enh if x[1]['max_change'] > 0.001]
    if not real_enh:
        print("No mutations significantly increased the score.")
    else:
        for p, info in real_enh[:5]:
            print(f"Pos {p:<5} | {info['best_mut']:<30} | Change: +{info['max_change']:.4f}")

def main():
    print(f"--- Loading Model (Optimized) ---")
    model = RisdiplamModel().to(device)
    try:
        model.load_state_dict(torch.load(MODEL_PATH, map_location=device))
    except Exception as e:
        print(f"Error: {e}")
        return
    
    top_candidates = get_top_responders(model, INPUT_FILE, n=TOP_N)
    results_map, boundaries_map = analyze_mutations(model, top_candidates)
    plot_exploratory(results_map, boundaries_map)

if __name__ == "__main__":
    main()