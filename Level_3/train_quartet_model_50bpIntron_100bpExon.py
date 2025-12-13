import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from scipy.stats import pearsonr
import random
import sys

# --- Config ---
DATA_FILE = "risdiplam_full_dataset.csv.gz"
BATCH_SIZE = 64
EPOCHS = 20
LEARNING_RATE = 0.0005 # Lower LR for regression stability
DPSI_THRESH = 0.10     # Threshold for "Interest" (for balancing only)
DECOY_RATIO = 5 

# Architecture Sizes
CTX_LEN = 50        
TARGET_LEN = 150    

# --- 1. The Quartet Architecture (Regression) ---
class QuartetNet(nn.Module):
    def __init__(self):
        super(QuartetNet, self).__init__()
        
        # Branch A/B: Flanking Contexts (50bp)
        self.context_conv = nn.Sequential(
            nn.Conv1d(4, 32, kernel_size=7, padding=3),
            nn.ReLU(),
            nn.MaxPool1d(2), 
            nn.Flatten()
        )
        
        # Branch C: Pseudoexon Acceptor
        self.acceptor_conv = nn.Sequential(
            nn.Conv1d(4, 64, kernel_size=9, padding=4),
            nn.ReLU(),
            nn.MaxPool1d(2), 
            nn.Conv1d(64, 128, kernel_size=5, padding=2),
            nn.ReLU(),
            nn.MaxPool1d(3), 
            nn.Flatten()
        )
        
        # Branch D: Pseudoexon Donor
        self.donor_conv = nn.Sequential(
            nn.Conv1d(4, 64, kernel_size=9, padding=4),
            nn.ReLU(),
            nn.MaxPool1d(2), 
            nn.Conv1d(64, 128, kernel_size=5, padding=2),
            nn.ReLU(),
            nn.MaxPool1d(3), 
            nn.Flatten()
        )
        
        # Fusion
        total_feats = 800 + 800 + 3200 + 3200 
        
        self.fc_fusion = nn.Sequential(
            nn.Linear(total_feats, 512),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(512, 128),
            nn.ReLU(),
            nn.Linear(128, 1)
            # NO SIGMOID: We want raw regression output (can be negative or >1)
        )

    def forward(self, xa, xb, xp_acc, xp_don):
        f_a = self.context_conv(xa)
        f_b = self.context_conv(xb)
        f_acc = self.acceptor_conv(xp_acc)
        f_don = self.donor_conv(xp_don)
        combined = torch.cat((f_a, f_b, f_acc, f_don), dim=1)
        return self.fc_fusion(combined)

# --- 2. Dataset ---
class QuartetDataset(Dataset):
    def __init__(self, seqs_a, seqs_b, seqs_p_acc, seqs_p_don, labels):
        self.seqs_a = seqs_a
        self.seqs_b = seqs_b
        self.seqs_p_acc = seqs_p_acc
        self.seqs_p_don = seqs_p_don
        self.labels = labels

    def __len__(self):
        return len(self.labels)

    def one_hot(self, seq):
        mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        target_len = len(seq) 
        arr = torch.zeros((4, target_len), dtype=torch.float32)
        for i, char in enumerate(seq):
            if char in mapping:
                arr[mapping[char], i] = 1.0
        return arr

    def __getitem__(self, idx):
        label = torch.tensor(self.labels[idx], dtype=torch.float32)
        return (
            self.one_hot(self.seqs_a[idx]),
            self.one_hot(self.seqs_b[idx]),
            self.one_hot(self.seqs_p_acc[idx]),
            self.one_hot(self.seqs_p_don[idx]),
            label
        )

# --- 3. Sequence Logic ---

def pad_sequence(seq, target_len, pad_side='right'):
    if len(seq) == target_len: return seq
    if len(seq) > target_len:
        if pad_side == 'right': return seq[:target_len]
        else: return seq[-target_len:]
    pad = "N" * (target_len - len(seq))
    if pad_side == 'right': return seq + pad
    else: return pad + seq

def score_donor_strength(seq_window):
    if len(seq_window) != 9: return 0
    score = 0
    if seq_window[0] in ['C', 'A']: score += 1
    if seq_window[1] == 'A': score += 1
    if seq_window[2] == 'G': score += 2 
    if seq_window[5] in ['A', 'G']: score += 1
    if seq_window[6] == 'A': score += 1
    if seq_window[7] == 'G': score += 1.5
    return score

def score_acceptor_strength(seq_window):
    if len(seq_window) != 15: return 0
    if seq_window[-2:] != 'AG': return 0
    py_tract = seq_window[:-2]
    py_count = py_tract.count('C') + py_tract.count('T')
    return py_count / len(py_tract)

def get_real_sample(row):
    try:
        ctx_a = row['seq_ExonA'][-25:] + row['seq_Intron1'][:25]
        ctx_b = row['seq_Intron2'][-25:] + row['seq_ExonB'][:25]
        
        p_exon = row['seq_PseudoExon']
        p_int1 = row['seq_Intron1']
        
        p_acc_intron = p_int1[-50:] 
        p_acc_exon = p_exon[:100] 
        p_acc = p_acc_intron + p_acc_exon
        p_acc = pad_sequence(p_acc, 150, pad_side='right')
        
        p_int2 = row['seq_Intron2']
        p_don_exon = p_exon[-100:] 
        p_don_intron = p_int2[:50]
        p_don = p_don_exon + p_don_intron
        p_don = pad_sequence(p_don, 150, pad_side='left') 
        
        # REGRESSION CHANGE: Return raw dPSI value
        label = float(row['dPSI'])
        
        if len(ctx_a)!=50 or len(ctx_b)!=50: return None
        return (ctx_a, ctx_b, p_acc, p_don, label)
    except: return None

def generate_adversarial_decoys(row, target_count=DECOY_RATIO):
    decoys = []
    try:
        source_seq = row['seq_Intron1'] if len(row['seq_Intron1']) > len(row['seq_Intron2']) else row['seq_Intron2']
        if len(source_seq) < 500: return []
        
        for _ in range(50):
            if len(decoys) >= target_count: break
            don_idx = random.randint(300, len(source_seq) - 100)
            
            donor_window = source_seq[don_idx-3 : don_idx+6]
            if source_seq[don_idx:don_idx+2] not in ['GT', 'GC']: continue
            if score_donor_strength(donor_window) < 3.0: continue
            
            exon_len = random.randint(60, 200)
            acc_idx = don_idx - exon_len
            
            if source_seq[acc_idx-2:acc_idx] != 'AG': continue
            acc_window = source_seq[acc_idx-15 : acc_idx] 
            if score_acceptor_strength(acc_window) < 0.5: continue
            
            fake_exon = source_seq[acc_idx : don_idx]
            fake_intron_before = source_seq[acc_idx-50 : acc_idx]
            fake_intron_after = source_seq[don_idx : don_idx+50]
            
            p_acc = fake_intron_before + fake_exon[:100]
            p_acc = pad_sequence(p_acc, 150, pad_side='right')
            p_don = fake_exon[-100:] + fake_intron_after
            p_don = pad_sequence(p_don, 150, pad_side='left')
            
            ctx_a = row['seq_ExonA'][-25:] + row['seq_Intron1'][:25]
            ctx_b = row['seq_Intron2'][-25:] + row['seq_ExonB'][:25]
            
            # REGRESSION CHANGE: Decoys have 0.0 dPSI
            decoys.append((ctx_a, ctx_b, p_acc, p_don, 0.0))
            
        return decoys
    except: return []

def smart_balance_dataset(data_arrays):
    """
    Regression Balancing.
    We classify samples into "Interesting" (abs(dPSI) > 0.1) and "Boring".
    We keep ALL Interesting samples and Decoys.
    We downsample the Boring ones.
    """
    X_a, X_b, X_acc, X_don, y, is_decoy = data_arrays
    
    # Define interest based on dPSI magnitude
    # Note: y is now a float array
    idx_interesting = np.where(np.abs(y) > DPSI_THRESH)[0]
    idx_decoy = np.where(is_decoy == 1)[0]
    idx_boring = np.where((np.abs(y) <= DPSI_THRESH) & (is_decoy == 0))[0]
    
    n_int = len(idx_interesting)
    n_decoy = len(idx_decoy)
    n_boring = len(idx_boring)
    
    print("\n--- Smart Balancing (Regression Mode) ---")
    print(f"  Interesting (dPSI > {DPSI_THRESH}): {n_int}")
    print(f"  Hard Decoys: {n_decoy}")
    print(f"  Boring (dPSI <= {DPSI_THRESH}): {n_boring}")
    
    target_boring = int((n_int + n_decoy) * 1.5)
    
    if target_boring < n_boring:
        print(f"  -> Downsampling Boring samples to {target_boring}...")
        idx_boring = np.random.choice(idx_boring, target_boring, replace=False)
    
    final_indices = np.concatenate([idx_interesting, idx_decoy, idx_boring])
    np.random.shuffle(final_indices)
    
    print(f"  Final Training Set Size: {len(final_indices)}")
    
    return (
        X_a[final_indices], X_b[final_indices], 
        X_acc[final_indices], X_don[final_indices], 
        y[final_indices]
    )

# --- 4. Main Workflow ---
def main():
    print("Loading Data...")
    if not pd.io.common.file_exists(DATA_FILE):
        print(f"Error: {DATA_FILE} not found.")
        sys.exit(1)

    df = pd.read_csv(DATA_FILE)
    df = df.dropna(subset=['seq_ExonA', 'seq_Intron1', 'seq_PseudoExon', 'seq_Intron2', 'seq_ExonB'])
    
    data_a, data_b, data_p_acc, data_p_don, data_y, data_is_decoy = [], [], [], [], [], []
    
    print(f"Generating Dataset (Regression Mode)...")
    
    for _, row in df.iterrows():
        res = get_real_sample(row)
        if res:
            data_a.append(res[0]); data_b.append(res[1])
            data_p_acc.append(res[2]); data_p_don.append(res[3])
            data_y.append(res[4]) # Appends FLOAT
            data_is_decoy.append(0)
            
            decoys = generate_adversarial_decoys(row, target_count=DECOY_RATIO)
            for decoy in decoys:
                data_a.append(decoy[0]); data_b.append(decoy[1])
                data_p_acc.append(decoy[2]); data_p_don.append(decoy[3])
                data_y.append(decoy[4]) # Appends 0.0
                data_is_decoy.append(1)

    arrays = (
        np.array(data_a), np.array(data_b), 
        np.array(data_p_acc), np.array(data_p_don), 
        np.array(data_y), np.array(data_is_decoy)
    )
    
    X_a, X_b, X_acc, X_don, y = smart_balance_dataset(arrays)
    
    indices = np.arange(len(y))
    train_idx, test_idx = train_test_split(indices, test_size=0.2, random_state=42)
    
    train_ds = QuartetDataset(X_a[train_idx], X_b[train_idx], X_acc[train_idx], X_don[train_idx], y[train_idx])
    test_ds = QuartetDataset(X_a[test_idx], X_b[test_idx], X_acc[test_idx], X_don[test_idx], y[test_idx])
    
    train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True)
    test_loader = DataLoader(test_ds, batch_size=BATCH_SIZE)
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Training on {device}...")
    
    model = QuartetNet().to(device)
    optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE)
    criterion = nn.MSELoss() # Regression Loss
    
    for epoch in range(EPOCHS):
        model.train()
        total_loss = 0
        for xa, xb, xpa, xpd, lbl in train_loader:
            xa, xb, xpa, xpd = xa.to(device), xb.to(device), xpa.to(device), xpd.to(device)
            lbl = lbl.to(device).unsqueeze(1)
            
            optimizer.zero_grad()
            out = model(xa, xb, xpa, xpd)
            loss = criterion(out, lbl)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
            
        print(f"Epoch {epoch+1} MSE Loss: {total_loss/len(train_loader):.5f}")

    model.eval()
    preds, targets = [], []
    with torch.no_grad():
        for xa, xb, xpa, xpd, lbl in test_loader:
            xa, xb, xpa, xpd = xa.to(device), xb.to(device), xpa.to(device), xpd.to(device)
            out = model(xa, xb, xpa, xpd)
            preds.extend(out.cpu().numpy().flatten())
            targets.extend(lbl.numpy().flatten())
            
    # Regression Metrics
    r_score, _ = pearsonr(targets, preds)
    print(f"Test Results:")
    print(f"  Pearson Correlation (R): {r_score:.4f}")
    print(f"  (Range -1 to 1. Higher is better.)")
    
    torch.save(model.state_dict(), "quartet_regression_net.pt")
    print("Model Saved.")

if __name__ == "__main__":
    main()