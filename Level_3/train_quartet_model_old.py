import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, average_precision_score
import random
import sys

# --- Config ---
DATA_FILE = "risdiplam_full_dataset.csv.gz"
BATCH_SIZE = 64
EPOCHS = 20
LEARNING_RATE = 0.001
DPSI_THRESH = 0.10

# Architecture Sizes
CTX_LEN = 50        # 25bp Exon + 25bp Intron (For A and B)
TARGET_LEN = 150    # 100bp Exon + 50bp Intron (For Pseudo Acceptor/Donor)

# --- 1. The Quartet Architecture ---
class QuartetNet(nn.Module):
    def __init__(self):
        super(QuartetNet, self).__init__()
        
        # Branch A/B: Flanking Contexts (50bp)
        self.context_conv = nn.Sequential(
            nn.Conv1d(4, 32, kernel_size=7, padding=3),
            nn.ReLU(),
            nn.MaxPool1d(2), # 50 -> 25
            nn.Flatten()
        )
        
        # Branch C: Pseudoexon Acceptor (The "Entry")
        # Input: 150bp (50 Intron + 100 Exon)
        self.acceptor_conv = nn.Sequential(
            nn.Conv1d(4, 64, kernel_size=9, padding=4),
            nn.ReLU(),
            nn.MaxPool1d(2), # 150 -> 75
            nn.Conv1d(64, 128, kernel_size=5, padding=2),
            nn.ReLU(),
            nn.MaxPool1d(3), # 75 -> 25
            nn.Flatten()
        )
        
        # Branch D: Pseudoexon Donor (The "Exit" & Drug Target)
        # Input: 150bp (100 Exon + 50 Intron)
        self.donor_conv = nn.Sequential(
            nn.Conv1d(4, 64, kernel_size=9, padding=4),
            nn.ReLU(),
            nn.MaxPool1d(2), # 150 -> 75
            nn.Conv1d(64, 128, kernel_size=5, padding=2),
            nn.ReLU(),
            nn.MaxPool1d(3), # 75 -> 25
            nn.Flatten()
        )
        
        # Fusion
        # Contexts: 25 * 32 = 800 features each
        # Targets: 25 * 128 = 3200 features each
        total_feats = 800 + 800 + 3200 + 3200 # 8000 inputs
        
        self.fc_fusion = nn.Sequential(
            nn.Linear(total_feats, 512),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(512, 128),
            nn.ReLU(),
            nn.Linear(128, 1),
            nn.Sigmoid()
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
        # Fixed size padding logic
        target_len = len(seq) 
        # Use torch.zeros directly
        arr = torch.zeros((4, target_len), dtype=torch.float32)
        for i, char in enumerate(seq):
            if char in mapping:
                arr[mapping[char], i] = 1.0
        return arr

    def __getitem__(self, idx):
        # one_hot now returns a tensor, so no need to wrap it
        # Convert label to tensor using torch.float32
        label = torch.tensor(self.labels[idx], dtype=torch.float32)
        
        return (
            self.one_hot(self.seqs_a[idx]),
            self.one_hot(self.seqs_b[idx]),
            self.one_hot(self.seqs_p_acc[idx]),
            self.one_hot(self.seqs_p_don[idx]),
            label
        )

# --- 3. Sequence Logic & Hard Negative Generation ---

def pad_sequence(seq, target_len, pad_side='right'):
    """Ensures sequence is exactly target_len by padding 'N' or truncating."""
    if len(seq) == target_len:
        return seq
    if len(seq) > target_len:
        if pad_side == 'right': return seq[:target_len]
        else: return seq[-target_len:]
    
    # Need padding
    pad = "N" * (target_len - len(seq))
    if pad_side == 'right': return seq + pad
    else: return pad + seq

def get_real_sample(row):
    try:
        # Context A: Last 25 Exon A, First 25 Intron 1
        ctx_a = row['seq_ExonA'][-25:] + row['seq_Intron1'][:25]
        
        # Context B: Last 25 Intron 2, First 25 Exon B
        ctx_b = row['seq_Intron2'][-25:] + row['seq_ExonB'][:25]
        
        # --- The Pseudoexon Overlap Strategy ---
        
        # P_Acc: 50bp Intron1 + 100bp PseudoExon
        p_exon = row['seq_PseudoExon']
        p_int1 = row['seq_Intron1']
        
        p_acc_intron = p_int1[-50:] 
        p_acc_exon = p_exon[:100] 
        p_acc = p_acc_intron + p_acc_exon
        p_acc = pad_sequence(p_acc, 150, pad_side='right')
        
        # P_Don: 100bp PseudoExon + 50bp Intron2
        p_int2 = row['seq_Intron2']
        
        p_don_exon = p_exon[-100:] 
        p_don_intron = p_int2[:50]
        p_don = p_don_exon + p_don_intron
        p_don = pad_sequence(p_don, 150, pad_side='left') 
        
        label = 1 if abs(row['dPSI']) > DPSI_THRESH else 0
        
        if len(ctx_a)!=50 or len(ctx_b)!=50: return None
        return (ctx_a, ctx_b, p_acc, p_don, label)
    except: return None

def generate_decoy_sample(row):
    try:
        source_seq = row['seq_Intron1'] if len(row['seq_Intron1']) > len(row['seq_Intron2']) else row['seq_Intron2']
        if len(source_seq) < 400: return None
        
        for _ in range(5):
            don_idx = random.randint(300, len(source_seq) - 100)
            if source_seq[don_idx:don_idx+2] not in ['GT', 'GC']: continue
            
            exon_len = random.randint(60, 200)
            acc_idx = don_idx - exon_len
            
            if source_seq[acc_idx-2:acc_idx] != 'AG': continue
            
            fake_exon = source_seq[acc_idx : don_idx]
            fake_intron_before = source_seq[acc_idx-50 : acc_idx]
            fake_intron_after = source_seq[don_idx : don_idx+50]
            
            p_acc = fake_intron_before + fake_exon[:100]
            p_acc = pad_sequence(p_acc, 150, pad_side='right')
            
            p_don = fake_exon[-100:] + fake_intron_after
            p_don = pad_sequence(p_don, 150, pad_side='left')
            
            ctx_a = row['seq_ExonA'][-25:] + row['seq_Intron1'][:25]
            ctx_b = row['seq_Intron2'][-25:] + row['seq_ExonB'][:25]
            
            return (ctx_a, ctx_b, p_acc, p_don, 0)
    except: return None

# --- 4. Main Workflow ---
def main():
    print("Loading Data...")
    if not pd.io.common.file_exists(DATA_FILE):
        print(f"Error: {DATA_FILE} not found. Run fetch_sequences.py first.")
        sys.exit(1)

    df = pd.read_csv(DATA_FILE)
    df = df.dropna(subset=['seq_ExonA', 'seq_Intron1', 'seq_PseudoExon', 'seq_Intron2', 'seq_ExonB'])
    
    data_a, data_b, data_p_acc, data_p_don, data_y = [], [], [], [], []
    
    print("Generating Dataset...")
    for _, row in df.iterrows():
        res = get_real_sample(row)
        if res:
            data_a.append(res[0]); data_b.append(res[1])
            data_p_acc.append(res[2]); data_p_don.append(res[3])
            data_y.append(res[4])
            
            decoy = generate_decoy_sample(row)
            if decoy:
                data_a.append(decoy[0]); data_b.append(decoy[1])
                data_p_acc.append(decoy[2]); data_p_don.append(decoy[3])
                data_y.append(decoy[4])

    print(f"Total Samples: {len(data_y)}")
    if len(data_y) == 0:
        print("Error: No valid samples generated. Check data quality.")
        sys.exit(1)
    
    X_a = np.array(data_a); X_b = np.array(data_b)
    X_acc = np.array(data_p_acc); X_don = np.array(data_p_don)
    y = np.array(data_y)
    
    # Split
    indices = np.arange(len(y))
    train_idx, test_idx = train_test_split(indices, test_size=0.2, random_state=42)
    
    train_ds = QuartetDataset(X_a[train_idx], X_b[train_idx], X_acc[train_idx], X_don[train_idx], y[train_idx])
    test_ds = QuartetDataset(X_a[test_idx], X_b[test_idx], X_acc[test_idx], X_don[test_idx], y[test_idx])
    
    train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True)
    test_loader = DataLoader(test_ds, batch_size=BATCH_SIZE)
    
    # Train
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Training on {device}...")
    
    model = QuartetNet().to(device)
    optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE)
    criterion = nn.BCELoss()
    
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
            
        print(f"Epoch {epoch+1} Loss: {total_loss/len(train_loader):.4f}")

    # Evaluate
    model.eval()
    preds, targets = [], []
    with torch.no_grad():
        for xa, xb, xpa, xpd, lbl in test_loader:
            xa, xb, xpa, xpd = xa.to(device), xb.to(device), xpa.to(device), xpd.to(device)
            out = model(xa, xb, xpa, xpd)
            preds.extend(out.cpu().numpy())
            targets.extend(lbl.numpy())
            
    print(f"AUROC: {roc_auc_score(targets, preds):.4f}")
    print(f"AUPRC: {average_precision_score(targets, preds):.4f}")
    
    torch.save(model.state_dict(), "quartet_full_net.pt")
    print("Model Saved.")

if __name__ == "__main__":
    main()