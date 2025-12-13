import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
import torch.utils.checkpoint as checkpoint
from torch.amp import autocast, GradScaler 
from torch.utils.data import DataLoader, Dataset, WeightedRandomSampler
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_score, recall_score, f1_score, confusion_matrix
from tqdm import tqdm
import sys
import os
import random

# --- 1. CONFIGURATION ---
VERSION = "v15_hybrid_multimodal"
INPUT_FILE = "dataset_components_v8.tsv"
MODEL_SAVE_PATH = f"deep_risdiplam_model_{VERSION}.pth"
LOG_FILE = f"training_log_{VERSION}.txt"

# VRAM OPTIMIZATION
MAX_SEQ_LEN = 1500      
MICRO_BATCH_SIZE = 32   
ACCUM_STEPS = 4         
EFFECTIVE_BATCH_SIZE = MICRO_BATCH_SIZE * ACCUM_STEPS
LEARNING_RATE = 1e-4
EPOCHS = 35 

# FEATURE CONFIG
# We use the top features identified by the Random Forest/LR exploration
NUM_SCALAR_FEATURES = 6 

torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.allow_tf32 = True
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# --- 2. HELPER FUNCTIONS ---
def get_gc(seq):
    if len(seq) == 0: return 0.0
    return (seq.count('G') + seq.count('C') + seq.count('g') + seq.count('c')) / len(seq)

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
    for base_val in range(4): one_hot[:, base_val, :] = (arr == base_val)
    return torch.tensor(one_hot)

# --- 3. DATA PIPELINE (With Explicit Features) ---
class ComponentDataset(Dataset):
    def __init__(self, df, is_training=False):
        self.df = df
        self.is_training = is_training
        
        # Franken-Bank (Keep V11 Augmentations)
        self.positive_intron1_middles = []
        self.positive_intron2_middles = []
        if is_training:
            positives = df[df['label'] == 1]
            for _, row in positives.iterrows():
                i1 = row['Intron1']; i2 = row['Intron2']
                if len(i1) > 100: self.positive_intron1_middles.append(i1[50:-50])
                if len(i2) > 100: self.positive_intron2_middles.append(i2[50:-50])

    def __len__(self):
        return len(self.df)

    def scramble_string(self, s):
        l = list(s)
        random.shuffle(l)
        return "".join(l)

    def apply_masking(self, seq, mask_prob=0.10):
        seq_list = list(seq)
        seq_len = len(seq_list)
        num_mask = int(seq_len * mask_prob)
        mask_indices = np.random.choice(seq_len, num_mask, replace=False)
        for idx in mask_indices: seq_list[idx] = 'N' 
        return "".join(seq_list)

    def get_scalar_features(self, ps, i1, i2):
        # NORMALIZATION IS CRITICAL FOR NEURAL NETS
        # We scale inputs to be roughly 0-1 range
        
        # 1. Pseudo Length (Mean ~100, Max ~500. Divide by 200)
        f_len_ps = len(ps) / 200.0
        
        # 2. Pseudo GC (Already 0-1)
        f_gc_ps = get_gc(ps)
        
        # 3. Intron 1 Length (Mean ~250-500. Divide by 500)
        f_len_i1 = len(i1) / 500.0
        
        # 4. Intron 1 GC
        f_gc_i1 = get_gc(i1)
        
        # 5. Canonical Donor (Binary 0 or 1)
        f_canon = 1.0 if i2.upper().startswith("GT") else 0.0
        
        # 6. Risdiplam Signature (End of Pseudo is G)
        f_g_minus1 = 1.0 if ps.upper().endswith("G") else 0.0
        
        return torch.tensor([f_len_ps, f_gc_ps, f_len_i1, f_gc_i1, f_canon, f_g_minus1], dtype=torch.float32)

    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        ea, i1, ps, i2, eb = row['ExonA'], row['Intron1'], row['Pseudoexon'], row['Intron2'], row['ExonB']
        label = float(row['label'])
        
        if self.is_training:
            # Franken-Swap
            if label == 1.0:
                if len(self.positive_intron1_middles) > 0 and len(i1) > 100 and np.random.rand() < 0.50:
                    new_mid = random.choice(self.positive_intron1_middles)
                    i1 = i1[:50] + new_mid + i1[-50:]
                if len(self.positive_intron2_middles) > 0 and len(i2) > 100 and np.random.rand() < 0.50:
                    new_mid = random.choice(self.positive_intron2_middles)
                    i2 = i2[:50] + new_mid + i2[-50:]
            
            # Jitter
            crop_a = np.random.randint(0, 16)
            if len(ea) > crop_a: ea = ea[crop_a:]
            crop_b = np.random.randint(0, 16)
            if len(eb) > crop_b: eb = eb[:-crop_b]

            # Decoys
            if label == 1.0:
                rand_val = np.random.rand()
                if rand_val < 0.15: 
                     i2 = "CC" + i2[2:] 
                     label = 0.0 
                elif rand_val < 0.30: 
                    if len(i1) > 2: i1 = self.scramble_string(i1[:-2]) + i1[-2:]
                    ps = self.scramble_string(ps)
                    if len(i2) > 2: i2 = i2[:2] + self.scramble_string(i2[2:])
                    label = 0.0

            # Extract Scalars BEFORE masking (we want true stats)
            scalars = self.get_scalar_features(ps, i1, i2)

            full_seq = ea + i1 + ps + i2 + eb
            full_seq = self.apply_masking(full_seq, mask_prob=0.10)
            return full_seq, scalars, torch.tensor([label], dtype=torch.float32)

        # Validation
        scalars = self.get_scalar_features(ps, i1, i2)
        full_seq = ea + i1 + ps + i2 + eb
        return full_seq, scalars, torch.tensor([label], dtype=torch.float32)

def collate_fn(batch):
    seqs, scalars, labels = zip(*batch)
    X_seq = vectorized_one_hot_encode(seqs)
    X_scal = torch.stack(scalars)
    y = torch.stack(labels)
    return X_seq, X_scal, y

# --- 4. MULTIMODAL ARCHITECTURE ---
class RisdiplamHybridMultimodal(nn.Module):
    def __init__(self):
        super(RisdiplamHybridMultimodal, self).__init__()
        
        # --- Deep Path (Sequence) ---
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
        
        # --- Fusion Layer ---
        # 128 (Seq Features) + 6 (Scalar Features) = 134
        self.fc1 = nn.Linear(128 + NUM_SCALAR_FEATURES, 64)
        self.dropout = nn.Dropout(0.65) 
        self.fc2 = nn.Linear(64, 1)

    def _attn_block(self, x):
        attn_out, _ = self.attention(x, x, x)
        return x + attn_out

    def forward(self, x_seq, x_scal):
        # 1. Process Sequence
        x = self.relu(self.bn1(self.conv1(x_seq)))
        x = self.pool(x)
        x = self.relu(self.bn2(self.conv2(x)))
        x = x.permute(0, 2, 1) 
        if self.training and x.requires_grad:
            x = checkpoint.checkpoint(self._attn_block, x, use_reentrant=False)
        else:
            x = self._attn_block(x)
        x = x.permute(0, 2, 1)
        x = self.global_pool(x)
        seq_feats = self.flatten(x) # [Batch, 128]
        
        # 2. Concatenate Scalars
        # [Batch, 128] + [Batch, 6] -> [Batch, 134]
        combined = torch.cat((seq_feats, x_scal), dim=1)
        
        # 3. Classify
        x = self.relu(self.fc1(combined))
        x = self.dropout(x)
        x = self.fc2(x)
        return x

# --- 5. TRAINING LOOP ---
def log_message(message):
    print(message)
    with open(LOG_FILE, "a") as f:
        f.write(message + "\n")

def evaluate(model, loader):
    model.eval()
    total_loss = 0
    all_preds = []
    all_labels = []
    criterion = nn.BCEWithLogitsLoss()
    
    with torch.inference_mode():
        for inputs_seq, inputs_scal, labels in loader:
            inputs_seq = inputs_seq.to(device)
            inputs_scal = inputs_scal.to(device)
            labels = labels.to(device)
            
            with autocast(device_type='cuda'):
                logits = model(inputs_seq, inputs_scal)
                loss = criterion(logits, labels)
            total_loss += loss.item()
            predicted = (torch.sigmoid(logits) > 0.5).float()
            all_preds.extend(predicted.cpu().numpy())
            all_labels.extend(labels.cpu().numpy())
            
    # Metrics
    prec = precision_score(all_labels, all_preds, zero_division=0)
    rec = recall_score(all_labels, all_preds, zero_division=0)
    f1 = f1_score(all_labels, all_preds, zero_division=0)
    cm = confusion_matrix(all_labels, all_preds, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()
    
    return total_loss / len(loader), prec, rec, f1, (tn, fp, fn, tp)

def train_model():
    log_message(f"--- Starting Multimodal V15 (Wide & Deep) Training on {torch.cuda.get_device_name(0)} ---")
    log_message(f"Feature: Merging CNN Sequence Features with 6 Explicit Biological Scalars.")
    
    df = pd.read_csv(INPUT_FILE, sep='\t')
    df = df.dropna(subset=['ExonA', 'Intron1', 'Pseudoexon', 'Intron2', 'ExonB'])
    train_df, val_df = train_test_split(df, test_size=0.1, random_state=42, stratify=df['label'])
    
    train_labels = train_df['label'].values.astype(int)
    class_counts = np.bincount(train_labels)
    weight_per_class = 1.0 / class_counts
    sample_weights = weight_per_class[train_labels]
    sampler = WeightedRandomSampler(weights=sample_weights, num_samples=len(sample_weights), replacement=True)
    
    train_ds = ComponentDataset(train_df, is_training=True)
    val_ds = ComponentDataset(val_df, is_training=False)
    
    train_loader = DataLoader(train_ds, batch_size=MICRO_BATCH_SIZE, sampler=sampler, shuffle=False, 
                              num_workers=4, collate_fn=collate_fn, pin_memory=True)
    val_loader = DataLoader(val_ds, batch_size=MICRO_BATCH_SIZE*2, shuffle=False, 
                            num_workers=4, collate_fn=collate_fn, pin_memory=True)

    model = RisdiplamHybridMultimodal().to(device)
    optimizer = optim.AdamW(model.parameters(), lr=LEARNING_RATE, weight_decay=1e-3)
    scaler = GradScaler('cuda')
    criterion = nn.BCEWithLogitsLoss()
    
    best_f1 = 0.0
    
    for epoch in range(EPOCHS):
        model.train()
        running_loss = 0.0
        optimizer.zero_grad()
        
        pbar = tqdm(train_loader, desc=f"Epoch {epoch+1}/{EPOCHS}")
        for i, (inputs_seq, inputs_scal, labels) in enumerate(pbar):
            inputs_seq = inputs_seq.to(device, non_blocking=True)
            inputs_scal = inputs_scal.to(device, non_blocking=True)
            labels = labels.to(device, non_blocking=True)
            
            with autocast(device_type='cuda'):
                logits = model(inputs_seq, inputs_scal)
                loss = criterion(logits, labels)
                loss = loss / ACCUM_STEPS
            scaler.scale(loss).backward()
            
            if (i + 1) % ACCUM_STEPS == 0:
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                scaler.step(optimizer)
                scaler.update()
                optimizer.zero_grad()
            
            running_loss += loss.item() * ACCUM_STEPS
            pbar.set_postfix({'loss': f"{running_loss / (i+1):.4f}"})
            
        val_loss, prec, rec, f1, (tn, fp, fn, tp) = evaluate(model, val_loader)
        
        log_msg = (f"Epoch {epoch+1}: Loss={val_loss:.4f} | Prec={prec:.4f} | Rec={rec:.4f} | F1={f1:.4f} "
                   f"| TP={tp} FN={fn}")
        
        if f1 > best_f1:
            best_f1 = f1
            torch.save(model.state_dict(), MODEL_SAVE_PATH)
            log_msg += " [SAVED BEST]"
        
        log_message(log_msg)

if __name__ == "__main__":
    train_model()