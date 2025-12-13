import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
import torch.utils.checkpoint as checkpoint
from torch.amp import autocast, GradScaler 
from torch.utils.data import DataLoader, Dataset, WeightedRandomSampler
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_score, recall_score
from tqdm import tqdm
import sys
import os
import random

# --- 1. CONFIGURATION ---
VERSION = "v13_hybrid_focal"
INPUT_FILE = "dataset_components_v8.tsv"
MODEL_SAVE_PATH = f"deep_risdiplam_model_{VERSION}.pth"
LOG_FILE = f"training_log_{VERSION}.txt"

# VRAM OPTIMIZATION (RTX 3070 Ti)
MAX_SEQ_LEN = 1500      
MICRO_BATCH_SIZE = 32   
ACCUM_STEPS = 4         
EFFECTIVE_BATCH_SIZE = MICRO_BATCH_SIZE * ACCUM_STEPS
LEARNING_RATE = 1e-4 # Low LR is crucial for Focal Loss stability
EPOCHS = 40          # Hard mining takes longer to converge

torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.allow_tf32 = True

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# --- 2. FOCAL LOSS IMPLEMENTATION ---
class FocalLoss(nn.Module):
    def __init__(self, alpha=1, gamma=2, reduction='mean'):
        super(FocalLoss, self).__init__()
        self.alpha = alpha
        self.gamma = gamma
        self.reduction = reduction
        self.bce = nn.BCEWithLogitsLoss(reduction='none')

    def forward(self, inputs, targets):
        bce_loss = self.bce(inputs, targets)
        pt = torch.exp(-bce_loss) # p_t is the probability of being right
        focal_loss = self.alpha * (1 - pt) ** self.gamma * bce_loss
        
        if self.reduction == 'mean':
            return focal_loss.mean()
        elif self.reduction == 'sum':
            return focal_loss.sum()
        else:
            return focal_loss

# --- 3. DATA PIPELINE ---

def log_message(message):
    print(message)
    with open(LOG_FILE, "a") as f:
        f.write(message + "\n")

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

class ComponentDataset(Dataset):
    def __init__(self, df, is_training=False):
        self.df = df
        self.is_training = is_training
        
        # Franken-Bank (From V11)
        self.positive_intron1_middles = []
        self.positive_intron2_middles = []
        
        if is_training:
            positives = df[df['label'] == 1]
            for _, row in positives.iterrows():
                i1 = row['Intron1']
                i2 = row['Intron2']
                if len(i1) > 100: self.positive_intron1_middles.append(i1[50:-50])
                if len(i2) > 100: self.positive_intron2_middles.append(i2[50:-50])
            print(f"Franken-Bank Loaded: {len(self.positive_intron1_middles)} I1, {len(self.positive_intron2_middles)} I2.")

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
        for idx in mask_indices:
            seq_list[idx] = 'N' 
        return "".join(seq_list)

    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        
        ea, i1, ps, i2, eb = row['ExonA'], row['Intron1'], row['Pseudoexon'], row['Intron2'], row['ExonB']
        label = float(row['label'])
        
        if self.is_training:
            # --- AUGMENTATION STACK (From V11) ---
            
            # 1. Franken-Swap (Positives Only)
            if label == 1.0:
                if len(self.positive_intron1_middles) > 0 and len(i1) > 100 and np.random.rand() < 0.50:
                    new_mid = random.choice(self.positive_intron1_middles)
                    i1 = i1[:50] + new_mid + i1[-50:]
                if len(self.positive_intron2_middles) > 0 and len(i2) > 100 and np.random.rand() < 0.50:
                    new_mid = random.choice(self.positive_intron2_middles)
                    i2 = i2[:50] + new_mid + i2[-50:]

            # 2. Jitter
            crop_a = np.random.randint(0, 16)
            if len(ea) > crop_a: ea = ea[crop_a:]
            crop_b = np.random.randint(0, 16)
            if len(eb) > crop_b: eb = eb[:-crop_b]

            # 3. Adversarial Decoys
            if label == 1.0:
                rand_val = np.random.rand()
                if rand_val < 0.15: # Destroy GT
                     i2 = "CC" + i2[2:] 
                     label = 0.0 
                elif rand_val < 0.30: # Scramble
                    if len(i1) > 2: i1 = self.scramble_string(i1[:-2]) + i1[-2:]
                    ps = self.scramble_string(ps)
                    if len(i2) > 2: i2 = i2[:2] + self.scramble_string(i2[2:])
                    label = 0.0

            # 4. Masking
            full_seq = ea + i1 + ps + i2 + eb
            full_seq = self.apply_masking(full_seq, mask_prob=0.10)
            return full_seq, torch.tensor([label], dtype=torch.float32)

        # Validation
        full_seq = ea + i1 + ps + i2 + eb
        return full_seq, torch.tensor([label], dtype=torch.float32)

def collate_fn(batch):
    seqs, labels = zip(*batch)
    X = vectorized_one_hot_encode(seqs)
    y = torch.stack(labels)
    return X, y

# --- 4. MODEL ARCHITECTURE (Standard V9/V10/V11) ---

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
        
        if self.training and x.requires_grad:
            x = checkpoint.checkpoint(self._attn_block, x, use_reentrant=False)
        else:
            x = self._attn_block(x)
            
        x = x.permute(0, 2, 1)
        x = self.global_pool(x)
        x = self.flatten(x)
        x = self.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.fc2(x)
        return x

# --- 5. TRAINING ENGINE ---

def evaluate(model, loader):
    model.eval()
    total_loss = 0
    all_preds = []
    all_labels = []
    criterion = nn.BCEWithLogitsLoss()
    
    with torch.inference_mode():
        for inputs, labels in loader:
            inputs, labels = inputs.to(device), labels.to(device)
            
            with autocast(device_type='cuda'):
                outputs = model(inputs)
                loss = criterion(outputs, labels)
                
            total_loss += loss.item()
            predicted = (torch.sigmoid(outputs) > 0.5).float()
            
            all_preds.extend(predicted.cpu().numpy())
            all_labels.extend(labels.cpu().numpy())
            
    precision = precision_score(all_labels, all_preds, zero_division=0)
    recall = recall_score(all_labels, all_preds, zero_division=0)
    
    all_preds = np.array(all_preds)
    all_labels = np.array(all_labels)
    acc = np.mean(all_preds == all_labels)
            
    return total_loss / len(loader), acc, precision, recall

def train_model():
    log_message(f"--- Starting Hybrid V13 (FOCAL LOSS) Training on {torch.cuda.get_device_name(0)} ---")
    log_message(f"Strategy: Focal Loss (Gamma=2) to mine hard positives.")
    
    # 1. Load Data
    log_message("Loading TSV...")
    df = pd.read_csv(INPUT_FILE, sep='\t')
    df = df.dropna(subset=['ExonA', 'Intron1', 'Pseudoexon', 'Intron2', 'ExonB'])
    
    train_df, val_df = train_test_split(df, test_size=0.1, random_state=42, stratify=df['label'])
    
    # Calculate Weights
    train_labels = train_df['label'].values.astype(int)
    class_counts = np.bincount(train_labels)
    weight_per_class = 1.0 / class_counts
    sample_weights = weight_per_class[train_labels]
    sample_weights = torch.DoubleTensor(sample_weights)
    
    sampler = WeightedRandomSampler(weights=sample_weights, num_samples=len(sample_weights), replacement=True)
    
    train_ds = ComponentDataset(train_df, is_training=True)
    val_ds = ComponentDataset(val_df, is_training=False)
    
    train_loader = DataLoader(train_ds, batch_size=MICRO_BATCH_SIZE, sampler=sampler, shuffle=False, 
                              num_workers=4, collate_fn=collate_fn, pin_memory=True)
    val_loader = DataLoader(val_ds, batch_size=MICRO_BATCH_SIZE*2, shuffle=False, 
                            num_workers=4, collate_fn=collate_fn, pin_memory=True)

    # 2. Setup Model & Optimizers
    model = RisdiplamHybrid().to(device)
    optimizer = optim.AdamW(model.parameters(), lr=LEARNING_RATE, weight_decay=1e-3)
    scaler = GradScaler('cuda') 
    
    # --- KEY CHANGE: Use Focal Loss ---
    criterion = FocalLoss(gamma=2.0)
    
    best_f1 = 0.0
    
    # 3. Training Loop
    for epoch in range(EPOCHS):
        model.train()
        running_loss = 0.0
        optimizer.zero_grad() 
        
        pbar = tqdm(train_loader, desc=f"Epoch {epoch+1}/{EPOCHS}")
        
        for i, (inputs, labels) in enumerate(pbar):
            inputs, labels = inputs.to(device, non_blocking=True), labels.to(device, non_blocking=True)
            
            with autocast(device_type='cuda'):
                outputs = model(inputs)
                loss = criterion(outputs, labels)
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
            
        # Validation
        val_loss, val_acc, val_prec, val_rec = evaluate(model, val_loader)
        
        if (val_prec + val_rec) > 0: f1 = 2 * (val_prec * val_rec) / (val_prec + val_rec)
        else: f1 = 0.0
        
        log_msg = f"Epoch {epoch+1}: Loss={val_loss:.4f} | Prec={val_prec:.4f} | Rec={val_rec:.4f} | F1={f1:.4f}"
        
        if f1 > best_f1:
            best_f1 = f1
            torch.save(model.state_dict(), MODEL_SAVE_PATH)
            log_msg += " [SAVED BEST]"
            
        log_message(log_msg)

if __name__ == "__main__":
    train_model()