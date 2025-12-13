import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
import torch.utils.checkpoint as checkpoint
from torch.cuda.amp import autocast, GradScaler
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import train_test_split
from tqdm import tqdm
import sys
import os

# --- 1. CONFIGURATION ---
VERSION = "v9_hybrid_optimized"
INPUT_FILE = "dataset_components_v8.tsv"
MODEL_SAVE_PATH = f"deep_risdiplam_model_{VERSION}.pth"
LOG_FILE = f"training_log_{VERSION}.txt"

# VRAM OPTIMIZATION SETTINGS
# RTX 3070 Ti (8GB) Profile
MAX_SEQ_LEN = 1500      # Max length after stitching
MICRO_BATCH_SIZE = 32   # Actual samples on GPU (Low to prevent OOM)
ACCUM_STEPS = 4         # How many batches to sum (32 * 4 = 128 effective)
EFFECTIVE_BATCH_SIZE = MICRO_BATCH_SIZE * ACCUM_STEPS
LEARNING_RATE = 1e-4
EPOCHS = 25

# Enable TensorFloat-32 (Free speedup on Ampere GPUs)
torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.allow_tf32 = True

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# --- 2. DATA PIPELINE ---

def log_message(message):
    print(message)
    with open(LOG_FILE, "a") as f:
        f.write(message + "\n")

def vectorized_one_hot_encode(seqs, max_len=MAX_SEQ_LEN):
    """
    Fast vectorized one-hot encoding.
    """
    batch_size = len(seqs)
    mapping = np.zeros(128, dtype=np.int8) + 4 # Default to 4 (padding)
    mapping[ord('A')] = 0; mapping[ord('a')] = 0
    mapping[ord('C')] = 1; mapping[ord('c')] = 1
    mapping[ord('G')] = 2; mapping[ord('g')] = 2
    mapping[ord('T')] = 3; mapping[ord('t')] = 3
    
    # Pre-allocate buffer
    arr = np.zeros((batch_size, max_len), dtype=np.int8) + 4
    
    for i, seq in enumerate(seqs):
        l = min(len(seq), max_len)
        # Convert string to ascii codes -> map to 0-3
        arr[i, :l] = mapping[np.frombuffer(seq[:l].encode('latin1'), dtype=np.uint8)]
        
    # Turn into (Batch, 4, Length) float tensor
    one_hot = np.zeros((batch_size, 4, max_len), dtype=np.float32)
    
    # Fancy indexing to fill the 1s
    # We only fill where the value is 0,1,2,3 (skip 4)
    for base_val in range(4):
        one_hot[:, base_val, :] = (arr == base_val)
        
    return torch.tensor(one_hot)

class ComponentDataset(Dataset):
    def __init__(self, df, is_training=False):
        self.df = df
        self.is_training = is_training
        
    def __len__(self):
        return len(self.df)
    
    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        
        # Stitch components
        ea, i1, ps, i2, eb = row['ExonA'], row['Intron1'], row['Pseudoexon'], row['Intron2'], row['ExonB']
        
        # Training Augmentation: Create Artificial Negatives
        # If label is 1, sometimes mutate the GT donor to simulate a dead site
        label = float(row['label'])
        
        if self.is_training and label == 1.0 and np.random.rand() < 0.15:
             # Artificial Negative: Destroy splice site
             i2 = "CC" + i2[2:] # Mutate GT -> CC
             label = 0.0 # Flip label
             
        full_seq = ea + i1 + ps + i2 + eb
        return full_seq, torch.tensor([label], dtype=torch.float32)

def collate_fn(batch):
    seqs, labels = zip(*batch)
    X = vectorized_one_hot_encode(seqs)
    y = torch.stack(labels)
    return X, y

# --- 3. HYBRID MODEL ARCHITECTURE ---

class RisdiplamHybrid(nn.Module):
    def __init__(self):
        super(RisdiplamHybrid, self).__init__()
        
        # --- CNN Feature Extractor ---
        # Reduces sequence length and finds local motifs
        self.conv1 = nn.Conv1d(4, 64, 11, padding='same')
        self.bn1 = nn.BatchNorm1d(64)
        self.relu = nn.ReLU()
        self.pool = nn.MaxPool1d(2) # Downsample 1500 -> 750
        
        self.conv2 = nn.Conv1d(64, 128, 11, padding='same', dilation=2)
        self.bn2 = nn.BatchNorm1d(128)
        
        # --- Transformer Context Engine ---
        self.d_model = 128
        self.nhead = 4
        # Batch First: (Batch, Seq, Feature)
        self.attention = nn.MultiheadAttention(embed_dim=self.d_model, num_heads=self.nhead, batch_first=True)
        
        # --- Classifier ---
        self.global_pool = nn.AdaptiveAvgPool1d(1)
        self.flatten = nn.Flatten()
        
        self.fc1 = nn.Linear(128, 64)
        self.dropout = nn.Dropout(0.5) 
        self.fc2 = nn.Linear(64, 1)

    def _attn_block(self, x):
        """Helper for checkpointing: self-attention + residual"""
        attn_out, _ = self.attention(x, x, x)
        return x + attn_out

    def forward(self, x):
        # 1. CNN (Local Features)
        x = self.relu(self.bn1(self.conv1(x)))
        x = self.pool(x) # 1500 -> 750 (Saves 4x Attention Memory)
        x = self.relu(self.bn2(self.conv2(x)))
        
        # 2. Transformer (Long Range Context)
        # Permute to (Batch, Seq, Channels) for Attention
        x = x.permute(0, 2, 1) 
        
        # CHECKPOINTING: The "Infinite VRAM" Trick
        # We don't save the massive attention map. We recompute it during backprop.
        if self.training and x.requires_grad:
            x = checkpoint.checkpoint(self._attn_block, x, use_reentrant=False)
        else:
            x = self._attn_block(x)
            
        # Permute back to (Batch, Channels, Seq) for Pooling
        x = x.permute(0, 2, 1)
        
        # 3. Classification
        x = self.global_pool(x)
        x = self.flatten(x)
        x = self.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.fc2(x)
        return x

# --- 4. TRAINING ENGINE ---

def evaluate(model, loader):
    model.eval()
    total_loss = 0
    correct = 0
    total = 0
    criterion = nn.BCEWithLogitsLoss()
    
    # Use inference_mode for max speed (no grad tracking)
    with torch.inference_mode():
        for inputs, labels in loader:
            inputs, labels = inputs.to(device), labels.to(device)
            
            # Use AMP for inference too
            with autocast():
                outputs = model(inputs)
                loss = criterion(outputs, labels)
                
            total_loss += loss.item()
            predicted = (torch.sigmoid(outputs) > 0.5).float()
            total += labels.size(0)
            correct += (predicted == labels).sum().item()
            
    return total_loss / len(loader), correct / total

def train_model():
    log_message(f"--- Starting Hybrid V9 Training on {torch.cuda.get_device_name(0)} ---")
    log_message(f"Config: MicroBatch={MICRO_BATCH_SIZE} | AccumSteps={ACCUM_STEPS} | EffectiveBatch={EFFECTIVE_BATCH_SIZE}")
    log_message(f"Optimizations: AMP + Checkpointing + Gradient Accumulation + TF32")
    
    # 1. Load Data
    log_message("Loading TSV...")
    df = pd.read_csv(INPUT_FILE, sep='\t')
    
    # Filter out empty strings if any
    df = df.dropna(subset=['ExonA', 'Intron1', 'Pseudoexon', 'Intron2', 'ExonB'])
    
    train_df, val_df = train_test_split(df, test_size=0.1, random_state=42, stratify=df['label'])
    
    log_message(f"Train: {len(train_df)} | Val: {len(val_df)}")
    
    train_ds = ComponentDataset(train_df, is_training=True)
    val_ds = ComponentDataset(val_df, is_training=False)
    
    train_loader = DataLoader(train_ds, batch_size=MICRO_BATCH_SIZE, shuffle=True, 
                              num_workers=4, collate_fn=collate_fn, pin_memory=True)
    val_loader = DataLoader(val_ds, batch_size=MICRO_BATCH_SIZE*2, shuffle=False, 
                            num_workers=4, collate_fn=collate_fn, pin_memory=True)

    # 2. Setup Model & Optimizers
    model = RisdiplamHybrid().to(device)
    optimizer = optim.AdamW(model.parameters(), lr=LEARNING_RATE, weight_decay=1e-4)
    criterion = nn.BCEWithLogitsLoss()
    scaler = GradScaler() # For AMP
    
    best_acc = 0.0
    
    # 3. Training Loop
    for epoch in range(EPOCHS):
        model.train()
        running_loss = 0.0
        optimizer.zero_grad() # Reset gradients at start of epoch
        
        pbar = tqdm(train_loader, desc=f"Epoch {epoch+1}/{EPOCHS}")
        
        for i, (inputs, labels) in enumerate(pbar):
            inputs, labels = inputs.to(device, non_blocking=True), labels.to(device, non_blocking=True)
            
            # --- THE OPTIMIZED FORWARD PASS ---
            with autocast():
                outputs = model(inputs)
                loss = criterion(outputs, labels)
                # Normalize loss because we sum gradients over ACCUM_STEPS
                loss = loss / ACCUM_STEPS 
            
            # --- BACKWARD PASS ---
            scaler.scale(loss).backward()
            
            # --- OPTIMIZER STEP (Delayed) ---
            if (i + 1) % ACCUM_STEPS == 0:
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                
                scaler.step(optimizer)
                scaler.update()
                optimizer.zero_grad()
                
            running_loss += loss.item() * ACCUM_STEPS # Scale back up for logging
            
            pbar.set_postfix({'loss': f"{running_loss / (i+1):.4f}"})
            
        # Validation
        val_loss, val_acc = evaluate(model, val_loader)
        
        log_msg = f"Epoch {epoch+1}: Loss={val_loss:.4f} | Acc={val_acc:.4f}"
        if val_acc > best_acc:
            best_acc = val_acc
            torch.save(model.state_dict(), MODEL_SAVE_PATH)
            log_msg += " [SAVED BEST]"
            
        log_message(log_msg)

if __name__ == "__main__":
    train_model()