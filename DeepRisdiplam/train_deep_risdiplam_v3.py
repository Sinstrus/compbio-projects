import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, IterableDataset
from sklearn.model_selection import train_test_split
import random
import sys
import os
from tqdm import tqdm
import datetime

# --- Version Control ---
VERSION = "v3"  # Updated Version

# --- Configuration ---
INPUT_X_FILE = "dataset_X_sequences.tsv"
INPUT_Y_FILE = "dataset_y_labels.tsv"

MODEL_SAVE_PATH = f"deep_risdiplam_model_{VERSION}.pth"
LOG_FILE = f"training_log_{VERSION}.txt"

# Structure Constants
OFFSET_END_FROM_RIGHT = 550 

# Training Params
BATCH_SIZE = 96 
EPOCHS = 50
MAX_SEQ_LEN = 1500 
LEARNING_RATE = 0.0005  # Lowered starting LR (was 0.001)
STEPS_PER_EPOCH = 100 
NUM_WORKERS = 8 

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def log_message(message):
    print(message)
    with open(LOG_FILE, "a") as f:
        f.write(message + "\n")

# --- Optimized Helper Functions ---

def vectorized_one_hot_encode(seqs, max_len=MAX_SEQ_LEN, random_pad=False):
    """Vectorized one-hot encoding."""
    batch_size = len(seqs)
    mapping = np.zeros(128, dtype=np.int8) + 4 
    mapping[ord('A')] = 0; mapping[ord('a')] = 0
    mapping[ord('C')] = 1; mapping[ord('c')] = 1
    mapping[ord('G')] = 2; mapping[ord('g')] = 2
    mapping[ord('T')] = 3; mapping[ord('t')] = 3
    
    tensor_batch = np.zeros((batch_size, 4, max_len), dtype=np.float32)
    
    for i, seq in enumerate(seqs):
        if len(seq) > max_len:
            start = (len(seq) - max_len) // 2
            seq = seq[start : start + max_len]
            
        indices = mapping[np.array([ord(c) for c in seq])]
        
        seq_len = len(indices)
        pad_total = max_len - seq_len
        pad_left = 0
        if pad_total > 0 and random_pad:
            pad_left = np.random.randint(0, pad_total + 1)
            
        for k in range(seq_len):
            idx = indices[k]
            if idx < 4:
                tensor_batch[i, idx, pad_left + k] = 1.0
                
    return torch.from_numpy(tensor_batch)

class BalancedIterableDataset(IterableDataset):
    """Produces balanced batches: 50% Pos, 50% Neg (Real/Dead/Scrambled)."""
    def __init__(self, x_pos, x_neg, batch_size):
        self.x_pos = x_pos
        self.x_neg = x_neg
        self.batch_size = batch_size
        self.n_pos = batch_size // 2
        remaining = batch_size - self.n_pos
        self.n_neg_real = remaining // 3
        self.n_neg_dead = remaining // 3
        self.n_neg_scram = remaining - (self.n_neg_real + self.n_neg_dead)

    def augment_seq(self, seq):
        if len(seq) < 200: return seq
        trim_start = random.randint(0, 15)
        trim_end = random.randint(0, 15)
        if len(seq) - (trim_start + trim_end) < 100: return seq
        return seq[trim_start : len(seq)-trim_end]

    def mutate_dead_donor(self, seq):
        idx_donor = len(seq) - OFFSET_END_FROM_RIGHT
        return seq[:idx_donor] + "CC" + seq[idx_donor+2:]

    def scramble_internal(self, seq):
        start_region = 550
        end_region = len(seq) - 550
        protect_start = start_region - 20
        protect_end_acceptor = start_region + 3
        protect_start_donor = end_region - 3
        s_idx = protect_end_acceptor
        e_idx = protect_start_donor
        if e_idx <= s_idx: return seq 
        subseq = list(seq[s_idx:e_idx])
        random.shuffle(subseq)
        return seq[:s_idx] + "".join(subseq) + seq[e_idx:]

    def __iter__(self):
        while True:
            batch_pos = np.random.choice(self.x_pos, self.n_pos)
            batch_neg_real = np.random.choice(self.x_neg, self.n_neg_real)
            batch_base_for_dead = np.random.choice(self.x_pos, self.n_neg_dead)
            batch_base_for_scram = np.random.choice(self.x_pos, self.n_neg_scram)
            
            batch_dead = [self.mutate_dead_donor(s) for s in batch_base_for_dead]
            batch_scram = [self.scramble_internal(s) for s in batch_base_for_scram]
            
            raw_seqs = np.concatenate([batch_pos, batch_neg_real, batch_dead, batch_scram])
            aug_seqs = [self.augment_seq(s) for s in raw_seqs]
            X = vectorized_one_hot_encode(aug_seqs, random_pad=True)
            
            y = torch.zeros(len(raw_seqs), 1)
            y[:self.n_pos] = 1.0
            yield X, y

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
        x = self.pool(self.relu(self.bn1(self.conv1(x))))
        x = self.pool(self.relu(self.bn2(self.conv2(x))))
        x = self.relu(self.bn3(self.conv3(x)))
        x = self.fc2(self.dropout(self.relu(self.fc1(self.flatten(self.global_pool(x))))))
        return x

def evaluate(model, val_loader, steps=20):
    model.eval()
    total_loss = 0
    correct = 0
    total = 0
    criterion = nn.BCEWithLogitsLoss()
    val_iter = iter(val_loader)
    
    with torch.no_grad():
        for _ in range(steps):
            try:
                inputs, labels = next(val_iter)
            except StopIteration:
                break
            inputs, labels = inputs.to(device), labels.to(device)
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            total_loss += loss.item()
            predicted = (torch.sigmoid(outputs) > 0.5).float()
            total += labels.size(0)
            correct += (predicted == labels).sum().item()
    return total_loss / steps, correct / total

def main():
    with open(LOG_FILE, "w") as f:
        f.write(f"DeepRisdiplam Training Log - {VERSION}\n")

    log_message(f"Starting V3 Training (Stable). Device: {device}")
    
    # 1. Load Data
    try:
        df_x = pd.read_csv(INPUT_X_FILE, sep='\t')
        df_y = pd.read_csv(INPUT_Y_FILE, sep='\t')
        df = pd.merge(df_x, df_y, on="event_id")
    except Exception as e:
        log_message(f"Error: {e}")
        return

    positives = df[df['label'] == 1]['minigene_sequence'].values
    negatives = df[df['label'] == 0]['minigene_sequence'].values
    
    pos_train, pos_val = train_test_split(positives, test_size=0.2, random_state=42)
    neg_train, neg_val = train_test_split(negatives, test_size=0.2, random_state=42)
    
    # 2. Datasets
    train_ds = BalancedIterableDataset(pos_train, neg_train, BATCH_SIZE)
    val_ds = BalancedIterableDataset(pos_val, neg_val, BATCH_SIZE)
    
    train_loader = DataLoader(train_ds, batch_size=None, num_workers=NUM_WORKERS, prefetch_factor=2)
    val_loader = DataLoader(val_ds, batch_size=None, num_workers=2)

    # 3. Model & Optimizer
    model = RisdiplamModel().to(device)
    criterion = nn.BCEWithLogitsLoss()
    
    # FIX 1: Weight Decay for regularization
    optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE, weight_decay=1e-4)
    
    # FIX 2: Scheduler
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='max', factor=0.5, patience=5, verbose=True)
    
    best_val_acc = 0.0
    train_iter = iter(train_loader)
    
    for epoch in range(EPOCHS):
        model.train()
        running_loss = 0.0
        correct = 0
        total = 0
        
        pbar = tqdm(range(STEPS_PER_EPOCH), desc=f"Epoch {epoch+1}/{EPOCHS}")
        
        for _ in pbar:
            inputs, labels = next(train_iter)
            inputs, labels = inputs.to(device, non_blocking=True), labels.to(device, non_blocking=True)
            
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            loss.backward()
            
            # FIX 3: Gradient Clipping (The most important fix for explosions)
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            
            optimizer.step()
            
            running_loss += loss.item()
            predicted = (torch.sigmoid(outputs) > 0.5).float()
            total += labels.size(0)
            correct += (predicted == labels).sum().item()
            
            pbar.set_postfix({'loss': f"{loss.item():.4f}", 'acc': f"{correct/total:.4f}"})
            
        val_loss, val_acc = evaluate(model, val_loader)
        
        # Scheduler Step
        scheduler.step(val_acc)
        
        log_msg = f"Epoch {epoch+1}: Val Loss={val_loss:.4f} | Val Acc={val_acc:.4f}"
        log_message(log_msg)
        
        if val_acc > best_val_acc:
            best_val_acc = val_acc
            torch.save(model.state_dict(), MODEL_SAVE_PATH)
            log_message(f"   Saved new best model.")
            
    log_message("\nDone!")

if __name__ == "__main__":
    main()