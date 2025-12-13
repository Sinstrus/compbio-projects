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

# --- Configuration ---
VERSION = "v8_components"
INPUT_FILE = "dataset_components_v8.tsv"
MODEL_SAVE_PATH = f"deep_risdiplam_model_{VERSION}.pth"
LOG_FILE = f"training_log_{VERSION}.txt"

BATCH_SIZE = 96
EPOCHS = 30
MAX_SEQ_LEN = 1500
LEARNING_RATE = 0.0001
STEPS_PER_EPOCH = 100
NUM_WORKERS = 8  # Safe to use workers now

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def log_message(message):
    print(message)
    with open(LOG_FILE, "a") as f:
        f.write(message + "\n")

def vectorized_one_hot_encode(seqs, max_len=MAX_SEQ_LEN, random_pad=False):
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

class ComponentDataset(IterableDataset):
    def __init__(self, data_df, batch_size):
        self.data = data_df
        self.positives = data_df[data_df['label'] == 1].to_dict('records')
        self.negatives = data_df[data_df['label'] == 0].to_dict('records')
        
        self.batch_size = batch_size
        self.n_pos = batch_size // 2
        remaining = batch_size - self.n_pos
        self.n_neg_real = remaining // 3
        self.n_neg_dead = remaining // 3
        self.n_neg_scram = remaining - (self.n_neg_real + self.n_neg_dead)

    def mutate_intron2(self, seq_i2):
        # Guaranteed to be at start because it's the component
        return "CC" + seq_i2[2:]

    def scramble_pseudo(self, seq_ps):
        if len(seq_ps) < 6: return seq_ps
        # Protect splice sites (first 3, last 3)
        mid = list(seq_ps[3:-3])
        random.shuffle(mid)
        return seq_ps[:3] + "".join(mid) + seq_ps[-3:]

    def stitch(self, row, mutate_type=None):
        ea = row['ExonA']
        i1 = row['Intron1']
        ps = row['Pseudoexon']
        i2 = row['Intron2']
        eb = row['ExonB']
        
        if mutate_type == "DEAD":
            i2 = self.mutate_intron2(i2)
        elif mutate_type == "SCRAM":
            ps = self.scramble_pseudo(ps)
            
        # Stitch
        full_seq = ea + i1 + ps + i2 + eb
        
        # Jitter (Trim start/end)
        if len(full_seq) > 200:
            trim_start = random.randint(0, 15)
            trim_end = random.randint(0, 15)
            if len(full_seq) - (trim_start + trim_end) > 100:
                full_seq = full_seq[trim_start : len(full_seq)-trim_end]
                
        return full_seq

    def __iter__(self):
        while True:
            # Sample Rows (Dictionaries)
            batch_pos_rows = np.random.choice(self.positives, self.n_pos)
            batch_neg_real_rows = np.random.choice(self.negatives, self.n_neg_real)
            
            # For synthetics, pick random positives to mutate
            batch_dead_rows = np.random.choice(self.positives, self.n_neg_dead)
            batch_scram_rows = np.random.choice(self.positives, self.n_neg_scram)
            
            seqs = []
            
            # 1. Positives
            for r in batch_pos_rows: seqs.append(self.stitch(r))
            # 2. Real Negatives
            for r in batch_neg_real_rows: seqs.append(self.stitch(r))
            # 3. Dead
            for r in batch_dead_rows: seqs.append(self.stitch(r, "DEAD"))
            # 4. Scrambled
            for r in batch_scram_rows: seqs.append(self.stitch(r, "SCRAM"))
            
            X = vectorized_one_hot_encode(seqs, random_pad=True)
            y = torch.zeros(len(seqs), 1)
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
        self.dropout = nn.Dropout(0.6)
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
    avg_pred = 0.0
    criterion = nn.BCEWithLogitsLoss()
    val_iter = iter(val_loader)
    with torch.no_grad():
        for _ in range(steps):
            try:
                inputs, labels = next(val_iter)
            except StopIteration: break
            inputs, labels = inputs.to(device), labels.to(device)
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            total_loss += loss.item()
            probs = torch.sigmoid(outputs)
            avg_pred += probs.mean().item()
            predicted = (probs > 0.5).float()
            total += labels.size(0)
            correct += (predicted == labels).sum().item()
    return total_loss / steps, correct / total, avg_pred / steps

def main():
    with open(LOG_FILE, "w") as f:
        f.write(f"DeepRisdiplam Component Training - {VERSION}\n")
    
    print("Loading Components...")
    try:
        df = pd.read_csv(INPUT_FILE, sep='\t')
        print(f"Loaded {len(df)} samples.")
    except Exception as e:
        print(f"Error: {e}")
        return

    # Split by ID to avoid data leakage
    # We create a column for splitting
    ids = df['event_id'].unique()
    train_ids, val_ids = train_test_split(ids, test_size=0.2, random_state=42)
    
    train_df = df[df['event_id'].isin(train_ids)]
    val_df = df[df['event_id'].isin(val_ids)]
    
    print(f"Train size: {len(train_df)} | Val size: {len(val_df)}")
    
    train_ds = ComponentDataset(train_df, BATCH_SIZE)
    val_ds = ComponentDataset(val_df, BATCH_SIZE)
    
    train_loader = DataLoader(train_ds, batch_size=None, num_workers=NUM_WORKERS, prefetch_factor=2)
    val_loader = DataLoader(val_ds, batch_size=None, num_workers=2)

    model = RisdiplamModel().to(device)
    criterion = nn.BCEWithLogitsLoss()
    optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE, weight_decay=1e-3)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='max', factor=0.5, patience=5, verbose=True)
    
    best_val_acc = 0.0
    train_iter = iter(train_loader)
    
    for epoch in range(EPOCHS):
        model.train()
        running_loss = 0.0
        correct = 0
        total = 0
        pbar = tqdm(range(STEPS_PER_EPOCH), desc=f"Epoch {epoch+1}")
        
        for _ in pbar:
            inputs, labels = next(train_iter)
            inputs, labels = inputs.to(device, non_blocking=True), labels.to(device, non_blocking=True)
            
            smooth_labels = labels * 0.9 + 0.05
            
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, smooth_labels)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=0.5)
            optimizer.step()
            
            running_loss += loss.item()
            predicted = (torch.sigmoid(outputs) > 0.5).float()
            total += labels.size(0)
            correct += (predicted == labels).sum().item()
            pbar.set_postfix({'loss': f"{loss.item():.4f}", 'acc': f"{correct/total:.4f}"})
            
        val_loss, val_acc, val_avg = evaluate(model, val_loader)
        scheduler.step(val_acc)
        log_message(f"Epoch {epoch+1}: Loss={val_loss:.4f} | Acc={val_acc:.4f} | AvgPred={val_avg:.2f}")
        
        if val_acc > best_val_acc:
            best_val_acc = val_acc
            torch.save(model.state_dict(), MODEL_SAVE_PATH)
            log_message("   Saved new best model.")

if __name__ == "__main__":
    main()