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
VERSION = "debug_synth_v2"
INPUT_X_FILE = "dataset_X_sequences_strict.tsv"
INPUT_Y_FILE = "dataset_y_labels_strict.tsv"
OFFSET_END_FROM_RIGHT = 550 
BATCH_SIZE = 96 
EPOCHS = 10 
MAX_SEQ_LEN = 1500 
LEARNING_RATE = 0.0001
STEPS_PER_EPOCH = 50

# FIX 1: Disable Multiprocessing to stop OSError/Crashing
NUM_WORKERS = 0 

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def vectorized_one_hot_encode(seqs, max_len=MAX_SEQ_LEN, random_pad=False):
    batch_size = len(seqs)
    mapping = np.zeros(128, dtype=np.int8) + 4 
    mapping[ord('A')] = 0; mapping[ord('a')] = 0
    mapping[ord('C')] = 1; mapping[ord('c')] = 1
    mapping[ord('G')] = 2; mapping[ord('g')] = 2
    mapping[ord('T')] = 3; mapping[ord('t')] = 3
    tensor_batch = np.zeros((batch_size, 4, max_len), dtype=np.float32)
    
    for i, seq in enumerate(seqs):
        # FIX 2: Restore Truncation Logic to prevent crashes
        if len(seq) > max_len:
            start = (len(seq) - max_len) // 2
            seq = seq[start : start + max_len]
            
        indices = mapping[np.array([ord(c) for c in seq])]
        
        # Calculate padding safely
        seq_len = len(indices)
        pad_total = max_len - seq_len
        pad_left = 0
        
        # FIX 3: Safe Random Int
        if pad_total > 0 and random_pad:
            pad_left = np.random.randint(0, pad_total + 1)
            
        for k in range(seq_len):
            idx = indices[k]
            if idx < 4: 
                pos = pad_left + k
                if pos < max_len:
                    tensor_batch[i, idx, pos] = 1.0
                    
    return torch.from_numpy(tensor_batch)

class DebugIterableDataset(IterableDataset):
    def __init__(self, x_pos, batch_size):
        self.x_pos = x_pos
        self.batch_size = batch_size
        self.n_pos = batch_size // 2
        remaining = batch_size - self.n_pos
        self.n_neg_dead = remaining // 2
        self.n_neg_scram = remaining - self.n_neg_dead

    def augment_seq(self, seq):
        if len(seq) < 200: return seq
        trim_start = random.randint(0, 15)
        trim_end = random.randint(0, 15)
        if len(seq) - (trim_start + trim_end) < 100: return seq
        return seq[trim_start : len(seq)-trim_end]

    def mutate_dead_donor(self, seq):
        idx_donor = len(seq) - OFFSET_END_FROM_RIGHT
        if idx_donor < 0 or idx_donor >= len(seq)-2: return seq
        return seq[:idx_donor] + "CC" + seq[idx_donor+2:]

    def scramble_internal(self, seq):
        start_region = 550
        end_region = len(seq) - 550
        if start_region >= end_region: return seq
        s_idx = start_region - 20 + 3
        e_idx = end_region - 3
        if s_idx >= e_idx: return seq
        subseq = list(seq[s_idx:e_idx])
        random.shuffle(subseq)
        return seq[:s_idx] + "".join(subseq) + seq[e_idx:]

    def __iter__(self):
        while True:
            batch_pos = np.random.choice(self.x_pos, self.n_pos)
            batch_base_dead = np.random.choice(self.x_pos, self.n_neg_dead)
            batch_base_scram = np.random.choice(self.x_pos, self.n_neg_scram)
            
            batch_dead = [self.mutate_dead_donor(s) for s in batch_base_dead]
            batch_scram = [self.scramble_internal(s) for s in batch_base_scram]
            
            raw_seqs = np.concatenate([batch_pos, batch_dead, batch_scram])
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
        self.fc2 = nn.Linear(64, 1) 
    def forward(self, x):
        x = self.pool(self.relu(self.bn1(self.conv1(x))))
        x = self.pool(self.relu(self.bn2(self.conv2(x))))
        x = self.relu(self.bn3(self.conv3(x)))
        x = self.fc2(self.relu(self.fc1(self.flatten(self.global_pool(x)))))
        return x

def evaluate(model, val_loader, steps=20):
    model.eval()
    correct = 0
    total = 0
    val_iter = iter(val_loader)
    with torch.no_grad():
        for _ in range(steps):
            try:
                inputs, labels = next(val_iter)
            except StopIteration:
                break
            inputs, labels = inputs.to(device), labels.to(device)
            outputs = model(inputs)
            predicted = (torch.sigmoid(outputs) > 0.5).float()
            correct += (predicted == labels).sum().item()
            total += labels.size(0)
    return correct / total

def main():
    print(f"--- Debug Training: Positives vs Synthetics ONLY ---")
    try:
        df_x = pd.read_csv(INPUT_X_FILE, sep='\t')
        df_y = pd.read_csv(INPUT_Y_FILE, sep='\t')
        df = pd.merge(df_x, df_y, on="event_id")
    except Exception: return

    positives = df[df['label'] == 1]['minigene_sequence'].values
    if len(positives) == 0:
        print("Error: No positive samples found.")
        return
        
    pos_train, pos_val = train_test_split(positives, test_size=0.2, random_state=42)
    
    # Debug Dataset uses ONLY positives to generate synthetic negatives
    train_ds = DebugIterableDataset(pos_train, BATCH_SIZE)
    val_ds = DebugIterableDataset(pos_val, BATCH_SIZE)
    train_loader = DataLoader(train_ds, batch_size=None, num_workers=NUM_WORKERS)
    val_loader = DataLoader(val_ds, batch_size=None, num_workers=0) # Fix workers here too

    model = RisdiplamModel().to(device)
    criterion = nn.BCEWithLogitsLoss()
    optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE)
    
    for epoch in range(EPOCHS):
        model.train()
        pbar = tqdm(range(STEPS_PER_EPOCH), desc=f"Epoch {epoch+1}")
        for _ in pbar:
            try:
                inputs, labels = next(iter(train_loader))
            except StopIteration:
                break
            inputs, labels = inputs.to(device), labels.to(device)
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            
        acc = evaluate(model, val_loader)
        print(f"Epoch {epoch+1} Val Accuracy: {acc:.4f}")

if __name__ == "__main__":
    main()