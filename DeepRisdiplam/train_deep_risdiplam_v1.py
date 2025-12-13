import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import train_test_split
import random
import sys
import os
from tqdm import tqdm
import datetime

# --- Version Control ---
VERSION = "v1"

# --- Configuration ---
INPUT_X_FILE = "dataset_X_sequences.tsv"
INPUT_Y_FILE = "dataset_y_labels.tsv"

# Versioned Outputs
MODEL_SAVE_PATH = f"deep_risdiplam_model_{VERSION}.pth"
LOG_FILE = f"training_log_{VERSION}.txt"

# Structure Constants (Must match digital_cloner.py)
EXON_A_LEN = 50
INTRON_1_PART_LEN = 250
INTRON_2_PART_LEN = 250
EXON_B_LEN = 50
OFFSET_START_PSEUDO = EXON_A_LEN + (INTRON_1_PART_LEN * 2) 
OFFSET_END_FROM_RIGHT = INTRON_2_PART_LEN * 2 + EXON_B_LEN 

# Training Params
BATCH_SIZE = 96 
EPOCHS = 50
MAX_SEQ_LEN = 1500 
LEARNING_RATE = 0.001
STEPS_PER_EPOCH = 100 

# Device Configuration
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def log_message(message):
    """Prints to console and appends to log file."""
    print(message)
    with open(LOG_FILE, "a") as f:
        f.write(message + "\n")

# --- Helper Functions ---

def one_hot_encode(seq, max_len=MAX_SEQ_LEN, random_pad=False):
    """
    One-hot encodes DNA sequence to tensor (C, L) format for PyTorch Conv1d.
    Shape: (4, MAX_SEQ_LEN)
    """
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3}
    
    # Truncate if too long
    if len(seq) > max_len:
        diff = len(seq) - max_len
        start_trim = diff // 2
        seq = seq[start_trim : start_trim + max_len]
        
    # Calculate padding
    pad_len = max_len - len(seq)
    pad_left = 0
    pad_right = pad_len
    
    if pad_len > 0:
        if random_pad:
            pad_left = random.randint(0, pad_len)
            pad_right = pad_len - pad_left
    
    tensor = torch.zeros((4, max_len), dtype=torch.float32)
    
    for i, char in enumerate(seq):
        if char in mapping:
            idx = mapping[char]
            tensor[idx, pad_left + i] = 1.0
            
    return tensor

class RisdiplamModel(nn.Module):
    def __init__(self):
        super(RisdiplamModel, self).__init__()
        
        # Block 1: Detect Motifs
        self.conv1 = nn.Conv1d(in_channels=4, out_channels=64, kernel_size=11, padding='same')
        self.bn1 = nn.BatchNorm1d(64)
        self.relu = nn.ReLU()
        self.pool = nn.MaxPool1d(kernel_size=2)
        
        # Block 2: Dilated Context (Dilation=2)
        self.conv2 = nn.Conv1d(in_channels=64, out_channels=128, kernel_size=11, padding='same', dilation=2)
        self.bn2 = nn.BatchNorm1d(128)
        
        # Block 3: Wide Context (Dilation=4)
        self.conv3 = nn.Conv1d(in_channels=128, out_channels=128, kernel_size=11, padding='same', dilation=4)
        self.bn3 = nn.BatchNorm1d(128)
        
        # Head
        self.global_pool = nn.AdaptiveAvgPool1d(1)
        self.flatten = nn.Flatten()
        self.fc1 = nn.Linear(128, 64)
        self.dropout = nn.Dropout(0.5)
        self.fc2 = nn.Linear(64, 1) 

    def forward(self, x):
        x = self.conv1(x)
        x = self.bn1(x)
        x = self.relu(x)
        x = self.pool(x)
        
        x = self.conv2(x)
        x = self.bn2(x)
        x = self.relu(x)
        x = self.pool(x)
        
        x = self.conv3(x)
        x = self.bn3(x)
        x = self.relu(x)
        
        x = self.global_pool(x)
        x = self.flatten(x)
        
        x = self.fc1(x)
        x = self.relu(x)
        x = self.dropout(x)
        
        x = self.fc2(x)
        return x

class BalancedBatchGenerator:
    """
    Generates balanced batches with on-the-fly augmentation.
    Batch Composition: 25% Pos, 25% Neg, 25% Dead Donor, 25% Scrambled.
    """
    def __init__(self, x_pos, x_neg, batch_size=96):
        self.x_pos = x_pos
        self.x_neg = x_neg
        self.batch_size = batch_size
        self.quarter_batch = batch_size // 4
        
    def augment_seq(self, seq):
        """Randomly trim start/end to simulate variable exon length (Jitter)."""
        if len(seq) < 200: return seq
        trim_start = random.randint(0, 15)
        trim_end = random.randint(0, 15)
        if len(seq) - (trim_start + trim_end) < 100: return seq
        return seq[trim_start : len(seq)-trim_end]

    def mutate_dead_donor(self, seq):
        """Kill the GT donor at Intron 2 start."""
        idx_donor = len(seq) - OFFSET_END_FROM_RIGHT
        seq_list = list(seq)
        seq_list[idx_donor] = 'C'
        seq_list[idx_donor + 1] = 'C'
        return "".join(seq_list)

    def scramble_internal(self, seq):
        """Shuffle internal sequence preserving splice sites and base composition."""
        start_region = 550
        end_region = len(seq) - 550
        
        protect_start = start_region - 20
        protect_end_acceptor = start_region + 3
        protect_start_donor = end_region - 3
        protect_end_donor = end_region + 6
        
        s_idx = protect_end_acceptor
        e_idx = protect_start_donor
        
        if e_idx <= s_idx: return seq 
            
        seq_list = list(seq)
        subseq = seq_list[s_idx:e_idx]
        random.shuffle(subseq) # Strict shuffling
        seq_list[s_idx:e_idx] = subseq
        return "".join(seq_list)

    def generate_batch(self):
        # 1. Sample Indices
        pos_indices = np.random.choice(len(self.x_pos), self.quarter_batch)
        neg_indices = np.random.choice(len(self.x_neg), self.quarter_batch)
        
        batch_pos = [self.x_pos[i] for i in pos_indices]
        batch_neg = [self.x_neg[i] for i in neg_indices]
        
        # 2. Create Synthetics (Before Jitter)
        batch_dead = [self.mutate_dead_donor(s) for s in batch_pos]
        batch_scrambled = [self.scramble_internal(s) for s in batch_pos]
        
        raw_seqs = batch_pos + batch_neg + batch_dead + batch_scrambled
        
        # 3. Augment (Jitter) & Encode (Random Pad)
        tensor_list = []
        for s in raw_seqs:
            aug_s = self.augment_seq(s)
            tensor_list.append(one_hot_encode(aug_s, random_pad=True))
            
        X = torch.stack(tensor_list)
        
        # 4. Labels (First quarter is 1, rest 0)
        y = torch.zeros(len(raw_seqs), 1)
        y[:self.quarter_batch] = 1.0
        
        return X, y

def evaluate(model, val_gen, steps=20):
    model.eval()
    total_loss = 0
    correct = 0
    total = 0
    criterion = nn.BCEWithLogitsLoss()
    
    with torch.no_grad():
        for _ in range(steps):
            inputs, labels = val_gen.generate_batch()
            inputs, labels = inputs.to(device), labels.to(device)
            
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            total_loss += loss.item()
            
            # Accuracy
            predicted = (torch.sigmoid(outputs) > 0.5).float()
            total += labels.size(0)
            correct += (predicted == labels).sum().item()
            
    return total_loss / steps, correct / total

def main():
    # Initialize Log
    with open(LOG_FILE, "w") as f:
        f.write(f"DeepRisdiplam Training Log - {VERSION}\n")
        f.write(f"Date: {datetime.datetime.now()}\n")
        f.write("-" * 40 + "\n")

    log_message(f"Starting Training Run: {VERSION}")
    log_message(f"Device: {device}")
    log_message("1. Loading Data...")
    
    if not os.path.exists(INPUT_X_FILE):
        log_message(f"Error: {INPUT_X_FILE} not found.")
        return

    try:
        df_x = pd.read_csv(INPUT_X_FILE, sep='\t')
        df_y = pd.read_csv(INPUT_Y_FILE, sep='\t')
        df = pd.merge(df_x, df_y, on="event_id")
    except Exception as e:
        log_message(f"Error loading data: {e}")
        return

    positives = df[df['label'] == 1]['minigene_sequence'].values
    negatives = df[df['label'] == 0]['minigene_sequence'].values
    
    log_message(f"   Positives: {len(positives)}")
    log_message(f"   Negatives: {len(negatives)}")
    
    if len(positives) < 10:
        log_message("Error: Too few positives.")
        return

    pos_train, pos_val = train_test_split(positives, test_size=0.2, random_state=42)
    neg_train, neg_val = train_test_split(negatives, test_size=0.2, random_state=42)
    
    train_gen = BalancedBatchGenerator(pos_train, neg_train, BATCH_SIZE)
    val_gen = BalancedBatchGenerator(pos_val, neg_val, BATCH_SIZE)

    log_message("2. Building PyTorch Model...")
    model = RisdiplamModel().to(device)
    criterion = nn.BCEWithLogitsLoss()
    optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE)
    
    log_message(f"3. Training (Batch Size: {BATCH_SIZE})...")
    
    best_val_acc = 0.0
    
    for epoch in range(EPOCHS):
        model.train()
        running_loss = 0.0
        correct = 0
        total = 0
        
        pbar = tqdm(range(STEPS_PER_EPOCH), desc=f"Epoch {epoch+1}/{EPOCHS}", unit="batch")
        
        for _ in pbar:
            # 1. Get Batch
            inputs, labels = train_gen.generate_batch()
            inputs, labels = inputs.to(device), labels.to(device)
            
            # 2. Forward
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, labels)
            
            # 3. Backward
            loss.backward()
            optimizer.step()
            
            # 4. Stats
            running_loss += loss.item()
            predicted = (torch.sigmoid(outputs) > 0.5).float()
            total += labels.size(0)
            correct += (predicted == labels).sum().item()
            
            pbar.set_postfix({'loss': f"{loss.item():.4f}", 'acc': f"{correct/total:.4f}"})
            
        # End of Epoch Validation
        val_loss, val_acc = evaluate(model, val_gen)
        
        # Log metrics
        log_msg = f"Epoch {epoch+1}: Val Loss={val_loss:.4f} | Val Acc={val_acc:.4f}"
        log_message(log_msg)
        
        # Checkpointing
        if val_acc > best_val_acc:
            best_val_acc = val_acc
            torch.save(model.state_dict(), MODEL_SAVE_PATH)
            log_message(f"   Saved new best model to {MODEL_SAVE_PATH}")
            
    log_message("\nDone!")

if __name__ == "__main__":
    main()clear
