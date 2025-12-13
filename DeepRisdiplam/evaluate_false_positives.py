import torch
import torch.nn as nn
import pandas as pd
import numpy as np
import csv
import sys
import matplotlib.pyplot as plt
from tqdm import tqdm

# --- Configuration ---
VERSION = "v8_components"
MODEL_PATH = f"deep_risdiplam_model_{VERSION}.pth"
INPUT_FILE = "dataset_components_v8.tsv"
MAX_SEQ_LEN = 3000  # Conservative max length for stitching
BATCH_SIZE = 64
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# --- Model Architecture (Reconstructed from snippets) ---
class RisdiplamModel(nn.Module):
    def __init__(self):
        super(RisdiplamModel, self).__init__()
        # Matches V8 architecture
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

# --- Helper: One Hot Encoding ---
def one_hot_encode(seq, max_len=MAX_SEQ_LEN):
    # Map A,C,G,T to 0,1,2,3
    mapping = {
        'A': 0, 'a': 0,
        'C': 1, 'c': 1,
        'G': 2, 'g': 2,
        'T': 3, 't': 3
    }
    
    seq_len = len(seq)
    # Pad or truncate
    if seq_len > max_len:
        seq = seq[:max_len]
        seq_len = max_len
        
    tensor = np.zeros((4, max_len), dtype=np.float32)
    
    for i, char in enumerate(seq):
        if i >= max_len: break
        if char in mapping:
            tensor[mapping[char], i] = 1.0
            
    return torch.tensor(tensor)

def main():
    print(f"--- Loading Model: {MODEL_PATH} ---")
    model = RisdiplamModel().to(device)
    try:
        model.load_state_dict(torch.load(MODEL_PATH, map_location=device))
    except Exception as e:
        print(f"Error loading model weights: {e}")
        print("Ensure the class architecture matches exactly what was saved.")
        return

    model.eval()
    
    print(f"--- Loading Negatives from {INPUT_FILE} ---")
    negatives = []
    
    # Read TSV and filter for label '0'
    with open(INPUT_FILE, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['label'] == '0':
                # Stitch the components exactly as the model expects
                full_seq = (
                    row['ExonA'] + 
                    row['Intron1'] + 
                    row['Pseudoexon'] + 
                    row['Intron2'] + 
                    row['ExonB']
                )
                negatives.append(full_seq)
    
    total_negatives = len(negatives)
    print(f"Total Negative Samples Found: {total_negatives}")
    
    if total_negatives == 0:
        print("No negatives found. Check your file path or labels.")
        return

    # --- Inference Loop ---
    scores = []
    print("Running Inference...")
    
    with torch.no_grad():
        for i in tqdm(range(0, total_negatives, BATCH_SIZE)):
            batch_seqs = negatives[i : i + BATCH_SIZE]
            
            # Encode batch
            tensors = [one_hot_encode(s) for s in batch_seqs]
            batch_tensor = torch.stack(tensors).to(device)
            
            # Forward pass
            outputs = model(batch_tensor)
            probs = torch.sigmoid(outputs).cpu().numpy().flatten()
            scores.extend(probs)

    scores = np.array(scores)
    
    # --- Statistics ---
    false_positives = np.sum(scores > 0.5)
    strict_false_positives = np.sum(scores > 0.9)
    
    fpr = (false_positives / total_negatives) * 100
    strict_fpr = (strict_false_positives / total_negatives) * 100
    
    print("\n" + "="*40)
    print(f"RESULTS for {total_negatives} NEGATIVE SAMPLES")
    print("="*40)
    print(f"False Positives (Score > 0.5): {false_positives} ({fpr:.2f}%)")
    print(f"Strict False Positives (Score > 0.9): {strict_false_positives} ({strict_fpr:.2f}%)")
    print(f"Mean Prediction Score: {np.mean(scores):.4f}")
    print(f"Median Prediction Score: {np.median(scores):.4f}")
    
    # --- Histogram ---
    plt.figure(figsize=(10, 6))
    plt.hist(scores, bins=50, color='skyblue', edgecolor='black', range=(0,1))
    plt.axvline(0.5, color='red', linestyle='--', label='Decision Threshold (0.5)')
    plt.yscale('log') # Log scale to see the outliers
    plt.title(f"Score Distribution for {total_negatives} Negative Samples (Log Scale)")
    plt.xlabel("Predicted Probability (1.0 = Drug Responder)")
    plt.ylabel("Count (Log Scale)")
    plt.legend()
    plt.savefig("false_positive_analysis.png")
    print("\nSaved histogram to 'false_positive_analysis.png'")

if __name__ == "__main__":
    main()