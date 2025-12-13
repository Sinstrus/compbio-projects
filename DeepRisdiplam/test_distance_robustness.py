import torch
import torch.nn as nn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import random

# --- Configuration ---
VERSION = "v8_components"
MODEL_PATH = f"deep_risdiplam_model_{VERSION}.pth"
INPUT_FILE = "dataset_components_v8.tsv"
MAX_SEQ_LEN = 3000
NUM_SAMPLES = 10  # Number of sequences to test
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# --- Model Architecture ---
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

def one_hot_encode(seq, max_len=MAX_SEQ_LEN):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    if len(seq) > max_len: seq = seq[:max_len]
    tensor = torch.zeros((4, max_len), dtype=torch.float32)
    for i, char in enumerate(seq):
        if char in mapping: tensor[mapping[char], i] = 1.0
    return tensor

def generate_random_bases(length):
    return "".join(random.choices("ACGT", k=length))

def get_base_score(model, row):
    seq = row['ExonA'] + row['Intron1'] + row['Pseudoexon'] + row['Intron2'] + row['ExonB']
    tensor = one_hot_encode(seq).unsqueeze(0).to(device)
    with torch.no_grad():
        return torch.sigmoid(model(tensor)).item()

def main():
    print(f"--- Running Multi-Sample Distance Robustness Test (n={NUM_SAMPLES}) ---")
    
    # 1. Load Model
    model = RisdiplamModel().to(device)
    try:
        model.load_state_dict(torch.load(MODEL_PATH))
    except FileNotFoundError:
        print("Model file not found.")
        return
    model.eval()
    
    # 2. Load Data & Select High Confidence Positives
    try:
        df = pd.read_csv(INPUT_FILE, sep='\t')
        positives = df[df['label'] == 1].copy()
    except FileNotFoundError: return

    print("Ranking positives by model confidence...")
    # Add a column for base prediction score
    positives['score'] = positives.apply(lambda r: get_base_score(model, r), axis=1)
    
    # Take top N highest scoring samples
    top_positives = positives.sort_values('score', ascending=False).head(NUM_SAMPLES)
    print(f"Selected top {NUM_SAMPLES} samples (Score range: {top_positives['score'].min():.4f} - {top_positives['score'].max():.4f})")

    # Define perturbations (Insert/Delete)
    # Include the specific small steps requested + larger steps
    deltas = sorted([-50, -25, -10, -5, -1, 0, 1, 5, 10, 25, 50, 100])
    
    all_scores = [] # List of lists
    
    plt.figure(figsize=(12, 7))
    
    # 3. Run Experiments
    for idx, row in top_positives.iterrows():
        ea, base_i1, ps, i2, eb = row['ExonA'], row['Intron1'], row['Pseudoexon'], row['Intron2'], row['ExonB']
        sample_scores = []
        mid_idx = len(base_i1) // 2
        
        for delta in deltas:
            if delta < 0:
                # Deletion
                num_delete = abs(delta)
                if num_delete >= len(base_i1): synth_i1 = "" 
                else:
                    start_cut = mid_idx - (num_delete // 2)
                    synth_i1 = base_i1[:start_cut] + base_i1[start_cut+num_delete:]
            elif delta > 0:
                # Insertion
                synth_i1 = base_i1[:mid_idx] + generate_random_bases(delta) + base_i1[mid_idx:]
            else:
                synth_i1 = base_i1
                
            full_seq = ea + synth_i1 + ps + i2 + eb
            tensor = one_hot_encode(full_seq).unsqueeze(0).to(device)
            with torch.no_grad():
                prob = torch.sigmoid(model(tensor)).item()
            sample_scores.append(prob)
            
        all_scores.append(sample_scores)
        # Plot individual faint line
        plt.plot(deltas, sample_scores, color='blue', alpha=0.15)

    # 4. Plot Mean Trend
    all_scores_np = np.array(all_scores)
    mean_scores = np.mean(all_scores_np, axis=0)
    
    plt.plot(deltas, mean_scores, color='red', linewidth=3, marker='o', label='Mean Prediction')
    
    plt.title(f"Model Robustness to Intron 1 Length Changes\n(Top {NUM_SAMPLES} High-Confidence Positives)")
    plt.xlabel("Change in Length (bp)\n(- = Deletion, + = Insertion)")
    plt.ylabel("Predicted Probability")
    plt.ylim(0, 1.05)
    plt.axvline(0, color='black', linestyle='--', alpha=0.5, label='Original')
    plt.axhline(0.5, color='gray', linestyle=':', label='Decision Threshold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.savefig("distance_robustness_multi.png")
    print("Saved plot to distance_robustness_multi.png")

if __name__ == "__main__":
    main()