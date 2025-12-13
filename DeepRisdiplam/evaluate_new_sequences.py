import torch
import torch.nn as nn
import pandas as pd
import numpy as np
import sys

# --- Configuration ---
VERSION = "v7_components"
MODEL_PATH = f"deep_risdiplam_model_{VERSION}.pth"
INPUT_CSV = "test_set_template.csv" # Your edited file
MAX_SEQ_LEN = 3000
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
    # Case Insensitive Mapping
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3} 
    
    seq = seq.upper() # FORCE UPPERCASE
    
    if len(seq) > max_len:
        # Center crop if too long
        start = (len(seq) - max_len) // 2
        seq = seq[start : start + max_len]
        
    tensor = torch.zeros((4, max_len), dtype=torch.float32)
    for i, char in enumerate(seq):
        if char in mapping:
            tensor[mapping[char], i] = 1.0
    return tensor

def main():
    print(f"--- Evaluating New Sequences (Model: {VERSION}) ---")
    
    # 1. Load Model
    model = RisdiplamModel().to(device)
    try:
        model.load_state_dict(torch.load(MODEL_PATH))
        model.eval()
    except FileNotFoundError:
        print(f"Error: Model {MODEL_PATH} not found.")
        return

    # 2. Load CSV
    try:
        df = pd.read_csv(INPUT_CSV)
    except FileNotFoundError:
        print(f"Error: {INPUT_CSV} not found.")
        return

    print(f"{'ID':<30} | {'Type':<15} | {'Score':<10} | {'Verdict'}")
    print("-" * 70)

    results = []

    for idx, row in df.iterrows():
        eid = str(row['event_id'])
        stype = str(row['type'])
        
        # Check for placeholder text
        if "Insert_" in str(row['ExonA']):
            print(f"{eid:<30} | {stype:<15} | SKIPPED    | (Placeholder Data)")
            continue
            
        # Stitch
        # Convert to string to handle potential NaN/float issues in CSV reading
        ea = str(row['ExonA']).strip()
        i1 = str(row['Intron1']).strip()
        ps = str(row['Pseudoexon']).strip()
        i2 = str(row['Intron2']).strip()
        eb = str(row['ExonB']).strip()
        
        full_seq = ea + i1 + ps + i2 + eb
        
        # Predict
        tensor = one_hot_encode(full_seq).unsqueeze(0).to(device)
        
        with torch.no_grad():
            logit = model(tensor)
            prob = torch.sigmoid(logit).item()
            
        verdict = "RESPONSIVE" if prob > 0.5 else "Silent"
        color = "\033[92m" if prob > 0.5 else "\033[91m" # Green/Red codes
        reset = "\033[0m"
        
        print(f"{eid:<30} | {stype:<15} | {prob:.4f}     | {color}{verdict}{reset}")
        
        results.append(prob)

    # Save results back to CSV
    df['Predicted_Score'] = pd.Series(results)
    df.to_csv("evaluated_results.csv", index=False)
    print("\nSaved detailed results to evaluated_results.csv")

if __name__ == "__main__":
    main()