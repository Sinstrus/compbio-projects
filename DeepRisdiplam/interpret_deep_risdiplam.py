import torch
import torch.nn as nn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# --- Configuration ---
VERSION = "v4"
MODEL_PATH = f"deep_risdiplam_model_{VERSION}.pth"
INPUT_X_FILE = "dataset_X_sequences.tsv"
INPUT_Y_FILE = "dataset_y_labels.tsv"

# Structure Constants (Must match training)
EXON_A_LEN = 50
INTRON_1_PART_LEN = 250
INTRON_2_PART_LEN = 250
EXON_B_LEN = 50
MAX_SEQ_LEN = 1500 

# Calculated Offsets
OFFSET_PSEUDO_START = EXON_A_LEN + (INTRON_1_PART_LEN * 2) # Start of Pseudoexon
OFFSET_PSEUDO_END_FROM_RIGHT = INTRON_2_PART_LEN * 2 + EXON_B_LEN # End of Pseudoexon

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# --- Re-define Model Architecture (Must match exactly) ---
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
    if len(seq) > max_len:
        start = (len(seq) - max_len) // 2
        seq = seq[start : start + max_len]
    
    pad_len = max_len - len(seq)
    # No random padding for interpretation - we want stable indices
    tensor = torch.zeros((4, max_len), dtype=torch.float32)
    for i, char in enumerate(seq):
        if char in mapping:
            tensor[mapping[char], i] = 1.0
    return tensor

def compute_saliency(model, seq_tensor):
    """
    Computes gradient of output score w.r.t input sequence.
    """
    model.eval()
    seq_tensor.requires_grad_()
    
    output = model(seq_tensor.unsqueeze(0).to(device))
    output.backward()
    
    # Saliency is the absolute value of the gradient
    # We take the max across the 4 channels (A,C,G,T) to get one score per position
    saliency = seq_tensor.grad.data.abs().cpu().numpy() # (4, 1500)
    saliency_per_pos = np.max(saliency, axis=0) # (1500,)
    
    return saliency_per_pos

def plot_importance(saliency_map, event_id, seq_len):
    """
    Plots the saliency map with annotations for the specific Minigene structure.
    """
    plt.figure(figsize=(15, 5))
    
    # Plot the raw importance
    # Use only the valid length of the sequence
    valid_map = saliency_map[:seq_len]
    plt.plot(valid_map, label='Gradient Importance', color='black', linewidth=0.8)
    
    # Shade Regions based on "Digital Minigene" Structure
    # 1. Exon A (First 50)
    plt.axvspan(0, EXON_A_LEN, color='red', alpha=0.1, label='Exon A')
    
    # 2. Intron 1 (Next 500)
    start_i1 = EXON_A_LEN
    end_i1 = start_i1 + (INTRON_1_PART_LEN * 2)
    plt.axvspan(start_i1, end_i1, color='gray', alpha=0.1, label='Intron 1')
    
    # 3. Pseudoexon (Variable Length)
    # We know the total length, and we know the end structure is fixed size.
    # So Pseudo Start = end_i1
    # Pseudo End = Total - Fixed_End_Size
    start_pseudo = end_i1
    end_pseudo = seq_len - OFFSET_PSEUDO_END_FROM_RIGHT
    plt.axvspan(start_pseudo, end_pseudo, color='blue', alpha=0.1, label='Pseudoexon')
    
    # 4. Intron 2 (Next 500)
    start_i2 = end_pseudo
    end_i2 = start_i2 + (INTRON_2_PART_LEN * 2)
    plt.axvspan(start_i2, end_i2, color='gray', alpha=0.1) # Intron 2
    
    # 5. Exon B (Last 50)
    plt.axvspan(end_i2, seq_len, color='red', alpha=0.1, label='Exon B')

    # Highlight Splice Sites (The Critical Check)
    # Pseudo Acceptor (Start of Pseudo)
    plt.axvline(start_pseudo, color='green', linestyle='--', alpha=0.5)
    # Pseudo Donor (End of Pseudo)
    plt.axvline(end_pseudo, color='green', linestyle='--', alpha=0.5)

    plt.title(f"Saliency Map for Event: {event_id}\n(Where is the model looking?)")
    plt.xlabel("Position (bp)")
    plt.ylabel("Importance Score (Gradient Norm)")
    plt.legend(loc='upper right')
    
    filename = f"saliency_{event_id.replace(':','_')}.png"
    plt.savefig(filename)
    print(f"   Saved plot to {filename}")
    plt.close()

def main():
    print(f"--- Interpreting Model {VERSION} ---")
    
    # 1. Load Model
    if not os.path.exists(MODEL_PATH):
        print("Model file not found.")
        return
        
    model = RisdiplamModel().to(device)
    model.load_state_dict(torch.load(MODEL_PATH))
    print("Model loaded.")

    # 2. Load Data (Positives Only)
    df_x = pd.read_csv(INPUT_X_FILE, sep='\t')
    df_y = pd.read_csv(INPUT_Y_FILE, sep='\t')
    df = pd.merge(df_x, df_y, on="event_id")
    
    positives = df[df['label'] == 1]
    
    if len(positives) == 0:
        print("No positives found to analyze.")
        return

    print(f"Found {len(positives)} positives. Analyzing top 5...")
    
    # 3. Analyze
    for idx, row in positives.head(5).iterrows():
        seq = row['minigene_sequence']
        eid = row['event_id']
        
        tensor = one_hot_encode(seq).to(device)
        saliency = compute_saliency(model, tensor)
        
        plot_importance(saliency, eid, len(seq))

    print("Done. Check the .png files.")

if __name__ == "__main__":
    main()