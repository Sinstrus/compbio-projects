import torch
import torch.nn as nn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
import sys

# --- Configuration ---
VERSION = "v8_components"
MODEL_PATH = f"deep_risdiplam_model_{VERSION}.pth"
INPUT_FILE = "dataset_components_v8.tsv"
MAX_SEQ_LEN = 1500
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# --- Model Architecture (Must Match V7) ---
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
    
    # Center the sequence if it fits, or truncate
    if len(seq) > max_len:
        start = (len(seq) - max_len) // 2
        seq = seq[start : start + max_len]
        pad_left = 0
    else:
        pad_left = 0 # No random padding for interpretation
        
    tensor = torch.zeros((4, max_len), dtype=torch.float32)
    for i, char in enumerate(seq):
        if char in mapping:
            tensor[mapping[char], pad_left + i] = 1.0
            
    return tensor, pad_left

def get_gradient_window(model, ea, i1, ps, i2, eb):
    """
    1. Stitches sequence.
    2. Runs Backward pass.
    3. Extracts 20bp window of GRADIENTS around the splice site.
    """
    # 1. Stitch
    full_seq = ea + i1 + ps + i2 + eb
    
    # Calculate exactly where the donor is (End of Pseudo + Start of Intron 2)
    # The junction is at index len(ea + i1 + ps)
    junction_idx = len(ea + i1 + ps)
    
    # 2. Encode
    tensor, pad_left = one_hot_encode(full_seq)
    tensor = tensor.unsqueeze(0).to(device)
    tensor.requires_grad_()
    
    # 3. Forward & Backward
    model.zero_grad()
    output = model(tensor)
    output.backward()
    
    # 4. Extract Gradients
    # Shape: (1, 4, 1500) -> (4, 1500)
    grads = tensor.grad.data.abs().cpu().squeeze(0).numpy()
    
    # We need the gradients corresponding to the specific bases present
    # Input was one-hot, so we take the gradient at the channel that was '1'
    # This gives us "Importance of the actual nucleotide"
    importance_seq = []
    
    # Define window indices relative to the FULL tensor
    # Note: If one_hot_encode truncated the sequence, we need to adjust junction_idx
    if len(full_seq) > MAX_SEQ_LEN:
        trim_start = (len(full_seq) - MAX_SEQ_LEN) // 2
        junction_idx = junction_idx - trim_start
    
    window_start = junction_idx - 10
    window_end = junction_idx + 10
    
    # Safety Check
    if window_start < 0 or window_end >= MAX_SEQ_LEN:
        return None, None

    # Extract score for each position in window
    scores = []
    bases = []
    
    for i in range(window_start, window_end):
        # Find which base is here
        col = tensor[0, :, i].cpu().detach().numpy()
        if col.sum() == 0: 
            # Padding/N
            scores.append(0)
            bases.append('N')
        else:
            base_idx = np.argmax(col)
            score = grads[base_idx, i]
            scores.append(score)
            bases.append("ACGT"[base_idx])
            
    return scores, bases

def main():
    print("--- Generating Model Importance Logo ---")
    
    # Load Model
    model = RisdiplamModel().to(device)
    try:
        model.load_state_dict(torch.load(MODEL_PATH))
    except FileNotFoundError:
        print("Model file not found. Train V7 first.")
        return
    model.eval()
    
    # Load Data
    try:
        df = pd.read_csv(INPUT_FILE, sep='\t')
    except FileNotFoundError:
        print("Data file not found.")
        return

    positives = df[df['label'] == 1]
    print(f"Analyzing gradients for {len(positives)} sequences...")
    
    # Accumulate
    aggregated_scores = np.zeros((20, 4)) # 20 positions, 4 bases
    count_valid = 0
    
    base_map = {'A':0, 'C':1, 'G':2, 'T':3}
    
    for idx, row in positives.iterrows():
        scores, bases = get_gradient_window(
            model, row['ExonA'], row['Intron1'], 
            row['Pseudoexon'], row['Intron2'], row['ExonB']
        )
        
        if scores is None: continue
        
        # Add to aggregate matrix
        for i in range(20):
            base = bases[i]
            if base in base_map:
                # We add the score to the specific base bucket
                aggregated_scores[i, base_map[base]] += scores[i]
                
        count_valid += 1
        if count_valid % 50 == 0: print(f"Processed {count_valid}...", end='\r')

    print(f"\nSuccessfully analyzed {count_valid} sequences.")
    
    # Normalize
    avg_scores = aggregated_scores / count_valid
    
    # Create DataFrame for Logomaker
    df_logo = pd.DataFrame(avg_scores, columns=['A','C','G','T'])
    
    # Plot
    plt.figure(figsize=(10, 4))
    logo = logomaker.Logo(df_logo, shade_below=.5, fade_below=.5, color_scheme='classic')
    
    logo.style_spines(visible=False)
    logo.style_spines(spines=['bottom', 'left'], visible=True)
    logo.ax.set_xticks(range(20))
    logo.ax.set_xticklabels(np.arange(-10, 10))
    logo.ax.set_ylabel("Importance (Gradient Magnitude)")
    logo.ax.set_title("What the Model Cares About\n(Gradient-Weighted Sequence Logo)")
    
    plt.axvline(9.5, color='black', linestyle='--', linewidth=1)
    
    plt.savefig("risdiplam_importance_logo.png")
    print("Saved logo to risdiplam_importance_logo.png")

if __name__ == "__main__":
    main()