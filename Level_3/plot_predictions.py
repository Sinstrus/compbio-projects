import pandas as pd
import numpy as np
import torch
import matplotlib.pyplot as plt
from torch.utils.data import DataLoader
from train_quartet_model import QuartetNet, QuartetDataset # Import your class structure
import sys

# --- Config ---
MODEL_PATH = "quartet_regression_net.pt"
DATA_FILE = "risdiplam_full_dataset.csv.gz"
BATCH_SIZE = 256

def get_real_sample(row):
    # Copy of the logic from train_quartet_model.py to ensure consistency
    try:
        # Helpers
        def pad_sequence(seq, target_len, pad_side='right'):
            if len(seq) == target_len: return seq
            if len(seq) > target_len:
                if pad_side == 'right': return seq[:target_len]
                else: return seq[-target_len:]
            pad = "N" * (target_len - len(seq))
            if pad_side == 'right': return seq + pad
            else: return pad + seq
            
        ctx_a = row['seq_ExonA'][-25:] + row['seq_Intron1'][:25]
        ctx_b = row['seq_Intron2'][-25:] + row['seq_ExonB'][:25]
        
        p_exon = row['seq_PseudoExon']
        p_int1 = row['seq_Intron1']
        p_acc = pad_sequence(p_int1[-50:] + p_exon[:100], 150, 'right')
        
        p_int2 = row['seq_Intron2']
        p_don = pad_sequence(p_exon[-100:] + p_int2[:50], 150, 'left')
        
        # Regression Label
        label = float(row['dPSI'])
        
        if len(ctx_a)!=50 or len(ctx_b)!=50: return None
        return (ctx_a, ctx_b, p_acc, p_don, label)
    except: return None

def main():
    print("Loading Data for Evaluation...")
    df = pd.read_csv(DATA_FILE)
    df = df.dropna(subset=['seq_ExonA', 'seq_Intron1', 'seq_PseudoExon', 'seq_Intron2', 'seq_ExonB'])
    
    # We only want to plot REAL data, not decoys, to see how well we predict actual biology
    # Filter for "Interesting" rows to keep the plot readable, or plot all real rows
    print(f"Processing {len(df)} real biological events...")
    
    data_a, data_b, data_p_acc, data_p_don, data_y = [], [], [], [], []
    
    for _, row in df.iterrows():
        res = get_real_sample(row)
        if res:
            data_a.append(res[0]); data_b.append(res[1])
            data_p_acc.append(res[2]); data_p_don.append(res[3])
            data_y.append(res[4])

    # To Tensor
    dataset = QuartetDataset(data_a, data_b, data_p_acc, data_p_don, data_y)
    loader = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=False)
    
    # Load Model
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Evaluating on {device}...")
    model = QuartetNet().to(device)
    model.load_state_dict(torch.load(MODEL_PATH, map_location=device))
    model.eval()
    
    preds = []
    targets = []
    
    with torch.no_grad():
        for xa, xb, xpa, xpd, lbl in loader:
            xa, xb, xpa, xpd = xa.to(device), xb.to(device), xpa.to(device), xpd.to(device)
            out = model(xa, xb, xpa, xpd)
            preds.extend(out.cpu().numpy().flatten())
            targets.extend(lbl.numpy().flatten())
            
    # Plot
    print("Generating Plot...")
    plt.figure(figsize=(10, 10))
    
    # Scatter plot with transparency
    plt.scatter(targets, preds, alpha=0.1, s=2, c='blue')
    
    # Perfect prediction line
    plt.plot([-1, 1], [-1, 1], color='red', linestyle='--', label='Perfect Prediction')
    
    plt.title(f"DeepRisdiplam: Predicted vs Actual dPSI\nPearson R = {np.corrcoef(targets, preds)[0,1]:.4f}")
    plt.xlabel("Actual dPSI (Biological Ground Truth)")
    plt.ylabel("Predicted dPSI (Model Output)")
    plt.xlim(-0.5, 1.0)
    plt.ylim(-0.5, 1.0)
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.savefig("results_scatter.png", dpi=300)
    print("Saved plot to results_scatter.png")

if __name__ == "__main__":
    main()