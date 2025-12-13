import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_score, recall_score, f1_score, confusion_matrix
from collections import Counter
import sys

# --- CONFIGURATION ---
INPUT_FILE = "dataset_components_v8.tsv"

# --- FEATURE ENGINEERING ---
def get_gc(seq):
    if len(seq) == 0: return 0
    return (seq.count('G') + seq.count('C')) / len(seq)

def get_kmers(seq, k=2):
    """Counts dinucleotides (AA, AC, AG, etc.) normalized by length"""
    if len(seq) < k: return {p:0 for p in ['AA']} # Dummy return
    
    counts = Counter([seq[i:i+k] for i in range(len(seq) - k + 1)])
    total = sum(counts.values())
    
    # Normalize
    freqs = {kmer: v/total for kmer, v in counts.items()}
    return freqs

def extract_features(row):
    """
    Turns a DNA sequence into a vector of biological statistics.
    No deep learning "black box" features.
    """
    ea, i1, ps, i2, eb = row['ExonA'], row['Intron1'], row['Pseudoexon'], row['Intron2'], row['ExonB']
    
    feats = {}
    
    # 1. Structural Features (Lengths)
    feats['len_ea'] = len(ea)
    feats['len_i1'] = len(i1)
    feats['len_ps'] = len(ps)
    feats['len_i2'] = len(i2)
    feats['len_eb'] = len(eb)
    
    # 2. Compositional Features (GC Content)
    feats['gc_ea'] = get_gc(ea)
    feats['gc_i1'] = get_gc(i1)
    feats['gc_ps'] = get_gc(ps)
    feats['gc_i2'] = get_gc(i2)
    feats['gc_eb'] = get_gc(eb)
    
    # 3. Motif Features (Specific Risdiplam Signals)
    # Splice Donor (Start of Intron 2): "GT..."
    # Risdiplam likes a 'G' at -1 (End of Pseudo) and 'T' at +1 (Start of Intron 2)
    feats['donor_is_canonical_gt'] = 1 if i2.startswith('GT') else 0
    feats['donor_minus1_is_g'] = 1 if ps.endswith('G') else 0
    
    # 4. Dinucleotide Counts (Texture) in the Pseudoexon
    # (High GA content often indicates Splicing Enhancers)
    dimers = get_kmers(ps, k=2)
    for dimer in ['GA', 'AG', 'GG', 'GT']:
        feats[f'ps_dimer_{dimer}'] = dimers.get(dimer, 0)

    return feats

def main():
    print("--- Training Classical Baseline Models ---")
    
    # 1. Load Data
    print("Loading data...")
    df = pd.read_csv(INPUT_FILE, sep='\t')
    df = df.dropna(subset=['ExonA', 'Intron1', 'Pseudoexon', 'Intron2', 'ExonB'])
    
    # 2. Extract Features
    print("Engineering biological features...")
    feature_list = []
    labels = []
    
    for _, row in df.iterrows():
        feature_list.append(extract_features(row))
        labels.append(int(row['label']))
        
    X = pd.DataFrame(feature_list)
    y = np.array(labels)
    
    print(f"Feature Matrix Shape: {X.shape}")
    
    # 3. Split (Stratified - SAME as Deep Learning)
    X_train, X_val, y_train, y_val = train_test_split(
        X, y, test_size=0.1, random_state=42, stratify=y
    )
    
    # 4. Train Random Forest (The "Hammer")
    print("\n--- Model 1: Random Forest (Non-Linear) ---")
    rf = RandomForestClassifier(n_estimators=200, max_depth=10, class_weight='balanced', random_state=42)
    rf.fit(X_train, y_train)
    
    y_pred_rf = rf.predict(X_val)
    
    # 5. Metrics RF
    cm = confusion_matrix(y_val, y_pred_rf, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()
    
    print(f"Confusion Matrix:\n{cm}")
    print(f"TP (Caught): {tp}")
    print(f"FN (Missed): {fn}")
    print(f"FP (False Pos): {fp}")
    print(f"Precision: {precision_score(y_val, y_pred_rf, zero_division=0):.4f}")
    print(f"Recall:    {recall_score(y_val, y_pred_rf, zero_division=0):.4f}")
    
    # 6. Feature Importance
    print("\nTop 5 Features driving the Random Forest:")
    importances = rf.feature_importances_
    indices = np.argsort(importances)[::-1]
    for i in range(5):
        print(f"{i+1}. {X.columns[indices[i]]} ({importances[indices[i]]:.4f})")

    # 7. Train Logistic Regression (Linear Baseline)
    print("\n--- Model 2: Logistic Regression (Linear) ---")
    # Increased max_iter to ensure convergence
    lr = LogisticRegression(class_weight='balanced', max_iter=2000)
    lr.fit(X_train, y_train)
    y_pred_lr = lr.predict(X_val)
    
    cm_lr = confusion_matrix(y_val, y_pred_lr, labels=[0, 1])
    tn_lr, fp_lr, fn_lr, tp_lr = cm_lr.ravel()
    
    print(f"Confusion Matrix:\n{cm_lr}")
    print(f"TP (Caught): {tp_lr}")
    print(f"FN (Missed): {fn_lr}")
    print(f"FP (False Pos): {fp_lr}")
    print(f"TN (True Neg): {tn_lr}")
    print(f"Precision: {precision_score(y_val, y_pred_lr, zero_division=0):.4f}")
    print(f"Recall:    {recall_score(y_val, y_pred_lr, zero_division=0):.4f}")
    print(f"F1 Score:  {f1_score(y_val, y_pred_lr, zero_division=0):.4f}")

    # 8. Conclusion Logic
    print("\n" + "="*40)
    print("VERDICT")
    print("="*40)
    
    print(f"Logistic Regression F1: {f1_score(y_val, y_pred_lr, zero_division=0):.4f}")
    print(f"Logistic Regression FP Count: {fp_lr}")
    
    if fp_lr > 1000:
        print("RESULT: Logistic Regression is spamming positives (Low Precision).")
        print("MEANING: The high TP count (Recall) is artificial. It's not a better model.")
    elif tp_lr > 9 and fp_lr < 500:
        print("RESULT: Logistic Regression genuinely found more hits with reasonable precision.")
        print("MEANING: Explicit scalar features are critical.")

if __name__ == "__main__":
    main()