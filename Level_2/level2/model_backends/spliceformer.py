# level2/model_backends/spliceformer.py
from __future__ import annotations
import torch
from transformers import AutoTokenizer, AutoModelForTokenClassification
from dataclasses import dataclass
from typing import List, Dict
from rich import print as rprint
import numpy as np

# --- CONFIGURATION ---
# Using a common SpliceFormer model ID. Verify this is the one you want.
MODEL_NAME = "ChatterjeeLab/SpliceFormer-common-human" 

@dataclass
class Site:
    pos: int
    prob: float
    kind: str

@dataclass
class OsaiResult:
    donors: List[Site]
    acceptors: List[Site]

# Global model loading
rprint(f"[bold blue]Loading SpliceFormer model: {MODEL_NAME}...[/bold blue]")
try:
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME, trust_remote_code=True)
    model = AutoModelForTokenClassification.from_pretrained(MODEL_NAME, trust_remote_code=True)
    
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model.to(device)
    model.eval()
    rprint(f"[green]Model loaded on {device}.[/green]")
except Exception as e:
    rprint(f"[bold red]Error loading model {MODEL_NAME}:[/bold red] {e}")
    model = None
    tokenizer = None

def score_batch(
    X_sequences: List[str],
    left_flank: str,
    right_flank: str,
    model_path: str = None, # Ignored, using global MODEL_NAME
    flanking_size: int = 0, # Ignored by Transformer context
) -> List[OsaiResult]:
    
    if model is None:
        raise RuntimeError("SpliceFormer model failed to load.")

    results = []
    L_left = len(left_flank)
    
    # Transformer batch size (smaller than CNN because of memory)
    BATCH_SIZE = 8 
    
    for i in range(0, len(X_sequences), BATCH_SIZE):
        batch_chunk = X_sequences[i : i + BATCH_SIZE]
        
        # 1. Prepare Inputs
        full_seqs = [left_flank + seq + right_flank for seq in batch_chunk]
        
        # 2. Tokenize
        # SpliceFormer typically handles sequences up to 2048 or 4096 tokens.
        # We rely on truncation if it exceeds limits.
        inputs = tokenizer(
            full_seqs, 
            return_tensors="pt", 
            padding=True, 
            truncation=True, 
            max_length=4096 
        )
        inputs = {k: v.to(device) for k, v in inputs.items()}
        
        # 3. Inference
        with torch.no_grad():
            outputs = model(**inputs)
            logits = outputs.logits # Shape: [Batch, SeqLen, NumClasses]
            probs = torch.softmax(logits, dim=-1)
            
        # 4. Parse Outputs
        # Class 0: None, Class 1: Acceptor, Class 2: Donor
        probs_np = probs.cpu().numpy()
        
        for j, seq_str in enumerate(batch_chunk):
            seq_len = len(seq_str)
            
            donors = []
            acceptors = []
            
            # Extract predictions for the "X" region
            # Assumes 1-to-1 token mapping + special tokens
            # Inputs usually start with [CLS], so index 0 is CLS.
            # Sequence starts at index 1.
            
            for k in range(seq_len):
                # Calculate index in the model output
                # 1 (CLS) + L_left (Flank) + k (Position in X)
                model_idx = 1 + L_left + k 
                
                if model_idx >= probs_np.shape[1] - 1: 
                    break
                
                p_acc = float(probs_np[j, model_idx, 1])
                p_don = float(probs_np[j, model_idx, 2])
                
                if p_acc > 0.01:
                    acceptors.append(Site(k, p_acc, "Acceptor"))
                if p_don > 0.01:
                    donors.append(Site(k, p_don, "Donor"))
            
            results.append(OsaiResult(donors, acceptors))
            
    return results