import numpy as np
import torch
import sys

# Constants
MAX_SEQ_LEN = 1500

def vectorized_one_hot_encode(seqs, max_len=MAX_SEQ_LEN, random_pad=False):
    batch_size = len(seqs)
    mapping = np.zeros(128, dtype=np.int8) + 4 
    mapping[ord('A')] = 0; mapping[ord('a')] = 0
    mapping[ord('C')] = 1; mapping[ord('c')] = 1
    mapping[ord('G')] = 2; mapping[ord('g')] = 2
    mapping[ord('T')] = 3; mapping[ord('t')] = 3
    
    tensor_batch = np.zeros((batch_size, 4, max_len), dtype=np.float32)
    
    for i, seq in enumerate(seqs):
        # Logic from train script
        if len(seq) > max_len:
            start = (len(seq) - max_len) // 2
            seq = seq[start : start + max_len]
            
        indices = mapping[np.array([ord(c) for c in seq])]
        
        seq_len = len(indices)
        pad_total = max_len - seq_len
        pad_left = 0
        if pad_total > 0 and random_pad:
            pad_left = np.random.randint(0, pad_total + 1)
            
        for k in range(seq_len):
            idx = indices[k]
            if idx < 4: 
                tensor_batch[i, idx, pad_left + k] = 1.0
                
    return torch.from_numpy(tensor_batch)

def main():
    print("--- Debugging Tensor Encoding ---")
    
    # Test Sequence
    seq = "ACGT" * 100 # 400bp
    print(f"Sequence Length: {len(seq)}")
    
    # Encode
    batch = [seq, seq]
    tensors = vectorized_one_hot_encode(batch, max_len=1500, random_pad=True)
    
    print(f"Tensor Shape: {tensors.shape}")
    print(f"Tensor Sum (Should be {len(seq)} * 2 = 800): {tensors.sum().item()}")
    
    if tensors.sum().item() == 0:
        print("CRITICAL ERROR: Encoding produced all zeros!")
        print("Check ASCII mapping.")
    else:
        print("Encoding seems operational.")
        
    # Check if 'N' breaks it
    seq_n = "ACGTN" * 10
    tensors_n = vectorized_one_hot_encode([seq_n], max_len=1500)
    print(f"Sequence with N sum (Should be 40): {tensors_n.sum().item()}")

if __name__ == "__main__":
    main()