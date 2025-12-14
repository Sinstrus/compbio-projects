#!/usr/bin/env python3
"""
Simple GPU test with ESM-2 (much smaller model) to verify setup
"""
import torch
from transformers import AutoTokenizer, EsmModel

print("=" * 80)
print("Simple GPU Test with ESM-2 (protein language model)")
print("=" * 80)

# Check GPU
print(f"\nCUDA available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"GPU: {torch.cuda.get_device_name(0)}")
    print(f"VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

# Test sequence
sequence = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSTADMGWFRQAPGKGRELVAAVSGSGFSTYSDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAKATISLYYAMDVWGQGTTVTVSS"

print(f"\nTest sequence length: {len(sequence)} amino acids")

# Load small ESM-2 model (only ~35MB)
print("\nLoading ESM-2 small model (fast, lightweight test)...")
tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
model = EsmModel.from_pretrained("facebook/esm2_t6_8M_UR50D")

# Move to GPU
model = model.cuda()
model.eval()

print("✓ Model loaded on GPU")
print(f"GPU Memory: {torch.cuda.memory_allocated() / 1e9:.2f} GB")

# Run inference
print("\nRunning inference on GPU...")
inputs = tokenizer(sequence, return_tensors="pt")
inputs = {k: v.cuda() for k, v in inputs.items()}

with torch.no_grad():
    outputs = model(**inputs)

embeddings = outputs.last_hidden_state

print(f"✓ Inference complete!")
print(f"Output shape: {embeddings.shape}")
print(f"GPU Memory: {torch.cuda.memory_allocated() / 1e9:.2f} GB")

print("\n" + "=" * 80)
print("SUCCESS! Your GPU setup works for protein models.")
print("=" * 80)
print("\nThe issue with ESMFold is that the full folding model is VERY large.")
print("For structure prediction, consider:")
print("  1. Use ColabFold locally (being installed)")
print("  2. Use ESMFold web server: https://esmatlas.com/resources?action=fold")
print("  3. Use AlphaFold Server: https://alphafoldserver.com/")
