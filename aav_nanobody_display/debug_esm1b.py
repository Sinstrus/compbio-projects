#!/usr/bin/env python3
"""
Debug script to identify where ESM-1b loading hangs.
Run with: python debug_esm1b.py
"""
import os
import sys
import time
import signal

# Set environment variables before importing torch
os.environ['TORCH_HOME'] = os.path.expanduser('~/.cache/torch')
os.environ['HF_HOME'] = os.path.expanduser('~/.cache/huggingface')

def log(msg):
    """Timestamped logging."""
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def timeout_handler(signum, frame):
    log("ERROR: Timeout reached! Process is stuck.")
    sys.exit(1)

def main():
    log("=" * 60)
    log("ESM-1b Debug Script")
    log("=" * 60)

    # Step 1: Import torch
    log("Step 1: Importing torch...")
    start = time.time()
    import torch
    log(f"  ✓ torch imported in {time.time() - start:.1f}s")

    # Step 2: Check CUDA
    log("Step 2: Checking CUDA...")
    device = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"  Device: {device}")
    if device == "cuda":
        log(f"  GPU: {torch.cuda.get_device_name(0)}")
        log(f"  VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

    # Step 3: Import esm
    log("Step 3: Importing esm library...")
    start = time.time()
    import esm
    log(f"  ✓ esm imported in {time.time() - start:.1f}s")

    # Step 4: Check model cache
    cache_dir = os.path.expanduser("~/.cache/torch/hub/checkpoints/")
    log(f"Step 4: Checking cache at {cache_dir}")
    if os.path.exists(cache_dir):
        files = os.listdir(cache_dir)
        log(f"  Cache files: {files if files else 'empty'}")
    else:
        log(f"  Cache directory doesn't exist")

    # Step 5: Load model (this is where it likely hangs)
    log("Step 5: Loading ESM-1b model...")
    log("  This will download ~1.5 GB on first run...")
    log("  Model: esm1b_t33_650M_UR50S")

    # Set a timeout (5 minutes for download + load)
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(300)  # 5 minute timeout

    start = time.time()
    try:
        # Try loading with explicit torch.hub parameters
        model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
        log(f"  ✓ Model loaded in {time.time() - start:.1f}s")
        signal.alarm(0)  # Cancel timeout
    except Exception as e:
        log(f"  ✗ Model load failed: {e}")
        import traceback
        traceback.print_exc()
        return

    # Step 6: Check model parameters
    log("Step 6: Model info")
    num_params = sum(p.numel() for p in model.parameters())
    log(f"  Parameters: {num_params:,} ({num_params/1e6:.1f}M)")

    # Step 7: Move to GPU
    log(f"Step 7: Moving model to {device}...")
    start = time.time()
    model = model.to(device)
    model.eval()
    log(f"  ✓ Model on {device} in {time.time() - start:.1f}s")

    if device == "cuda":
        log(f"  GPU Memory used: {torch.cuda.memory_allocated() / 1e9:.2f} GB")

    # Step 8: Quick inference test
    log("Step 8: Testing inference...")
    batch_converter = alphabet.get_batch_converter()

    test_seq = "EVQLVESGGGLVQPGGSLRL"  # Short test
    data = [("test", test_seq)]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_tokens = batch_tokens.to(device)

    start = time.time()
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=True)
    log(f"  ✓ Inference completed in {time.time() - start:.1f}s")

    contacts = results["contacts"][0]
    log(f"  Contact map shape: {contacts.shape}")

    if device == "cuda":
        log(f"  Peak GPU Memory: {torch.cuda.max_memory_allocated() / 1e9:.2f} GB")

    log("=" * 60)
    log("SUCCESS! ESM-1b is working correctly.")
    log("=" * 60)

    # Check cache after download
    if os.path.exists(cache_dir):
        for f in os.listdir(cache_dir):
            path = os.path.join(cache_dir, f)
            size = os.path.getsize(path) / 1e9
            log(f"  Cached: {f} ({size:.2f} GB)")

if __name__ == "__main__":
    main()
