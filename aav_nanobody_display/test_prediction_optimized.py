#!/usr/bin/env python3
"""
Optimized VHH structure prediction using ESMFold from Hugging Face
With explicit GPU monitoring and memory optimization
"""
import logging
import torch
from pathlib import Path
from transformers import AutoTokenizer, EsmForProteinFolding
import sys

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def print_gpu_memory():
    """Print current GPU memory usage."""
    if torch.cuda.is_available():
        allocated = torch.cuda.memory_allocated() / 1e9
        reserved = torch.cuda.memory_reserved() / 1e9
        print(f"  GPU Memory: {allocated:.2f} GB allocated, {reserved:.2f} GB reserved")

def predict_vhh_structure(sequence: str, output_path: str = "predictions/test_vhh_optimized.pdb"):
    """
    Predict VHH structure using ESMFold with optimizations for RTX 3070 Ti.
    """
    print("=" * 80)
    print("VHH Structure Prediction (Optimized for RTX 3070 Ti)")
    print("=" * 80)
    print(f"\nSequence length: {len(sequence)} amino acids")
    print(f"Sequence: {sequence[:50]}...")

    # Check GPU
    if not torch.cuda.is_available():
        print("\nERROR: CUDA not available!")
        return None

    device = "cuda"
    gpu_name = torch.cuda.get_device_name(0)
    gpu_mem = torch.cuda.get_device_properties(0).total_memory / 1e9

    print(f"\nDevice: {device}")
    print(f"GPU: {gpu_name}")
    print(f"Total VRAM: {gpu_mem:.1f} GB")

    # Test GPU with simple tensor
    print("\nTesting GPU allocation...")
    test_tensor = torch.randn(100, 100).cuda()
    print(f"✓ GPU test successful - tensor on {test_tensor.device}")
    del test_tensor
    torch.cuda.empty_cache()

    print("\n" + "=" * 80)
    print("Loading ESMFold model (optimized)...")
    print("=" * 80)

    # Load with memory optimizations
    print("\nStep 1/3: Loading tokenizer...")
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    print("✓ Tokenizer loaded")

    print("\nStep 2/3: Loading model (this may take a few minutes on first run)...")
    print("  - Using torch_dtype=torch.float16 for memory efficiency")
    print("  - Loading directly to GPU")

    try:
        # Load model with FP16 for memory efficiency and direct GPU loading
        model = EsmForProteinFolding.from_pretrained(
            "facebook/esmfold_v1",
            torch_dtype=torch.float16,  # Use half precision to save memory
            low_cpu_mem_usage=True,
            device_map="cuda:0"  # Load directly to GPU
        )
        print("✓ Model loaded")
        print_gpu_memory()

    except Exception as e:
        print(f"\nERROR loading model: {e}")
        return None

    model.eval()
    print("\nStep 3/3: Model ready for prediction")

    print("\n" + "=" * 80)
    print("Running structure prediction on GPU...")
    print("=" * 80)

    # Tokenize
    print("\nTokenizing sequence...")
    tokenized_input = tokenizer([sequence], return_tensors="pt", add_special_tokens=False)
    tokenized_input = {k: v.cuda() for k, v in tokenized_input.items()}
    print("✓ Input prepared")
    print_gpu_memory()

    # Predict with explicit GPU usage monitoring
    print("\nRunning inference on GPU (you should see GPU usage in nvidia-smi now)...")
    with torch.no_grad():
        output = model(**tokenized_input)

    print("✓ Prediction complete!")
    print_gpu_memory()

    # Convert to PDB
    print("\nConverting to PDB format...")
    pdb_string = model.output_to_pdb(output)[0]

    # Save
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as f:
        f.write(pdb_string)

    print("\n" + "=" * 80)
    print("SUCCESS!")
    print("=" * 80)
    print(f"Structure saved to: {output_file.absolute()}")
    print(f"\nVisualization options:")
    print(f"  - PyMOL: pymol {output_file}")
    print(f"  - ChimeraX: chimerax {output_file}")
    print(f"  - Online: https://www.rcsb.org/3d-view")

    # Cleanup
    torch.cuda.empty_cache()

    return str(output_file.absolute())

if __name__ == "__main__":
    vhh_sequence = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSTADMGWFRQAPGKGRELVAAVSGSGFSTYSDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAKATISLYYAMDVWGQGTTVTVSS"

    try:
        pdb_path = predict_vhh_structure(vhh_sequence)
        if pdb_path:
            print(f"\n✓ All done! Check {pdb_path}")
        else:
            print("\n✗ Prediction failed")
            sys.exit(1)
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
