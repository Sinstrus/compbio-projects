#!/usr/bin/env python3
"""
Test VHH structure prediction using ESMFold from Hugging Face Transformers
"""
import logging
import torch
from pathlib import Path
from transformers import AutoTokenizer, EsmForProteinFolding

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def predict_vhh_structure(sequence: str, output_path: str = "predictions/test_vhh_hf.pdb"):
    """
    Predict VHH structure using ESMFold from Hugging Face.

    Args:
        sequence: Amino acid sequence
        output_path: Where to save the PDB file
    """
    print("=" * 80)
    print("VHH Structure Prediction Test (Hugging Face ESMFold)")
    print("=" * 80)
    print(f"\nSequence length: {len(sequence)} amino acids")
    print(f"Sequence: {sequence[:50]}...")

    # Check GPU availability
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"\nDevice: {device}")
    if device == "cuda":
        gpu_name = torch.cuda.get_device_name(0)
        gpu_mem = torch.cuda.get_device_properties(0).total_memory / 1e9
        print(f"GPU: {gpu_name}")
        print(f"VRAM: {gpu_mem:.1f} GB")

    print("\n" + "=" * 80)
    print("Loading ESMFold model from Hugging Face...")
    print("(This will download ~3GB on first run)")
    print("=" * 80)

    # Load model and tokenizer
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)

    # Move to GPU if available
    if device == "cuda":
        model = model.cuda()
        print("\nModel loaded on GPU")
    else:
        print("\nModel loaded on CPU (this will be slow)")

    model.eval()

    print("\n" + "=" * 80)
    print("Running structure prediction...")
    print("=" * 80)

    # Tokenize sequence
    tokenized_input = tokenizer([sequence], return_tensors="pt", add_special_tokens=False)

    # Move input to same device as model
    if device == "cuda":
        tokenized_input = {k: v.cuda() for k, v in tokenized_input.items()}

    # Predict structure
    with torch.no_grad():
        output = model(**tokenized_input)

    # Extract PDB string
    pdb_string = model.output_to_pdb(output)[0]

    # Save PDB file
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as f:
        f.write(pdb_string)

    print("\n" + "=" * 80)
    print("Prediction Complete!")
    print("=" * 80)
    print(f"Structure saved to: {output_file}")
    print(f"\nYou can visualize this structure using:")
    print(f"  - PyMOL: pymol {output_file}")
    print(f"  - ChimeraX: chimerax {output_file}")
    print(f"  - Online: Upload to https://www.rcsb.org/3d-view")

    return str(output_file)

if __name__ == "__main__":
    # User's VHH sequence
    vhh_sequence = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSTADMGWFRQAPGKGRELVAAVSGSGFSTYSDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAKATISLYYAMDVWGQGTTVTVSS"

    try:
        pdb_path = predict_vhh_structure(vhh_sequence)
        print(f"\nSuccess! Structure prediction completed.")
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        exit(1)
