#!/usr/bin/env python3
"""
Test VHH structure prediction using ESMFold
"""
import logging
import sys
from pathlib import Path

# Add the modeling directory to Python path
sys.path.insert(0, str(Path(__file__).parent))

from modeling.nanobody import ESMFoldPredictor, NanobodySequence

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

def main():
    # User's VHH sequence
    vhh_sequence = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSTADMGWFRQAPGKGRELVAAVSGSGFSTYSDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAKATISLYYAMDVWGQGTTVTVSS"

    print("=" * 80)
    print("VHH Structure Prediction Test")
    print("=" * 80)
    print(f"\nSequence length: {len(vhh_sequence)} amino acids")
    print(f"Sequence: {vhh_sequence[:50]}...")

    # Create NanobodySequence object
    nanobody = NanobodySequence(
        name="test_vhh",
        sequence=vhh_sequence
    )

    # Create output directory
    output_dir = Path("predictions")
    output_dir.mkdir(exist_ok=True)

    # Initialize ESMFold predictor with auto device detection
    print("\n" + "=" * 80)
    print("Initializing ESMFold predictor...")
    print("=" * 80)
    predictor = ESMFoldPredictor(device="auto", chunk_size=64)

    # Run prediction
    print("\n" + "=" * 80)
    print("Running structure prediction...")
    print("=" * 80)
    pdb_path = predictor.predict(nanobody, output_dir=output_dir)

    print("\n" + "=" * 80)
    print("Prediction Complete!")
    print("=" * 80)
    print(f"Structure saved to: {pdb_path}")
    print(f"\nYou can visualize this structure using:")
    print(f"  - PyMOL: pymol {pdb_path}")
    print(f"  - ChimeraX: chimerax {pdb_path}")
    print(f"  - Online: https://www.rcsb.org/3d-view")

    return pdb_path

if __name__ == "__main__":
    try:
        pdb_path = main()
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
