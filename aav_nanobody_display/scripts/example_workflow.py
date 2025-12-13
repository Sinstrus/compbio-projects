#!/usr/bin/env python3
"""
Example: Complete workflow for AAV9 capsid with nanobody at VR-VIII

This script demonstrates the full pipeline:
1. Fetch AAV9 capsid structure
2. Analyze VR loops
3. (Optional) Predict nanobody structure
4. Generate visualizations

Run with:
    python scripts/example_workflow.py
"""

import logging
from pathlib import Path
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)
logger = logging.getLogger(__name__)


def main():
    """Run example workflow."""
    
    # Create output directory
    output_dir = Path("example_output")
    output_dir.mkdir(exist_ok=True)
    
    # =========================================
    # Step 1: Fetch and Parse AAV9 Capsid
    # =========================================
    logger.info("=" * 50)
    logger.info("Step 1: Fetching AAV9 capsid structure")
    logger.info("=" * 50)
    
    from structures import fetch_capsid, parse_capsid
    
    capsid_path = fetch_capsid("AAV9")
    logger.info(f"Downloaded to: {capsid_path}")
    
    parser, monomer = parse_capsid(capsid_path, "AAV9")
    
    logger.info(f"\nMonomer information:")
    logger.info(f"  Chain: {monomer.chain_id}")
    logger.info(f"  Residues: {monomer.residue_start} - {monomer.residue_end}")
    logger.info(f"  Sequence length: {monomer.length}")
    
    logger.info(f"\nVR Loops identified:")
    for name, loop in monomer.vr_loops.items():
        logger.info(f"  {name}: residues {loop.start_residue}-{loop.end_residue} "
                   f"({loop.length} residues)")
    
    # =========================================
    # Step 2: Visualize Capsid with VR Loops
    # =========================================
    logger.info("\n" + "=" * 50)
    logger.info("Step 2: Creating capsid visualization")
    logger.info("=" * 50)
    
    from visualization import CapsidVisualizer, VisualizationConfig
    
    config = VisualizationConfig(
        width=1000,
        height=800,
        color_scheme="publication",
        show_labels=True,
    )
    
    viz = CapsidVisualizer(config)
    
    # Highlight all VR loops
    vr_loops = {
        name: (loop.start_residue, loop.end_residue)
        for name, loop in monomer.vr_loops.items()
    }
    
    html_path = output_dir / "aav9_monomer_vr_loops.html"
    viz.visualize_monomer(
        capsid_path,
        vr_loops=vr_loops,
        output_path=html_path
    )
    logger.info(f"Saved visualization to: {html_path}")
    
    # =========================================
    # Step 3: Get VR-VIII Insertion Site Details
    # =========================================
    logger.info("\n" + "=" * 50)
    logger.info("Step 3: Analyzing VR-VIII insertion site")
    logger.info("=" * 50)
    
    vr8_info = parser.get_vr_loop_coordinates("VR-VIII")
    logger.info(f"VR-VIII details:")
    logger.info(f"  Start residue: {vr8_info['start_residue']}")
    logger.info(f"  End residue: {vr8_info['end_residue']}")
    logger.info(f"  Midpoint residue: {vr8_info['midpoint_residue']}")
    if 'n_anchor_coord' in vr8_info:
        logger.info(f"  N-anchor coord: {vr8_info['n_anchor_coord']}")
    if 'c_anchor_coord' in vr8_info:
        logger.info(f"  C-anchor coord: {vr8_info['c_anchor_coord']}")
    
    # =========================================
    # Step 4: (Optional) Predict Nanobody Structure
    # =========================================
    logger.info("\n" + "=" * 50)
    logger.info("Step 4: Nanobody structure prediction")
    logger.info("=" * 50)
    
    # Example anti-GFP nanobody sequence
    example_nanobody = """
    QVQLVESGGGLVQAGGSLRLSCAASGFPVNRYSMRWYRQAPGKEREWVAGMSSAGDRSSYEDSVKGRFTISRDDARNTVYLQMNSLKPEDTAVYYCNVNVGFEYWGQGTQVTVSS
    """.strip().replace("\n", "").replace(" ", "")
    
    logger.info(f"Example nanobody sequence ({len(example_nanobody)} residues)")
    
    # Check if ESMFold is available
    try:
        import esm
        esmfold_available = True
        logger.info("ESMFold is available!")
    except ImportError:
        esmfold_available = False
        logger.info("ESMFold not installed - skipping prediction")
        logger.info("Install with: pip install fair-esm torch")
    
    if esmfold_available:
        from modeling import predict_nanobody_structure
        
        logger.info("Predicting nanobody structure (this may take ~30 seconds)...")
        try:
            nb_path = predict_nanobody_structure(
                example_nanobody,
                name="anti_GFP_nanobody",
                method="esmfold",
                output_dir=output_dir
            )
            logger.info(f"Predicted structure saved to: {nb_path}")
            
            # Visualize the nanobody
            nb_html = output_dir / "nanobody_structure.html"
            viz.visualize_monomer(nb_path, output_path=nb_html)
            logger.info(f"Nanobody visualization saved to: {nb_html}")
            
        except Exception as e:
            logger.error(f"Prediction failed: {e}")
            logger.info("This might be due to GPU memory - try with device='cpu'")
    
    # =========================================
    # Step 5: Create Fusion Construct (Sequence Only)
    # =========================================
    logger.info("\n" + "=" * 50)
    logger.info("Step 5: Creating VP-nanobody fusion construct")
    logger.info("=" * 50)
    
    from modeling import create_vr_loop_fusion, FusionPredictor
    
    # Create fusion construct at VR-VIII
    fusion = create_vr_loop_fusion(
        vp_sequence=monomer.sequence,
        nanobody_sequence=example_nanobody,
        loop_start=585,  # VR-VIII start
        loop_end=596,    # VR-VIII end
        name="AAV9_VR8_antiGFP",
        n_linker="GSG",
        c_linker="GSG"
    )
    
    logger.info(f"Fusion construct: {fusion.name}")
    logger.info(f"Total length: {fusion.length} residues")
    
    # Check memory requirements
    predictor = FusionPredictor()
    mem_info = predictor.estimate_memory_requirement(fusion.length)
    logger.info(f"Memory estimate: {mem_info['estimated_vram_gb']} GB VRAM")
    logger.info(f"Feasible on 8GB GPU: {mem_info['feasible_on_8gb_gpu']}")
    
    # Save fusion sequence
    fusion_fasta = output_dir / f"{fusion.name}.fasta"
    fusion.to_fasta(fusion_fasta)
    logger.info(f"Saved fusion sequence to: {fusion_fasta}")
    
    # =========================================
    # Summary
    # =========================================
    logger.info("\n" + "=" * 50)
    logger.info("Workflow Complete!")
    logger.info("=" * 50)
    logger.info(f"\nOutput files in: {output_dir.absolute()}")
    for f in output_dir.glob("*"):
        logger.info(f"  - {f.name}")
    
    logger.info("\nNext steps:")
    logger.info("  1. Open HTML files in browser for interactive 3D views")
    logger.info("  2. Use the fusion FASTA for structure prediction")
    logger.info("  3. Modify this script for your specific nanobody sequences")


if __name__ == "__main__":
    main()
