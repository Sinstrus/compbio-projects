#!/usr/bin/env python3
"""
Command-line interface for AAV nanobody display pipeline.

Usage:
    python -m aav_nanobody_display.cli <command> [options]
    
Commands:
    fetch       - Download capsid structures
    predict     - Predict nanobody structures
    fuse        - Create VP-nanobody fusion models
    visualize   - Generate visualizations
    pipeline    - Run full pipeline
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def setup_logging(verbose: bool = False):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%H:%M:%S"
    )


def cmd_fetch(args):
    """Fetch capsid structures from RCSB."""
    from .structures import fetch_capsid
    
    logger.info(f"Fetching {args.serotype} capsid structure...")
    path = fetch_capsid(args.serotype)
    print(f"Downloaded to: {path}")
    
    if args.output:
        import shutil
        shutil.copy(path, args.output)
        print(f"Copied to: {args.output}")


def cmd_predict(args):
    """Predict nanobody structure."""
    from .modeling.nanobody import (
        NanobodySequence,
        ESMFoldPredictor,
        ColabFoldPredictor
    )

    # Read sequence from file or argument
    if args.sequence_file:
        with open(args.sequence_file) as f:
            content = f.read()
            # Handle FASTA format
            if content.startswith(">"):
                lines = content.strip().split("\n")
                sequence = "".join(lines[1:])
            else:
                sequence = content.strip()
    else:
        sequence = args.sequence

    if not sequence:
        logger.error("No sequence provided. Use --sequence or --sequence-file")
        sys.exit(1)

    output_dir = Path(args.output_dir) if args.output_dir else Path("./predictions")
    output_dir.mkdir(parents=True, exist_ok=True)

    nanobody = NanobodySequence(name=args.name, sequence=sequence)

    logger.info(f"╔════════════════════════════════════════════════════════════╗")
    logger.info(f"║ Nanobody Structure Prediction                             ║")
    logger.info(f"╚════════════════════════════════════════════════════════════╝")
    logger.info(f"Sequence: {args.name} ({nanobody.length} residues)")
    logger.info(f"Method: {args.method}")
    logger.info(f"Device: {args.device}")

    # Auto-select method if needed
    method = args.method
    if method == "auto":
        colabfold = ColabFoldPredictor(msa_mode=args.msa_mode)
        if colabfold.is_available():
            method = "colabfold"
            logger.info("Auto-detected: Using ColabFold (higher quality)")
        else:
            method = "esmfold"
            logger.info("Auto-detected: Using ESMFold (ColabFold not found)")

    # Run prediction
    if method == "colabfold":
        predictor = ColabFoldPredictor(msa_mode=args.msa_mode)
        path = predictor.predict(nanobody, output_dir)
    else:
        predictor = ESMFoldPredictor(device=args.device, chunk_size=args.chunk_size)
        path = predictor.predict(nanobody, output_dir)

    print(f"\n{'='*60}")
    print(f"✓ Predicted structure saved to: {path}")
    print(f"{'='*60}")


def cmd_fuse(args):
    """Create and predict VP-nanobody fusion."""
    from .modeling import predict_fusion_structure
    import yaml
    
    # Load VP sequence from config or file
    if args.vp_sequence_file:
        with open(args.vp_sequence_file) as f:
            vp_sequence = f.read().strip()
    else:
        # Use default AAV9 VP3 sequence
        from .structures import fetch_capsid, parse_capsid
        pdb_path = fetch_capsid(args.serotype)
        _, monomer = parse_capsid(pdb_path, args.serotype)
        vp_sequence = monomer.sequence
    
    # Load nanobody sequence
    with open(args.nanobody_file) as f:
        content = f.read()
        if content.startswith(">"):
            lines = content.strip().split("\n")
            nb_sequence = "".join(lines[1:])
        else:
            nb_sequence = content.strip()
    
    output_dir = Path(args.output_dir) if args.output_dir else Path("./fusions")
    
    logger.info(f"Creating fusion: {args.insertion_site}")
    logger.info(f"VP length: {len(vp_sequence)}, Nb length: {len(nb_sequence)}")
    
    path, construct = predict_fusion_structure(
        vp_sequence=vp_sequence,
        nanobody_sequence=nb_sequence,
        insertion_site=args.insertion_site,
        name=args.name,
        output_dir=output_dir,
        method=args.method,
        linker=args.linker
    )
    
    print(f"Fusion structure saved to: {path}")
    print(f"Total fusion length: {construct.length} residues")


def cmd_visualize(args):
    """Generate visualization."""
    from .visualization import CapsidVisualizer, VisualizationConfig
    
    config = VisualizationConfig(
        width=args.width,
        height=args.height,
        color_scheme=args.color_scheme,
        show_surface=args.surface,
        show_labels=not args.no_labels,
    )
    
    viz = CapsidVisualizer(config)
    
    pdb_path = Path(args.pdb)
    output_path = Path(args.output) if args.output else pdb_path.with_suffix(".html")
    
    # Parse VR loops if provided
    vr_loops = None
    if args.highlight_loops:
        import yaml
        config_path = Path(__file__).parent / "config" / "serotypes.yaml"
        with open(config_path) as f:
            serotype_config = yaml.safe_load(f)
        
        if args.serotype in serotype_config["serotypes"]:
            all_loops = serotype_config["serotypes"][args.serotype].get("vr_loops", {})
            vr_loops = {
                name: (info["start"], info["end"])
                for name, info in all_loops.items()
                if name in args.highlight_loops or "all" in args.highlight_loops
            }
    
    logger.info(f"Creating visualization: {pdb_path}")
    viz.visualize_monomer(pdb_path, vr_loops=vr_loops, output_path=output_path)
    
    print(f"Visualization saved to: {output_path}")
    print(f"Open in browser to view interactive 3D structure")


def cmd_pipeline(args):
    """Run full visualization pipeline."""
    from .structures import fetch_capsid, parse_capsid
    from .modeling import predict_nanobody_structure, place_nanobody_at_site
    from .visualization import CapsidVisualizer
    import yaml
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Fetch capsid
    logger.info(f"Step 1: Fetching {args.serotype} capsid...")
    capsid_path = fetch_capsid(args.serotype)
    
    # Step 2: Parse and get VR loop info
    logger.info("Step 2: Parsing capsid structure...")
    config_path = Path(__file__).parent / "config" / "serotypes.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)
    serotype_config = config["serotypes"][args.serotype]
    
    parser, monomer = parse_capsid(capsid_path, args.serotype)
    
    # Step 3: Predict nanobody structure (if sequence provided)
    if args.nanobody_file:
        logger.info("Step 3: Predicting nanobody structure...")
        with open(args.nanobody_file) as f:
            content = f.read()
            if content.startswith(">"):
                lines = content.strip().split("\n")
                nb_sequence = "".join(lines[1:])
            else:
                nb_sequence = content.strip()
        
        nb_path = predict_nanobody_structure(
            sequence=nb_sequence,
            name="nanobody",
            method=args.method,
            output_dir=output_dir
        )
        
        # Step 4: Place nanobody at site
        logger.info(f"Step 4: Placing nanobody at {args.insertion_site}...")
        
        # Get insertion residue
        if args.insertion_site.startswith("VR"):
            loop_info = serotype_config["vr_loops"].get(args.insertion_site)
            if loop_info:
                insertion_res = (loop_info["start"] + loop_info["end"]) // 2
            else:
                logger.error(f"Unknown loop: {args.insertion_site}")
                sys.exit(1)
        else:
            insertion_res = int(args.insertion_site)
        
        combined_path = place_nanobody_at_site(
            nanobody_pdb=nb_path,
            capsid_pdb=capsid_path,
            insertion_residue=insertion_res,
            output_path=output_dir / "capsid_with_nanobody.pdb"
        )
    else:
        combined_path = None
    
    # Step 5: Generate visualizations
    logger.info("Step 5: Generating visualizations...")
    viz = CapsidVisualizer()
    
    # Capsid only
    vr_loops = {
        name: (info["start"], info["end"])
        for name, info in serotype_config.get("vr_loops", {}).items()
    }
    viz.visualize_monomer(
        capsid_path,
        vr_loops=vr_loops,
        output_path=output_dir / "capsid_only.html"
    )
    
    # With nanobody (if created)
    if combined_path:
        viz.visualize_capsid_with_nanobody(
            combined_path,
            output_path=output_dir / "capsid_with_nanobody.html"
        )
        
        # Comparison view
        viz.compare_with_without_nanobody(
            capsid_path,
            combined_path,
            output_path=output_dir / "comparison.html"
        )
    
    print(f"\n{'='*50}")
    print("Pipeline complete!")
    print(f"{'='*50}")
    print(f"Output directory: {output_dir}")
    print(f"Files generated:")
    for f in output_dir.glob("*.html"):
        print(f"  - {f.name}")
    for f in output_dir.glob("*.pdb"):
        print(f"  - {f.name}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="AAV Nanobody Display Visualization Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Fetch AAV9 capsid structure
  python -m aav_nanobody_display.cli fetch --serotype AAV9
  
  # Predict nanobody structure
  python -m aav_nanobody_display.cli predict --sequence-file nb.fasta --name my_nanobody
  
  # Visualize capsid with VR loop highlighting
  python -m aav_nanobody_display.cli visualize --pdb capsid.pdb --highlight-loops VR-VIII
  
  # Run full pipeline
  python -m aav_nanobody_display.cli pipeline --serotype AAV9 --nanobody-file nb.fasta --insertion-site VR-VIII
"""
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # Fetch command
    fetch_parser = subparsers.add_parser("fetch", help="Download capsid structures")
    fetch_parser.add_argument("--serotype", default="AAV9", help="AAV serotype (default: AAV9)")
    fetch_parser.add_argument("--output", "-o", help="Output path")
    
    # Predict command
    predict_parser = subparsers.add_parser("predict", help="Predict nanobody structure")
    predict_parser.add_argument("--sequence", "-s", help="Amino acid sequence")
    predict_parser.add_argument("--sequence-file", "-f", help="FASTA file with sequence")
    predict_parser.add_argument("--name", "-n", default="nanobody", help="Output name")
    predict_parser.add_argument("--method", "-m", default="auto",
                               choices=["auto", "esmfold", "colabfold"],
                               help="Prediction method (auto=detect best available)")
    predict_parser.add_argument("--output-dir", "-o", help="Output directory")
    predict_parser.add_argument("--device", default="auto", choices=["auto", "cuda", "cpu"],
                               help="Device for ESMFold (auto=detect GPU, cuda=require GPU, cpu=force CPU)")
    predict_parser.add_argument("--chunk-size", type=int, default=64,
                               help="ESMFold chunk size (default: 64 for 8GB VRAM, reduce to 32 if OOM)")
    predict_parser.add_argument("--msa-mode", default="single_sequence",
                               choices=["mmseqs2+prefilter", "single_sequence"],
                               help="ColabFold MSA mode (single_sequence recommended for 8GB VRAM)")
    
    # Fuse command  
    fuse_parser = subparsers.add_parser("fuse", help="Create VP-nanobody fusion")
    fuse_parser.add_argument("--nanobody-file", "-nb", required=True, help="Nanobody FASTA")
    fuse_parser.add_argument("--serotype", default="AAV9", help="AAV serotype")
    fuse_parser.add_argument("--vp-sequence-file", help="VP sequence file (optional)")
    fuse_parser.add_argument("--insertion-site", "-i", required=True,
                            help="Insertion site: N-term, C-term, or residue number")
    fuse_parser.add_argument("--linker", default="GSG", help="Linker sequence")
    fuse_parser.add_argument("--name", "-n", default="fusion", help="Output name")
    fuse_parser.add_argument("--method", "-m", default="auto",
                            choices=["auto", "esmfold", "colabfold"])
    fuse_parser.add_argument("--output-dir", "-o", help="Output directory")
    
    # Visualize command
    viz_parser = subparsers.add_parser("visualize", help="Generate visualization")
    viz_parser.add_argument("--pdb", "-p", required=True, help="PDB file to visualize")
    viz_parser.add_argument("--output", "-o", help="Output HTML path")
    viz_parser.add_argument("--serotype", default="AAV9", help="Serotype for VR loop info")
    viz_parser.add_argument("--highlight-loops", nargs="+", help="VR loops to highlight")
    viz_parser.add_argument("--width", type=int, default=1000, help="Viewer width")
    viz_parser.add_argument("--height", type=int, default=800, help="Viewer height")
    viz_parser.add_argument("--color-scheme", default="capsid_default",
                           choices=["capsid_default", "publication", "colorblind_safe"])
    viz_parser.add_argument("--surface", action="store_true", help="Show surface")
    viz_parser.add_argument("--no-labels", action="store_true", help="Hide labels")
    
    # Pipeline command
    pipe_parser = subparsers.add_parser("pipeline", help="Run full pipeline")
    pipe_parser.add_argument("--serotype", default="AAV9", help="AAV serotype")
    pipe_parser.add_argument("--nanobody-file", "-nb", help="Nanobody FASTA (optional)")
    pipe_parser.add_argument("--insertion-site", "-i", default="VR-VIII",
                            help="Insertion site")
    pipe_parser.add_argument("--method", "-m", default="auto",
                            choices=["auto", "esmfold", "colabfold"])
    pipe_parser.add_argument("--output-dir", "-o", default="./pipeline_output",
                            help="Output directory")
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    setup_logging(args.verbose)
    
    # Dispatch to command handler
    commands = {
        "fetch": cmd_fetch,
        "predict": cmd_predict,
        "fuse": cmd_fuse,
        "visualize": cmd_visualize,
        "pipeline": cmd_pipeline,
    }
    
    try:
        commands[args.command](args)
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
