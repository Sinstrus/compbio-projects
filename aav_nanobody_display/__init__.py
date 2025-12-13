"""
AAV Nanobody Display Visualization Pipeline

A toolkit for visualizing AAV capsids with nanobody insertions at various sites:
- VP2 N-terminus
- VR loops (VR-I through VR-IX, especially VR-VIII)
- Full 60-mer capsid assembly

Features:
- Fetch capsid structures from RCSB PDB
- Predict nanobody structures (ESMFold, ColabFold)
- Model VP-nanobody fusions
- Generate interactive 3D visualizations (py3Dmol)

Quick Start:
    from aav_nanobody_display import (
        fetch_capsid,
        predict_nanobody_structure,
        visualize_capsid,
    )
    
    # Fetch AAV9 capsid
    capsid_path = fetch_capsid("AAV9")
    
    # Predict nanobody structure
    nb_path = predict_nanobody_structure(
        "QVQLVESGGGLVQAGGSLRLSCAAS...",
        name="my_nanobody"
    )
    
    # Visualize
    visualize_capsid(capsid_path, output_html="view.html")

Command-line usage:
    python -m aav_nanobody_display.cli --help
"""

__version__ = "0.1.0"
__author__ = "AAV Nanobody Display Pipeline"

# Convenience imports
from .structures import fetch_capsid, parse_capsid, CapsidParser
from .modeling import (
    predict_nanobody_structure,
    predict_fusion_structure,
    NanobodySequence,
    FusionConstruct,
)
from .visualization import visualize_capsid, CapsidVisualizer

__all__ = [
    # Structure handling
    "fetch_capsid",
    "parse_capsid",
    "CapsidParser",
    # Modeling
    "predict_nanobody_structure",
    "predict_fusion_structure",
    "NanobodySequence",
    "FusionConstruct",
    # Visualization
    "visualize_capsid",
    "CapsidVisualizer",
]
