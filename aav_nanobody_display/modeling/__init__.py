"""
Modeling module for AAV nanobody display pipeline.

Submodules:
- nanobody: Predict nanobody structures (ESMFold, ColabFold)
- fusion: Create and model VP-nanobody fusion constructs
- dock: Quick placement for visualization (no modeling)
"""

from .nanobody import (
    NanobodySequence,
    ESMFoldPredictor,
    ColabFoldPredictor,
    predict_nanobody_structure,
)

from .fusion import (
    FusionConstruct,
    FusionPredictor,
    create_vr_loop_fusion,
    create_vp2_nterm_fusion,
    predict_fusion_structure,
)

from .dock import (
    place_nanobody_at_site,
    place_multiple_nanobodies,
)

__all__ = [
    # Nanobody prediction
    "NanobodySequence",
    "ESMFoldPredictor", 
    "ColabFoldPredictor",
    "predict_nanobody_structure",
    # Fusion modeling
    "FusionConstruct",
    "FusionPredictor",
    "create_vr_loop_fusion",
    "create_vp2_nterm_fusion",
    "predict_fusion_structure",
    # Quick docking
    "place_nanobody_at_site",
    "place_multiple_nanobodies",
]
