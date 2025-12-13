"""
VP-Nanobody fusion modeling module.

Handles creating fusion constructs and predicting their structures using
AlphaFold Multimer or similar approaches.
"""

import logging
from pathlib import Path
from typing import Optional, Literal, Dict, Any, Tuple
from dataclasses import dataclass
import subprocess
import tempfile

logger = logging.getLogger(__name__)


@dataclass
class FusionConstruct:
    """Represents a VP-nanobody fusion construct."""
    name: str
    vp_sequence: str
    nanobody_sequence: str
    insertion_site: str  # "N-term", "C-term", or residue number for loop insertion
    n_linker: str = "GSG"
    c_linker: str = "GSG"
    
    @property
    def fusion_sequence(self) -> str:
        """Generate the full fusion sequence."""
        if self.insertion_site == "N-term":
            return f"{self.nanobody_sequence}{self.c_linker}{self.vp_sequence}"
        elif self.insertion_site == "C-term":
            return f"{self.vp_sequence}{self.n_linker}{self.nanobody_sequence}"
        else:
            # Loop insertion - need to split VP sequence
            try:
                insert_pos = int(self.insertion_site)
                vp_n = self.vp_sequence[:insert_pos]
                vp_c = self.vp_sequence[insert_pos:]
                return f"{vp_n}{self.n_linker}{self.nanobody_sequence}{self.c_linker}{vp_c}"
            except ValueError:
                raise ValueError(f"Invalid insertion site: {self.insertion_site}")
    
    @property
    def length(self) -> int:
        return len(self.fusion_sequence)
    
    def to_fasta(self, output_path: Optional[Path] = None) -> str:
        """Export as FASTA."""
        fasta = f">{self.name}\n{self.fusion_sequence}"
        if output_path:
            with open(output_path, "w") as f:
                f.write(fasta)
        return fasta


def create_vr_loop_fusion(
    vp_sequence: str,
    nanobody_sequence: str,
    loop_start: int,
    loop_end: int,
    name: str = "fusion",
    n_linker: str = "GSG",
    c_linker: str = "GSG",
    replace_loop: bool = True
) -> FusionConstruct:
    """
    Create a fusion construct with nanobody inserted at a VR loop.
    
    Args:
        vp_sequence: Full VP sequence (VP3 region)
        nanobody_sequence: Nanobody sequence to insert
        loop_start: Start residue of VR loop (VP1 numbering)
        loop_end: End residue of VR loop
        name: Name for the construct
        n_linker: Linker sequence on N-terminal side of nanobody
        c_linker: Linker sequence on C-terminal side of nanobody
        replace_loop: If True, replace the loop; if False, insert at midpoint
        
    Returns:
        FusionConstruct object
    """
    # Convert to 0-indexed for sequence slicing
    # Assuming VP3 starts at residue 203 in most numbering schemes
    vp3_start = 203
    local_start = loop_start - vp3_start
    local_end = loop_end - vp3_start
    
    if replace_loop:
        # Replace the entire loop with nanobody
        vp_n = vp_sequence[:local_start]
        vp_c = vp_sequence[local_end:]
    else:
        # Insert at midpoint of loop
        midpoint = (local_start + local_end) // 2
        vp_n = vp_sequence[:midpoint]
        vp_c = vp_sequence[midpoint:]
    
    # Build fusion sequence manually
    fusion_seq = f"{vp_n}{n_linker}{nanobody_sequence}{c_linker}{vp_c}"
    
    return FusionConstruct(
        name=name,
        vp_sequence=vp_sequence,
        nanobody_sequence=nanobody_sequence,
        insertion_site=str(local_start),
        n_linker=n_linker,
        c_linker=c_linker
    )


def create_vp2_nterm_fusion(
    vp2_nterm_sequence: str,
    nanobody_sequence: str,
    name: str = "VP2_Nb_fusion",
    linker: str = "GGGGS"
) -> FusionConstruct:
    """
    Create a VP2 N-terminal fusion construct.
    
    The nanobody is fused to the N-terminus of VP2, which sits inside
    the capsid and can be externalized.
    
    Args:
        vp2_nterm_sequence: VP2 N-terminal sequence (before VP3 start)
        nanobody_sequence: Nanobody sequence
        name: Construct name
        linker: Flexible linker between nanobody and VP2
        
    Returns:
        FusionConstruct object
    """
    return FusionConstruct(
        name=name,
        vp_sequence=vp2_nterm_sequence,
        nanobody_sequence=nanobody_sequence,
        insertion_site="N-term",
        n_linker="",
        c_linker=linker
    )


class FusionPredictor:
    """
    Predict structures for VP-nanobody fusion constructs.
    
    Uses AlphaFold Multimer or ColabFold for structure prediction.
    """
    
    def __init__(
        self,
        colabfold_path: Optional[Path] = None,
        use_templates: bool = True,
        num_recycle: int = 3
    ):
        """
        Initialize the fusion predictor.
        
        Args:
            colabfold_path: Path to colabfold_batch
            use_templates: Whether to use template structures
            num_recycle: Number of recycling iterations
        """
        self.colabfold_path = colabfold_path
        self.use_templates = use_templates
        self.num_recycle = num_recycle
        
        # Find colabfold if not provided
        if self.colabfold_path is None:
            import shutil
            path = shutil.which("colabfold_batch")
            if path:
                self.colabfold_path = Path(path)
    
    def is_available(self) -> bool:
        """Check if prediction is available."""
        return self.colabfold_path is not None and self.colabfold_path.exists()
    
    def estimate_memory_requirement(self, sequence_length: int) -> Dict[str, Any]:
        """
        Estimate VRAM requirement for prediction.
        
        Based on empirical observations with AlphaFold/ColabFold.
        """
        # Rough estimates - actual may vary
        if sequence_length < 400:
            vram_gb = 8
            feasible_on_3070ti = True
        elif sequence_length < 600:
            vram_gb = 12
            feasible_on_3070ti = True  # Might need to reduce MSA depth
        elif sequence_length < 800:
            vram_gb = 16
            feasible_on_3070ti = False  # Would need --msa-mode single_sequence
        else:
            vram_gb = 24
            feasible_on_3070ti = False
        
        return {
            "sequence_length": sequence_length,
            "estimated_vram_gb": vram_gb,
            "feasible_on_8gb_gpu": feasible_on_3070ti,
            "recommendation": (
                "Use ESMFold for initial model" if not feasible_on_3070ti 
                else "Should work with ColabFold"
            )
        }
    
    def predict_with_colabfold(
        self,
        construct: FusionConstruct,
        output_dir: Path,
        num_models: int = 1,
        msa_mode: Literal["mmseqs2", "single_sequence"] = "mmseqs2"
    ) -> Path:
        """
        Predict fusion structure using ColabFold.
        
        Args:
            construct: FusionConstruct to predict
            output_dir: Output directory
            num_models: Number of models to generate
            msa_mode: MSA mode ("single_sequence" for faster/less memory)
            
        Returns:
            Path to predicted structure
        """
        if not self.is_available():
            raise RuntimeError("ColabFold not available")
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Check memory requirements
        mem_info = self.estimate_memory_requirement(construct.length)
        logger.info(f"Sequence length: {construct.length}, "
                   f"Estimated VRAM: {mem_info['estimated_vram_gb']} GB")
        
        if not mem_info["feasible_on_8gb_gpu"] and msa_mode != "single_sequence":
            logger.warning("Sequence may be too long for 8GB GPU. "
                         "Consider using msa_mode='single_sequence'")
        
        # Write FASTA
        fasta_path = output_dir / f"{construct.name}.fasta"
        construct.to_fasta(fasta_path)
        
        # Build command
        cmd = [
            str(self.colabfold_path),
            str(fasta_path),
            str(output_dir),
            "--num-models", str(num_models),
            "--num-recycle", str(self.num_recycle),
            "--msa-mode", msa_mode,
        ]
        
        if not self.use_templates:
            cmd.append("--templates")
            cmd.append("false")
        
        logger.info(f"Running ColabFold for {construct.name}...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"ColabFold failed: {result.stderr}")
            raise RuntimeError(f"Prediction failed: {result.stderr}")
        
        # Find output
        output_files = list(output_dir.glob(f"{construct.name}*rank_001*.pdb"))
        if not output_files:
            output_files = list(output_dir.glob(f"{construct.name}*.pdb"))
        
        if not output_files:
            raise RuntimeError(f"No output found for {construct.name}")
        
        return output_files[0]
    
    def predict_with_esmfold(
        self,
        construct: FusionConstruct,
        output_dir: Path,
        device: str = "cuda"
    ) -> Path:
        """
        Predict fusion structure using ESMFold.
        
        Faster than ColabFold but may be less accurate for multi-domain proteins.
        """
        from .nanobody import ESMFoldPredictor, NanobodySequence
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # ESMFold just needs a sequence
        predictor = ESMFoldPredictor(device=device)
        
        # Create a NanobodySequence-like object (ESMFold doesn't care about domains)
        seq_obj = NanobodySequence(
            name=construct.name,
            sequence=construct.fusion_sequence
        )
        
        return predictor.predict(seq_obj, output_dir)


def predict_fusion_structure(
    vp_sequence: str,
    nanobody_sequence: str,
    insertion_site: str,
    name: str = "VP_Nb_fusion",
    output_dir: Optional[Path] = None,
    method: Literal["auto", "colabfold", "esmfold"] = "auto",
    linker: str = "GSG"
) -> Tuple[Path, FusionConstruct]:
    """
    Convenience function to create and predict a fusion structure.
    
    Args:
        vp_sequence: VP protein sequence
        nanobody_sequence: Nanobody sequence
        insertion_site: Where to insert ("N-term", "C-term", or residue number)
        name: Construct name
        output_dir: Output directory
        method: Prediction method
        linker: Linker sequence
        
    Returns:
        Tuple of (predicted PDB path, FusionConstruct object)
    """
    if output_dir is None:
        output_dir = Path(tempfile.mkdtemp())
    
    construct = FusionConstruct(
        name=name,
        vp_sequence=vp_sequence,
        nanobody_sequence=nanobody_sequence,
        insertion_site=insertion_site,
        n_linker=linker,
        c_linker=linker
    )
    
    logger.info(f"Created fusion construct: {construct.length} residues")
    
    predictor = FusionPredictor()
    mem_info = predictor.estimate_memory_requirement(construct.length)
    
    if method == "auto":
        if mem_info["feasible_on_8gb_gpu"] and predictor.is_available():
            method = "colabfold"
        else:
            method = "esmfold"
    
    if method == "colabfold":
        msa_mode = "mmseqs2" if mem_info["feasible_on_8gb_gpu"] else "single_sequence"
        pdb_path = predictor.predict_with_colabfold(construct, output_dir, msa_mode=msa_mode)
    else:
        pdb_path = predictor.predict_with_esmfold(construct, output_dir)
    
    return pdb_path, construct


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    # Example VP3 sequence snippet (first 100 residues of AAV9 VP3)
    vp3_snippet = "MAADGYLPDWLEDNLSEGIREWWALKPGAPQPKANQQHQDNARGLVLPGYKYLGPFNGLDKGEPVNAADAAALEHDKAYDQQLKAGDNPYLKYNHADAEFQERLKEDTSF"
    
    # Example nanobody sequence
    nanobody = "QVQLVESGGGLVQAGGSLRLSCAASGFPVNRYSMRWYRQAPGKEREWVAGMSSAGDRSSYEDSVKGRFTISRDDARNTVYLQMNSLKPEDTAVYYCNVNVGFEYWGQGTQVTVSS"
    
    # Create N-terminal fusion
    construct = create_vp2_nterm_fusion(
        vp2_nterm_sequence=vp3_snippet,
        nanobody_sequence=nanobody,
        name="VP2_anti_GFP_Nb"
    )
    
    print(f"Fusion construct: {construct.name}")
    print(f"Total length: {construct.length} residues")
    
    # Check memory requirements
    predictor = FusionPredictor()
    mem_info = predictor.estimate_memory_requirement(construct.length)
    print(f"Memory estimate: {mem_info}")
