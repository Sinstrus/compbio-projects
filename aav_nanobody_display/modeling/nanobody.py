"""
Nanobody structure prediction module.

Supports prediction via:
- ESMFold (fast, local, good for initial models)
- ColabFold/LocalColabFold (higher quality, uses MSA)
- AlphaFold API (if available)
"""

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Literal, Dict, Any
from dataclasses import dataclass
import json

logger = logging.getLogger(__name__)


@dataclass
class NanobodySequence:
    """Represents a nanobody sequence for structure prediction."""
    name: str
    sequence: str
    cdr1: Optional[tuple] = None  # (start, end) residue indices
    cdr2: Optional[tuple] = None
    cdr3: Optional[tuple] = None
    
    def __post_init__(self):
        # Clean sequence
        self.sequence = self.sequence.strip().upper().replace(" ", "").replace("\n", "")
    
    @property
    def length(self) -> int:
        return len(self.sequence)
    
    def to_fasta(self, output_path: Optional[Path] = None) -> str:
        """Export as FASTA format."""
        fasta = f">{self.name}\n{self.sequence}"
        if output_path:
            with open(output_path, "w") as f:
                f.write(fasta)
        return fasta


class ESMFoldPredictor:
    """
    Predict structures using ESMFold.
    
    ESMFold is fast and requires no MSA, making it ideal for quick predictions.
    Works well for single-domain proteins like nanobodies.
    
    Requirements:
        pip install torch fair-esm
    """
    
    def __init__(self, device: str = "auto", chunk_size: int = 64):
        """
        Initialize ESMFold predictor.

        Args:
            device: "cuda" for GPU, "cpu" for CPU (slow), "auto" to auto-detect (default)
            chunk_size: Chunk size for attention operations (default: 64 for RTX 3070 Ti 8GB)
                       Reduce to 32 if OOM, increase to 128 for GPUs with >12GB VRAM
        """
        self.device = device
        self.model = None
        self.default_chunk_size = chunk_size
        self._cuda_checked = False
        
    def _check_cuda(self):
        """Check CUDA availability and provide helpful guidance."""
        if self._cuda_checked:
            return

        try:
            import torch
            cuda_available = torch.cuda.is_available()

            if self.device == "auto":
                self.device = "cuda" if cuda_available else "cpu"
                if cuda_available:
                    logger.info("Auto-detected CUDA - using GPU acceleration")
                else:
                    logger.warning("No CUDA detected - falling back to CPU (50-100x slower)")

            if self.device == "cuda" and not cuda_available:
                raise RuntimeError(
                    "╔════════════════════════════════════════════════════════════╗\n"
                    "║ CUDA REQUESTED BUT NOT AVAILABLE                          ║\n"
                    "╚════════════════════════════════════════════════════════════╝\n"
                    "You are on your Work Laptop (no GPU).\n"
                    "\n"
                    "Options:\n"
                    "  1. Transfer this command to your RTX 3070 Ti home laptop\n"
                    "  2. Use --device cpu (WARNING: 50-100x slower)\n"
                    "  3. Use --device auto to automatically select available device\n"
                    f"\n"
                    f"Current system: CUDA available = {cuda_available}"
                )

            self._cuda_checked = True

        except ImportError:
            logger.warning("PyTorch not installed - cannot check CUDA")

    def _load_model(self):
        """Lazy-load the ESMFold model."""
        if self.model is not None:
            return

        # Check CUDA before loading
        self._check_cuda()

        try:
            import torch
            import esm

            logger.info("Loading ESMFold model (this may take a minute)...")
            self.model = esm.pretrained.esmfold_v1()
            self.model = self.model.eval()

            if self.device == "cuda":
                self.model = self.model.cuda()
                gpu_name = torch.cuda.get_device_name(0)
                gpu_mem = torch.cuda.get_device_properties(0).total_memory / 1e9
                logger.info(f"╔════════════════════════════════════════════════════════════╗")
                logger.info(f"║ ESMFold loaded on GPU: {gpu_name:<30} ║")
                logger.info(f"║ VRAM: {gpu_mem:.1f} GB | chunk_size={self.default_chunk_size:<28} ║")
                logger.info(f"╚════════════════════════════════════════════════════════════╝")
            else:
                logger.info("ESMFold loaded on CPU")
                logger.warning("⚠️  CPU mode is 50-100x slower than GPU")

        except ImportError as e:
            raise ImportError(
                "ESMFold not installed. Install with:\n"
                "  pip install torch fair-esm\n"
                "\n"
                "For CUDA support (RTX 3070 Ti):\n"
                "  pip install torch --index-url https://download.pytorch.org/whl/cu118\n"
                "  pip install fair-esm\n"
                f"\nOriginal error: {e}"
            )
    
    def predict(
        self,
        nanobody: NanobodySequence,
        output_dir: Optional[Path] = None,
        chunk_size: Optional[int] = None
    ) -> Path:
        """
        Predict structure for a nanobody sequence.

        Args:
            nanobody: NanobodySequence object
            output_dir: Directory to save output PDB
            chunk_size: Chunk size for attention (default: use value from __init__)
                       RTX 3070 Ti (8GB): 64 recommended
                       If OOM: reduce to 32

        Returns:
            Path to predicted PDB file
        """
        self._load_model()

        import torch

        if output_dir is None:
            output_dir = Path(tempfile.mkdtemp())
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Use default chunk_size from __init__ if not specified
        if chunk_size is None:
            chunk_size = self.default_chunk_size

        logger.info(f"Predicting structure for {nanobody.name} ({nanobody.length} residues)")
        logger.info(f"Memory optimization: chunk_size={chunk_size}")

        with torch.no_grad():
            # Set chunk size to manage memory
            self.model.set_chunk_size(chunk_size)

            # Run prediction
            output = self.model.infer_pdb(nanobody.sequence)

        # Save PDB
        output_path = output_dir / f"{nanobody.name}_esmfold.pdb"
        with open(output_path, "w") as f:
            f.write(output)

        logger.info(f"Saved prediction to: {output_path}")
        return output_path
    
    def predict_batch(
        self,
        nanobodies: list,
        output_dir: Path,
        chunk_size: Optional[int] = None
    ) -> Dict[str, Path]:
        """Predict structures for multiple nanobodies."""
        results = {}
        for nb in nanobodies:
            try:
                path = self.predict(nb, output_dir, chunk_size)
                results[nb.name] = path
            except Exception as e:
                logger.error(f"Failed to predict {nb.name}: {e}")
                results[nb.name] = None
        return results


class ColabFoldPredictor:
    """
    Predict structures using LocalColabFold.

    Higher quality than ESMFold due to MSA usage, but slower.

    Requirements:
        Install LocalColabFold: https://github.com/YoshitakaMo/localcolabfold
    """

    def __init__(
        self,
        colabfold_path: Optional[Path] = None,
        use_amber: bool = True,
        num_recycle: int = 3,
        msa_mode: str = "mmseqs2+prefilter"
    ):
        """
        Initialize ColabFold predictor.

        Args:
            colabfold_path: Path to colabfold_batch executable
            use_amber: Whether to use Amber relaxation
            num_recycle: Number of recycling iterations
            msa_mode: MSA generation mode
                     "mmseqs2+prefilter" (default, high quality)
                     "single_sequence" (fast, low memory, recommended for 8GB VRAM)
        """
        self.colabfold_path = colabfold_path or self._find_colabfold()
        self.use_amber = use_amber
        self.num_recycle = num_recycle
        self.msa_mode = msa_mode
        
    def _find_colabfold(self) -> Optional[Path]:
        """Try to find colabfold_batch in PATH."""
        import shutil
        path = shutil.which("colabfold_batch")
        if path:
            return Path(path)
        
        # Check common installation locations
        common_paths = [
            Path.home() / "localcolabfold" / "colabfold-conda" / "bin" / "colabfold_batch",
            Path("/usr/local/bin/colabfold_batch"),
        ]
        for p in common_paths:
            if p.exists():
                return p
        
        return None
    
    def is_available(self) -> bool:
        """Check if ColabFold is available."""
        return self.colabfold_path is not None and self.colabfold_path.exists()
    
    def predict(
        self,
        nanobody: NanobodySequence,
        output_dir: Path,
        num_models: int = 1,
        msa_mode: Optional[str] = None
    ) -> Path:
        """
        Predict structure using ColabFold.

        Args:
            nanobody: NanobodySequence object
            output_dir: Directory for output files
            num_models: Number of models to generate (1-5)
            msa_mode: Override default MSA mode
                     "single_sequence" recommended for RTX 3070 Ti (8GB VRAM)

        Returns:
            Path to best predicted structure
        """
        if not self.is_available():
            raise RuntimeError(
                "ColabFold not found. Install LocalColabFold:\n"
                "https://github.com/YoshitakaMo/localcolabfold"
            )

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Write FASTA
        fasta_path = output_dir / f"{nanobody.name}.fasta"
        nanobody.to_fasta(fasta_path)

        # Use instance msa_mode if not overridden
        msa_mode = msa_mode or self.msa_mode

        # Build command
        cmd = [
            str(self.colabfold_path),
            str(fasta_path),
            str(output_dir),
            "--num-models", str(num_models),
            "--num-recycle", str(self.num_recycle),
            "--msa-mode", msa_mode,
        ]

        if self.use_amber:
            cmd.append("--amber")

        logger.info(f"Running ColabFold for {nanobody.name}...")
        logger.info(f"MSA mode: {msa_mode} ({'low memory' if msa_mode == 'single_sequence' else 'high quality'})")
        logger.debug(f"Command: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"ColabFold failed: {result.stderr}")
            raise RuntimeError(f"ColabFold prediction failed: {result.stderr}")

        # Find output PDB (ColabFold names it with _relaxed or _unrelaxed suffix)
        output_files = list(output_dir.glob(f"{nanobody.name}*rank_001*.pdb"))
        if not output_files:
            output_files = list(output_dir.glob(f"{nanobody.name}*.pdb"))

        if not output_files:
            raise RuntimeError(f"No output PDB found for {nanobody.name}")

        best_model = output_files[0]
        logger.info(f"Best model: {best_model}")
        return best_model


def predict_nanobody_structure(
    sequence: str,
    name: str = "nanobody",
    method: Literal["esmfold", "colabfold", "auto"] = "auto",
    output_dir: Optional[Path] = None,
    device: str = "cuda"
) -> Path:
    """
    Convenience function to predict a nanobody structure.
    
    Args:
        sequence: Amino acid sequence
        name: Name for the output files
        method: Prediction method ("esmfold", "colabfold", or "auto")
        output_dir: Where to save results
        device: Device for ESMFold ("cuda" or "cpu")
        
    Returns:
        Path to predicted PDB
    """
    nanobody = NanobodySequence(name=name, sequence=sequence)
    
    if output_dir is None:
        output_dir = Path(tempfile.mkdtemp())
    
    if method == "auto":
        # Try ColabFold first (better quality), fall back to ESMFold
        colabfold = ColabFoldPredictor()
        if colabfold.is_available():
            method = "colabfold"
            logger.info("Using ColabFold (auto-detected)")
        else:
            method = "esmfold"
            logger.info("Using ESMFold (ColabFold not found)")
    
    if method == "colabfold":
        predictor = ColabFoldPredictor()
        return predictor.predict(nanobody, output_dir)
    else:
        predictor = ESMFoldPredictor(device=device)
        return predictor.predict(nanobody, output_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    # Example nanobody sequence (anti-GFP nanobody)
    test_sequence = """
    QVQLVESGGGLVQAGGSLRLSCAASGFPVNRYSMRWYRQAPGKEREWVAGMSSAGDRSSYEDSVKGRFTISRDDARNTVYLQMNSLKPEDTAVYYCNVNVGFEYWGQGTQVTVSS
    """.strip()
    
    print("Testing nanobody structure prediction...")
    print(f"Sequence length: {len(test_sequence.replace(' ', '').replace(chr(10), ''))}")
    
    # Check what's available
    colabfold = ColabFoldPredictor()
    print(f"ColabFold available: {colabfold.is_available()}")
    
    try:
        import esm
        print("ESMFold available: True")
    except ImportError:
        print("ESMFold available: False")
