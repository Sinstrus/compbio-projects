"""
Structure fetching utilities for AAV capsid visualization pipeline.

Downloads PDB/mmCIF files from RCSB and caches them locally.
"""

import os
import urllib.request
import hashlib
from pathlib import Path
from typing import Optional, Literal
import logging

logger = logging.getLogger(__name__)

# Default cache directory
DEFAULT_CACHE_DIR = Path.home() / ".cache" / "aav_nanobody_display" / "structures"


class StructureFetcher:
    """Fetch and cache protein structures from RCSB PDB."""
    
    RCSB_BASE_URL = "https://files.rcsb.org/download"
    
    def __init__(self, cache_dir: Optional[Path] = None):
        """
        Initialize the structure fetcher.
        
        Args:
            cache_dir: Directory to cache downloaded structures.
                      Defaults to ~/.cache/aav_nanobody_display/structures
        """
        self.cache_dir = Path(cache_dir) if cache_dir else DEFAULT_CACHE_DIR
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Structure cache directory: {self.cache_dir}")
    
    def get_structure_path(
        self, 
        pdb_id: str, 
        format: Literal["pdb", "cif"] = "pdb",
        force_download: bool = False
    ) -> Path:
        """
        Get path to a structure file, downloading if necessary.
        
        Args:
            pdb_id: 4-character PDB ID (e.g., "3UX1")
            format: File format - "pdb" or "cif" (mmCIF)
            force_download: If True, re-download even if cached
            
        Returns:
            Path to the structure file
        """
        pdb_id = pdb_id.upper()
        
        if format == "pdb":
            filename = f"{pdb_id}.pdb"
            url = f"{self.RCSB_BASE_URL}/{pdb_id}.pdb"
        else:
            filename = f"{pdb_id}.cif"
            url = f"{self.RCSB_BASE_URL}/{pdb_id}.cif"
        
        local_path = self.cache_dir / filename
        
        if local_path.exists() and not force_download:
            logger.info(f"Using cached structure: {local_path}")
            return local_path
        
        logger.info(f"Downloading {pdb_id} from RCSB...")
        try:
            urllib.request.urlretrieve(url, local_path)
            logger.info(f"Downloaded to: {local_path}")
        except Exception as e:
            logger.error(f"Failed to download {pdb_id}: {e}")
            raise RuntimeError(f"Could not download PDB {pdb_id}") from e
        
        return local_path
    
    def get_alphafold_structure(
        self,
        uniprot_id: str,
        force_download: bool = False
    ) -> Path:
        """
        Fetch a predicted structure from AlphaFold DB.
        
        Args:
            uniprot_id: UniProt accession (e.g., "Q6JC40" for AAV9)
            force_download: If True, re-download even if cached
            
        Returns:
            Path to the structure file
        """
        filename = f"AF-{uniprot_id}-F1-model_v4.pdb"
        url = f"https://alphafold.ebi.ac.uk/files/{filename}"
        local_path = self.cache_dir / filename
        
        if local_path.exists() and not force_download:
            logger.info(f"Using cached AlphaFold structure: {local_path}")
            return local_path
        
        logger.info(f"Downloading AlphaFold model for {uniprot_id}...")
        try:
            urllib.request.urlretrieve(url, local_path)
            logger.info(f"Downloaded to: {local_path}")
        except Exception as e:
            logger.error(f"Failed to download AlphaFold model: {e}")
            raise RuntimeError(f"Could not download AlphaFold model for {uniprot_id}") from e
        
        return local_path


def fetch_capsid(
    serotype: str = "AAV9",
    format: Literal["pdb", "cif"] = "pdb",
    cache_dir: Optional[Path] = None
) -> Path:
    """
    Convenience function to fetch an AAV capsid structure.
    
    Args:
        serotype: AAV serotype name (e.g., "AAV9", "AAV2")
        format: File format
        cache_dir: Optional cache directory
        
    Returns:
        Path to the structure file
    """
    import yaml
    
    # Load serotype config
    config_path = Path(__file__).parent.parent / "config" / "serotypes.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)
    
    if serotype not in config["serotypes"]:
        available = list(config["serotypes"].keys())
        raise ValueError(f"Unknown serotype '{serotype}'. Available: {available}")
    
    pdb_id = config["serotypes"][serotype]["pdb_id"]
    
    fetcher = StructureFetcher(cache_dir)
    return fetcher.get_structure_path(pdb_id, format)


if __name__ == "__main__":
    # Test the fetcher
    logging.basicConfig(level=logging.INFO)
    
    fetcher = StructureFetcher()
    
    # Fetch AAV9 structure
    path = fetcher.get_structure_path("3UX1")
    print(f"AAV9 structure: {path}")
    
    # Fetch AAV2 structure
    path = fetcher.get_structure_path("1LP3")
    print(f"AAV2 structure: {path}")
