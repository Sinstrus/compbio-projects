"""
Structure parsing utilities for AAV capsid analysis.

Extracts chains, identifies VR loops, and parses capsid-specific information.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
import numpy as np

logger = logging.getLogger(__name__)

# Try to import BioPython
try:
    from Bio.PDB import PDBParser, MMCIFParser, Structure, Chain, Residue
    from Bio.PDB.Polypeptide import is_aa, three_to_one
    from Bio.PDB import NeighborSearch
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    logger.warning("BioPython not installed. Some features will be unavailable.")


@dataclass
class VRLoop:
    """Represents a variable region loop on the AAV capsid."""
    name: str
    start_residue: int
    end_residue: int
    chain_id: str
    sequence: str = ""
    center_of_mass: Optional[np.ndarray] = None
    surface_exposed: bool = True
    
    @property
    def length(self) -> int:
        return self.end_residue - self.start_residue + 1
    
    @property
    def midpoint_residue(self) -> int:
        return (self.start_residue + self.end_residue) // 2


@dataclass
class CapsidMonomer:
    """Represents a single VP monomer from the capsid."""
    chain_id: str
    sequence: str
    residue_start: int
    residue_end: int
    vr_loops: Dict[str, VRLoop] = field(default_factory=dict)
    coordinates: Optional[np.ndarray] = None  # CA coordinates
    
    @property
    def length(self) -> int:
        return len(self.sequence)


class CapsidParser:
    """Parse AAV capsid structures and extract relevant information."""
    
    def __init__(self, structure_path: Path, serotype_config: Dict[str, Any]):
        """
        Initialize the parser.
        
        Args:
            structure_path: Path to PDB or mmCIF file
            serotype_config: Configuration dict for this serotype from serotypes.yaml
        """
        if not HAS_BIOPYTHON:
            raise ImportError("BioPython is required for structure parsing. "
                            "Install with: pip install biopython")
        
        self.structure_path = Path(structure_path)
        self.config = serotype_config
        self.structure = self._load_structure()
        
    def _load_structure(self) -> "Structure":
        """Load the structure from file."""
        if self.structure_path.suffix.lower() == ".cif":
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
        
        structure = parser.get_structure(
            self.structure_path.stem, 
            str(self.structure_path)
        )
        logger.info(f"Loaded structure with {len(list(structure.get_chains()))} chains")
        return structure
    
    def get_chains(self) -> List[str]:
        """Get list of chain IDs in the structure."""
        model = self.structure[0]  # First model
        return [chain.id for chain in model]
    
    def extract_monomer(self, chain_id: str = "A") -> CapsidMonomer:
        """
        Extract a single VP monomer from the structure.
        
        Args:
            chain_id: Chain identifier
            
        Returns:
            CapsidMonomer object with sequence and VR loop info
        """
        model = self.structure[0]
        chain = model[chain_id]
        
        # Extract sequence and coordinates
        sequence = []
        ca_coords = []
        residue_ids = []
        
        for residue in chain:
            if is_aa(residue, standard=True):
                try:
                    one_letter = three_to_one(residue.get_resname())
                    sequence.append(one_letter)
                    residue_ids.append(residue.id[1])
                    
                    # Get CA coordinates
                    if "CA" in residue:
                        ca_coords.append(residue["CA"].get_coord())
                except KeyError:
                    continue
        
        monomer = CapsidMonomer(
            chain_id=chain_id,
            sequence="".join(sequence),
            residue_start=min(residue_ids) if residue_ids else 0,
            residue_end=max(residue_ids) if residue_ids else 0,
            coordinates=np.array(ca_coords) if ca_coords else None
        )
        
        # Extract VR loops
        if "vr_loops" in self.config:
            for loop_name, loop_info in self.config["vr_loops"].items():
                vr_loop = self._extract_vr_loop(
                    chain, loop_name, 
                    loop_info["start"], 
                    loop_info["end"]
                )
                if vr_loop:
                    monomer.vr_loops[loop_name] = vr_loop
        
        logger.info(f"Extracted monomer from chain {chain_id}: "
                   f"{monomer.length} residues, {len(monomer.vr_loops)} VR loops")
        return monomer
    
    def _extract_vr_loop(
        self, 
        chain: "Chain", 
        loop_name: str,
        start: int, 
        end: int
    ) -> Optional[VRLoop]:
        """Extract a single VR loop from a chain."""
        sequence = []
        coords = []
        
        for residue in chain:
            res_id = residue.id[1]
            if start <= res_id <= end and is_aa(residue, standard=True):
                try:
                    sequence.append(three_to_one(residue.get_resname()))
                    if "CA" in residue:
                        coords.append(residue["CA"].get_coord())
                except KeyError:
                    continue
        
        if not sequence:
            logger.warning(f"Could not extract {loop_name} (residues {start}-{end})")
            return None
        
        center_of_mass = np.mean(coords, axis=0) if coords else None
        
        return VRLoop(
            name=loop_name,
            start_residue=start,
            end_residue=end,
            chain_id=chain.id,
            sequence="".join(sequence),
            center_of_mass=center_of_mass
        )
    
    def get_vr_loop_coordinates(
        self, 
        loop_name: str, 
        chain_id: str = "A"
    ) -> Dict[str, Any]:
        """
        Get detailed coordinates for a VR loop insertion site.
        
        Returns N-terminal anchor, C-terminal anchor, and suggested
        insertion point coordinates.
        """
        model = self.structure[0]
        chain = model[chain_id]
        
        loop_config = self.config["vr_loops"].get(loop_name)
        if not loop_config:
            raise ValueError(f"Unknown VR loop: {loop_name}")
        
        start = loop_config["start"]
        end = loop_config["end"]
        midpoint = (start + end) // 2
        
        result = {
            "loop_name": loop_name,
            "chain_id": chain_id,
            "start_residue": start,
            "end_residue": end,
            "midpoint_residue": midpoint,
        }
        
        # Get anchor coordinates
        for residue in chain:
            res_id = residue.id[1]
            if res_id == start and "CA" in residue:
                result["n_anchor_coord"] = residue["CA"].get_coord().tolist()
                result["n_anchor_resname"] = residue.get_resname()
            elif res_id == end and "CA" in residue:
                result["c_anchor_coord"] = residue["CA"].get_coord().tolist()
                result["c_anchor_resname"] = residue.get_resname()
            elif res_id == midpoint and "CA" in residue:
                result["midpoint_coord"] = residue["CA"].get_coord().tolist()
        
        return result
    
    def get_icosahedral_operators(self) -> Optional[List[np.ndarray]]:
        """
        Extract icosahedral symmetry operators from the structure.
        
        Returns list of 4x4 transformation matrices, or None if not present.
        """
        # This would parse BIOMT records from PDB or equiv in mmCIF
        # For now, return None - full implementation would read these from header
        logger.warning("Icosahedral operator extraction not yet implemented. "
                      "Use assemble.py for capsid assembly.")
        return None


def parse_capsid(
    structure_path: Path,
    serotype: str = "AAV9"
) -> Tuple[CapsidParser, CapsidMonomer]:
    """
    Convenience function to parse a capsid structure.
    
    Args:
        structure_path: Path to structure file
        serotype: Serotype name for configuration
        
    Returns:
        Tuple of (parser, monomer) for further analysis
    """
    import yaml
    
    config_path = Path(__file__).parent.parent / "config" / "serotypes.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)
    
    serotype_config = config["serotypes"][serotype]
    parser = CapsidParser(structure_path, serotype_config)
    monomer = parser.extract_monomer()
    
    return parser, monomer


if __name__ == "__main__":
    # Test the parser
    logging.basicConfig(level=logging.INFO)
    
    from fetch import fetch_capsid
    
    # Fetch and parse AAV9
    structure_path = fetch_capsid("AAV9")
    parser, monomer = parse_capsid(structure_path, "AAV9")
    
    print(f"\nMonomer info:")
    print(f"  Chain: {monomer.chain_id}")
    print(f"  Residues: {monomer.residue_start}-{monomer.residue_end}")
    print(f"  Sequence length: {monomer.length}")
    print(f"\nVR Loops:")
    for name, loop in monomer.vr_loops.items():
        print(f"  {name}: residues {loop.start_residue}-{loop.end_residue}, "
              f"sequence: {loop.sequence[:20]}...")
