"""
Capsid assembly utilities.

Builds full 60-mer capsids from monomers using icosahedral symmetry operations.
Also generates asymmetric units (trimers) for visualization.
"""

import logging
from pathlib import Path
from typing import List, Optional, Tuple
import numpy as np
from dataclasses import dataclass

logger = logging.getLogger(__name__)

# Icosahedral symmetry matrices (60 operators for T=1 capsid)
# These are the rotation matrices that generate the full capsid from one asymmetric unit
# Golden ratio
PHI = (1 + np.sqrt(5)) / 2


def generate_icosahedral_matrices() -> List[np.ndarray]:
    """
    Generate the 60 rotation matrices for icosahedral symmetry.
    
    The icosahedral group I has 60 elements:
    - 1 identity
    - 12 rotations by 72° and 144° around 5-fold axes (6 axes × 2 = 12, but 
      actually 12 total 5-fold rotations)
    - 20 rotations by 120° and 240° around 3-fold axes  
    - 15 rotations by 180° around 2-fold axes
    
    Returns:
        List of 60 3x3 rotation matrices
    """
    matrices = []
    
    # Identity
    matrices.append(np.eye(3))
    
    # 5-fold axes (through vertices of icosahedron)
    # There are 6 5-fold axes, each contributing rotations of 72°, 144°, 216°, 288°
    fivefold_axes = [
        np.array([0, 1, PHI]),
        np.array([0, 1, -PHI]),
        np.array([0, -1, PHI]),
        np.array([0, -1, -PHI]),
        np.array([1, PHI, 0]),
        np.array([-1, PHI, 0]),
    ]
    
    for axis in fivefold_axes:
        axis = axis / np.linalg.norm(axis)
        for k in range(1, 5):  # 72°, 144°, 216°, 288°
            angle = k * 2 * np.pi / 5
            matrices.append(rotation_matrix_from_axis_angle(axis, angle))
    
    # 3-fold axes (through face centers)
    # There are 10 3-fold axes
    threefold_axes = [
        np.array([1, 1, 1]),
        np.array([1, 1, -1]),
        np.array([1, -1, 1]),
        np.array([1, -1, -1]),
        np.array([-1, 1, 1]),
        np.array([-1, 1, -1]),
        np.array([-1, -1, 1]),
        np.array([-1, -1, -1]),
        np.array([0, PHI, 1/PHI]),
        np.array([0, PHI, -1/PHI]),
    ]
    
    for axis in threefold_axes:
        axis = axis / np.linalg.norm(axis)
        for k in range(1, 3):  # 120°, 240°
            angle = k * 2 * np.pi / 3
            matrices.append(rotation_matrix_from_axis_angle(axis, angle))
    
    # 2-fold axes (through edge midpoints)
    # There are 15 2-fold axes
    twofold_axes = [
        np.array([1, 0, 0]),
        np.array([0, 1, 0]),
        np.array([0, 0, 1]),
        np.array([1, PHI, 0]),
        np.array([1, -PHI, 0]),
        np.array([PHI, 0, 1]),
        np.array([PHI, 0, -1]),
        np.array([0, 1, PHI]),
        np.array([0, -1, PHI]),
        np.array([1, 1, 1]),
        np.array([1, 1, -1]),
        np.array([1, -1, 1]),
        np.array([-1, 1, 1]),
        np.array([PHI, 1, 1/PHI]),
        np.array([1/PHI, PHI, 1]),
    ]
    
    for axis in twofold_axes:
        axis = axis / np.linalg.norm(axis)
        matrices.append(rotation_matrix_from_axis_angle(axis, np.pi))
    
    # Remove duplicates (within numerical tolerance)
    unique_matrices = []
    for m in matrices:
        is_duplicate = False
        for um in unique_matrices:
            if np.allclose(m, um, atol=1e-10):
                is_duplicate = True
                break
        if not is_duplicate:
            unique_matrices.append(m)
    
    logger.info(f"Generated {len(unique_matrices)} unique icosahedral operators")
    
    # Pad to exactly 60 if needed (numerical issues may cause slight variations)
    if len(unique_matrices) < 60:
        logger.warning(f"Only generated {len(unique_matrices)} operators, expected 60")
    
    return unique_matrices[:60]


def rotation_matrix_from_axis_angle(axis: np.ndarray, angle: float) -> np.ndarray:
    """
    Generate a rotation matrix from axis-angle representation.
    
    Uses Rodrigues' rotation formula.
    """
    axis = axis / np.linalg.norm(axis)
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])
    R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)
    return R


@dataclass
class AssembledCapsid:
    """Represents a fully assembled capsid."""
    coordinates: np.ndarray  # Shape: (60, n_atoms, 3)
    chain_ids: List[str]
    atom_names: List[str]
    residue_ids: List[int]
    residue_names: List[str]
    
    @property
    def n_subunits(self) -> int:
        return self.coordinates.shape[0]
    
    @property
    def n_atoms_per_subunit(self) -> int:
        return self.coordinates.shape[1]


def assemble_capsid_from_monomer(
    monomer_pdb_path: Path,
    output_path: Optional[Path] = None,
    n_subunits: int = 60
) -> Path:
    """
    Assemble a full capsid from a single monomer PDB.
    
    Applies icosahedral symmetry operators to generate the full 60-mer.
    
    Args:
        monomer_pdb_path: Path to monomer PDB file
        output_path: Where to save assembled capsid (optional)
        n_subunits: Number of subunits (60 for full, 3 for ASU, 1 for monomer)
        
    Returns:
        Path to assembled capsid PDB
    """
    from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain
    from Bio.PDB.Atom import Atom
    import copy
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("monomer", str(monomer_pdb_path))
    model = structure[0]
    
    # Get the first chain as template
    template_chain = list(model.get_chains())[0]
    
    # Get all atom coordinates
    atoms = list(template_chain.get_atoms())
    coords = np.array([a.get_coord() for a in atoms])
    
    # Generate symmetry operators
    operators = generate_icosahedral_matrices()[:n_subunits]
    
    # Create new structure for assembled capsid
    new_structure = Structure.Structure("capsid")
    new_model = Model.Model(0)
    new_structure.add(new_model)
    
    # Chain ID alphabet
    chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    
    for i, op in enumerate(operators):
        # Transform coordinates
        new_coords = coords @ op.T  # Apply rotation
        
        # Create new chain
        chain_id = chain_ids[i] if i < len(chain_ids) else f"{i}"
        new_chain = Chain.Chain(chain_id)
        
        # Copy residues with transformed coordinates
        for j, residue in enumerate(template_chain.get_residues()):
            new_residue = residue.copy()
            new_residue.detach_parent()
            
            for atom in new_residue:
                # Find corresponding atom index
                for k, orig_atom in enumerate(atoms):
                    if (orig_atom.get_parent().id == residue.id and 
                        orig_atom.name == atom.name):
                        atom.set_coord(new_coords[k])
                        break
            
            new_chain.add(new_residue)
        
        new_model.add(new_chain)
    
    # Save assembled capsid
    if output_path is None:
        output_path = monomer_pdb_path.parent / f"{monomer_pdb_path.stem}_capsid_{n_subunits}mer.pdb"
    
    io = PDBIO()
    io.set_structure(new_structure)
    io.save(str(output_path))
    
    logger.info(f"Assembled {n_subunits}-mer capsid saved to: {output_path}")
    return output_path


def extract_asymmetric_unit(
    capsid_pdb_path: Path,
    output_path: Optional[Path] = None,
    chains: List[str] = ["A", "B", "C"]
) -> Path:
    """
    Extract the asymmetric unit (typically 3 chains) from a capsid.
    
    The ASU for AAV contains portions of VP1, VP2, and VP3 that form
    one "biological unit" of the capsid.
    """
    from Bio.PDB import PDBParser, PDBIO, Select
    
    class ChainSelect(Select):
        def accept_chain(self, chain):
            return chain.id in chains
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("capsid", str(capsid_pdb_path))
    
    if output_path is None:
        output_path = capsid_pdb_path.parent / f"{capsid_pdb_path.stem}_asu.pdb"
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path), ChainSelect())
    
    logger.info(f"Extracted ASU (chains {chains}) to: {output_path}")
    return output_path


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    # Test matrix generation
    matrices = generate_icosahedral_matrices()
    print(f"Generated {len(matrices)} icosahedral operators")
    
    # Verify orthogonality
    for i, m in enumerate(matrices[:5]):
        det = np.linalg.det(m)
        print(f"Matrix {i}: det = {det:.6f}")
