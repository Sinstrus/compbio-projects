"""
Quick nanobody docking/placement module.

For rough visualization without full structural modeling.
Simply positions a nanobody structure at an insertion site using
geometric transformations.
"""

import logging
from pathlib import Path
from typing import Optional, Tuple, List
import numpy as np

logger = logging.getLogger(__name__)

try:
    from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain
    from Bio.PDB.vectors import Vector, rotmat
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False


def get_terminal_vectors(pdb_path: Path, chain_id: str = "A") -> dict:
    """
    Get N-terminal and C-terminal positions and directions from a structure.
    
    Returns dict with:
        - n_term_coord: N-terminus CA position
        - c_term_coord: C-terminus CA position  
        - n_term_direction: Vector pointing into protein from N-term
        - c_term_direction: Vector pointing into protein from C-term
    """
    if not HAS_BIOPYTHON:
        raise ImportError("BioPython required")
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", str(pdb_path))
    chain = structure[0][chain_id]
    
    # Get all CA atoms
    ca_atoms = []
    for residue in chain:
        if "CA" in residue:
            ca_atoms.append(residue["CA"])
    
    if len(ca_atoms) < 3:
        raise ValueError("Structure has fewer than 3 CA atoms")
    
    # N-terminus
    n_term_coord = ca_atoms[0].get_coord()
    n_term_next = ca_atoms[1].get_coord()
    n_term_direction = n_term_next - n_term_coord
    n_term_direction = n_term_direction / np.linalg.norm(n_term_direction)
    
    # C-terminus
    c_term_coord = ca_atoms[-1].get_coord()
    c_term_prev = ca_atoms[-2].get_coord()
    c_term_direction = c_term_prev - c_term_coord
    c_term_direction = c_term_direction / np.linalg.norm(c_term_direction)
    
    return {
        "n_term_coord": n_term_coord,
        "c_term_coord": c_term_coord,
        "n_term_direction": n_term_direction,
        "c_term_direction": c_term_direction,
    }


def calculate_transformation(
    source_point: np.ndarray,
    source_direction: np.ndarray,
    target_point: np.ndarray,
    target_direction: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate rotation matrix and translation to align source to target.
    
    Args:
        source_point: Point to move from
        source_direction: Direction at source (will be aligned to -target_direction)
        target_point: Point to move to
        target_direction: Direction at target
        
    Returns:
        (rotation_matrix, translation_vector)
    """
    # We want source_direction to point opposite to target_direction
    # (so the nanobody points away from the insertion site)
    target_dir_reversed = -target_direction
    
    # Calculate rotation matrix using Rodrigues formula
    v = np.cross(source_direction, target_dir_reversed)
    s = np.linalg.norm(v)
    c = np.dot(source_direction, target_dir_reversed)
    
    if s < 1e-10:  # Vectors are parallel
        if c > 0:  # Same direction
            R = np.eye(3)
        else:  # Opposite direction
            # Rotate 180 degrees around any perpendicular axis
            perp = np.array([1, 0, 0]) if abs(source_direction[0]) < 0.9 else np.array([0, 1, 0])
            perp = perp - np.dot(perp, source_direction) * source_direction
            perp = perp / np.linalg.norm(perp)
            R = 2 * np.outer(perp, perp) - np.eye(3)
    else:
        vx = np.array([
            [0, -v[2], v[1]],
            [v[2], 0, -v[0]],
            [-v[1], v[0], 0]
        ])
        R = np.eye(3) + vx + vx @ vx * (1 - c) / (s * s)
    
    # Translation: move rotated source_point to target_point
    rotated_source = R @ source_point
    translation = target_point - rotated_source
    
    return R, translation


def place_nanobody_at_site(
    nanobody_pdb: Path,
    capsid_pdb: Path,
    insertion_residue: int,
    chain_id: str = "A",
    output_path: Optional[Path] = None,
    offset_distance: float = 5.0,  # Angstroms away from surface
    attach_terminus: str = "C"  # Which terminus of nanobody attaches to capsid
) -> Path:
    """
    Place a nanobody structure at a capsid insertion site.
    
    This is a ROUGH placement for visualization purposes only.
    For accurate structures, use the fusion modeling module.
    
    Args:
        nanobody_pdb: Path to nanobody structure
        capsid_pdb: Path to capsid monomer structure
        insertion_residue: Residue number at insertion site
        chain_id: Chain ID in capsid
        output_path: Where to save combined structure
        offset_distance: Distance to offset nanobody from insertion point
        attach_terminus: Which nanobody terminus attaches ("N" or "C")
        
    Returns:
        Path to combined PDB with placed nanobody
    """
    if not HAS_BIOPYTHON:
        raise ImportError("BioPython required")
    
    parser = PDBParser(QUIET=True)
    
    # Load structures
    capsid = parser.get_structure("capsid", str(capsid_pdb))
    nanobody = parser.get_structure("nanobody", str(nanobody_pdb))
    
    capsid_chain = capsid[0][chain_id]
    nb_chain = list(nanobody[0].get_chains())[0]
    
    # Find insertion site coordinates
    target_residue = None
    for residue in capsid_chain:
        if residue.id[1] == insertion_residue:
            target_residue = residue
            break
    
    if target_residue is None:
        raise ValueError(f"Residue {insertion_residue} not found in chain {chain_id}")
    
    if "CA" not in target_residue:
        raise ValueError(f"No CA atom in residue {insertion_residue}")
    
    target_coord = target_residue["CA"].get_coord()
    
    # Get direction for placement (pointing outward from capsid center)
    # Approximate: just use vector from origin to insertion point
    capsid_center = np.zeros(3)  # Assume capsid is centered
    outward_direction = target_coord - capsid_center
    outward_direction = outward_direction / np.linalg.norm(outward_direction)
    
    # Get nanobody terminal info
    nb_vectors = get_terminal_vectors(nanobody_pdb, nb_chain.id)
    
    if attach_terminus == "C":
        source_coord = nb_vectors["c_term_coord"]
        source_direction = nb_vectors["c_term_direction"]
    else:
        source_coord = nb_vectors["n_term_coord"]
        source_direction = nb_vectors["n_term_direction"]
    
    # Calculate transformation
    R, T = calculate_transformation(
        source_coord, source_direction,
        target_coord + outward_direction * offset_distance,
        -outward_direction
    )
    
    # Apply transformation to nanobody
    for atom in nanobody.get_atoms():
        new_coord = R @ atom.get_coord() + T
        atom.set_coord(new_coord)
    
    # Combine structures
    # Add nanobody as new chain to capsid
    new_chain_id = chr(ord(chain_id) + 1)  # Next letter
    nb_chain.id = new_chain_id
    capsid[0].add(nb_chain)
    
    # Save
    if output_path is None:
        output_path = capsid_pdb.parent / f"{capsid_pdb.stem}_with_nanobody.pdb"
    
    io = PDBIO()
    io.set_structure(capsid)
    io.save(str(output_path))
    
    logger.info(f"Saved capsid with placed nanobody to: {output_path}")
    return output_path


def place_multiple_nanobodies(
    nanobody_pdb: Path,
    capsid_pdb: Path,
    insertion_residues: List[int],
    chain_ids: List[str],
    output_path: Optional[Path] = None
) -> Path:
    """
    Place nanobodies at multiple sites (e.g., all 60 subunits).
    
    Warning: This creates a very large structure file.
    """
    if not HAS_BIOPYTHON:
        raise ImportError("BioPython required")
    
    parser = PDBParser(QUIET=True)
    capsid = parser.get_structure("capsid", str(capsid_pdb))
    nb_template = parser.get_structure("nanobody", str(nanobody_pdb))
    
    nb_vectors = get_terminal_vectors(nanobody_pdb)
    
    for i, (res_num, chain_id) in enumerate(zip(insertion_residues, chain_ids)):
        try:
            capsid_chain = capsid[0][chain_id]
            
            # Find insertion site
            target_residue = None
            for residue in capsid_chain:
                if residue.id[1] == res_num:
                    target_residue = residue
                    break
            
            if target_residue is None or "CA" not in target_residue:
                logger.warning(f"Skipping {chain_id}:{res_num} - residue not found")
                continue
            
            target_coord = target_residue["CA"].get_coord()
            outward_direction = target_coord / np.linalg.norm(target_coord)
            
            # Clone nanobody
            import copy
            nb_copy = copy.deepcopy(nb_template)
            
            # Transform
            R, T = calculate_transformation(
                nb_vectors["c_term_coord"],
                nb_vectors["c_term_direction"],
                target_coord + outward_direction * 5.0,
                -outward_direction
            )
            
            for atom in nb_copy.get_atoms():
                new_coord = R @ atom.get_coord() + T
                atom.set_coord(new_coord)
            
            # Add as new chain
            new_chain_id = f"{i:02d}"[-2:]  # Two-character chain ID
            nb_chain = list(nb_copy[0].get_chains())[0]
            nb_chain.id = new_chain_id
            capsid[0].add(nb_chain)
            
        except Exception as e:
            logger.warning(f"Failed to place nanobody at {chain_id}:{res_num}: {e}")
    
    if output_path is None:
        output_path = capsid_pdb.parent / f"{capsid_pdb.stem}_with_nanobodies.pdb"
    
    io = PDBIO()
    io.set_structure(capsid)
    io.save(str(output_path))
    
    logger.info(f"Saved capsid with {len(insertion_residues)} nanobodies to: {output_path}")
    return output_path


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    print("Quick docking module loaded.")
    print("Use place_nanobody_at_site() for rough visualization.")
    print("For accurate modeling, use the fusion module instead.")
