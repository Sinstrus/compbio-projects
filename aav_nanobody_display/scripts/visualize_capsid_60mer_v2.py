#!/usr/bin/env python3
"""
Visualize full 60-mer AAV9 capsid - CORRECTED VERSION

Uses either:
1. Pre-built biological assembly from RCSB (recommended)
2. BIOMT records from PDB header

Usage:
    python scripts/visualize_capsid_60mer_v2.py
"""

import sys
import gzip
import urllib.request
from pathlib import Path
from io import StringIO

# Check for BioPython
try:
    import Bio
    from Bio.PDB import PDBParser, MMCIFParser, PDBIO
except ImportError:
    print("Error: BioPython is required. Please install it: pip install biopython")
    sys.exit(1)

import numpy as np


def download_biological_assembly(pdb_id: str = "3UX1", output_dir: Path = None) -> Path:
    """Download the pre-built biological assembly from RCSB."""
    if output_dir is None:
        output_dir = Path.home() / ".cache" / "aav_nanobody_display" / "structures"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Try CIF format first (more reliable for large assemblies)
    cif_path = output_dir / f"{pdb_id}_assembly1.cif"
    
    if cif_path.exists():
        print(f"Using cached assembly: {cif_path}")
        return cif_path
    
    url = f"https://files.rcsb.org/download/{pdb_id}-assembly1.cif.gz"
    gz_path = output_dir / f"{pdb_id}-assembly1.cif.gz"
    
    print(f"Downloading biological assembly from RCSB...")
    print(f"URL: {url}")
    
    try:
        urllib.request.urlretrieve(url, gz_path)
        
        # Decompress
        with gzip.open(gz_path, 'rb') as f_in:
            with open(cif_path, 'wb') as f_out:
                f_out.write(f_in.read())
        
        gz_path.unlink()  # Remove compressed file
        print(f"Downloaded: {cif_path} ({cif_path.stat().st_size / 1e6:.1f} MB)")
        return cif_path
        
    except Exception as e:
        print(f"CIF download failed: {e}")
        print("Trying PDB format...")
        
        # Fallback to PDB format
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb1.gz"
        pdb_path = output_dir / f"{pdb_id}_assembly1.pdb"
        gz_path = output_dir / f"{pdb_id}.pdb1.gz"
        
        urllib.request.urlretrieve(pdb_url, gz_path)
        
        with gzip.open(gz_path, 'rb') as f_in:
            with open(pdb_path, 'wb') as f_out:
                f_out.write(f_in.read())
        
        gz_path.unlink()
        print(f"Downloaded: {pdb_path}")
        return pdb_path


def parse_biomt_from_pdb(pdb_path: Path) -> list:
    """
    Parse BIOMT (biological assembly) transformation matrices from PDB header.
    """
    matrices = []
    current_matrix = {}
    
    print(f"Parsing BIOMT records from {pdb_path}...")
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith("REMARK 350") and "BIOMT" in line:
                    parts = line.split()
                    # Robust parsing to handle spacing variations
                    # Format: REMARK 350   BIOMT1   1  1.000000  0.000000...
                    try:
                        # Find the index where 'BIOMT' appears
                        biomt_idx = -1
                        for i, part in enumerate(parts):
                            if "BIOMT" in part:
                                biomt_idx = i
                                break
                        
                        if biomt_idx == -1 or len(parts) < biomt_idx + 5:
                            continue
                            
                        biomt_type = parts[biomt_idx]  # BIOMT1, BIOMT2, BIOMT3
                        matrix_num = int(parts[biomt_idx + 1])
                        values = [float(x) for x in parts[biomt_idx + 2 : biomt_idx + 6]]
                        
                        if matrix_num not in current_matrix:
                            current_matrix[matrix_num] = {'rot': [], 'trans': []}
                        
                        current_matrix[matrix_num]['rot'].append(values[:3])
                        current_matrix[matrix_num]['trans'].append(values[3])
                        
                        # When we have all 3 rows, create the matrix
                        if len(current_matrix[matrix_num]['rot']) == 3:
                            rot = np.array(current_matrix[matrix_num]['rot'])
                            trans = np.array(current_matrix[matrix_num]['trans'])
                            matrices.append((rot, trans))
                            
                    except (ValueError, IndexError):
                        continue
                        
    except Exception as e:
        print(f"Warning: Failed to parse BIOMT records: {e}")
        return []
    
    print(f"Parsed {len(matrices)} BIOMT transformation matrices")
    return matrices


def create_60mer_from_biomt(monomer_pdb: Path, matrices: list) -> str:
    """Apply BIOMT matrices to create full assembly."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("monomer", str(monomer_pdb))
    chain = list(structure[0].get_chains())[0]
    
    # Collect atom data
    atoms = []
    for residue in chain:
        for atom in residue:
            atoms.append({
                'coord': atom.get_coord(),
                'name': atom.get_name(),
                'resname': residue.get_resname(),
                'resid': residue.id[1],
                'element': atom.element if atom.element else atom.get_name()[0]
            })
    
    coords = np.array([a['coord'] for a in atoms])
    
    # Generate PDB with all subunits
    pdb_lines = []
    atom_num = 1
    # Extended chain ID list for 60 subunits
    chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    
    for i, (rot, trans) in enumerate(matrices):
        chain_id = chain_ids[i] if i < len(chain_ids) else "0"
        transformed = (coords @ rot.T) + trans
        
        for j, atom in enumerate(atoms):
            x, y, z = transformed[j]
            pdb_lines.append(
                f"ATOM  {atom_num:5d}  {atom['name']:<3s} {atom['resname']:3s} "
                f"{chain_id}{atom['resid']:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  "
                f"1.00  0.00           {atom['element']:>2s}"
            )
            atom_num += 1
        pdb_lines.append(f"TER")
        atom_num += 1
    
    pdb_lines.append("END")
    return "\n".join(pdb_lines)


def load_cif_assembly(cif_path: Path) -> str:
    """
    Load biological assembly from mmCIF file and convert to PDB format.
    CRITICAL FIX: Remaps extended chain IDs (e.g., 'A-2') to single chars for PDB compatibility.
    """
    print("Parsing mmCIF file (this may take a moment for 60 chains)...")
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("capsid", str(cif_path))
    
    chains = list(structure[0].get_chains())
    print(f"Loaded {len(chains)} chains from CIF structure")
    
    # FIX: Remap chain IDs to single characters to satisfy PDB format limits
    # AAV has 60 subunits; standard PDB allows ~62 chars (A-Z, a-z, 0-9)
    valid_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    
    count = 0
    for chain in chains:
        if count < len(valid_ids):
            chain.id = valid_ids[count]
        else:
            chain.id = "0" # Overflow fallback
        count += 1
        
    print("Remapped chain IDs to single characters for PDB export.")

    # Write to PDB format string
    io = PDBIO()
    io.set_structure(structure)
    
    output = StringIO()
    io.save(output)
    return output.getvalue()


def create_visualization_html(pdb_string: str, output_path: Path):
    """Create interactive HTML visualization."""
    
    # Escape backticks and backslashes for JavaScript
    pdb_escaped = pdb_string.replace('\\', '\\\\').replace('`', '\\`').replace('\n', '\\n')
    
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>AAV9 Full Capsid (60-mer) - Biological Assembly</title>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{ margin: 0; font-family: Arial, sans-serif; }}
        #viewer {{ width: 100vw; height: 100vh; position: absolute; top: 0; left: 0; }}
        .legend {{
            position: absolute;
            top: 10px;
            left: 10px;
            background: rgba(255,255,255,0.95);
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.2);
            z-index: 100;
            max-width: 280px;
        }}
        .legend h3 {{ margin-top: 0; }}
        .legend-item {{
            display: flex;
            align-items: center;
            margin: 8px 0;
        }}
        .color-box {{
            width: 20px;
            height: 20px;
            margin-right: 10px;
            border-radius: 3px;
            flex-shrink: 0;
        }}
        .controls {{
            margin-top: 15px;
            padding-top: 10px;
            border-top: 1px solid #ddd;
        }}
        button {{
            margin: 3px;
            padding: 5px 10px;
            cursor: pointer;
        }}
    </style>
</head>
<body>
    <div id="viewer"></div>
    <div class="legend">
        <h3>AAV9 Capsid (60-mer)</h3>
        <p style="font-size: 12px; color: #666; margin: 5px 0;">
            Biological assembly from RCSB PDB
        </p>
        <div class="legend-item">
            <div class="color-box" style="background: #E8E8E8;"></div>
            <span>VP3 core (60 subunits)</span>
        </div>
        <div class="legend-item">
            <div class="color-box" style="background: #DC143C;"></div>
            <span>VR-VIII (585-596) - primary Nb site</span>
        </div>
        <div class="legend-item">
            <div class="color-box" style="background: #4682B4;"></div>
            <span>VR-IV (451-474) - alternative site</span>
        </div>
        <div class="controls">
            <strong>View:</strong><br>
            <button onclick="viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray'}}}}); highlightLoops(); viewer.render();">Cartoon</button>
            <button onclick="viewer.setStyle({{}}, {{sphere: {{radius: 1.5, color: 'lightgray'}}}}); viewer.render();">Surface</button>
            <button onclick="viewer.spin(true);">Spin</button>
            <button onclick="viewer.spin(false);">Stop</button>
        </div>
        <p style="font-size: 11px; color: #888; margin-top: 10px; margin-bottom: 0;">
            Drag to rotate | Scroll to zoom | Right-drag to pan
        </p>
    </div>
    <script>
        let viewer = $3Dmol.createViewer('viewer', {{
            backgroundColor: 'white'
        }});
        
        // Load the PDB data
        viewer.addModel(`{pdb_escaped}`, 'pdb');
        
        function highlightLoops() {{
            // VR-VIII in red
            viewer.setStyle(
                {{resi: [585,586,587,588,589,590,591,592,593,594,595,596]}}, 
                {{cartoon: {{color: '#DC143C'}}}}
            );
            // VR-IV in blue
            for (let i = 451; i <= 474; i++) {{
                viewer.setStyle({{resi: i}}, {{cartoon: {{color: '#4682B4'}}}});
            }}
        }}
        
        viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray'}}}});
        highlightLoops();
        viewer.zoomTo();
        viewer.render();
    </script>
</body>
</html>"""
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Saved: {output_path}")


def main():
    print("=" * 60)
    print("AAV9 Full Capsid Visualization (Corrected)")
    print("=" * 60)
    
    output_dir = Path(".")
    
    # Method 1: Download pre-built biological assembly (RECOMMENDED)
    print("\nMethod: Downloading biological assembly from RCSB...")
    
    pdb_string = None
    
    try:
        assembly_path = download_biological_assembly("3UX1")
        
        if assembly_path.suffix == ".cif":
            # This is the path that requires the fix
            pdb_string = load_cif_assembly(assembly_path)
        else:
            with open(assembly_path, 'r') as f:
                pdb_string = f.read()
        
        # Also save as PDB for other uses
        pdb_output = output_dir / "aav9_60mer_biological.pdb"
        with open(pdb_output, 'w') as f:
            f.write(pdb_string)
        print(f"Saved remapped PDB: {pdb_output}")
        
    except Exception as e:
        print(f"Download/Parsing failed: {e}")
        import traceback
        traceback.print_exc()
        print("\nFalling back to BIOMT method...")
        
        # Method 2: Use BIOMT records from original PDB
        # We need to fetch the monomer PDB first
        from structures import fetch_capsid
        try:
            monomer_path = fetch_capsid("AAV9")
            matrices = parse_biomt_from_pdb(monomer_path)
            
            if not matrices:
                print("ERROR: No BIOMT records found and download failed.")
                sys.exit(1)
            
            pdb_string = create_60mer_from_biomt(monomer_path, matrices)
            
            pdb_output = output_dir / "aav9_60mer_biomt.pdb"
            with open(pdb_output, 'w') as f:
                f.write(pdb_string)
            print(f"Saved PDB: {pdb_output}")
            
        except Exception as e2:
            print(f"Fallback failed: {e2}")
            sys.exit(1)
    
    if pdb_string:
        # Create visualization
        html_output = output_dir / "aav9_60mer_complete.html"
        create_visualization_html(pdb_string, html_output)
        
        print("\n" + "=" * 60)
        print("Done!")
        print(f"Open with: explorer.exe {html_output}")
        print("=" * 60)


if __name__ == "__main__":
    main()