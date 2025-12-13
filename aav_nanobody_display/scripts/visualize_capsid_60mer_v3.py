#!/usr/bin/env python3
"""
Visualize full 60-mer AAV9 capsid - CIF NATIVE VERSION

Bypasses PDB atom limit (99,999) by using mmCIF format directly.

Usage:
    python scripts/visualize_capsid_60mer_v2.py
"""

import sys
import gzip
import urllib.request
from pathlib import Path

# Fix module search path to include project root
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# Check for BioPython (still needed for fallback)
try:
    import Bio
    from Bio.PDB import PDBParser
except ImportError:
    print("Error: BioPython is required. Please install it: pip install biopython")
    sys.exit(1)

import numpy as np

def download_biological_assembly(pdb_id: str = "3UX1", output_dir: Path = None) -> Path:
    """Download the pre-built biological assembly from RCSB."""
    if output_dir is None:
        output_dir = Path.home() / ".cache" / "aav_nanobody_display" / "structures"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Try CIF format first (Preferred for large structures)
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
        
        gz_path.unlink()
        print(f"Downloaded: {cif_path}")
        return cif_path
        
    except Exception as e:
        print(f"CIF download failed: {e}")
        # Fallback logic omitted for brevity as CIF is standard for 3UX1
        return None

def parse_biomt_from_pdb(pdb_path: Path) -> list:
    """Parse BIOMT transformation matrices from PDB header (Fallback method)."""
    matrices = []
    current_matrix = {}
    
    print(f"Parsing BIOMT records from {pdb_path}...")
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith("REMARK 350") and "BIOMT" in line:
                    parts = line.split()
                    try:
                        # Locate BIOMT keyword index
                        biomt_idx = -1
                        for i, part in enumerate(parts):
                            if "BIOMT" in part:
                                biomt_idx = i
                                break
                        
                        if biomt_idx == -1 or len(parts) < biomt_idx + 5:
                            continue
                            
                        matrix_num = int(parts[biomt_idx + 1])
                        values = [float(x) for x in parts[biomt_idx + 2 : biomt_idx + 6]]
                        
                        if matrix_num not in current_matrix:
                            current_matrix[matrix_num] = {'rot': [], 'trans': []}
                        
                        current_matrix[matrix_num]['rot'].append(values[:3])
                        current_matrix[matrix_num]['trans'].append(values[3])
                        
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
    """Apply BIOMT matrices to create full assembly (Returns PDB format string)."""
    # NOTE: This fallback will still fail if >99,999 atoms are generated.
    # It is kept only for smaller test cases or partial assemblies.
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("monomer", str(monomer_pdb))
    chain = list(structure[0].get_chains())[0]
    
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
    pdb_lines = []
    atom_num = 1
    chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    
    for i, (rot, trans) in enumerate(matrices):
        chain_id = chain_ids[i] if i < len(chain_ids) else "0"
        transformed = (coords @ rot.T) + trans
        
        for j, atom in enumerate(atoms):
            # Reset atom_num if it exceeds limit to prevent crash, 
            # though some viewers might glitch.
            if atom_num > 99999: atom_num = 1
            
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

def create_visualization_html(data_string: str, output_path: Path, format_type: str = 'cif'):
    """Create interactive HTML visualization supporting both CIF and PDB."""
    
    # Escape for JS inclusion
    data_escaped = data_string.replace('\\', '\\\\').replace('`', '\\`').replace('\n', '\\n')
    
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>AAV9 Full Capsid - {format_type.upper()}</title>
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
        .legend-item {{ display: flex; align-items: center; margin: 8px 0; }}
        .color-box {{ width: 20px; height: 20px; margin-right: 10px; border-radius: 3px; }}
        .controls {{ margin-top: 15px; padding-top: 10px; border-top: 1px solid #ddd; }}
        button {{ margin: 3px; padding: 5px 10px; cursor: pointer; }}
    </style>
</head>
<body>
    <div id="viewer"></div>
    <div class="legend">
        <h3>AAV9 Capsid (60-mer)</h3>
        <p style="font-size: 12px; color: #666;">Format: {format_type.upper()}</p>
        <div class="legend-item">
            <div class="color-box" style="background: #E8E8E8;"></div>
            <span>VP3 core</span>
        </div>
        <div class="legend-item">
            <div class="color-box" style="background: #DC143C;"></div>
            <span>VR-VIII (585-596)</span>
        </div>
        <div class="legend-item">
            <div class="color-box" style="background: #4682B4;"></div>
            <span>VR-IV (451-474)</span>
        </div>
        <div class="controls">
            <button onclick="resetView()">Reset View</button>
            <button onclick="viewer.spin(true);">Spin</button>
            <button onclick="viewer.spin(false);">Stop</button>
        </div>
    </div>
    <script>
        let viewer = $3Dmol.createViewer('viewer', {{ backgroundColor: 'white' }});
        let data = `{data_escaped}`;
        
        // Load data explicitly as CIF or PDB
        viewer.addModel(data, '{format_type}');
        
        function highlightLoops() {{
            // VR-VIII in red
            viewer.setStyle(
                {{resi: [585,586,587,588,589,590,591,592,593,594,595,596]}}, 
                {{cartoon: {{color: '#DC143C'}}}}
            );
            // VR-IV in blue
            let vr_iv = [];
            for (let i = 451; i <= 474; i++) vr_iv.push(i);
            viewer.setStyle({{resi: vr_iv}}, {{cartoon: {{color: '#4682B4'}}}});
        }}
        
        function resetView() {{
            viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray'}}}});
            highlightLoops();
            viewer.zoomTo();
            viewer.render();
        }}
        
        resetView();
    </script>
</body>
</html>"""
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Saved: {output_path}")

def main():
    print("=" * 60)
    print("AAV9 Full Capsid Visualization (CIF Native)")
    print("=" * 60)
    
    output_dir = Path(".")
    
    # Method 1: Download pre-built biological assembly (CIF)
    print("\nMethod: Downloading biological assembly from RCSB...")
    
    try:
        assembly_path = download_biological_assembly("3UX1")
        
        if assembly_path and assembly_path.suffix == ".cif":
            print("Reading CIF data directly (Bypassing PDB atom limit)...")
            with open(assembly_path, 'r') as f:
                cif_string = f.read()
            
            html_output = output_dir / "aav9_60mer_complete.html"
            create_visualization_html(cif_string, html_output, format_type='cif')
            
            print("\nSUCCESS: Native CIF visualization created.")
            print(f"Open with: explorer.exe {html_output}")
            return

    except Exception as e:
        print(f"Primary method failed: {e}")
        import traceback
        traceback.print_exc()

    print("\nAttempting fallback to BIOMT method...")
    try:
        # Fallback requires 'structures' module
        from structures import fetch_capsid
        monomer_path = fetch_capsid("AAV9")
        matrices = parse_biomt_from_pdb(monomer_path)
        
        if not matrices:
            print("ERROR: No BIOMT records found.")
            sys.exit(1)
            
        pdb_string = create_60mer_from_biomt(monomer_path, matrices)
        
        html_output = output_dir / "aav9_60mer_biomt.html"
        create_visualization_html(pdb_string, html_output, format_type='pdb')
        print(f"Saved Fallback: {html_output}")
        
    except Exception as e2:
        print(f"Fallback failed: {e2}")
        print("Ensure you are running from the project root and 'structures' module exists.")
        sys.exit(1)

if __name__ == "__main__":
    main()