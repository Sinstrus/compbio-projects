#!/usr/bin/env python3
"""
Visualize full 60-mer AAV9 capsid with VR loops highlighted.

Usage:
    python scripts/visualize_capsid_60mer.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from structures import fetch_capsid
from structures.assemble import generate_icosahedral_matrices
import numpy as np

def create_60mer_visualization():
    """Create full capsid visualization using py3Dmol."""
    import py3Dmol
    from Bio.PDB import PDBParser
    
    print("1. Fetching AAV9 monomer...")
    monomer_path = fetch_capsid("AAV9")
    
    print("2. Loading structure...")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("aav9", str(monomer_path))
    model = structure[0]
    chain = list(model.get_chains())[0]
    
    # Get all atom coordinates
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
    
    print("3. Generating 60 icosahedral copies...")
    matrices = generate_icosahedral_matrices()
    
    # Build PDB string for all 60 subunits
    pdb_lines = []
    atom_num = 1
    chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz01234567"
    
    for i, rot_matrix in enumerate(matrices[:60]):
        chain_id = chain_ids[i] if i < len(chain_ids) else "X"
        rotated_coords = coords @ rot_matrix.T
        
        for j, atom in enumerate(atoms):
            x, y, z = rotated_coords[j]
            pdb_lines.append(
                f"ATOM  {atom_num:5d}  {atom['name']:<3s} {atom['resname']:3s} "
                f"{chain_id}{atom['resid']:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  "
                f"1.00  0.00           {atom['element']:>2s}"
            )
            atom_num += 1
        pdb_lines.append(f"TER   {atom_num:5d}      {atoms[-1]['resname']:3s} {chain_id}")
        atom_num += 1
    
    pdb_lines.append("END")
    pdb_string = "\n".join(pdb_lines)
    
    print(f"4. Created {len(matrices)} subunits, {atom_num} atoms total")
    
    # Save full capsid PDB (optional, useful for other tools)
    capsid_pdb = Path("aav9_60mer.pdb")
    with open(capsid_pdb, 'w') as f:
        f.write(pdb_string)
    print(f"   Saved: {capsid_pdb}")
    
    print("5. Creating visualization...")
    view = py3Dmol.view(width=1200, height=900)
    view.addModel(pdb_string, "pdb")
    
    # Base style - light gray cartoon
    view.setStyle({}, {"cartoon": {"color": "lightgray"}})
    
    # Highlight VR-VIII (585-596) on ALL chains - red
    view.setStyle(
        {"resi": "585-596"}, 
        {"cartoon": {"color": "firebrick"}}
    )
    
    # Highlight VR-IV (451-474) on ALL chains - blue  
    view.setStyle(
        {"resi": "451-474"},
        {"cartoon": {"color": "steelblue"}}
    )
    
    view.setBackgroundColor("white")
    view.zoomTo()
    
    # Generate HTML
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>AAV9 Full Capsid (60-mer)</title>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{ margin: 0; font-family: Arial, sans-serif; }}
        #viewer {{ width: 100vw; height: 100vh; }}
        .legend {{
            position: absolute;
            top: 10px;
            left: 10px;
            background: rgba(255,255,255,0.9);
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.2);
            z-index: 100;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            margin: 5px 0;
        }}
        .color-box {{
            width: 20px;
            height: 20px;
            margin-right: 10px;
            border-radius: 3px;
        }}
    </style>
</head>
<body>
    <div class="legend">
        <h3 style="margin-top:0">AAV9 Capsid (60-mer)</h3>
        <div class="legend-item">
            <div class="color-box" style="background: lightgray;"></div>
            <span>VP3 core</span>
        </div>
        <div class="legend-item">
            <div class="color-box" style="background: firebrick;"></div>
            <span>VR-VIII (585-596) - 60 sites</span>
        </div>
        <div class="legend-item">
            <div class="color-box" style="background: steelblue;"></div>
            <span>VR-IV (451-474) - 60 sites</span>
        </div>
        <p style="font-size: 12px; color: #666; margin-bottom: 0;">
            Drag to rotate | Scroll to zoom
        </p>
    </div>
    <div id="viewer"></div>
    <script>
        let viewer = $3Dmol.createViewer('viewer', {{backgroundColor: 'white'}});
        let pdbData = `{pdb_string}`;
        viewer.addModel(pdbData, 'pdb');
        viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray'}}}});
        viewer.setStyle({{resi: ['585','586','587','588','589','590','591','592','593','594','595','596']}}, 
                       {{cartoon: {{color: 'firebrick'}}}});
        viewer.setStyle({{resi: ['451','452','453','454','455','456','457','458','459','460','461','462','463','464','465','466','467','468','469','470','471','472','473','474']}},
                       {{cartoon: {{color: 'steelblue'}}}});
        viewer.zoomTo();
        viewer.render();
    </script>
</body>
</html>"""
    
    output_html = Path("aav9_60mer_capsid.html")
    with open(output_html, 'w') as f:
        f.write(html_content)
    
    print(f"\n{'='*50}")
    print(f"Done! Output files:")
    print(f"  - {capsid_pdb} (PDB structure)")
    print(f"  - {output_html} (interactive visualization)")
    print(f"\nOpen with: explorer.exe {output_html}")
    print(f"{'='*50}")
    
    return output_html


if __name__ == "__main__":
    create_60mer_visualization()