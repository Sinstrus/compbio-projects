#!/usr/bin/env python3
"""
Local VHH structure prediction using ESM contact prediction + simple folding
This uses models that actually fit in your GPU memory!
"""
import torch
import numpy as np
from pathlib import Path
from Bio.PDB import PDBIO, Structure, Model, Chain, Residue, Atom
from Bio.PDB.vectors import Vector
import esm

def predict_contacts_gpu(sequence):
    """
    Use ESM-1b to predict residue contacts on GPU.
    This model is much smaller (~650MB) and includes contact prediction.
    """
    print("=" * 80)
    print("Local Structure Prediction using ESM Contact Prediction")
    print("=" * 80)

    print(f"\nSequence length: {len(sequence)} amino acids")
    print(f"Sequence: {sequence[:50]}...")

    # Check GPU
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"\nDevice: {device}")
    if device == "cuda":
        print(f"GPU: {torch.cuda.get_device_name(0)}")
        print(f"VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

    print("\n" + "=" * 80)
    print("Loading ESM-1b model (contact prediction)...")
    print("=" * 80)

    # Load ESM-1b model (smaller, has contact head)
    model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
    batch_converter = alphabet.get_batch_converter()
    model = model.to(device)
    model.eval()

    print(f"✓ Model loaded on {device}")
    if device == "cuda":
        print(f"GPU Memory: {torch.cuda.memory_allocated() / 1e9:.2f} GB")

    print("\n" + "=" * 80)
    print("Predicting contacts...")
    print("=" * 80)

    # Prepare data
    data = [("protein", sequence)]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_tokens = batch_tokens.to(device)

    # Predict
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=True)

    # Get contact map
    contacts = results["contacts"][0].cpu().numpy()

    print(f"✓ Contact map computed: {contacts.shape}")
    if device == "cuda":
        print(f"GPU Memory: {torch.cuda.memory_allocated() / 1e9:.2f} GB")

    # Get representations
    token_representations = results["representations"][33][0].cpu().numpy()

    return contacts, token_representations

def build_simple_structure(sequence, contacts, output_pdb="predictions/vhh_local.pdb"):
    """
    Build a simple 3D structure from contact predictions.
    This is a simplified approach - not as accurate as AlphaFold but works locally!
    """
    print("\n" + "=" * 80)
    print("Building 3D structure from contacts...")
    print("=" * 80)

    n = len(sequence)

    # Initialize with extended chain (simple backbone)
    coords = np.zeros((n, 3))

    # Place residues along a helix-like path as starting point
    for i in range(n):
        # Simple helical geometry
        angle = i * (2 * np.pi / 3.6)  # ~3.6 residues per turn
        coords[i] = [
            3.8 * np.cos(angle),
            3.8 * np.sin(angle),
            i * 1.5  # 1.5 Å rise per residue
        ]

    # Refine using contact constraints (simple distance geometry)
    # This is a very simplified approach - real folding is much more complex!
    print("Refining structure using contact constraints...")

    # Get top contacts
    threshold = np.percentile(contacts, 95)  # Top 5% contacts
    strong_contacts = contacts > threshold

    # Simple iterative refinement (gradient descent on distances)
    learning_rate = 0.1
    for iteration in range(100):
        forces = np.zeros_like(coords)

        for i in range(n):
            for j in range(i+5, n):  # Only consider i,j pairs separated by 5+
                if strong_contacts[i, j]:
                    # Ideal distance for a contact
                    ideal_dist = 8.0  # Ångströms

                    vec = coords[j] - coords[i]
                    dist = np.linalg.norm(vec)

                    if dist > 0:
                        # Force proportional to distance error
                        force = (dist - ideal_dist) * vec / dist
                        forces[i] += force
                        forces[j] -= force

        # Update positions
        coords -= learning_rate * forces

        if iteration % 20 == 0:
            avg_force = np.mean(np.linalg.norm(forces, axis=1))
            print(f"  Iteration {iteration}: avg force = {avg_force:.3f}")

    print("✓ Structure refined")

    # Create PDB structure
    print("\nCreating PDB file...")
    structure = Structure.Structure("predicted")
    model = Model.Model(0)
    chain = Chain.Chain("A")

    # Standard residue names
    three_letter = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
        'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
        'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
    }

    for i, (aa, coord) in enumerate(zip(sequence, coords)):
        resname = three_letter.get(aa, 'ALA')
        residue = Residue.Residue((' ', i+1, ' '), resname, '')

        # Add CA atom (simplified - only backbone alpha carbon)
        atom = Atom.Atom(
            'CA',
            Vector(coord),
            1.0,  # b-factor
            1.0,  # occupancy
            ' ',
            'CA',
            i+1,
            'C'
        )
        residue.add(atom)
        chain.add(residue)

    model.add(chain)
    structure.add(model)

    # Save PDB
    output_path = Path(output_pdb)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path))

    print(f"✓ Structure saved: {output_path.absolute()}")

    return str(output_path.absolute())

def visualize_structure(pdb_path):
    """Create interactive 3D visualization."""
    import py3Dmol

    print("\n" + "=" * 80)
    print("Creating visualization...")
    print("=" * 80)

    # Load PDB
    with open(pdb_path) as f:
        pdb_data = f.read()

    # Create view
    view = py3Dmol.view(width=1000, height=600)
    view.addModel(pdb_data, 'pdb')

    # Style: cartoon + spectrum coloring (N-term=blue, C-term=red)
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.setBackgroundColor('white')
    view.zoomTo()

    # Save HTML
    html_path = Path(pdb_path).with_suffix('.html')
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>VHH Structure Prediction</title>
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <style>
        body {{ margin: 0; padding: 20px; font-family: Arial, sans-serif; }}
        h2 {{ color: #333; }}
        #info {{ margin-bottom: 20px; padding: 15px; background: #f0f0f0; border-radius: 5px; }}
        #viewer {{ width: 100%; height: 600px; border: 2px solid #333; }}
    </style>
</head>
<body>
    <h2>VHH Structure Prediction (Local GPU)</h2>
    <div id="info">
        <strong>Method:</strong> ESM-1b Contact Prediction + Distance Geometry<br>
        <strong>Model:</strong> RTX 3070 Ti GPU<br>
        <strong>Note:</strong> This is a simplified structure based on predicted contacts.
        For publication-quality structures, use AlphaFold or ESMFold.
    </div>
    <div id="viewer"></div>
    <script>
        let viewer = $3Dmol.createViewer("viewer");
        viewer.addModel(`{pdb_data}`, "pdb");
        viewer.setStyle({{}}, {{cartoon: {{color: 'spectrum'}}}});
        viewer.setBackgroundColor('white');
        viewer.zoomTo();
        viewer.render();
    </script>
</body>
</html>
"""

    with open(html_path, 'w') as f:
        f.write(html_content)

    print(f"✓ Visualization saved: {html_path.absolute()}")
    print(f"\nOpen {html_path.absolute()} in your browser!")

    return str(html_path.absolute())

if __name__ == "__main__":
    # Your VHH sequence
    sequence = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSTADMGWFRQAPGKGRELVAAVSGSGFSTYSDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAKATISLYYAMDVWGQGTTVTVSS"

    try:
        # Step 1: Predict contacts on GPU
        contacts, representations = predict_contacts_gpu(sequence)

        # Step 2: Build structure from contacts
        pdb_path = build_simple_structure(sequence, contacts)

        # Step 3: Visualize
        html_path = visualize_structure(pdb_path)

        print("\n" + "=" * 80)
        print("SUCCESS!")
        print("=" * 80)
        print(f"\n✓ Structure: {pdb_path}")
        print(f"✓ Visualization: {html_path}")
        print("\nNOTE: This is a simplified prediction. For higher accuracy, use:")
        print("  - AlphaFold Server: https://alphafoldserver.com/")
        print("  - ESMFold Online: https://esmatlas.com/resources?action=fold")

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
