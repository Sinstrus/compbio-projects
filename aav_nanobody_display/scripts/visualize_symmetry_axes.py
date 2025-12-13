#!/usr/bin/env python3
"""
AAV9 Symmetry Axis Explorer for VHH Binder Design

Generates an interactive HTML viewer to inspect symmetry-related interfaces:
- 5-fold axis (pentamer/pore): 1 reference chain + 4 nearest neighbors
- 3-fold axis (trimer/spike): 1 reference chain + 2 nearest neighbors
- 2-fold axis (dimer): 1 reference chain + 1 nearest neighbor

Uses geometric center-of-mass (COM) calculations to identify neighbors,
avoiding hardcoded chain ID assumptions for biological assemblies.

Usage:
    python scripts/visualize_symmetry_axes.py

Output:
    aav9_symmetry_explorer.html - Interactive viewer with symmetry toggles
"""

import sys
import gzip
import urllib.request
from pathlib import Path
from io import StringIO
from collections import defaultdict

# Check for BioPython
try:
    import Bio
    from Bio.PDB import MMCIFParser, PDBIO
except ImportError:
    print("Error: BioPython is required. Please install it: pip install biopython")
    sys.exit(1)

import numpy as np


def download_biological_assembly(pdb_id: str = "3UX1", output_dir: Path = None) -> Path:
    """Download the pre-built biological assembly from RCSB."""
    if output_dir is None:
        output_dir = Path.home() / ".cache" / "aav_nanobody_display" / "structures"
    output_dir.mkdir(parents=True, exist_ok=True)

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
        print(f"Download failed: {e}")
        raise


def load_cif_assembly(cif_path: Path) -> tuple:
    """
    Load biological assembly from mmCIF file.
    Returns (cif_string, structure_object) for downstream COM analysis.

    Note: We keep the CIF format to avoid PDB's 99,999 atom limit.
    py3dmol can render CIF files directly.
    """
    print("Parsing mmCIF file (this may take a moment for 60 chains)...")
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("capsid", str(cif_path))

    chains = list(structure[0].get_chains())
    print(f"Loaded {len(chains)} chains from CIF structure")

    # Read the raw CIF file
    with open(cif_path, 'r') as f:
        cif_string = f.read()

    print(f"CIF file size: {len(cif_string) / 1e6:.2f} MB")

    return cif_string, structure


def get_interface_residues(chain, interface_type: str) -> list:
    """
    Get residues that define specific symmetry interfaces.

    Args:
        chain: BioPython chain object
        interface_type: '5fold', '3fold', or '2fold'

    Returns:
        List of residues at the interface
    """
    # Expanded interface definitions based on AAV structural biology
    interface_regions = {
        # 5-fold: Pore-forming regions (VR-I + VR-II + N-term shoulder)
        '5fold': list(range(260, 275)) + list(range(326, 338)) + list(range(380, 395)),
        # 3-fold: Spike regions (VR-IV + VR-VIII interdigitation)
        '3fold': list(range(450, 475)) + list(range(580, 596)),
        # 2-fold: Dimer interface (HI loop + VR-IX)
        '2fold': list(range(500, 515)) + list(range(660, 675)),
    }

    if interface_type not in interface_regions:
        raise ValueError(f"Unknown interface type: {interface_type}")

    residue_range = interface_regions[interface_type]
    interface_residues = []

    for residue in chain:
        if residue.id[1] in residue_range:
            interface_residues.append(residue)

    return interface_residues


def find_interface_contacts(structure, ref_chain_id: str, interface_type: str,
                            cutoff: float = 8.0) -> list:
    """
    Find chains that contact the reference chain at a specific interface.

    Uses C-alpha atoms for efficiency. A chain is considered a contact if ANY
    of its C-alpha atoms are within 'cutoff' Angstroms of the interface residues.

    Args:
        structure: BioPython structure object
        ref_chain_id: Reference chain ID (e.g., 'A')
        interface_type: '5fold', '3fold', or '2fold'
        cutoff: Distance cutoff in Angstroms

    Returns:
        List of chain IDs that contact this interface (including reference)
    """
    from Bio.PDB.NeighborSearch import NeighborSearch

    chains = list(structure[0].get_chains())

    # Find reference chain
    ref_chain = None
    for chain in chains:
        if chain.id == ref_chain_id:
            ref_chain = chain
            break

    if ref_chain is None:
        print(f"Warning: Reference chain '{ref_chain_id}' not found. Using first chain.")
        ref_chain = chains[0]
        ref_chain_id = ref_chain.id

    # Get interface residues from reference chain
    interface_residues = get_interface_residues(ref_chain, interface_type)

    if not interface_residues:
        print(f"Warning: No interface residues found for {interface_type}")
        return [ref_chain_id]

    # Collect C-alpha atoms from interface residues
    interface_atoms = []
    for residue in interface_residues:
        if 'CA' in residue:
            interface_atoms.append(residue['CA'])

    if not interface_atoms:
        print(f"Warning: No CA atoms found in interface residues")
        return [ref_chain_id]

    # Collect all C-alpha atoms from all other chains
    all_ca_atoms = []
    for chain in chains:
        if chain.id != ref_chain_id:
            for residue in chain:
                if 'CA' in residue:
                    all_ca_atoms.append(residue['CA'])

    # Use NeighborSearch for efficient spatial queries
    ns = NeighborSearch(all_ca_atoms)

    # Find chains with atoms near the interface
    contacting_chains = set([ref_chain_id])

    for interface_atom in interface_atoms:
        # Search for neighbors within cutoff
        neighbors = ns.search(interface_atom.coord, cutoff, level='A')
        for neighbor_atom in neighbors:
            # Get the chain ID of this atom
            parent_chain = neighbor_atom.get_parent().get_parent()
            contacting_chains.add(parent_chain.id)

    contacting_list = sorted(list(contacting_chains))

    print(f"  {interface_type.upper()} interface ({cutoff}Å cutoff): {len(contacting_list)} chains")
    print(f"    Reference: {ref_chain_id}")
    print(f"    Contacts: {[c for c in contacting_list if c != ref_chain_id][:10]}")

    return contacting_list


def find_symmetry_neighbors(structure, ref_chain_id: str, symmetry_type: str) -> list:
    """
    Find chains that form a specific symmetry cluster with the reference chain.

    Args:
        structure: BioPython structure object
        ref_chain_id: Reference chain ID
        symmetry_type: 'pentamer' (5-fold), 'trimer' (3-fold), or 'dimer' (2-fold)

    Returns:
        List of chain IDs forming the symmetry cluster
    """
    interface_map = {
        'pentamer': '5fold',
        'trimer': '3fold',
        'dimer': '2fold',
    }

    if symmetry_type not in interface_map:
        raise ValueError(f"Unknown symmetry type: {symmetry_type}")

    interface_type = interface_map[symmetry_type]

    # Different cutoffs for different interfaces (pentamers may need more distance)
    cutoffs = {
        'pentamer': 15.0,  # 5-fold pore has wider spacing
        'trimer': 10.0,    # 3-fold spike is tighter
        'dimer': 8.0,      # 2-fold interface (reduced to get true dimer)
    }

    cutoff = cutoffs[symmetry_type]
    contacts = find_interface_contacts(structure, ref_chain_id, interface_type, cutoff=cutoff)

    # Expected cluster sizes (with tolerance)
    expected_ranges = {
        'pentamer': (4, 6),  # 5 ± 1 chains
        'trimer': (3, 4),    # 3-4 chains
        'dimer': (2, 5),     # 2-5 chains (broader context okay)
    }

    expected_min, expected_max = expected_ranges[symmetry_type]
    actual = len(contacts)

    if actual < expected_min:
        print(f"    WARNING: Found {actual} chains, expected {expected_min}-{expected_max}")
        print(f"             May need to increase cutoff for {symmetry_type}")
    elif actual > expected_max:
        print(f"    NOTE: Found {actual} chains (expected {expected_min}-{expected_max})")
        print(f"          Showing broader interface context")

    return contacts


def create_symmetry_explorer_html(cif_string: str, structure, output_path: Path,
                                   ref_chain: str = 'A'):
    """
    Create interactive HTML with symmetry axis exploration controls.

    Pre-calculates neighbor sets for 2-fold, 3-fold, and 5-fold axes using
    interface contact detection.
    """
    print("\nIdentifying symmetry interfaces using contact detection...")
    print("=" * 70)

    # Get all chain IDs from structure
    all_chain_ids = [chain.id for chain in structure[0].get_chains()]

    # Use interface-based contact detection
    dimer_chains = find_symmetry_neighbors(structure, ref_chain, 'dimer')       # 2-fold
    trimer_chains = find_symmetry_neighbors(structure, ref_chain, 'trimer')     # 3-fold
    pentamer_chains = find_symmetry_neighbors(structure, ref_chain, 'pentamer') # 5-fold

    print("=" * 70)

    # Convert to JavaScript arrays
    dimer_js = str(dimer_chains).replace("'", '"')
    trimer_js = str(trimer_chains).replace("'", '"')
    pentamer_js = str(pentamer_chains).replace("'", '"')

    # Escape CIF string for JavaScript
    cif_escaped = cif_string.replace('\\', '\\\\').replace('`', '\\`')

    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>AAV9 Symmetry Axis Explorer - VHH Binder Design</title>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{
            margin: 0;
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: #f5f5f5;
        }}
        #viewer {{
            width: 100vw;
            height: 100vh;
            position: absolute;
            top: 0;
            left: 0;
        }}
        .control-panel {{
            position: absolute;
            top: 15px;
            left: 15px;
            background: rgba(255, 255, 255, 0.97);
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.15);
            z-index: 100;
            max-width: 320px;
            font-size: 14px;
        }}
        .control-panel h2 {{
            margin: 0 0 8px 0;
            font-size: 18px;
            color: #2c3e50;
        }}
        .control-panel h3 {{
            margin: 15px 0 10px 0;
            font-size: 14px;
            color: #34495e;
            border-bottom: 2px solid #3498db;
            padding-bottom: 5px;
        }}
        .subtitle {{
            font-size: 11px;
            color: #7f8c8d;
            margin: 0 0 15px 0;
        }}
        .btn-group {{
            display: flex;
            flex-direction: column;
            gap: 8px;
            margin-bottom: 15px;
        }}
        button {{
            padding: 10px 15px;
            cursor: pointer;
            border: 2px solid #3498db;
            background: white;
            color: #3498db;
            border-radius: 6px;
            font-size: 13px;
            font-weight: 600;
            transition: all 0.2s;
        }}
        button:hover {{
            background: #3498db;
            color: white;
            transform: translateY(-1px);
            box-shadow: 0 2px 8px rgba(52, 152, 219, 0.3);
        }}
        button.active {{
            background: #3498db;
            color: white;
        }}
        .info-box {{
            background: #ecf0f1;
            padding: 10px;
            border-radius: 5px;
            font-size: 12px;
            margin-top: 10px;
            color: #2c3e50;
        }}
        .info-box strong {{
            color: #e74c3c;
        }}
        .chain-selector {{
            margin-top: 10px;
        }}
        .chain-selector label {{
            display: block;
            margin-bottom: 5px;
            font-weight: 600;
            color: #34495e;
        }}
        .chain-selector select {{
            width: 100%;
            padding: 8px;
            border: 2px solid #bdc3c7;
            border-radius: 5px;
            font-size: 13px;
        }}
        .rendering-controls {{
            display: flex;
            gap: 5px;
            margin-top: 10px;
        }}
        .rendering-controls button {{
            flex: 1;
            padding: 6px 8px;
            font-size: 11px;
        }}
        input[type="number"] {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        }}
        input[type="number"]:focus {{
            border-color: #3498db;
            outline: none;
        }}
    </style>
</head>
<body>
    <div id="viewer"></div>
    <div class="control-panel">
        <h2>🔬 AAV9 Symmetry Explorer</h2>
        <p class="subtitle">Steric constraint analysis for VHH insertion</p>

        <h3>Symmetry Axis Views</h3>
        <div class="btn-group">
            <button onclick="show5Fold()" id="btn5fold">
                🔵 5-Fold Axis (Pentamer Pore)
            </button>
            <button onclick="show3Fold()" id="btn3fold">
                🟢 3-Fold Axis (Trimer Spike)
            </button>
            <button onclick="show2Fold()" id="btn2fold">
                🟡 2-Fold Axis (Dimer Interface)
            </button>
            <button onclick="showAll()" id="btnAll">
                ⚪ All Subunits (Reset)
            </button>
        </div>

        <div class="chain-selector">
            <label for="refChain">Reference Chain:</label>
            <select id="refChain" onchange="updateReference()">
                {_generate_chain_options(all_chain_ids, ref_chain)}
            </select>
        </div>

        <div class="info-box" id="infoBox">
            Click a symmetry button to isolate specific interfaces.
            <br><br>
            <strong>Focused chains:</strong> Opaque surface<br>
            <strong>Context chains:</strong> Transparent cartoon
        </div>

        <h3>Residue Highlighter</h3>
        <div style="margin-bottom: 10px;">
            <label for="residueInput" style="font-size: 12px; color: #34495e; display: block; margin-bottom: 5px;">
                Enter Residue Number (200-736):
            </label>
            <input type="number" id="residueInput" min="200" max="736" placeholder="e.g., 585"
                   style="width: 100%; padding: 8px; font-size: 14px; border: 2px solid #bdc3c7; border-radius: 5px;"
                   onchange="highlightResidue(this.value)">
            <div style="font-size: 10px; color: #7f8c8d; margin-top: 5px;">
                Type a number and press Enter
            </div>
        </div>
        <div style="display: flex; gap: 5px; margin-bottom: 8px;">
            <button onclick="jumpToResidue(451)" style="flex: 1; padding: 6px; font-size: 11px;">
                VR-IV (451)
            </button>
            <button onclick="jumpToResidue(585)" style="flex: 1; padding: 6px; font-size: 11px;">
                VR-VIII (585)
            </button>
        </div>
        <div style="margin-bottom: 15px;">
            <button onclick="clearResidueHighlight()" style="width: 100%; padding: 6px; font-size: 11px;">
                Clear Highlight
            </button>
        </div>

        <h3>View Controls</h3>
        <div class="rendering-controls">
            <button onclick="viewer.spin('y', 1); viewer.render();">🔄 Spin</button>
            <button onclick="viewer.spin(false); viewer.render();">⏸️ Stop</button>
            <button onclick="viewer.zoomTo(); viewer.render();">🎯 Reset</button>
        </div>

        <p style="font-size: 10px; color: #95a5a6; margin-top: 15px; margin-bottom: 0;">
            Drag: rotate | Scroll: zoom | Right-drag: pan
        </p>
    </div>

    <script>
        let viewer = $3Dmol.createViewer('viewer', {{
            backgroundColor: '#1a1a1a',
            antialias: true
        }});

        // Load CIF data (avoids PDB's 99,999 atom limit)
        viewer.addModel(`{cif_escaped}`, 'cif');

        // Fix clipping issues for large structures - disable slab/clipping
        viewer.enableFog(false);
        viewer.setSlab(-1000, 1000);  // Very wide clipping planes

        // Pre-calculated neighbor sets (from Python)
        const DIMER_CHAINS = {dimer_js};
        const TRIMER_CHAINS = {trimer_js};
        const PENTAMER_CHAINS = {pentamer_js};

        let currentView = 'all';

        function clearActiveButtons() {{
            document.querySelectorAll('.btn-group button').forEach(btn => {{
                btn.classList.remove('active');
            }});
        }}

        function updateInfoBox(title, focused, total) {{
            const infoBox = document.getElementById('infoBox');
            infoBox.innerHTML = `
                <strong>${{title}}</strong><br>
                Focused: ${{focused.join(', ')}}<br>
                Context: ${{total - focused.length}} transparent chains
            `;
        }}

        function show5Fold() {{
            clearActiveButtons();
            document.getElementById('btn5fold').classList.add('active');
            currentView = '5fold';

            // Reset all to transparent cartoon
            viewer.setStyle({{}}, {{
                cartoon: {{
                    color: '#95a5a6',
                    opacity: 0.3
                }}
            }});

            // Highlight pentamer with opaque cartoon only (no surface/balls)
            PENTAMER_CHAINS.forEach(chain => {{
                viewer.setStyle(
                    {{chain: chain}},
                    {{
                        cartoon: {{
                            color: '#3498db',
                            opacity: 1.0,
                            thickness: 0.8
                        }}
                    }}
                );
            }});

            // Highlight VR-VIII and VR-IV loops in distinct colors
            highlightLoops(PENTAMER_CHAINS);

            viewer.render();
            updateInfoBox('5-Fold Axis (Pore)', PENTAMER_CHAINS, 60);
        }}

        function show3Fold() {{
            clearActiveButtons();
            document.getElementById('btn3fold').classList.add('active');
            currentView = '3fold';

            // Reset all to transparent cartoon
            viewer.setStyle({{}}, {{
                cartoon: {{
                    color: '#95a5a6',
                    opacity: 0.3
                }}
            }});

            // Highlight trimer with opaque cartoon only (no surface/balls)
            TRIMER_CHAINS.forEach(chain => {{
                viewer.setStyle(
                    {{chain: chain}},
                    {{
                        cartoon: {{
                            color: '#2ecc71',
                            opacity: 1.0,
                            thickness: 0.8
                        }}
                    }}
                );
            }});

            highlightLoops(TRIMER_CHAINS);

            viewer.render();
            updateInfoBox('3-Fold Axis (Spike)', TRIMER_CHAINS, 60);
        }}

        function show2Fold() {{
            clearActiveButtons();
            document.getElementById('btn2fold').classList.add('active');
            currentView = '2fold';

            // Reset all to transparent cartoon
            viewer.setStyle({{}}, {{
                cartoon: {{
                    color: '#95a5a6',
                    opacity: 0.3
                }}
            }});

            // Highlight dimer with opaque cartoon only (no surface/balls)
            DIMER_CHAINS.forEach(chain => {{
                viewer.setStyle(
                    {{chain: chain}},
                    {{
                        cartoon: {{
                            color: '#f39c12',
                            opacity: 1.0,
                            thickness: 0.8
                        }}
                    }}
                );
            }});

            highlightLoops(DIMER_CHAINS);

            viewer.render();
            updateInfoBox('2-Fold Axis (Dimer)', DIMER_CHAINS, 60);
        }}

        function showAll() {{
            clearActiveButtons();
            document.getElementById('btnAll').classList.add('active');
            currentView = 'all';

            // Show all chains as cartoon
            viewer.setStyle({{}}, {{
                cartoon: {{
                    color: '#ecf0f1',
                    opacity: 0.8
                }}
            }});

            // Highlight VR-VIII and VR-IV on all chains
            highlightLoops(null);

            viewer.render();

            const infoBox = document.getElementById('infoBox');
            infoBox.innerHTML = 'All 60 subunits displayed as cartoon.<br><br>Select a symmetry axis to focus.';
        }}

        function highlightLoops(chainList) {{
            // VR-VIII (585-596) in red - primary Nb binding site
            const vrVIII = {{resi: [585,586,587,588,589,590,591,592,593,594,595,596]}};

            // VR-IV (451-474) in blue - alternative site
            const vrIV = {{resi: Array.from({{length: 24}}, (_, i) => 451 + i)}};

            if (chainList) {{
                // Highlight only on focused chains
                chainList.forEach(chain => {{
                    viewer.setStyle(
                        {{chain: chain, ...vrVIII}},
                        {{
                            cartoon: {{
                                color: '#e74c3c',
                                opacity: 1.0
                            }}
                        }}
                    );
                    viewer.setStyle(
                        {{chain: chain, ...vrIV}},
                        {{
                            cartoon: {{
                                color: '#3498db',
                                opacity: 1.0
                            }}
                        }}
                    );
                }});
            }} else {{
                // Highlight on all chains
                viewer.setStyle(
                    vrVIII,
                    {{
                        cartoon: {{
                            color: '#e74c3c',
                            opacity: 1.0
                        }}
                    }}
                );
                viewer.setStyle(
                    vrIV,
                    {{
                        cartoon: {{
                            color: '#3498db',
                            opacity: 1.0
                        }}
                    }}
                );
            }}
        }}

        function updateReference() {{
            const newRef = document.getElementById('refChain').value;
            alert(`Reference chain changed to: ${{newRef}}\\n\\nNote: This requires re-running the Python script with ref_chain='${{newRef}}' to recalculate neighbors.`);
        }}

        // Residue highlighting functions
        let currentHighlightedResidue = null;

        function highlightResidue(residueNum) {{
            const resNum = parseInt(residueNum);

            // Validate input
            if (isNaN(resNum) || resNum < 200 || resNum > 736) {{
                alert('Please enter a valid residue number between 200 and 736');
                return;
            }}

            currentHighlightedResidue = resNum;

            // Re-render current view to reset styles
            switch(currentView) {{
                case '5fold':
                    show5Fold();
                    break;
                case '3fold':
                    show3Fold();
                    break;
                case '2fold':
                    show2Fold();
                    break;
                default:
                    showAll();
            }}

            // Highlight ONLY the C-alpha atom of selected residue (one atom per chain)
            viewer.setStyle(
                {{resi: resNum, atom: 'CA'}},
                {{
                    sphere: {{
                        color: '#ff00ff',  // Magenta
                        radius: 3.0
                    }}
                }}
            );

            viewer.render();
        }}

        function clearResidueHighlight() {{
            currentHighlightedResidue = null;
            document.getElementById('residueInput').value = '';

            // Re-render current view
            switch(currentView) {{
                case '5fold':
                    show5Fold();
                    break;
                case '3fold':
                    show3Fold();
                    break;
                case '2fold':
                    show2Fold();
                    break;
                default:
                    showAll();
            }}
        }}

        function jumpToResidue(residueNum) {{
            document.getElementById('residueInput').value = residueNum;
            highlightResidue(residueNum);
        }}

        // Initialize with all subunits view
        showAll();
    </script>
</body>
</html>"""

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"\n✓ Saved interactive viewer: {output_path}")


def _generate_chain_options(chain_ids: list, current_ref: str) -> str:
    """Generate <option> elements for chain selector."""
    options = []
    for chain_id in chain_ids:
        selected = ' selected' if chain_id == current_ref else ''
        options.append(f'<option value="{chain_id}"{selected}>Chain {chain_id}</option>')
    return '\n                '.join(options)


def main():
    print("=" * 70)
    print("AAV9 SYMMETRY AXIS EXPLORER - VHH Binder Design Toolkit")
    print("=" * 70)
    print("\nPurpose: Isolate specific symmetry axes to inspect steric constraints")
    print("         for nanobody insertion at VR-IV/VR-VIII interfaces.\n")

    output_dir = Path(".")

    # Download and load biological assembly
    try:
        assembly_path = download_biological_assembly("3UX1")
        cif_string, structure = load_cif_assembly(assembly_path)
        print(f"✓ Assembly loaded successfully ({len(structure[0])} chains)")

    except Exception as e:
        print(f"ERROR: Failed to load assembly: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    # Generate interactive HTML
    print("\n" + "-" * 70)
    print("Generating interactive symmetry explorer...")
    print("-" * 70)

    html_output = output_dir / "aav9_symmetry_explorer.html"
    create_symmetry_explorer_html(cif_string, structure, html_output, ref_chain='A')

    print("\n" + "=" * 70)
    print("SUCCESS!")
    print("=" * 70)
    print(f"\nOpen in browser: {html_output.absolute()}")
    print("\nInteractive Features:")
    print("  - Symmetry axis selection: 5-fold (blue), 3-fold (green), 2-fold (yellow)")
    print("  - VR-VIII (585-596) and VR-IV (451-474) loops in red/blue")
    print("  - NEW! Residue highlighter (residues 200-736)")
    print("       → Type a residue number and press Enter")
    print("       → Magenta sphere appears on C-alpha of that residue (all 60 chains)")
    print("       → Quick-jump buttons for VR-IV (451) and VR-VIII (585)")
    print("       → Identify most exposed positions for binder insertion")
    print("\nTechnical Details:")
    print("  - Symmetry interfaces identified by residue contact detection")
    print("  - 5-fold: VR-I + VR-II + N-term shoulder (15Å cutoff)")
    print("  - 3-fold: VR-IV + VR-VIII interdigitation (10Å cutoff)")
    print("  - 2-fold: HI loop + VR-IX (8Å cutoff)")
    print("  - Rendering: Cartoon ribbons only (no surfaces)")
    print("  - Context chains: 30% opacity grey")
    print("  - Focused chains: 100% opacity colored")
    print("  - Wide clipping planes to prevent disappearing chains")
    print("=" * 70)


if __name__ == "__main__":
    main()
