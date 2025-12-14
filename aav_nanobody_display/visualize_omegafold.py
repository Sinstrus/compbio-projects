#!/usr/bin/env python3
"""
Visualize OmegaFold predicted VHH structure
Creates an interactive 3D HTML visualization using py3Dmol
"""
import py3Dmol
from pathlib import Path

def visualize_structure(pdb_path, output_html=None):
    """
    Create interactive 3D visualization of predicted structure.

    Args:
        pdb_path: Path to PDB file
        output_html: Optional output HTML path (defaults to same name as PDB)
    """
    pdb_file = Path(pdb_path)

    if not pdb_file.exists():
        print(f"Error: PDB file not found: {pdb_path}")
        return

    print(f"Visualizing: {pdb_file}")

    # Read PDB file
    with open(pdb_file) as f:
        pdb_data = f.read()

    # Create 3D viewer
    view = py3Dmol.view(width=1000, height=800)
    view.addModel(pdb_data, 'pdb')

    # Style: cartoon representation with spectrum coloring (N-term=blue, C-term=red)
    view.setStyle({'cartoon': {'color': 'spectrum'}})

    # Add surface representation (semi-transparent)
    view.addStyle({'cartoon': {'color': 'spectrum', 'opacity': 0.8}})

    # Center and zoom
    view.zoomTo()

    # Determine output path
    if output_html is None:
        output_html = pdb_file.with_suffix('.html')
    else:
        output_html = Path(output_html)

    # Save HTML
    html_content = view._make_html()
    with open(output_html, 'w') as f:
        f.write(html_content)

    print(f"✓ Interactive visualization saved: {output_html}")
    print(f"\nTo view:")
    print(f"  Open {output_html} in your web browser")
    print(f"\nVisualization features:")
    print(f"  • Drag to rotate")
    print(f"  • Scroll to zoom")
    print(f"  • Right-click drag to pan")
    print(f"  • Spectrum coloring: N-terminus (blue) → C-terminus (red)")

    return output_html

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Visualize OmegaFold predicted structure')
    parser.add_argument('pdb_file', nargs='?',
                       default='predictions_omegafold/test_vhh.pdb',
                       help='Path to PDB file (default: predictions_omegafold/test_vhh.pdb)')
    parser.add_argument('-o', '--output',
                       help='Output HTML file path (default: same as PDB with .html extension)')

    args = parser.parse_args()

    print("=" * 80)
    print("VHH Structure Visualization")
    print("=" * 80)
    print()

    visualize_structure(args.pdb_file, args.output)

    print()
    print("=" * 80)
