#!/usr/bin/env python3
"""
Quick Start: Visualize AAV9 capsid in under 30 seconds

This minimal script:
1. Downloads AAV9 structure
2. Creates an interactive HTML visualization
3. Opens it in your default browser

Usage:
    python scripts/quick_start.py
"""

from pathlib import Path
import sys
import webbrowser

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

def main():
    print("AAV9 Quick Visualization")
    print("=" * 40)
    
    # Fetch structure
    print("\n1. Downloading AAV9 capsid structure...")
    from structures import fetch_capsid
    capsid_path = fetch_capsid("AAV9")
    print(f"   Downloaded: {capsid_path}")
    
    # Create visualization
    print("\n2. Creating visualization...")
    from visualization import visualize_capsid
    
    output_html = Path("aav9_quick_view.html")
    
    # Highlight VR-VIII (common insertion site)
    vr_loops = {
        "VR-VIII": (585, 596),  # Most common for nanobody display
        "VR-IV": (451, 474),    # Alternative site
    }
    
    visualize_capsid(
        capsid_path,
        output_html=output_html,
        vr_loops=vr_loops,
        show_surface=False,
        show_labels=True
    )
    print(f"   Saved: {output_html.absolute()}")
    
    # Open in browser
    print("\n3. Opening in browser...")
    webbrowser.open(f"file://{output_html.absolute()}")
    
    print("\n" + "=" * 40)
    print("Done! You can rotate the structure with your mouse.")
    print("VR-VIII (blue) is the most common nanobody insertion site.")


if __name__ == "__main__":
    main()
