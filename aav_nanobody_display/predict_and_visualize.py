#!/usr/bin/env python3
"""
VHH Structure Prediction & Visualization
Uses ESMFold API (their GPU) + local visualization (your browser)
"""
import requests
import sys
from pathlib import Path

def predict_structure_online(sequence, output_pdb="predictions/vhh_predicted.pdb"):
    """
    Predict structure using ESMFold API (runs on their servers).
    Fast, reliable, no local setup needed!
    """
    print("=" * 80)
    print("VHH Structure Prediction via ESMFold API")
    print("=" * 80)

    print(f"\nSequence length: {len(sequence)} amino acids")
    print(f"Sequence: {sequence[:50]}...")

    print("\nSending request to ESMFold server...")
    print("(This uses META's GPU servers - completely free!)")

    # ESMAtlas API endpoint
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"

    try:
        response = requests.post(url, data=sequence, timeout=60)
        response.raise_for_status()

        pdb_data = response.text

        # Save PDB
        output_path = Path(output_pdb)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, 'w') as f:
            f.write(pdb_data)

        print(f"\n✓ Structure predicted successfully!")
        print(f"✓ Saved to: {output_path.absolute()}")

        return str(output_path.absolute()), pdb_data

    except requests.exceptions.Timeout:
        print("\n✗ Request timed out. The server might be busy.")
        print("Try again in a moment, or use the web interface:")
        print("https://esmatlas.com/resources?action=fold")
        return None, None
    except Exception as e:
        print(f"\n✗ Error: {e}")
        print("\nAlternative: Use the web interface:")
        print("https://esmatlas.com/resources?action=fold")
        return None, None

def visualize_structure(pdb_path, pdb_data=None):
    """
    Create interactive 3D visualization (works in any browser).
    Similar to the AAV capsid visualization!
    """
    print("\n" + "=" * 80)
    print("Creating Interactive 3D Visualization")
    print("=" * 80)

    if pdb_data is None:
        with open(pdb_path) as f:
            pdb_data = f.read()

    # Create beautiful interactive HTML visualization
    html_path = Path(pdb_path).with_suffix('.html')

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>VHH Structure - 3D Viewer</title>
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', 'Roboto', sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
            min-height: 100vh;
        }}

        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 20px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
            overflow: hidden;
        }}

        header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
        }}

        h1 {{
            font-size: 32px;
            margin-bottom: 10px;
        }}

        .subtitle {{
            opacity: 0.9;
            font-size: 16px;
        }}

        .info-panel {{
            padding: 20px 30px;
            background: #f8f9fa;
            border-bottom: 1px solid #e9ecef;
        }}

        .info-row {{
            display: flex;
            gap: 40px;
            margin: 10px 0;
            flex-wrap: wrap;
        }}

        .info-item {{
            display: flex;
            align-items: center;
            gap: 8px;
        }}

        .info-label {{
            font-weight: 600;
            color: #495057;
        }}

        .info-value {{
            color: #6c757d;
        }}

        .viewer-container {{
            padding: 30px;
        }}

        #viewer {{
            width: 100%;
            height: 600px;
            border-radius: 15px;
            background: #000;
            box-shadow: inset 0 2px 10px rgba(0,0,0,0.5);
        }}

        .controls {{
            padding: 20px 30px 30px 30px;
            display: flex;
            gap: 15px;
            flex-wrap: wrap;
            justify-content: center;
        }}

        button {{
            padding: 12px 24px;
            border: none;
            border-radius: 8px;
            font-size: 14px;
            font-weight: 600;
            cursor: pointer;
            transition: all 0.2s;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}

        button:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(0,0,0,0.15);
        }}

        .btn-primary {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
        }}

        .btn-secondary {{
            background: white;
            color: #667eea;
            border: 2px solid #667eea;
        }}

        .sequence-box {{
            margin-top: 20px;
            padding: 15px;
            background: #f1f3f5;
            border-radius: 8px;
            font-family: 'Courier New', monospace;
            font-size: 12px;
            word-break: break-all;
            line-height: 1.6;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>🧬 VHH Nanobody Structure</h1>
            <div class="subtitle">3D Structure Predicted by ESMFold</div>
        </header>

        <div class="info-panel">
            <div class="info-row">
                <div class="info-item">
                    <span class="info-label">Sequence Length:</span>
                    <span class="info-value">{len(sequence)} amino acids</span>
                </div>
                <div class="info-item">
                    <span class="info-label">Method:</span>
                    <span class="info-value">ESMFold (META AI)</span>
                </div>
                <div class="info-item">
                    <span class="info-label">Model:</span>
                    <span class="info-value">RTX 3070 Ti (Visualization)</span>
                </div>
            </div>
            <div class="sequence-box">
                <strong>Sequence:</strong><br>
                {sequence}
            </div>
        </div>

        <div class="viewer-container">
            <div id="viewer"></div>
        </div>

        <div class="controls">
            <button class="btn-primary" onclick="resetView()">🔄 Reset View</button>
            <button class="btn-secondary" onclick="toggleStyle('cartoon')">📊 Cartoon</button>
            <button class="btn-secondary" onclick="toggleStyle('stick')">🔗 Stick</button>
            <button class="btn-secondary" onclick="toggleStyle('sphere')">⚫ Sphere</button>
            <button class="btn-primary" onclick="toggleSpin()">🌀 Toggle Spin</button>
        </div>
    </div>

    <script>
        let viewer = $3Dmol.createViewer("viewer", {{backgroundColor: 'black'}});
        let spinning = false;
        let currentStyle = 'cartoon';

        // Load structure
        viewer.addModel(`{pdb_data}`, "pdb");
        viewer.setStyle({{}}, {{cartoon: {{color: 'spectrum'}}}});
        viewer.setBackgroundColor('#0a0a0a');
        viewer.zoomTo();
        viewer.render();

        function resetView() {{
            viewer.zoomTo();
            viewer.render();
        }}

        function toggleStyle(style) {{
            currentStyle = style;
            viewer.setStyle({{}}, {{}});

            if (style === 'cartoon') {{
                viewer.setStyle({{}}, {{cartoon: {{color: 'spectrum'}}}});
            }} else if (style === 'stick') {{
                viewer.setStyle({{}}, {{stick: {{color: 'spectrum'}}}});
            }} else if (style === 'sphere') {{
                viewer.setStyle({{}}, {{sphere: {{color: 'spectrum'}}}});
            }}

            viewer.render();
        }}

        function toggleSpin() {{
            spinning = !spinning;
            if (spinning) {{
                viewer.spin(true);
            }} else {{
                viewer.spin(false);
            }}
        }}

        // Mouse controls info
        console.log('Controls:');
        console.log('- Left mouse: Rotate');
        console.log('- Right mouse or Ctrl+Left: Zoom');
        console.log('- Middle mouse or Shift+Left: Pan');
    </script>
</body>
</html>
"""

    with open(html_path, 'w') as f:
        f.write(html_content)

    print(f"\n✓ Visualization created: {html_path.absolute()}")
    print(f"\n{'='*80}")
    print("SUCCESS! 🎉")
    print('='*80)
    print(f"\nOpen this file in your browser:")
    print(f"  {html_path.absolute()}")
    print(f"\nControls:")
    print(f"  • Left mouse: Rotate")
    print(f"  • Scroll: Zoom")
    print(f"  • Use buttons to change visualization style")

    return str(html_path.absolute())

if __name__ == "__main__":
    # Your VHH sequence
    sequence = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSTADMGWFRQAPGKGRELVAAVSGSGFSTYSDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAKATISLYYAMDVWGQGTTVTVSS"

    print("\n" + "=" * 80)
    print("  VHH STRUCTURE PREDICTION & VISUALIZATION")
    print("=" * 80)

    # Step 1: Predict structure (using ESMFold API)
    pdb_path, pdb_data = predict_structure_online(sequence)

    if pdb_path:
        # Step 2: Create visualization
        html_path = visualize_structure(pdb_path, pdb_data)
    else:
        sys.exit(1)
