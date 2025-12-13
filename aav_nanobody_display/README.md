# AAV Nanobody Display Visualization Pipeline

A Python toolkit for visualizing AAV capsids with nanobody insertions. Supports structure prediction, fusion modeling, and interactive 3D visualization.

## Features

- **Capsid Structure Handling**
  - Fetch AAV structures from RCSB PDB (AAV9, AAV2, AAV8, AAV5)
  - Parse VP monomers and identify VR loops
  - Assemble full 60-mer capsids from monomers

- **Nanobody Structure Prediction**
  - ESMFold: Fast local prediction (~30 sec on GPU)
  - ColabFold: Higher quality with MSA (~5-10 min)
  
- **VP-Nanobody Fusion Modeling**
  - VR loop insertions (VR-I through VR-IX)
  - VP2 N-terminal fusions
  - Automatic memory estimation for your GPU

- **Visualization**
  - Interactive 3D views with py3Dmol
  - Standalone HTML files (no server needed)
  - Side-by-side comparison views
  - Publication-ready color schemes

## Installation

### 1. Basic Installation

```bash
# Clone or copy the pipeline
cd aav_nanobody_display

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or: venv\Scripts\activate  # Windows

# Install core dependencies
pip install numpy biopython pyyaml py3dmol
```

### 2. Structure Prediction Setup

You have two options for predicting nanobody structures:

#### Option A: ESMFold (Recommended for Quick Predictions)

ESMFold is fast and runs locally. Good for nanobodies (~120 residues).

```bash
# Install PyTorch first (with CUDA support for GPU)
# See: https://pytorch.org/get-started/locally/
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

# Install ESMFold
pip install fair-esm

# Test installation
python -c "import esm; print('ESMFold ready!')"
```

**Memory Requirements:**
- Nanobody alone (~120 residues): ~4 GB VRAM
- VP-nanobody fusion (~650 residues): ~8 GB VRAM ✓ (fits on RTX 3070 Ti)

#### Option B: LocalColabFold (Higher Quality)

ColabFold uses MSA for better predictions but is slower.

```bash
# Install LocalColabFold
# Full instructions: https://github.com/YoshitakaMo/localcolabfold

# Quick install on Linux:
wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabfold_linux.sh
bash install_colabfold_linux.sh

# Add to PATH
export PATH="/path/to/localcolabfold/colabfold-conda/bin:$PATH"

# Test
colabfold_batch --help
```

**Memory Requirements:**
- Nanobody alone: ~6 GB VRAM
- VP-nanobody fusion (~650 residues): ~12 GB VRAM (may need `--msa-mode single_sequence`)
- Full VP3 monomer (~530 residues) + nanobody: ~10-12 GB VRAM

### 3. Verify Installation

```bash
cd aav_nanobody_display
python -c "
from structures import fetch_capsid
from visualization import visualize_capsid

# Fetch AAV9
path = fetch_capsid('AAV9')
print(f'Downloaded: {path}')

# Create visualization
visualize_capsid(path, output_html='test_view.html')
print('Open test_view.html in your browser!')
"
```

## Quick Start

### Command Line

```bash
# Fetch capsid structure
python -m aav_nanobody_display.cli fetch --serotype AAV9

# Predict nanobody structure
python -m aav_nanobody_display.cli predict \
    --sequence-file my_nanobody.fasta \
    --name anti_target_nb \
    --method esmfold

# Visualize capsid with VR loop highlighting
python -m aav_nanobody_display.cli visualize \
    --pdb 3UX1.pdb \
    --highlight-loops VR-VIII VR-IV \
    --output capsid_view.html

# Run full pipeline
python -m aav_nanobody_display.cli pipeline \
    --serotype AAV9 \
    --nanobody-file my_nanobody.fasta \
    --insertion-site VR-VIII \
    --output-dir ./results
```

### Python API

```python
from aav_nanobody_display import (
    fetch_capsid,
    parse_capsid,
    predict_nanobody_structure,
    visualize_capsid,
    CapsidVisualizer,
)

# 1. Fetch AAV9 capsid
capsid_path = fetch_capsid("AAV9")

# 2. Parse and examine VR loops
parser, monomer = parse_capsid(capsid_path, "AAV9")
print(f"VR loops: {list(monomer.vr_loops.keys())}")

# 3. Predict nanobody structure
nb_sequence = "QVQLVESGGGLVQAGGSLRLSCAASGFPVN..."
nb_path = predict_nanobody_structure(
    nb_sequence,
    name="my_nanobody",
    method="esmfold"  # or "colabfold" for better quality
)

# 4. Visualize
vr_loops = {
    "VR-VIII": (585, 596),
    "VR-IV": (451, 474),
}
visualize_capsid(
    capsid_path,
    output_html="aav9_capsid.html",
    vr_loops=vr_loops
)
```

## Project Structure

```
aav_nanobody_display/
├── config/
│   ├── serotypes.yaml       # VR loop positions per serotype
│   └── insertion_sites.yaml # Insertion site configurations
├── structures/
│   ├── fetch.py             # Download from RCSB/AlphaFold
│   ├── parse.py             # Parse capsid structures
│   └── assemble.py          # Build 60-mer from monomer
├── modeling/
│   ├── nanobody.py          # ESMFold/ColabFold prediction
│   ├── fusion.py            # VP-nanobody fusion modeling
│   └── dock.py              # Quick placement (no modeling)
├── visualization/
│   └── render.py            # py3Dmol visualization
├── cli.py                   # Command-line interface
└── requirements.txt
```

## Working with Claude Code

This pipeline is designed to be developed iteratively with Claude Code. Here's how to work with it:

### Typical Workflow

1. **Open terminal** in the pipeline directory
2. **Ask Claude Code** to modify/extend functionality
3. **Test changes** by running scripts or CLI commands
4. **View results** by opening HTML files in browser

### Example Claude Code Requests

```
"Add support for AAV-DJ serotype to the config"

"Modify the visualization to show surface electrostatics"

"Create a function to calculate steric clashes between nanobody and neighboring subunits"

"Add a batch processing mode for multiple nanobody sequences"
```

### Files Claude Code Can Edit

- `config/*.yaml` - Add new serotypes or modify VR loop positions
- `structures/*.py` - Modify structure handling
- `modeling/*.py` - Add new prediction methods or modify fusion logic
- `visualization/*.py` - Customize rendering styles
- `cli.py` - Add new commands

## Supported Serotypes

| Serotype | PDB ID | Common Insertion Sites |
|----------|--------|----------------------|
| AAV9     | 3UX1   | VR-VIII, VR-IV, VP2 N-term |
| AAV2     | 1LP3   | 588 (VR-VIII), VR-IV |
| AAV8     | 2QA0   | VR-VIII, VR-IV |
| AAV5     | 3NTT   | VR-VIII |

## Memory/Hardware Guidelines

For your RTX 3070 Ti (8 GB VRAM):

| Task | Method | Feasible? | Notes |
|------|--------|-----------|-------|
| Nanobody prediction | ESMFold | ✓ Yes | ~30 sec |
| Nanobody prediction | ColabFold | ✓ Yes | ~5 min |
| VP3 + Nb fusion (~650 aa) | ESMFold | ✓ Yes | ~2 min |
| VP3 + Nb fusion | ColabFold | ⚠️ Tight | Use `--msa-mode single_sequence` |
| Full 60-mer modeling | Any | ✗ No | Model monomer, apply symmetry |

## Troubleshooting

### ESMFold CUDA Out of Memory
```python
# Reduce chunk size
predictor = ESMFoldPredictor(device="cuda")
predictor.predict(nanobody, chunk_size=64)  # Default is 128
```

### ColabFold Too Slow / OOM
```bash
# Use single sequence mode (faster, less memory)
colabfold_batch input.fasta output/ --msa-mode single_sequence
```

### py3Dmol Not Rendering
- Make sure you're opening the HTML file directly in a browser
- Some browsers block local file JavaScript - try Firefox or disable security

## License

MIT License - feel free to modify and distribute.

## References

- [AAV capsid structures](https://www.rcsb.org/)
- [ESMFold](https://github.com/facebookresearch/esm)
- [LocalColabFold](https://github.com/YoshitakaMo/localcolabfold)
- [py3Dmol](https://3dmol.csb.pitt.edu/)
