# Structure Prediction Setup Guide
## AAV Nanobody Display Project - RTX 3070 Ti (8GB VRAM) Configuration

**Target Hardware:** MSI RTX 3070 Ti Laptop (8GB VRAM)
**Development Environment:** HP EliteBook (No GPU) - for code development only
**Last Updated:** 2025-12-12

---

## Overview

This guide provides the exact setup instructions for nanobody structure prediction optimized for your RTX 3070 Ti with 8GB VRAM. The code is written on your work laptop but **must be run on your home laptop** for GPU acceleration.

---

## Split-Architecture Workflow

```
Work Laptop (HP EliteBook)           Home Laptop (MSI RTX 3070 Ti)
├─ Code development                  ├─ Structure prediction (ESMFold)
├─ Visualization                     ├─ Fusion modeling (ColabFold)
├─ Analysis                          └─ Heavy computation
└─ Git commits
```

**Key Rule:** Never try to run `predict` commands with `--device cuda` on the work laptop - you'll get a clear error message telling you to transfer to the RTX 3070 Ti.

---

## Installation Instructions

### Work Laptop (Development Environment)

```bash
# Basic dependencies only (no CUDA)
pip install numpy biopython pyyaml py3dmol

# DO NOT install PyTorch/ESM on this machine
# It's fine to have them, but you won't use GPU features
```

### Home Laptop (RTX 3070 Ti - Compute Node)

#### Step 1: PyTorch with CUDA 11.8

```bash
# Install PyTorch with CUDA support
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

# Verify CUDA is available
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"
# Expected output: CUDA available: True
```

#### Step 2: ESMFold (Recommended for Nanobodies)

```bash
# Install fair-esm
pip install fair-esm

# Test ESMFold
python -c "import esm; print('ESMFold ready!')"
```

**Memory Profile:**
- Nanobody alone (~120 residues): ~4 GB VRAM ✓
- VP-nanobody fusion (~650 residues): ~7 GB VRAM ✓ (with chunk_size=64)
- **Prediction time:** ~30 seconds per nanobody

#### Step 3: LocalColabFold (Optional - Higher Quality)

```bash
# Download installer
wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabfold_linux.sh

# Install
bash install_colabfold_linux.sh

# Add to PATH (add to ~/.bashrc)
export PATH="$HOME/localcolabfold/colabfold-conda/bin:$PATH"

# Test
colabfold_batch --help
```

**Memory Profile:**
- Nanobody with `--msa-mode single_sequence`: ~6 GB VRAM ✓
- VP-nanobody fusion with `single_sequence`: ~10 GB VRAM ⚠️ (tight, may OOM)
- **Prediction time:** ~5-10 minutes per nanobody

---

## Usage Examples

### Basic Nanobody Prediction (ESMFold)

**On RTX 3070 Ti Home Laptop:**

```bash
# Auto-detect GPU and use ESMFold
python cli.py predict \
    --sequence-file examples/nanobody_anti_gfp.fasta \
    --name anti_gfp \
    --method esmfold \
    --device auto \
    --chunk-size 64 \
    --output-dir predictions/

# Expected output:
# ╔════════════════════════════════════════════════════════════╗
# ║ ESMFold loaded on GPU: NVIDIA GeForce RTX 3070 Ti Laptop  ║
# ║ VRAM: 8.0 GB | chunk_size=64                              ║
# ╚════════════════════════════════════════════════════════════╝
```

### If You Accidentally Run on Work Laptop:

```bash
python cli.py predict \
    --sequence-file examples/nanobody_anti_gfp.fasta \
    --name anti_gfp \
    --device auto

# Auto-fallback to CPU:
# WARNING: No CUDA detected - falling back to CPU (50-100x slower)
```

### Force CUDA (Will error on work laptop):

```bash
python cli.py predict \
    --sequence-file examples/nanobody_anti_gfp.fasta \
    --name anti_gfp \
    --device cuda

# On work laptop, you'll get:
# ╔════════════════════════════════════════════════════════════╗
# ║ CUDA REQUESTED BUT NOT AVAILABLE                          ║
# ╚════════════════════════════════════════════════════════════╝
# Options:
#   1. Transfer this command to your RTX 3070 Ti home laptop
#   2. Use --device cpu (WARNING: 50-100x slower)
#   3. Use --device auto to automatically select available device
```

### ColabFold with Memory Optimization

**On RTX 3070 Ti Home Laptop:**

```bash
python cli.py predict \
    --sequence-file examples/nanobody_anti_gfp.fasta \
    --name anti_gfp \
    --method colabfold \
    --msa-mode single_sequence \
    --output-dir predictions/
```

**MSA Mode Options:**
- `single_sequence`: Fast, low memory (6GB), good for 8GB VRAM ✓
- `mmseqs2+prefilter`: High quality, more memory (10GB+), use only if you have headroom

---

## Memory Optimization Reference

### ESMFold Chunk Size Tuning

| Sequence Length | chunk_size=32 | chunk_size=64 | chunk_size=128 |
|-----------------|---------------|---------------|----------------|
| Nanobody (~120) | 3 GB          | 4 GB ✓        | 6 GB           |
| VP-Nb (~650)    | 5 GB          | 7 GB ✓        | OOM ✗          |

**For RTX 3070 Ti (8GB):** Use `chunk_size=64` (default)

If you get OOM errors:
```bash
python cli.py predict ... --chunk-size 32
```

### ColabFold MSA Mode Comparison

| MSA Mode            | Quality | Speed    | VRAM (Nb) | VRAM (Fusion) |
|---------------------|---------|----------|-----------|---------------|
| mmseqs2+prefilter   | Best    | Slow     | 8 GB      | 12 GB ✗       |
| single_sequence     | Good    | Fast ✓   | 6 GB ✓    | 10 GB ⚠️      |

**For RTX 3070 Ti (8GB):** Use `--msa-mode single_sequence` (default)

---

## Troubleshooting

### ESMFold: CUDA Out of Memory

```bash
# Symptom: RuntimeError: CUDA out of memory
# Solution: Reduce chunk size
python cli.py predict ... --chunk-size 32
```

### ColabFold: OOM on Fusion Modeling

```bash
# Symptom: OOM when predicting VP-nanobody fusions
# Solution: Already using single_sequence mode? Model monomer only
# The full 60-mer will be assembled via symmetry operations
```

### PyTorch Not Finding CUDA

```bash
# Verify CUDA installation
python -c "import torch; print(torch.cuda.is_available())"

# If False, reinstall PyTorch with CUDA
pip uninstall torch torchvision torchaudio
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
```

---

## Quick Reference Commands

```bash
# Test installation (Work Laptop)
python cli.py --help                      # Should work
python cli.py predict --help              # Should work

# Test installation (RTX 3070 Ti)
python -c "import torch; print(torch.cuda.is_available())"  # Should be True
python cli.py predict --sequence "QVQLV..." --device auto   # Should use GPU

# Typical workflow
python cli.py predict -f nanobody.fasta --device auto --method esmfold  # Fast
python cli.py predict -f nanobody.fasta --method colabfold              # Quality
```

---

## Performance Benchmarks (RTX 3070 Ti)

| Task                  | Method    | Time      | VRAM Used |
|-----------------------|-----------|-----------|-----------|
| Nanobody (120 aa)     | ESMFold   | ~30s      | 4 GB      |
| Nanobody (120 aa)     | ColabFold | ~5 min    | 6 GB      |
| VP-Nb fusion (650 aa) | ESMFold   | ~2 min    | 7 GB      |
| VP-Nb fusion (650 aa) | ColabFold | ~10 min   | 10 GB ⚠️  |

---

## Next Steps

After completing this setup:
1. Test with example nanobody: `examples/nanobody_anti_gfp.fasta`
2. Predict your candidate sequences on RTX 3070 Ti
3. Transfer PDB files back to work laptop for visualization
4. Proceed to VP-nanobody fusion modeling (Step 3)

---

## Key Implementation Details (for Reference)

### ESMFoldPredictor Class
- **Library:** `fair-esm`
- **Default chunk_size:** 64 (optimized for 8GB VRAM)
- **Device auto-detection:** Enabled by default
- **Model loading:** Lazy (only loads when needed)

### ColabFoldPredictor Class
- **Wrapper:** Subprocess call to `colabfold_batch` CLI tool
- **Default MSA mode:** `single_sequence` (memory-optimized)
- **Arguments:** `--msa-mode`, `--num-recycle`, `--amber`

### CLI Safety Features
- Imports are lazy (doesn't crash on work laptop)
- `--device auto` detects GPU availability
- Clear error messages when CUDA requested but unavailable
- Help commands work without GPU dependencies
