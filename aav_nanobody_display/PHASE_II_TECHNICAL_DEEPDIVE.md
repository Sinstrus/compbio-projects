# Phase II: AAV Nanobody Display - Technical Deep Dive

## Executive Summary

This pipeline is a **computational protein engineering toolkit** for designing AAV capsids that display nanobodies (VHH domains) on their surface. The goal is to engineer AAV vectors with enhanced targeting, shielding, or functional capabilities by genetically fusing nanobodies at strategic surface locations.

---

## Table of Contents

1. [System Architecture](#system-architecture)
2. [Data Flow](#data-flow)
3. [Core Components Deep Dive](#core-components-deep-dive)
4. [Design Philosophy & Trade-offs](#design-philosophy--trade-offs)
5. [Hardware Considerations](#hardware-considerations)
6. [Visualization Strategy](#visualization-strategy)
7. [Current Limitations](#current-limitations)

---

## 1. System Architecture

### Component Hierarchy

```
aav_nanobody_display/
├── structures/          # Capsid structure handling
│   ├── fetch.py        # Download PDB files from RCSB
│   ├── parse.py        # Parse VP sequences, extract VR loops
│   └── assemble.py     # Biological assembly construction
│
├── modeling/           # Structure prediction & fusion
│   ├── nanobody.py    # Nanobody structure prediction (ESMFold/ColabFold)
│   ├── fusion.py      # VP-Nanobody fusion constructs
│   └── dock.py        # Positioning nanobodies on capsid surface
│
├── visualization/      # Interactive 3D visualization
│   └── render.py      # py3Dmol-based HTML generation
│
├── config/            # Serotype-specific data
│   ├── serotypes.yaml # VR loop positions, PDB IDs
│   └── insertion_sites.yaml
│
└── cli.py             # Command-line interface
```

### Technology Stack

| Component | Technology | Why This Choice |
|-----------|------------|-----------------|
| **Structure Prediction** | ESMFold, ColabFold | ESMFold: Fast, no MSA (2-5 min on RTX 3070 Ti). ColabFold: Higher quality for multi-domain proteins |
| **Visualization** | py3Dmol (3Dmol.js) | Client-side JavaScript rendering - works on any laptop, no GPU needed. Generates standalone HTML files |
| **Structure Parsing** | BioPython | Industry-standard for PDB/mmCIF parsing, well-maintained |
| **GPU Management** | PyTorch (CUDA) | ESMFold uses PyTorch. Auto-detection handles work laptop (CPU) vs home laptop (RTX 3070 Ti) |

---

## 2. Data Flow

### Full Pipeline Workflow

```
┌─────────────────┐
│ INPUT           │
│ - Nanobody seq  │
│ - AAV serotype  │
│ - Insertion site│
└────────┬────────┘
         │
         v
┌─────────────────────────┐
│ 1. FETCH CAPSID         │  fetch.py
│ - Download from RCSB    │  → Downloads 3UX1.pdb (AAV9)
│ - Cache locally         │
└────────┬────────────────┘
         │
         v
┌─────────────────────────┐
│ 2. PARSE STRUCTURE      │  parse.py
│ - Extract VP3 sequence  │  → MAADGYLPDWLED...
│ - Identify VR loops     │  → VR-VIII: 585-596
│ - Get coordinates       │
└────────┬────────────────┘
         │
         v
┌─────────────────────────┐
│ 3. PREDICT NANOBODY     │  nanobody.py
│ - ESMFold/ColabFold     │  → nanobody.pdb
│ - GPU if available      │  (3D structure)
│ - ~3-5 min on 8GB GPU   │
└────────┬────────────────┘
         │
         v
┌─────────────────────────┐
│ 4. CREATE FUSION        │  fusion.py
│ - Build fusion sequence │  → VP_N + linker + Nb + linker + VP_C
│ - Predict fold          │  → fusion.pdb
│ - Memory estimation     │
└────────┬────────────────┘
         │
         v
┌─────────────────────────┐
│ 5. VISUALIZE            │  render.py
│ - Generate HTML         │  → interactive_view.html
│ - Highlight VR loops    │  (28 KB - 28 MB depending on structure)
│ - Color-code domains    │
└────────┬────────────────┘
         │
         v
┌─────────────────┐
│ OUTPUT          │
│ - Fusion PDB    │
│ - HTML viewer   │
│ - Metadata JSON │
└─────────────────┘
```

---

## 3. Core Components Deep Dive

### 3.1 Structure Prediction (`modeling/nanobody.py`)

#### **Why Two Prediction Methods?**

**ESMFold:**
- **Pros:**
  - Fast (2-5 min for 120-residue nanobody)
  - No MSA required (doesn't need sequence databases)
  - Works well for single-domain proteins
  - Low memory footprint
- **Cons:**
  - Less accurate for multi-domain proteins
  - No confidence metrics like pLDDT from MSA
- **When to use:** Initial models, quick iterations, single nanobodies

**ColabFold (AlphaFold2):**
- **Pros:**
  - Higher quality (uses MSA for evolutionary information)
  - Better for multi-domain proteins (VP-Nb fusions)
  - Provides pLDDT confidence scores
- **Cons:**
  - Slower (10-30 min with MSA)
  - Higher memory usage
  - Requires external MSA database access
- **When to use:** Final models, fusion constructs, publication-quality structures

#### **GPU Memory Management**

The code implements adaptive chunking:

```python
# RTX 3070 Ti (8GB VRAM) optimization
chunk_size = 64  # Default for 8GB
# If OOM error: reduce to 32
# If 12GB+ VRAM: increase to 128 for speed
```

**How chunking works:**
- ESMFold processes attention operations in chunks
- Larger chunks = faster but more memory
- Smaller chunks = slower but fits in limited VRAM
- The code auto-detects GPU and warns if insufficient memory

#### **Device Auto-Detection Logic**

```python
if device == "auto":
    device = "cuda" if torch.cuda.is_available() else "cpu"

if device == "cuda" and not torch.cuda.is_available():
    raise RuntimeError("CUDA requested but not available...")
```

**Why this matters:**
- Work laptop: No GPU → Falls back to CPU (slow but works)
- Home laptop: RTX 3070 Ti → Uses GPU (50-100x faster)
- Prevents silent failures where user expects GPU speed but runs on CPU

---

### 3.2 Fusion Modeling (`modeling/fusion.py`)

#### **Fusion Construct Architecture**

Three fusion modes supported:

**1. N-Terminal Fusion (VP2 Strategy)**
```
[Nanobody] - [Linker] - [VP2 N-term] - [VP3]
```
- Nanobody sits inside capsid initially
- Can be externalized via conformational changes
- Less likely to disrupt capsid assembly

**2. C-Terminal Fusion**
```
[VP3] - [Linker] - [Nanobody]
```
- Nanobody extends from C-terminus
- Maximum surface exposure
- Risk of steric clashes during assembly

**3. Loop Insertion (Most Common)**
```
[VP3_N] - [Linker] - [Nanobody] - [Linker] - [VP3_C]
```
- Insert at VR-IV (451-474) or VR-VIII (585-596)
- Replaces loop or inserts at midpoint
- Balances exposure with assembly compatibility

#### **Linker Design**

Default linker: **`GSG`** (Glycine-Serine-Glycine)
- **Glycine:** Flexible, allows conformational freedom
- **Serine:** Hydrophilic, prevents aggregation
- Short (3 residues) to minimize steric bulk

Alternative: **`GGGGS`** (5 residues) for more flexibility

#### **Memory Estimation for Fusion Prediction**

```python
# Empirical VRAM requirements
if length < 400:   # ~120 Nb + 200 VP = 320
    vram_gb = 8   # Fits on RTX 3070 Ti
elif length < 600:
    vram_gb = 12  # Tight on 8GB, use --msa-mode single_sequence
else:
    vram_gb = 16  # Won't fit on 8GB GPU
```

**Why VP-Nb fusions are challenging:**
- Full VP3: ~533 residues
- Nanobody: ~120 residues
- Fusion: ~660 residues total
- AlphaFold2 memory scales as O(N²) with sequence length
- For 8GB GPU: Must use `--msa-mode single_sequence` (faster, less accurate MSA)

---

### 3.3 Visualization (`visualization/render.py`)

#### **py3Dmol Architecture**

**Key Design Choice:** Generate standalone HTML files, not interactive Python sessions

**Why py3Dmol over PyMOL/ChimeraX?**

| Feature | py3Dmol | PyMOL | ChimeraX |
|---------|---------|-------|----------|
| **Installation** | `pip install` | Complex install | Complex install |
| **Output** | HTML file | Screenshots | Screenshots |
| **Interactivity** | Full (JavaScript) | Session-only | Session-only |
| **Sharing** | Email HTML | Send session file | Send session file |
| **GPU Required** | No | No | Yes (preferred) |
| **File Size** | 28 KB - 28 MB | N/A | N/A |

**py3Dmol strengths for this project:**
1. **Portability:** HTML works on any device (laptop, phone, tablet)
2. **No installation:** Recipient just opens in browser
3. **Embedded structures:** Everything in one file
4. **Fast rendering:** 3Dmol.js uses WebGL (GPU-accelerated in browser)

#### **Color Scheme Strategy**

```python
COLOR_SCHEMES = {
    "capsid_default": {
        "vp3_core": "lightgray",     # Neutral background
        "vr_loops": "steelblue",     # Highlights surface loops
        "nanobody": "firebrick",     # Stands out (red)
        "linker": "gold",            # Marks fusion point
    },
    "publication": { ... },          # Higher contrast for figures
    "colorblind_safe": { ... },      # Deuteranopia/protanopia safe
}
```

**Why these colors:**
- **Gray capsid:** Doesn't compete with highlighted regions
- **Red nanobody:** High visual salience (danger/attention color)
- **Gold linker:** Distinct from both capsid and nanobody
- **Colorblind-safe palette:** Uses orange (#D55E00) instead of red for 8% of male population

#### **HTML File Size Analysis**

| Structure | Size | Why |
|-----------|------|-----|
| Monomer | 28 KB | ~700 atoms, minimal PDB data |
| 60-mer | 28 MB | ~42,000 atoms, full biological assembly |
| With Nb | +5 KB | Additional ~100 atoms per nanobody |

**60-mer Size Breakdown:**
- PDB coordinates: ~25 MB (1000x increase from monomer)
- HTML/CSS/JS: ~50 KB
- Embedded 3Dmol.js: Via CDN (not in file)

**Trade-off:** Large files but self-contained and shareable

---

## 4. Design Philosophy & Trade-offs

### 4.1 Why BioPython Instead of MDAnalysis?

| Feature | BioPython | MDAnalysis |
|---------|-----------|------------|
| **PDB Parsing** | Robust, mature | More powerful |
| **mmCIF Support** | Excellent | Good |
| **Sequence Extraction** | Built-in | Manual |
| **Community** | Larger | Niche (MD sims) |
| **Dependencies** | Minimal | NumPy-heavy |

**Decision:** BioPython
- **Reason:** This is structural biology, not molecular dynamics
- PDB/mmCIF parsing is the core need
- Sequence extraction for fusion is straightforward
- Smaller dependency footprint

### 4.2 Why Store VR Loop Positions in YAML?

**Alternative approaches:**
1. **Hardcode in Python:** Difficult to maintain, version control issues
2. **Database:** Overkill for ~4 serotypes
3. **YAML config:** Easy to edit, version control friendly

**YAML Benefits:**
```yaml
VR-VIII:
  start: 585
  end: 596
  surface_exposed: true
  common_insertion: true
  notes: "Most commonly used for large insertions"
```
- Human-readable
- Comments allowed
- Easy to add new serotypes
- Separates data from code

### 4.3 Local vs Cloud Prediction

**Why not just use AlphaFold Server?**

| Aspect | Local (ESMFold/ColabFold) | AlphaFold Server |
|--------|---------------------------|------------------|
| **Speed** | 2-5 min (GPU) | Hours (queue) |
| **Privacy** | Local data | Uploaded to Google |
| **Iterations** | Unlimited | Rate-limited |
| **Customization** | Full control | Black box |
| **Cost** | GPU power (~$0.50/hr) | Free (limited) |

**Decision:** Hybrid approach
- **Initial screening:** Local ESMFold (fast)
- **Final models:** ColabFold (quality) or AlphaFold Server (if no GPU)
- **Proprietary sequences:** Local only (IP protection)

---

## 5. Hardware Considerations

### 5.1 Two-Laptop Strategy

**Work Laptop (No GPU):**
- Structure parsing ✓
- Visualization generation ✓
- ESMFold prediction ✗ (50-100x slower)
- ColabFold prediction ✗ (impractical)

**Home Laptop (RTX 3070 Ti, 8GB VRAM):**
- All tasks ✓
- ESMFold: 2-5 min per nanobody
- ColabFold: 10-30 min per fusion (with single_sequence MSA)

### 5.2 Memory Bottlenecks

**ESMFold on 8GB VRAM:**
- Max sequence length: ~600 residues (with chunk_size=32)
- Typical nanobody (120 res): Comfortable at chunk_size=64
- VP-Nb fusion (660 res): Tight, may need chunk_size=32 or CPU

**ColabFold on 8GB VRAM:**
- **Full MSA mode (`mmseqs2+prefilter`):**
  - Max ~400 residues
  - VP-Nb fusion will OOM
- **Single-sequence mode:**
  - Max ~800 residues
  - VP-Nb fusion works but less accurate
  - Recommended for 8GB GPU

**Memory Usage Equation:**
```
VRAM(GB) ≈ (sequence_length² / 50000) * MSA_depth_factor

Where MSA_depth_factor:
  - mmseqs2+prefilter: 3-5x
  - single_sequence: 1x
```

---

## 6. Visualization Strategy

### 6.1 Monomer vs 60-mer Visualization

**Monomer View (Small, Fast):**
- Use case: Examine fusion construct details
- Shows: VP + Nanobody with labeled domains
- File size: ~30 KB
- Loads instantly

**60-mer View (Large, Context):**
- Use case: Understand surface distribution, symmetry
- Shows: All 60 subunits, icosahedral symmetry
- File size: ~28 MB
- Loads in 2-5 seconds on modern browser

**Hybrid Approach:**
- Generate both
- Monomer for detailed inspection
- 60-mer for overall capsid architecture

### 6.2 Surface vs Cartoon Representation

**Cartoon:**
- Shows protein backbone (alpha helices, beta sheets)
- Good for understanding fold
- Fast rendering
- **Use:** Fusion construct examination

**Surface:**
- Shows van der Waals surface (solvent-accessible area)
- Good for understanding steric clashes
- Slower rendering (especially for 60-mer)
- **Use:** Checking if nanobody fits without clashing

### 6.3 Interactive Features

**Implemented:**
- Rotation/zoom/pan (3Dmol.js native)
- Chain highlighting (color-coded)
- Loop labeling (VR-IV, VR-VIII, etc.)
- Side-by-side comparison (with/without nanobody)

**Possible Future Additions:**
- Distance measurements between chains
- Hydrogen bond visualization
- Electrostatic surface coloring
- Clash detection highlighting

---

## 7. Current Limitations

### 7.1 Structural Prediction Limitations

**ESMFold:**
- No confidence metrics (unlike AlphaFold's pLDDT)
- Less accurate for:
  - Multi-domain proteins
  - Protein-protein interfaces
  - Disordered regions
- Cannot predict quaternary structure (capsid assembly)

**ColabFold:**
- Memory-intensive for long sequences
- MSA quality depends on internet access
- Single-sequence mode sacrifices accuracy for speed

**Fundamental Issue:**
Neither method predicts **capsid assembly with nanobody insertions**
- AlphaFold Multimer: Limited to ~20 chains (we have 60)
- Full capsid folding is beyond current AI capacity
- Assumption: Nanobody doesn't disrupt assembly (needs validation)

### 7.2 Lack of Experimental Validation

**What's Missing:**
- No MD simulations (capsid stability with nanobody)
- No experimental structures (would need cryo-EM)
- No functional assays (does it still infect?)

**Why Not Included:**
- MD simulations: Require weeks of compute (beyond scope)
- Cryo-EM: Requires wet-lab experiments
- This is a **design tool**, not validation pipeline

### 7.3 Linker Optimization

**Current Approach:**
- Fixed linker sequences (GSG, GGGGS)
- No optimization for:
  - Flexibility vs rigidity trade-off
  - Proteolytic stability
  - Immunogenicity

**Better Approach (Future):**
- Test multiple linker lengths (3, 5, 10, 15 residues)
- Use Rosetta or similar for linker design
- Predict linker flexibility with MD

### 7.4 VR Loop Numbering Ambiguity

**Problem:**
- Different papers use different numbering schemes
  - VP1 numbering (includes VP1/VP2 N-term)
  - VP3 numbering (starts at VP3 start codon)
  - Some use AAV2 convention, others use serotype-specific

**Our Solution:**
- YAML config explicitly states numbering scheme
- Conversion functions in `parse.py`
- Comments in YAML explain assumptions

**Remaining Risk:**
- User must verify loop positions against their reference
- Off-by-one errors can place nanobody incorrectly

---

## 8. Why This Design Works for AAV VHH Project v1

### Strengths Aligned with Project Goals

1. **Rapid Iteration:**
   - ESMFold: 2-5 min per design
   - Can test 10+ insertion sites per hour

2. **Visual Feedback:**
   - HTML visualizations aid decision-making
   - Easy to share with collaborators

3. **Flexible:**
   - Supports all AAV serotypes (via YAML config)
   - Works on both work laptop (CPU) and home laptop (GPU)

4. **Reproducible:**
   - All predictions are deterministic (given same sequence)
   - Config files version-controlled

5. **Extensible:**
   - Easy to add new serotypes (edit YAML)
   - Modular design allows swapping prediction methods

### Acknowledged Gaps

1. **No Assembly Prediction:**
   - Cannot predict if nanobody disrupts capsid formation
   - Mitigation: Use well-characterized insertion sites (VR-VIII)

2. **No Binding Affinity Prediction:**
   - If nanobody is a binder, doesn't predict KD
   - Mitigation: Use experimentally validated nanobodies

3. **No In Silico Immunogenicity:**
   - Doesn't predict immune response to fusion
   - Mitigation: Use non-immunogenic nanobody scaffolds (e.g., humanized)

---

## 9. Recommended Workflow

### For Initial Design

```bash
# Step 1: Get capsid structure
python -m aav_nanobody_display.cli fetch --serotype AAV9

# Step 2: Predict nanobody structure (fast screening)
python -m aav_nanobody_display.cli predict \
  --sequence-file nanobody.fasta \
  --method esmfold \
  --device cuda \
  --name my_nanobody

# Step 3: Create fusion at VR-VIII
python -m aav_nanobody_display.cli fuse \
  --nanobody-file nanobody.fasta \
  --insertion-site VR-VIII \
  --method esmfold \
  --name AAV9_VRIII_Nb

# Step 4: Visualize
python -m aav_nanobody_display.cli visualize \
  --pdb AAV9_VRIII_Nb_esmfold.pdb \
  --highlight-loops VR-VIII \
  --output fusion_view.html
```

### For Final Publication Model

```bash
# Use ColabFold for highest quality
python -m aav_nanobody_display.cli fuse \
  --nanobody-file nanobody.fasta \
  --insertion-site VR-VIII \
  --method colabfold \
  --name AAV9_VRIII_Nb_final

# Generate publication-quality visualization
python -m aav_nanobody_display.cli visualize \
  --pdb AAV9_VRIII_Nb_final.pdb \
  --color-scheme publication \
  --surface \
  --output figure_1.html
```

---

## 10. Future Enhancements (Phase III+)

1. **Rosetta Integration:**
   - Dock nanobody onto capsid surface
   - Optimize linker geometry
   - Minimize steric clashes

2. **MD Simulations:**
   - GROMACS/AMBER for 100 ns simulations
   - Assess capsid stability with nanobody
   - Predict conformational changes

3. **High-Throughput Screening:**
   - Batch process 100+ nanobody sequences
   - Automated scoring (clash detection, solvent exposure)
   - Rank-order by predicted success

4. **Web Interface:**
   - Upload nanobody sequence via browser
   - Select insertion site from dropdown
   - Download prediction + visualization
   - No local installation needed

5. **Experimental Feedback Loop:**
   - Parse cryo-EM density maps
   - Compare prediction vs experiment
   - Refine prediction models based on accuracy

---

## Conclusion

This pipeline is a **computational microscope** for AAV capsid engineering. It accelerates the design-build-test cycle by providing rapid, visual feedback on nanobody insertion strategies.

**Core Philosophy:**
- Fast over perfect (ESMFold for screening)
- Visual over numerical (HTML views over RMSD tables)
- Flexible over rigid (YAML configs, modular code)
- Reproducible over convenient (version-controlled configs)

**Best Used For:**
- Screening insertion sites
- Comparing nanobody orientations
- Generating hypotheses for wet-lab testing

**Not a Substitute For:**
- Experimental validation (cryo-EM, infectivity assays)
- Full MD simulations (capsid assembly dynamics)
- Immunogenicity prediction

**Success Metric:**
If a computational biologist can go from nanobody sequence to publishable visualization in < 30 minutes, the system succeeds.

---

**Document Version:** 1.0
**Last Updated:** December 2025
**Author:** AAV VHH Project v1 Team
