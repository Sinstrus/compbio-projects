# AAV9 Symmetry Axis Explorer - Technical Notes

## Problem Solved: "Geometric Neighbor Confusion"

### Original Issue
The initial implementation used **center-of-mass (COM)** distance to identify symmetry-related chains. This failed because:
- AAV capsid is a hollow spherical shell
- COM of each VP3 monomer is deep in the protein core
- COM distances don't correlate with biological surface interfaces
- Result: Wrong chains selected (e.g., trimer + 2 random chains instead of pentamer)

### Solution: Interface Contact Detection
Switched to **residue-based contact detection** using BioPython's `NeighborSearch`:

#### Interface Definitions
1. **5-Fold Axis (Pentamer/Pore)**
   - Interface residues: VR-I (260-275) + VR-II (326-338) + N-term shoulder (380-395)
   - Cutoff: 15.0 Å
   - Result: 6 chains (A, A-2, A-23, A-3, A-4, A-5)
   - Biological interpretation: Captures the pore-forming pentamer

2. **3-Fold Axis (Trimer/Spike)**
   - Interface residues: VR-IV (450-475) + VR-VIII (580-596)
   - Cutoff: 10.0 Å
   - Result: 3 chains (A, A-10, A-23) ✓ Perfect!
   - Biological interpretation: The interdigitating spike residues

3. **2-Fold Axis (Dimer Interface)**
   - Interface residues: HI loop (500-515) + VR-IX (660-675)
   - Cutoff: 8.0 Å
   - Result: 5 chains (A, A-2, A-23, A-42, A-5)
   - Biological interpretation: Dimer + proximal neighbors

### Visualization Improvements
- Increased context chain opacity from 0.15 → 0.3 (fixes "giant hole" artifact)
- Increased focused chain opacity from 0.9 → 0.95
- VR-VIII (585-596) highlighted in red
- VR-IV (451-474) highlighted in blue

### Performance
- Uses C-alpha atoms only for efficiency
- BioPython NeighborSearch for O(log N) spatial queries
- Runs in ~5-10 seconds on standard laptop (no GPU)

### Key Files
- Script: `scripts/visualize_symmetry_axes.py`
- Output: `aav9_symmetry_explorer.html` (28 MB)
- Structure: 3UX1 biological assembly 1 (mmCIF format)

### Usage
```bash
python scripts/visualize_symmetry_axes.py
# Open aav9_symmetry_explorer.html in browser
```

### Interactive Features

#### Symmetry Axis Selection
- **5-Fold (Blue)**: Pentameric pore formed by VR-I + VR-II
- **3-Fold (Green)**: Trimeric spike with VR-IV/VR-VIII interdigitation
- **2-Fold (Yellow/Orange)**: Dimer interface at HI loop

#### Residue Highlighter
NEW! Identify exposed positions by highlighting individual residues:
- **Input Box**: Type residue number (200-736) and press Enter
- **Quick Jump Buttons**: VR-IV (451), VR-VIII (585)
- **Highlight Color**: Magenta sphere (3.0Å radius) on C-alpha atom only
- **Display**: Shows C-alpha across all 60 chains simultaneously (60 spheres)
- **Use Case**: Find the most extended/exposed residues for binder insertion
- **Performance**: Instant rendering (no lag)

**How to use:**
1. Select a symmetry view (e.g., 3-Fold)
2. Type a residue number (e.g., 585) and press Enter
3. A single magenta sphere appears on the C-alpha of that residue on all chains
4. Visually identify which residues "stick out" from the capsid surface
5. Click "Clear Highlight" to remove

### References
- PDB: 3UX1 (AAV9 capsid structure)
- Interface detection based on AAV structural virology literature
- VR (Variable Region) nomenclature from Govindasamy et al.
