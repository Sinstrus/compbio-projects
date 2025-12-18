# DNA Engineer Agent — Design Tools

This directory contains specialized tools for DNA sequence design and analysis.

## Available Tools

### 1. `silent_sites.py` — Silent Restriction Site Finder

**Purpose:** Find restriction enzyme sites that can be introduced into a DNA sequence while preserving the encoded protein (silent mutations only).

**Capability:**
- Finds existing sites (0 mutations required)
- Finds "one-out" sites (1 mutation required)  
- Finds "two-out" sites (2 mutations required)
- Configurable up to any number of mutations (with performance tradeoff)

**When to Use:**
- Adding "handles" for future modifications (e.g., flanking BsaI sites for Golden Gate)
- Creating diagnostic restriction sites for clone verification
- Planning cloning strategies that require unique restriction sites
- Preparing sequences for gene synthesis services (e.g., Genscript FLASH)

**Risk Policy:**
Silent mutations may change codon usage, which could theoretically affect expression. However:
- The marginal utility of adding a restriction site (enables fast modification) far outweighs this small risk
- Codon changes are documented in the output
- The tool only proposes changes — the agent/user decides whether to apply them

**Usage (CLI):**
```bash
# Find sites requiring up to 2 silent mutations
python scripts/tools/silent_sites.py \
    --dna "ATGGCTAGCGATATCGAATTCGGATCC..." \
    --protein "MASDIGGSKLF..." \
    --mutations 2 \
    --min-length 6

# For large sequences, use ROI to limit search area
python scripts/tools/silent_sites.py \
    --dna "FULL_PLASMID_SEQUENCE..." \
    --protein "YOUR_PROTEIN..." \
    --roi "CODING_REGION_ONLY..." \
    --mutations 2
```

**Usage (Programmatic — for agent scripts):**
```python
import sys
sys.path.insert(0, 'scripts/tools')
from silent_sites import find_candidates, find_protein_in_dna

# Find the coding region
dna = "ATGGCTAGCGATATCGAATTCGGATCC..."
protein = "MASDIGGSKLF..."
start, end, frame = find_protein_in_dna(dna, protein)

# Find candidate restriction sites (up to 2 mutations)
candidates = find_candidates(
    dna_seq=dna,
    protein_seq=protein,
    max_mutations=2,
    min_length=6,
    roi_seq=None  # or specify a region of interest for large sequences
)

# Filter for silent mutations only
silent_candidates = [c for c in candidates if c.mutation_type == "Silent"]

# Filter for unique sites only
unique_sites = [c for c in silent_candidates if c.uniqueness_dna == "Unique"]

# Prioritize by number of edits (prefer 0 > 1 > 2)
unique_sites.sort(key=lambda c: c.edits_required)
```

**Output Fields:**
| Field | Description |
|-------|-------------|
| position | 1-indexed start position in DNA |
| enzyme | Restriction enzyme name |
| site_sequence | Recognition sequence (may contain IUPAC codes) |
| edits_required | Number of mutations needed (0 = already present) |
| mutations | List of specific changes (e.g., ["A125T", "G128C"]) |
| mutation_type | "Silent", "Non-Silent (Missense)", etc. |
| uniqueness_dna | "Unique" or "Non-unique (X sites)" in full sequence |
| uniqueness_roi | Same but within region of interest |
| original_window | Original DNA at that position |

**Key Concepts:**
- **"Zero-out" sites:** Sites already present (`edits_required == 0`)
- **"One-out" sites:** Sites requiring exactly 1 mutation
- **"Two-out" sites:** Sites requiring exactly 2 mutations
- **Unique sites:** Sites that appear only once after mutation (ideal for cloning)
- **ROI (Region of Interest):** Limit search to a specific region for performance

**Performance Notes:**
- `--mutations 1`: Fast, good for quick searches
- `--mutations 2`: Slower but finds many more options (recommended default)
- `--mutations 3+`: Use with ROI to avoid long runtimes
- See `TIPS.md` for performance optimization

---

## Integration with Agent Workflow

### When Modifying DNA:

After planning an insertion (e.g., VHH at VR-IV), consider:

1. **Do we need handles for future modifications?**
   - If yes, run `silent_sites.py` on the insertion region
   - Look for BsaI/BbsI sites (Golden Gate compatible)
   - Prefer unique sites with fewer mutations

2. **Prioritize sites:**
   - 0 mutations (already present) > 1 mutation > 2 mutations
   - Unique in full construct > unique in ROI only
   - Golden Gate enzymes (BsaI, BbsI, BsmBI) for swappable parts

3. **Document any codon changes:**
   - Record in the modification report
   - Note original and new codons
   - Flag if changing from high-usage to low-usage codon

4. **Risk assessment:**
   - Silent mutations are LOW risk
   - Multiple mutations = slightly higher risk but still acceptable
   - Document the tradeoff: "2 silent mutations enable Golden Gate swapping"

### Example Workflow:

```
User: "Insert VHH-X at VR-IV with flanking BsaI sites for future swapping"

Agent:
1. Plan the insertion (linker + VHH + linker)
2. Run silent_sites.py on the insertion DNA with --mutations 2
3. Find BsaI sites, sorted by edits_required
4. Select unique, silent sites at 5' and 3' ends
5. Apply mutations to create BsaI sites
6. Document all codon changes in modification report
7. Generate final construct
```

---

## Adding New Tools

To add a new design tool:

1. Place the script in `scripts/tools/`
2. Document it in this file
3. If it should be a mandatory step, add to `AGENT_INSTRUCTIONS_v2.md`
4. Create example usage in `test_data/` if appropriate
