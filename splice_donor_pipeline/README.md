# Splice Donor H-Bond Prediction Pipeline

A Python pipeline for predicting splice donor strength based on hydrogen bonding between the 5' splice site (5'SS) and U1 snRNA.

## Overview

This tool:
1. Scores H-bonds for any 11-mer 5'SS sequence (positions -3 to +8)
2. Generates synthetic sequences with target H-bond counts
3. Exports results as CSV with visualization

## Background: 5'SS–U1 snRNA Base Pairing

The 5' splice site base-pairs with U1 snRNA in an antiparallel manner:

```
5'SS:   -3 -2 -1 | +1 +2 +3 +4 +5 +6 +7 +8
U1:     11 10  9 |  8  7  6  5  4  3  2  1
U1nt:    G  U  C |  C  A  Ψ  Ψ  A  C  U  G
```

### H-bond Rules

| Base Pair | H-bonds | Notes |
|-----------|---------|-------|
| G-C, C-G | 3 | Watson-Crick |
| A-U, U-A | 2 | Watson-Crick |
| G-U, U-G | 2 | Wobble pair |
| G-Ψ, A-Ψ | 2 | Pseudouridine pairs like U |
| Mismatch | 0 | Any other combination |

### Optimal Sequence

| Position | U1 Partner | Optimal 5'SS (DNA) | Max H-bonds | Constraint |
|----------|------------|-------------------|-------------|------------|
| -3 | G | C | 3 | None |
| -2 | U | A | 2 | None |
| -1 | C | G | 3 | None |
| +1 | C | G | 3 | **Invariant** |
| +2 | A | T | 2 | **Invariant** |
| +3 | Ψ | A or G | 2 | None |
| +4 | Ψ | A or G | 2 | None |
| +5 | A | T | 2 | None |
| +6 | C | G | 3 | None |
| +7 | U | A | 2 | None |
| +8 | G | C | 3 | None |

**Theoretical maximum: 27 H-bonds** (all positions optimal)

## Installation

The pipeline requires Python 3.7+ and uses these packages:
- `biopython` (optional, for Seq operations)
- `matplotlib` (optional, for distribution plot)
- `logomaker` (optional, for sequence logos)

```bash
pip install matplotlib logomaker pandas
```

## Usage

### Score a Single Sequence

```bash
python scripts/tools/splice_donor_hbonds.py --score CAGGTAAGTAT
```

Output:
```
Sequence: CAGGTAAGTAT
Total H-bonds: 21
Per position:
  -3: 3
  -2: 2
  -1: 3
  +1: 3
  +2: 2
  +3: 2
  +4: 2
  +5: 0
  +6: 0
  +7: 2
  +8: 2
Longest continuous stretch: 7 (starting at -3)
Number of base-paired positions: 9/11
Valid splice site (GT at +1/+2): Yes
```

### Validate Against Known Sequences

```bash
python scripts/tools/splice_donor_hbonds.py --validate --output-dir projects/splice_donor_pipeline
```

### Generate All Sequences

```bash
# Generate and export selected variants (100 per H-bond tier)
python scripts/tools/splice_donor_hbonds.py --generate --output-dir projects/splice_donor_pipeline/output

# Export ALL 262,144 sequences
python scripts/tools/splice_donor_hbonds.py --generate --export-full --output-dir projects/splice_donor_pipeline/output

# With visualization
python scripts/tools/splice_donor_hbonds.py --generate --plot --output-dir projects/splice_donor_pipeline/output

# With sequence logos (requires logomaker)
python scripts/tools/splice_donor_hbonds.py --generate --logos --output-dir projects/splice_donor_pipeline/output
```

### Full Pipeline

```bash
python scripts/tools/splice_donor_hbonds.py \
    --generate \
    --summary \
    --plot \
    --output-dir projects/splice_donor_pipeline/output
```

## Output Files

| File | Description |
|------|-------------|
| `output/synthetic_5ss_sequences.csv` | Selected variants (default: 100 per tier) |
| `output/all_5ss_sequences.csv` | All 262,144 sequences (with `--export-full`) |
| `output/hbond_distribution.png` | Histogram of sequences by H-bond count |
| `output/sequence_logos/logo_*hbonds.png` | Sequence logos per tier |
| `validation/known_sequences_check.txt` | Validation results |

### CSV Format

| Column | Description | Example |
|--------|-------------|---------|
| sequence | 11-mer DNA sequence | CAGGTAAGTAC |
| total_hbonds | Sum of H-bonds | 21 |
| hbonds_pos_m3 | H-bonds at -3 | 3 |
| hbonds_pos_m2 | H-bonds at -2 | 2 |
| hbonds_pos_m1 | H-bonds at -1 | 3 |
| hbonds_pos_p1 | H-bonds at +1 | 3 |
| ... | (all 11 positions) | ... |
| longest_stretch | Consecutive base pairs | 7 |
| stretch_start_pos | Position of stretch start | -3 |
| num_base_pairs | Positions with >0 H-bonds | 9 |
| notes | Annotations | |
| generation_method | How sequence was generated | enumeration |

## Testing

Run unit tests:

```bash
cd scripts/tools
pytest tests/test_splice_donor_hbonds.py -v
```

## API Usage

```python
from scripts.tools.splice_donor_hbonds import (
    calculate_hbonds,
    generate_all_valid_sequences,
    is_valid_splice_site,
)

# Score a sequence
result = calculate_hbonds("CAGGTAAGTAT")
print(f"Total H-bonds: {result.total_hbonds}")
print(f"Per position: {result.per_position}")
print(f"Longest stretch: {result.longest_stretch}")

# Generate all sequences grouped by H-bond count
all_seqs = generate_all_valid_sequences()
print(f"Sequences with 20 H-bonds: {len(all_seqs[20])}")

# Check if valid splice site
print(is_valid_splice_site("CAGGTAAGTAT"))  # True
```

## References

- Freund M, et al. (2005) A novel approach to describe a U1 snRNA binding site. *Nucleic Acids Res.* 33(3):1125-1135
- Zhuang Y, Weiner AM (1986) A compensatory base change in U1 snRNA suppresses a 5' splice site mutation. *Cell* 46(6):827-835

## Directory Structure

```
projects/splice_donor_pipeline/
├── README.md                           # This file
├── output/
│   ├── synthetic_5ss_sequences.csv     # Main output table
│   ├── hbond_distribution.png          # Visualization
│   └── sequence_logos/                 # Per-tier logos
└── validation/
    └── known_sequences_check.txt       # Literature validation results
```
