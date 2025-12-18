# Performance Tips for Large Sequences

## Issue: Long Runtime on Large Sequences

For sequences >1kb with `--mutations 2`, the tool may take several minutes to complete.

## Quick Fixes

### 1. Use Region of Interest (ROI) - BEST OPTION
```bash
# Only search for mutations within your coding region
# Example: protein at positions 2726-2911 in a 4821bp sequence
python silent_sites.py \
    --dna "YOUR_FULL_SEQUENCE..." \
    --protein "YOURPROTEIN..." \
    --roi "EXTRACT_JUST_THE_CODING_REGION_DNA..." \
    --mutations 2
```
**Impact**: 10-20x faster! Still tracks uniqueness in full sequence.
**Best when**: You only care about mutating a specific region (e.g., just the coding sequence, not the UTRs/flanking regions).

### 2. Reduce Mutations
```bash
# Instead of --mutations 2, use --mutations 1 or 0
python silent_sites.py --dna "..." --protein "..." --mutations 1
```
**Impact**: 2-5x faster. Most useful sites require 0-1 mutations anyway.

### 3. Increase Minimum Site Length
```bash
# Only search for longer sites (8+ bp)
python silent_sites.py --dna "..." --protein "..." --mutations 2 --min-length 8
```
**Impact**: ~50% fewer enzymes to search.

### 4. Run in Background
```bash
# Run and go get coffee
python silent_sites.py --dna "..." --protein "..." --mutations 2 &
# Check progress
tail -f silent_sites_results_*.csv
```

## Performance Expectations

| Sequence Size | Mutations | Enzymes | Time |
|--------------|-----------|---------|------|
| 500 bp | 2 | 250 | ~5 seconds |
| 1000 bp | 2 | 250 | ~15 seconds |
| 5000 bp | 2 | 250 | ~2-3 minutes |
| 5000 bp | 1 | 250 | ~30 seconds |
| 5000 bp | 2 | 100 (min-length=8) | ~1 minute |

## For Your Specific Case

Your sequence: 4821 bp, 62 aa protein
- Protein found at: position 2726-2911 (frame 1)

### Recommended Command:
```bash
python silent_sites.py \
    --dna "YOUR_DNA..." \
    --protein "GGGGSGGGGSGGGGSGGGGSGGGGSAPGKKRPVEQSPQEPDSSAGIGKSGAQPAKKRLNFGQ" \
    --mutations 1 \
    --min-length 6
```

This should complete in ~30-60 seconds and will find most useful sites.
