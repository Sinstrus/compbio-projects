# Silent Restriction Site Finder

A high-performance Python tool for synthetic biology that identifies restriction enzyme recognition sites that can be introduced into DNA sequences with minimal mutations while maintaining the encoded protein sequence (silent mutations).

## Features

- **IUPAC Ambiguity Support**: Handles all IUPAC nucleotide codes (R, Y, N, W, etc.) in restriction enzyme recognition sequences
- **Edit Distance Algorithm**: Efficient Levenshtein distance calculation with IUPAC pattern matching
- **Silent Mutation Detection**: Automatically classifies mutations as silent, missense, frameshift, or premature stop
- **Comprehensive Enzyme Database**: 250+ restriction enzymes hardcoded (no external files needed)
- **Fast Performance**: Processes 10kb sequences in <5 seconds
- **Frame Detection**: Automatically finds protein coding region in all reading frames

## Installation

No installation required! Just needs Python 3.7+

```bash
chmod +x silent_sites.py
```

## Usage

### Basic Usage

```bash
python silent_sites.py \
    --dna "ATGGCTAGCGATATCGAATTCGGATCCAAGCTT..." \
    --protein "MASDIGGSKLF..." \
    --mutations 2
```

### Arguments

- `--dna` (required): DNA sequence (case insensitive)
- `--protein` (required): Protein sequence in single-letter codes
- `--mutations` (default: 2): Maximum number of nucleotide edits allowed
- `--min-length` (default: 6): Minimum length of restriction sites to consider
- `--roi` (optional): Region of interest - DNA subsequence to search for mutations (must be found within `--dna`)
- `--max-mutations-cap` (default: 3): Safety limit to prevent slow execution
- `--show-all`: Show both silent and non-silent mutations (default: silent only)

### Example with Real Sequences

```bash
# GFP example (partial)
python silent_sites.py \
    --dna "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA" \
    --protein "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK" \
    --mutations 1 \
    --min-length 6
```

### Region of Interest (ROI) - For Large Sequences

When analyzing large sequences (>1kb), you can dramatically improve performance by limiting the mutation search to a specific region while still tracking uniqueness across the entire sequence.

```bash
# Example: 5kb plasmid with GFP insert - only mutate within GFP coding region
python silent_sites.py \
    --dna "AGCTTGCATGC...FULL_5KB_PLASMID..." \
    --protein "MVSKGEELFTGVVPILVELD..." \
    --roi "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGG...GFP_CODING_REGION..." \
    --mutations 2
```

**Benefits:**
- **10-20x faster** - Only searches for mutations within the ROI
- **Dual uniqueness tracking** - Shows if site is unique in ROI AND in full DNA
- **Practical for large constructs** - Makes `--mutations 2` feasible for plasmids/constructs

**Output with ROI:**
```
Position Enzyme    Site       Edits  Mutations  Type    Unique(DNA)  Unique(ROI)
---------------------------------------------------------------------------------
1234     EcoRI     GAATTC     1      A1235T     Silent  Unique       Unique
5678     BamHI     GGATCC     0      None       Silent  Non-unique   Unique
```

The dual uniqueness columns tell you:
- **Unique(DNA)**: Is this site unique in the entire plasmid? (Important for diagnostic digests)
- **Unique(ROI)**: Is this site unique within the coding region? (Important for local cloning)

## Output Format

The tool outputs results in two formats:
1. **Console table** - Human-readable formatted table
2. **CSV file** - Machine-readable data file (auto-generated in the same directory)

### Console Table

The formatted table includes these columns:

| Position | Enzyme | Site Sequence | Edits | Mutations | Type | Uniqueness |
|----------|--------|---------------|-------|-----------|------|------------|
| 123 | EcoRI | GAATTC | 1 | A125T | Silent | Unique |
| 456 | BamHI | GGATCC | 2 | C457G, T458A | Silent | Non-unique (3 sites) |

**Columns:**
- **Position**: 1-indexed start position of the restriction site in DNA
- **Enzyme**: Name of the restriction enzyme
- **Site Sequence**: Recognition sequence (may contain IUPAC codes)
- **Edits**: Number of mutations required (0 = already present, no changes needed)
- **Mutations**: Specific changes needed (e.g., "A123T" = change A at position 123 to T)
- **Type**: Impact classification
- **Uniqueness**: Whether the site appears only once after mutation

**Mutation Types:**
- `Silent`: No amino acid change
- `Non-Silent (Missense, X aa changed)`: Amino acid substitutions
- `Non-Silent (Frameshift)`: Insertion/deletion causing frameshift
- `Non-Silent (Premature Stop)`: Creates stop codon

**Uniqueness Values:**
- `Unique`: The restriction site appears only once in the modified sequence (ideal for diagnostic digestion)
- `Non-unique (X sites)`: The site appears multiple times (X = total count including the new site)

### CSV Output

Results are automatically exported to a CSV file with timestamp:
```
silent_sites_results_YYYYMMDD_HHMMSS.csv
```

**CSV Columns:**
- Position
- Enzyme
- Site_Sequence
- Edits_Required
- Mutations
- Type
- Uniqueness
- Original_Window (original DNA sequence at that position)

**Example:**
```csv
Position,Enzyme,Site_Sequence,Edits_Required,Mutations,Type,Uniqueness,Original_Window
34,EcoRI,GAATTC,0,None (exact match),Silent,Unique,GAATTC
40,BamHI,GGATCC,0,None (exact match),Silent,Unique,GGATCC
```

The CSV file can be opened in Excel, imported into databases, or processed programmatically.

## Algorithm Details

### 1. Protein Location
- Translates DNA in all 3 forward reading frames
- Finds exact match for protein sequence
- Identifies coding region boundaries

### 2. Edit Distance Calculation
- Uses dynamic programming (Levenshtein distance)
- Custom matching function for IUPAC ambiguity codes
- Example: 'A' matches 'N' (no edit cost), 'A' does not match 'C' (1 edit)

### 3. IUPAC Code Matching
Supports all standard codes:
- `R` = A or G (puRine)
- `Y` = C or T (pYrimidine)
- `N` = any base
- `W` = A or T (Weak)
- `S` = G or C (Strong)
- etc.

### 4. Mutation Impact Analysis
For each candidate site:
1. Apply proposed mutations to DNA
2. Translate modified coding region
3. Compare protein sequences
4. Classify impact (silent, missense, frameshift, stop)

## Performance

- **10kb sequence**: <5 seconds
- **Complexity**: O(n × m × k) where:
  - n = DNA sequence length
  - m = number of enzymes
  - k = average enzyme site length

Optimizations:
- Early termination for edit distances > threshold
- Efficient IUPAC matching (set operations)
- Minimal sequence copying

## Restriction Enzyme Database

Includes 250+ commonly used enzymes:
- Type II enzymes (EcoRI, BamHI, HindIII, etc.)
- Rare cutters (NotI, AscI, PacI, etc.)
- IUPAC degenerate sites (BsaI, BbsI, etc.)
- Nicking enzymes (Nt.BstNBI, Nb.BsmI, etc.)

## Example Use Cases

### 1. Golden Gate Assembly
Find BsaI/BbsI sites with 1-2 mutations:
```bash
python silent_sites.py --dna "..." --protein "..." --mutations 2 --min-length 6
```

### 2. Gibson Assembly
Identify unique restriction sites for verification (look for "Unique" in output):
```bash
python silent_sites.py --dna "..." --protein "..." --mutations 1 --min-length 8
```

### 3. Diagnostic Digestion
Find unique rare cutters (NotI, PacI) for clone verification:
```bash
python silent_sites.py --dna "..." --protein "..." --mutations 1 --min-length 8
# Look for enzymes marked "Unique" in the Uniqueness column
```

## Limitations

- **Codons only**: Does not optimize codon usage or avoid RNA secondary structures
- **Frame detection**: Only finds exact protein matches (no partial matches or mismatches)
- **Forward strand**: Does not check reverse complement
- **Standard genetic code**: Uses universal code only

## Advanced Options

### Show All Mutations (Including Non-Silent)
```bash
python silent_sites.py --dna "..." --protein "..." --show-all
```

### Increase Mutation Threshold
```bash
python silent_sites.py --dna "..." --protein "..." --mutations 3 --max-mutations-cap 5
```

### Filter for Longer Sites Only
```bash
python silent_sites.py --dna "..." --protein "..." --min-length 8
```

## FAQ

**Q: Why aren't any sites found?**
- Try increasing `--mutations` threshold
- Reduce `--min-length` to include shorter sites
- Use `--show-all` to see non-silent options

**Q: What if my protein isn't found?**
- Check that DNA encodes the protein exactly
- Verify you're using the correct genetic code
- Ensure DNA is in 5' to 3' orientation (forward strand)

**Q: How do I interpret IUPAC codes?**
- N = any base (most flexible)
- R/Y/S/W/K/M = 2 bases
- B/D/H/V = 3 bases
- See IUPAC_CODES dictionary in source

**Q: Can I add my own enzymes?**
- Yes! Edit the `RESTRICTION_ENZYMES` list in the script
- Format: `("Name", "SEQUENCE")`
- Supports any IUPAC codes

## License

MIT License - free for academic and commercial use

## Citation

If you use this tool in published research, please cite:
```
Silent Restriction Site Finder (2025)
https://github.com/[your-repo]/silent_sites
```

## Contributing

Contributions welcome! Areas for improvement:
- Reverse strand search
- Alternative genetic codes (mitochondrial, etc.)
- Codon optimization integration
- RNA structure prediction
- Export to GenBank/SnapGene formats
