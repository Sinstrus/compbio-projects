# Changelog

## [1.5.0] - 2025-12-11

### Fixed
- **Critical: Double-Strand Uniqueness Counting**: Fixed uniqueness calculation to properly count restriction sites on both DNA strands
  - **Palindromic sites** (e.g., EcoRI - GAATTC): Appear on both strands at same position, counted once
  - **Non-palindromic sites** (e.g., BbsI - GAAGAC): Can appear independently on either strand, both counted
  - New functions: `reverse_complement_iupac()`, `is_palindromic()`
  - Updated `count_pattern_occurrences()` to check both forward and reverse complement strands
  - Updated `count_pattern_occurrences_in_range()` with same logic for ROI support
  - **Impact**: Previously, non-palindromic sites on the reverse strand were not counted, leading to incorrect uniqueness reporting

### Why This Matters
- **Biological Accuracy**: Restriction enzymes cut double-stranded DNA, so a site can exist on either strand
- **Example Before Fix**:
  ```
  DNA: ...GAAGAC...GTCTTC...
  BbsI reported as "Unique" (only counted GAAGAC on forward strand)
  ```
- **Example After Fix**:
  ```
  DNA: ...GAAGAC...GTCTTC...
  BbsI reported as "Non-unique (2 sites)" (counts both GAAGAC and its reverse complement)
  ```
- **Critical for Cloning**: Ensures you know all positions where an enzyme will actually cut

## [1.4.0] - 2025-12-11

### Added
- **Region of Interest (ROI) Support**: New `--roi` parameter to limit mutation search to a specific DNA subsequence
  - Only searches for mutations within the ROI (dramatically improves performance for large sequences)
  - Tracks uniqueness both within the full DNA sequence AND within the ROI
  - ROI must be found within the provided DNA sequence (validated at runtime)
  - Dual uniqueness columns in output: `Unique(DNA)` and `Unique(ROI)`
  - CSV export includes both `Uniqueness_DNA` and `Uniqueness_ROI` columns when ROI is used
- New helper function: `count_pattern_occurrences_in_range()` for range-constrained pattern counting
  - Counts sites that START within a range (sites can extend beyond boundary)
  - Fixes edge case where sites near ROI boundary weren't counted correctly

### Changed
- Output table now shows dual uniqueness columns when ROI is specified
- Summary statistics now include separate counts for DNA-wide and ROI-specific unique sites
- Performance: Large sequences (>1kb) can now be analyzed by restricting search to coding regions

### Why This Matters
- **Performance**: Analyze only the region you care about (e.g., 300 bp coding region within a 5 kb plasmid)
  - 10-20x speedup for typical use cases
  - Makes `--mutations 2` practical for large sequences
- **Context Awareness**: Know if a site is unique locally (ROI) vs globally (full DNA)
  - Important for cloning strategies where local uniqueness matters for the feature of interest
  - Global uniqueness still tracked for diagnostic purposes
- **Example Use Case**:
  ```bash
  # Search only within the GFP coding sequence but check uniqueness in full plasmid
  python silent_sites.py \
      --dna "FULL_5KB_PLASMID..." \
      --protein "MSKGEEL..." \
      --roi "ATGAGTAAAGGAGAAGAACT..." \
      --mutations 1
  ```

## [1.3.0] - 2025-12-11

### Fixed
- **Stop Codon Handling**: Fixed critical bug where protein search would stop at the first stop codon
  - Added `stop_at_terminator` parameter to `translate()` function
  - `find_protein_in_dna()` now translates through stop codons to find proteins anywhere in the sequence
  - This allows finding proteins in sequences with UTRs, intergenic regions, or multiple ORFs
  - **Impact**: Proteins that appear after stop codons in the DNA sequence can now be found

### Why This Matters
- Many cloning constructs have multiple elements (promoters, UTRs, tags, etc.) that may contain stop codons
- The tool can now handle realistic plasmid sequences, not just clean coding sequences
- Enables analysis of fusion proteins, tagged constructs, and polycistronic sequences

## [1.2.0] - 2025-12-11

### Added
- **CSV Export**: Results are now automatically exported to a CSV file with timestamp
  - Format: `silent_sites_results_YYYYMMDD_HHMMSS.csv`
  - Includes all columns: Position, Enzyme, Site_Sequence, Edits_Required, Mutations, Type, Uniqueness, Original_Window
  - Can be imported into Excel, databases, or processed programmatically
- New function: `export_to_csv()` for structured data export

### Why This Matters
- **Data Analysis**: Import results into R, Python, or Excel for further analysis
- **Record Keeping**: Timestamped files provide permanent records of design sessions
- **Automation**: CSV format enables integration with downstream workflows
- **Collaboration**: Easy to share results with colleagues who need structured data

## [1.1.0] - 2025-12-11

### Added
- **Uniqueness Column**: New output column showing whether each restriction site is unique in the modified sequence
  - `Unique`: Site appears only once (ideal for diagnostic restriction digestion)
  - `Non-unique (X sites)`: Site appears multiple times in the sequence
- Enhanced summary statistics now include unique vs non-unique counts
- New helper functions:
  - `matches_iupac_pattern()`: Check if DNA window matches IUPAC pattern
  - `count_pattern_occurrences()`: Count how many times a pattern appears in sequence

### Changed
- Output table width increased from 100 to 120 characters to accommodate uniqueness column
- Mutation column width reduced from 30 to 25 characters to fit new column

### Why This Matters
When designing constructs for cloning or verification:
- **Unique sites** are perfect for diagnostic digestion (single, predictable cut)
- **Non-unique sites** may still be useful but require careful planning
- Helps prioritize which restriction sites to introduce for downstream applications

## [1.0.0] - 2025-12-11

### Initial Release
- IUPAC ambiguity code support for 250+ restriction enzymes
- Edit distance algorithm with dynamic programming
- Silent vs non-silent mutation classification
- Automatic protein finding in DNA sequence
- Support for frameshifts, missense, and premature stop detection
