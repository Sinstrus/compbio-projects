# AVD VHH Display Platform - Phase I: ALPL Target

**Project:** AAV VHH Display Platform
**Target:** ALPL (Alkaline Phosphatase)
**Date:** 2026-01-14
**Status:** Ready for GenScript ordering

---

## Overview

This project creates VP1-only AAV capsid constructs with anti-ALPL VHH3 nanobody displayed at Variable Region IV (VR-IV). The constructs use trans-complementation strategy where VP1-VHH is mixed with VP2+3 helper for capsid assembly.

## Constructs

| Construct | Parent | Description | Size |
|-----------|--------|-------------|------|
| **AVD005** | AVD003 | EF1a-VP1-VHH3-ALPL-bGH | 6,705 bp |
| **AVD006** | AVD002 | Rep2Mut2Cap9-VP1-VHH3-ALPL | 7,536 bp |

## Design Features

- **VHH insertion site:** VR-IV at amino acid 456 of VP1 (optimal for surface display)
- **Linker design:** D2 asymmetric flexible (from Biogen patent/literature)
  - **N-terminal:** (GGGGS)×4 = 20 aa = 60 bp (long, flexible "flagpole")
  - **C-terminal:** (GGGGS)×1 = 5 aa = 15 bp (SHORT anchor, NOT direct fusion)
  - **Rationale:** 4:1 asymmetry provides >10× enhancement vs symmetric designs
  - **See:** `docs/DESIGN_PATTERN_Asymmetric_Linkers.md` for full biophysical analysis
- **VP2 knockout:** ACG→ACC at codon 138 (silent)
- **VP3 knockout:** ATG→CTG at codon 203 (M→L)
- **Result:** VP1-only expression with VHH at 3-fold spike apex

### v2.0 Optimization (dnachisel)

The VHH+linker insert was optimized to address synthesis issues:

| Metric | Original | Optimized |
|--------|----------|-----------|
| Linker GC | 93.3% | 66.7% |
| VHH GC | 70.3% | 59.9% |
| Insert GC | 73.2% | 60.9% |
| Direct repeats | 718 | 8 |

Linker uses varied codons (GGT/GGA/GGC for Gly, TCT/TCA/AGT for Ser) instead of repetitive high-GC codons.

## Synthetic Fragments for Ordering

| Construct | Boundaries | Length | Cost |
|-----------|-----------|--------|------|
| AVD005 | AvrII (1676) → BsrGI (3373) | 1,698 bp | $424-$594 |
| AVD006 | SmaI (2519) → BsrGI (4231) | 1,713 bp | $428-$599 |
| **TOTAL** | | **3,411 bp** | **$852-$1,193** |

## Directory Structure

```
AVD_VHH_Display_ALPL/
├── plasmids/           # Source and final plasmid files (.dna, .gb)
├── synthetic_fragments/  # FASTA files for ordering, CSV spreadsheet
├── docs/               # Design verification, ordering instructions
├── scripts/            # Build and analysis scripts
└── README.md           # This file
```

## Files

### Plasmids
- `AVD001-004.dna` - Parent plasmids (SnapGene format)
- `AVD005-EF1A-VP1-VHH3-ALPL-bGH.gb` - Final construct (GenBank)
- `AVD006-Rep2Mut2Cap9-VP1-VHH3-ALPL.gb` - Final construct (GenBank)

### Synthetic Fragments
- `AVD005_synthetic_fragment_FINAL.fasta` - Order this (AvrII to BsrGI)
- `AVD006_synthetic_fragment_FINAL.fasta` - Order this (SmaI to BsrGI)
- `SYNTHETIC_FRAGMENTS_FINAL.csv` - Full ordering spreadsheet

### Documentation
- `FINAL_ORDERING_INSTRUCTIONS.md` - GenScript order form info
- `GENSCRIPT_ORDERING_SUMMARY.md` - Complete ordering summary
- `DESIGN_VERIFICATION_AVD005_AVD006.md` - Design verification report
- `DESIGN_PATTERN_Asymmetric_Linkers.md` - **Comprehensive analysis of D2 linker strategy**

### Scripts
- `build_avd005_006_genbank.py` - Builds final GenBank files
- `design_avd005_006.py` - Initial design analysis

## Timeline

- **Target:** Mouse experiments by end of March 2026
- **Order synthesis:** Week 1
- **Receive fragments:** Week 3-4
- **Clone and verify:** Week 5-7
- **Ready for transfection:** Week 8

## Notes

- AVD006 uses SmaI (blunt end) - optimize ligation conditions
- Both fragments verified to contain ONLY intentional mutations
- See `Lessons_learned.md` in parent directory for design lessons
