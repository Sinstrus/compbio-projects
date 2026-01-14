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
| **AVD005** | AVD003 | EF1a-VP1-VHH3-ALPL-bGH | 6,690 bp |
| **AVD006** | AVD002 | Rep2Mut2Cap9-VP1-VHH3-ALPL | 7,521 bp |

## Design Features

- **VHH insertion site:** VR-IV at amino acid 456 of VP1
- **Linker design:** D2 asymmetric - N-terminal (GGGGS)x4, C-terminal direct fusion
- **VP2 knockout:** ACG→ACC at codon 138 (silent)
- **VP3 knockout:** ATG→CTG at codon 203 (M→L)

## Synthetic Fragments for Ordering

| Construct | Boundaries | Length | Cost |
|-----------|-----------|--------|------|
| AVD005 | AvrII (1676) → BsrGI (3358) | 1,683 bp | $420-$589 |
| AVD006 | SmaI (2519) → BsrGI (4216) | 1,698 bp | $424-$594 |
| **TOTAL** | | **3,381 bp** | **$844-$1,183** |

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
