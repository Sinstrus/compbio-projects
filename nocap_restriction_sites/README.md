# NoCap Restriction Sites Project

**Status:** Completed
**Date:** 2025-12-17
**Plasmid:** BASE-DRAFT-AAV9-RepCap-NOCAP.gb

## Overview

This project designed and validated a set of six unique restriction enzyme sites across the AAV9 RepCap coding sequence, enabling modular cloning strategies for capsid engineering. All sites require only single silent mutations and avoid critical functional boundaries.

## Objective

Identify restriction enzyme sites in six strategic regions of the AAV9 capsid coding sequence for:
- Modular capsid engineering
- Peptide insertion sites
- Golden Gate assembly compatibility
- Minimal sequence changes (single silent mutations preferred)

## Key Accomplishments

### Final Site Selection

Six unique restriction sites identified, each requiring one silent mutation:

| Region | Enzyme | Position | Recognition | Mutation | Verified |
|--------|--------|----------|-------------|----------|----------|
| Region 1 | SmaI | 2505 | CCCGGG | T→C silent | ✅ |
| Region 2 | BbvCI | 2828 | CCTCAGC | C→A silent | ✅ |
| Region 3 | AgeI | 3583 | ACCGGT | G→C silent | ✅ |
| Region 4 | BsrGI | 3780 | TGTACA | C→A silent | ✅ |
| Region 5 | BmtI | 3937 | GCTAGC | C→T silent | ✅ |
| Region 6 | BstZ17I | 4163 | GTATAC | A→T silent | ✅ |

### Critical Bug Fixes

Discovered and corrected two major bugs in the `silent_sites.py` tool:

1. **Frame offset calculation** - Used incorrect formula `(start - vp1_start) % 3`
   - Fixed: `frame_offset = (3 - position_in_codon) % 3 if position_in_codon != 0 else 0`

2. **Impact classification** - Hardcoded `frame=0` instead of detected frame
   - Fixed: Line 903 now uses `frame=frame`

### Tools Created

- **silent_sites_curated.py** - Curated list of 61 reliable restriction enzymes
  - Excludes nicking enzymes, double-cutters, and fussy enzymes
  - Includes common 6-bp/8-bp cutters and reliable Type IIS enzymes

## Project Structure

```
nocap_restriction_sites/
├── README.md                          # This file
├── analysis/                          # Analysis scripts
│   ├── analyze_six_regions_CORRECTED.py    # Final corrected analysis
│   ├── analyze_region1_revised.py
│   ├── analyze_region6_narrowed.py
│   ├── analyze_restriction_sites.py
│   ├── analyze_six_regions.py
│   ├── apply_mutations_to_plasmid.py
│   ├── check_region1_alternatives.py
│   ├── generate_detailed_report.py
│   ├── generate_final_report.py
│   └── generate_revised_six_regions_report.py
├── verification/                      # Verification scripts
│   ├── verify_all_mutation_positions.py
│   ├── verify_agei.py
│   ├── verify_avrii_mutation.py
│   ├── verify_bstz17i_region6.py
│   ├── verify_features_preserved.py
│   ├── verify_region3_alternatives.py
│   ├── verify_region4_alternatives.py
│   ├── verify_region5.py
│   ├── verify_smai_mutation.py
│   └── verify_xbai_mutation.py
└── output/                            # Final summaries
    ├── FINAL_six_regions_summary.md   # Comprehensive summary
    └── FINAL_PLASMID_SUMMARY.md       # Original summary
```

## Key Output Files

- **FINAL_six_regions_summary.md** - Complete analysis with all six sites, bug fixes, and verification
- **FINAL_PLASMID_SUMMARY.md** - Original analysis summary

## Region Specifications

1. **Region 1:** Rep68-stop to VP2-start (2415-2775)
   - Narrowed to avoid splice acceptor in VP1 N-terminus

2. **Region 2:** VP2-AAP Intergenic (2779-2890)
   - BbvCI for Golden Gate assembly compatibility

3. **Region 3:** AAP-stop to VR4 (3485-3717)
   - AgeI (reliable workhorse enzyme)

4. **Region 4:** VR4 to VR5 (3745-3825)
   - BsrGI (reliable single-site enzyme)

5. **Region 5:** VR5 to VR8 (3880-4104)
   - BmtI (NheI compatible)

6. **Region 6:** Post-VR8 (4144-4343)
   - Narrowed to 200bp window after VR8
   - BstZ17I at position 4163

## Verification Summary

- ✅ All 6 sites are unique in the 7,078 bp plasmid
- ✅ All mutations are silent (codon synonyms)
- ✅ No boundary overlaps with critical features
- ✅ Region 6 constraint satisfied (within 200bp of VR8 end)
- ✅ Frame offsets calculated correctly for all regions

## How to Re-run Analysis

The main corrected analysis script:

```bash
python3 analysis/analyze_six_regions_CORRECTED.py
```

To verify all mutation positions:

```bash
python3 verification/verify_all_mutation_positions.py
```

## Dependencies

- Python 3.8+
- Biopython
- silent_sites.py tool (in scripts/tools/)

## Next Steps

To apply these mutations to the plasmid:

1. Use `analysis/apply_mutations_to_plasmid.py`
2. Verify with restriction digest
3. Sequence confirm all six sites

## Related Documentation

See the main repository `scripts/tools/` directory for:
- `silent_sites.py` - Core tool for silent mutation identification
- `silent_sites_curated.py` - Curated enzyme list
