# Quick Reference: 6 Restriction Sites for AAV9 Rep-Cap Engineering

**Plasmid:** BASE-DRAFT-AAV9-RepCap-NOCAP.gb
**Date:** 2025-12-17
**Total Mutations Required:** 5 silent mutations (+ 1 site naturally present)

---

## Sites Summary

| Region | Enzyme | Position | Site | Mutations | Codon Change |
|--------|--------|----------|------|-----------|--------------|
| **1. VP1 Unique** | EagI | 2369 | CGGCCG | 1 | GCT→GCG (A→A) |
| **2. VP2-AAP** | BbvCI | 2828 | CCTCAGC | 1 | TCC→TCA (S→S) |
| **3. AAP-VR4** | EcoNI | 3487 | CCTNNNNNAGG | 1 | CGG→AGG (R→R) |
| **4. VR4-VR5** | FseI | 3759 | GGCCGGCC | 1 | GCC→GCG (A→A) |
| **5. VR5-VR8** | BmtI | 3937 | GCTAGC | 1 | GCC→GCT (A→A) |
| **6. Post-VR8** | BaeI | 4318 | ACNNNNGTAYC | 0 | None (exists!) |

**Notes:**
- All sites are **unique** in the 7,078 bp plasmid
- All sites **do not overlap** critical boundaries (start codons, VRs)
- All mutations are **silent** (no amino acid changes) except Region 5 (see below)

---

## Wait - Region 5 Issue!

⚠️ **IMPORTANT:** The codon change for Region 5 (BmtI @ 3937) shows GCG→GTG (A→V), which is **NOT SILENT**. This changes Alanine to Valine.

**Action Required:** Verify this site or select an alternative from Region 5 candidates:
- AgeI @ 3959 (1 mutation)
- EcoNI @ 3970 (1 mutation)
- BspEI @ 3976 (1 mutation)

I need to double-check Region 5's analysis.

---

## Region Boundaries

| Region | Start | End | Length | Description |
|--------|-------|-----|--------|-------------|
| 1 | 2365 | 2775 | 411 bp | VP1 unique (VP1 to VP2 start) |
| 2 | 2779 | 2890 | 112 bp | VP2-AAP intergenic |
| 3 | 3485 | 3717 | 233 bp | AAP stop to VR4 |
| 4 | 3745 | 3825 | 81 bp | VR4 to VR5 |
| 5 | 3880 | 4104 | 225 bp | VR5 to VR8 |
| 6 | 4144 | 4575 | 432 bp | Post-VR8 |

---

## Critical Boundaries Protected

✓ VP1 start (2365)
✓ VP2 start (2776)
✓ AAP start (2891)
✓ VP3 start (2971)
✓ AAP stop (3484)
✓ VR-IV: 3718-3744
✓ VR-V: 3826-3879
✓ VR-VIII: 4105-4143

---

## Primers for Site-Directed Mutagenesis

### Region 1: EagI (T2370G)
```
Forward: 5'-GCCGATGGTTATCTTCCAGATTGGCGGAGGACAACCTTAGTGAAGG-3'
                                      ^
                                    T→G
```

### Region 2: BbvCI (C2832A)
```
Forward: 5'-CAGTCTCCTCAGGAACCGGACTCCTCAGCGGGTATTGGCAAATCG-3'
                                        ^
                                      C→A
```

### Region 3: EcoNI (C3495A)
```
Forward: 5'-CATGATTCCTCAGTAAGGTATCTGACGCTTAATGATGGAAGC-3'
                          ^
                        C→A
```

### Region 4: FseI (A3765C)
```
Forward: 5'-CCAGCAACATGGCCGGCCGTCCAGGGAAGAAACTACATACC-3'
                          ^
                        A→C
```

### Region 5: BmtI (C3939T) - **VERIFY FIRST**
```
Forward: 5'-GGAGTCAAGACCATCGCTAGCAACCTTACCAGCACGGTCCAG-3'
                            ^
                          C→T
```

---

## Implementation Options

**Option 1: Gene Synthesis** (Recommended)
- Order Cap gene (2365-4575) with all 5 mutations
- Cost: ~$220-330 (2,211 bp × $0.10-0.15/bp)
- Guaranteed accuracy

**Option 2: Site-Directed Mutagenesis**
- Use Q5 or QuikChange for individual regions
- 5 separate reactions needed
- More flexible for testing individual sites

---

## Files Generated

- `reports/AAV9_RepCap_SixRegions_Analysis.md` - Full detailed report
- `reports/AAV9_RepCap_QuickReference.md` - This file
- `analyze_six_regions.py` - Analysis script
- `generate_detailed_report.py` - Report generator

---

**Status:** ⚠️ Need to verify Region 5 silent mutation claim before proceeding!
