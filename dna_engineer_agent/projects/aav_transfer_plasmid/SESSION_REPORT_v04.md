# AAV Transfer Plasmid v4.0 - Session Report

**Date:** 2025-12-18
**Status:** ✅ COMPLETE - ssAAV plasmid generated with cloning strategy

---

## Executive Summary

Generated v04 transfer plasmid using a **cloning strategy** instead of stuffer replacement. The EF1a-VP1-polyA cassette is designed to be synthesized and cloned into HindIII/XbaI sites in the backbone.

**Key Changes from v03:**
1. Replaced XbaI (Region 1) with AvrII to avoid conflict with XbaI cloning site
2. Use full 1194 bp EF1α promoter (v03 was truncated to 1172 bp)
3. Preserve both HindIII and XbaI cloning sites in final plasmid

## Important Corrections Made

**Initial v04 (before corrections):**
- ❌ EF1α was 1172 bp (missing 22 bp)
- ❌ XbaI site was not present in final plasmid

**Corrected v04 (current):**
- ✅ EF1α is 1194 bp (full version from v01)
- ✅ XbaI site preserved at position 3732
- ✅ HindIII site at position 189
- ✅ Proper cloning simulation with both sites regenerated

---

## Motivation for v04

### User Requirements
1. **ITRs must be annotated** - New v02 backbone has ITR128 annotations
2. **Use cloning strategy** - Cassette synthesized and cloned into restriction sites, not stuffer replacement
3. **Cloning sites:** HindIII (upstream) and XbaI (downstream)
4. **Only ssAAV version** needed

### Critical Issue Identified
- **XbaI conflict:** v03 used XbaI as a unique site in VP1 Region 1
- XbaI is also the downstream cloning site in the backbone
- After cloning, XbaI would NOT be unique (one in VP1, one at cloning junction)
- **Solution:** Replace XbaI with AvrII in Region 1

---

## Final Site Selection (v04)

| Region | Enzyme | Type | Mutations | RepCap Position | ssAAV v04 Position |
|--------|--------|------|-----------|----------------|-------------------|
| **1** | **AvrII** | Standard | 1 silent | 2522 | 1522 |
| **2** | **BspEI** | Standard | 2 silent | 2850 | 1850 |
| **3** | **BsmBI** | Type IIS | 1 silent | 3574 | 2574 |
| **4** | **BsrGI** | Standard | 1 silent | 3782 | 2782 |
| **5** | **BmtI** | Standard | 1 silent | 3943 | 2943 |
| **6** | **BstZ17I** | Standard | 1 silent | 4167 | 3167 |

**Total mutations:** 7 silent mutations
**All sites verified:** ✅ Unique in RepCap and ssAAV v04

---

## Mutations Applied

### Region 1: AvrII (NEW - replaces XbaI)
- **Position 2523:** T→A | CTT→CTA (L→L)
- Creates recognition site: CCTAGG
- Only 1 mutation needed (vs 2 for XbaI)
- VP1 amino acid 53

### Region 2: BspEI (kept from v03)
- **Position 2850:** G→C | TCG→TCC (S→S)
- **Position 2853:** T→A | GGT→GGA (G→G)
- VP1 amino acids 162-163

### Region 3: BsmBI (kept from v03)
- **Position 3570:** G→T | TCG→TCT (S→S)
- Type IIS enzyme for Golden Gate assembly
- VP1 amino acid 402

### Regions 4-6: (kept from v03)
- **BsrGI** - Position 3783: C→A (VP1 aa 473)
- **BmtI** - Position 3939: C→T (VP1 aa 525)
- **BstZ17I** - Position 4164: A→T (VP1 aa 600)

---

## Assembly Strategy

### Previous Approach (v01-v03)
```
Replace stuffer region (positions 159-2350) end-to-end with cassette
```

### New Approach (v04)
```
1. Synthesize EF1a-VP1-polyA cassette with HindIII/XbaI flanking sites
2. Digest backbone with HindIII + XbaI
3. Ligate cassette into backbone
4. Result: seamless integration with proper ITR annotations preserved
```

### Cloning Sites in Backbone
- **HindIII:** Position 187 (upstream)
- **XbaI:** Position 2278 (downstream)
- **Removed region:** 187-2278 (2091 bp)
- **Inserted cassette:** 3515 bp (EF1a + VP1 + polyA, without flanking sites)
- **Net change:** +1424 bp

### Backbone Conflict Resolution
- v02 backbone had AvrII at position 2259
- This AvrII was **between** HindIII (187) and XbaI (2278)
- During cloning, this region is removed
- Result: Only the VP1 AvrII remains → **unique site** ✅

---

## Files Generated

### RepCap Plasmid
```
test_data/AAV9-RepCap-NOCAP-v04.gb (7,078 bp)
```
- 7 silent mutations applied to BASE-DRAFT
- All 6 restriction sites verified unique
- AvrII replaces XbaI in Region 1

### Transfer Plasmid
```
test_data/pGS-ssAAV-EF1A-VP1-rBG_v04.gb (6,462 bp)
```
- Assembled using cloning strategy
- EF1α promoter: 1194 bp (full version with intron)
- ITR128 annotations: positions 27-154 and 3758-3885
- Cloning sites preserved: HindIII (189) and XbaI (3732)
- All 6 VP1 sites verified unique
- Only ssAAV version generated

### Analysis Scripts
```
projects/aav_transfer_plasmid/analysis/
├── find_avrii_region1.py           # Identify AvrII silent mutation
├── generate_v04_repcap.py          # Generate RepCap v04
└── assemble_v04_cloning.py         # Assemble transfer plasmid (cloning strategy)
```

---

## Validation Results

### RepCap v04 Verification
```
✅ AvrII      1 site @ position 2522
✅ BspEI      1 site @ position 2850
✅ BsmBI      1 site @ position 3574
✅ BsrGI      1 site @ position 3782
✅ BmtI       1 site @ position 3943
✅ BstZ17I    1 site @ position 4167
```

### Transfer Plasmid v04 Verification

**Cloning sites:**
```
✅ HindIII    1 site @ position 189
✅ XbaI       1 site @ position 3732
```

**VP1 unique sites:**
```
✅ AvrII      1 site @ position 1544
✅ BspEI      1 site @ position 1872
✅ BsmBI      1 site @ position 2596
✅ BsrGI      1 site @ position 2804
✅ BmtI       1 site @ position 2965
✅ BstZ17I    1 site @ position 3189
```

### ITR Annotations
```
✅ ITR128 (left):  27-154
✅ ITR128 (right): 3758-3885
```

### EF1α Promoter
```
✅ Length: 1194 bp (positions 193-1386)
✅ Full version with intron
```

**✅ ALL REQUIREMENTS MET**

---

## Comparison: v03 vs v04

| Aspect | v03 | v04 |
|--------|-----|-----|
| **Region 1 enzyme** | XbaI (2 mutations) | AvrII (1 mutation) ✅ |
| **Assembly strategy** | Stuffer replacement | Cloning into sites ✅ |
| **ITR annotations** | Present but shifted | Annotated correctly ✅ |
| **Backbone** | v01 (no ITRs) | v02 (with ITRs) ✅ |
| **Cloning sites** | N/A | HindIII/XbaI ✅ |
| **EF1α length** | 1172 bp (truncated) | 1194 bp (full) ✅ |
| **XbaI site preserved** | N/A | Yes (position 3732) ✅ |
| **Total mutations** | 8 | 7 ✅ |
| **Plasmid size** | 6,436 bp | 6,462 bp |

---

## Synthesis Strategy

### For Wet Lab Implementation

1. **Synthesize cassette with overhangs:**
   ```
   5'-AAGCTT-[EF1a]-[VP1]-[polyA]-TCTAGA-3'
   ```
   - Length: ~3,527 bp (including HindIII and XbaI sites)
   - 7 silent mutations in VP1 (see mutations table above)

2. **Digest backbone:**
   ```
   pGS-ssAAV-ITR128-Amp-empty_v02.gb
   Digest with: HindIII + XbaI
   ```

3. **Ligate:**
   ```
   Insert cassette between HindIII and XbaI
   Transform into E. coli
   Select on ampicillin
   ```

4. **Verify:**
   ```
   - Sanger sequencing across VP1
   - Restriction digest with all 6 enzymes
   - Confirm single bands for each enzyme
   ```

---

## Applications

### Restriction Sites Enable:

**Standard Cloning:**
- AvrII, BspEI, BsrGI, BmtI, BstZ17I - all widely available

**Golden Gate Assembly:**
- BsmBI (Type IIS) - scarless peptide insertion

**Variable Region Engineering:**
- Sites distributed across VP1 for targeting different VR regions
- Positions: aa 53, 162, 402, 473, 525, 600

**Diagnostic Mapping:**
- Quick QC after cloning
- Verify plasmid identity

---

## Technical Notes

### Why AvrII Over XbaI?

| Criterion | XbaI | AvrII |
|-----------|------|-------|
| Mutations needed | 2 | 1 ✅ |
| Enzyme availability | Very common | Very common |
| Recognition site | TCTAGA | CCTAGG |
| Conflict with cloning | YES ❌ | NO ✅ |

**Decision:** AvrII is superior - fewer mutations, no cloning conflicts

### ITR Preservation

The cloning strategy preserves ITR annotations because:
1. ITRs are outside the HindIII-XbaI cloning region
2. Left ITR: 27-154 (before HindIII at 187)
3. Right ITR: 3736-3863 (after inserted cassette)
4. Features after insertion are shifted by +1424 bp

---

## Session Timeline

1. **User request:** Reassemble with ITRs, use cloning strategy, only ssAAV
2. **Identified XbaI conflict** with backbone cloning site
3. **Searched for AvrII alternative** in Region 1
4. **Found AvrII site** requiring only 1 silent mutation (position 2523)
5. **Generated v04 RepCap** with AvrII instead of XbaI
6. **Assembled v04 transfer plasmid** using cloning strategy
7. **Verified all sites unique** and ITRs annotated

---

## Next Steps (If Needed)

1. **Order synthesis** of EF1a-VP1-polyA cassette with v04 mutations
2. **Clone into backbone** using HindIII/XbaI
3. **Sequence verify** all mutations present and no errors
4. **Functional testing** of VP1 capsid assembly
5. **Use BsmBI site** for Golden Gate insertion of VHH nanobodies

---

## Repository State

```
test_data/
├── AAV9-RepCap-NOCAP-v04.gb              # v04 RepCap
├── pGS-ssAAV-EF1A-VP1-rBG_v04.gb         # v04 transfer plasmid
├── pGS-ssAAV-ITR128-Amp-empty_v02.gb     # Backbone with ITRs
└── [v03 files preserved for comparison]

projects/aav_transfer_plasmid/
├── analysis/
│   ├── find_avrii_region1.py
│   ├── generate_v04_repcap.py
│   └── assemble_v04_cloning.py
├── output/
│   └── SESSION_REPORT_v04.md              # This file
└── [previous analysis files preserved]
```

---

**Session Status:** ✅ COMPLETE
**Final Design:** v4.0
**Transfer Plasmid:** Generated with cloning strategy
**ITRs:** Annotated
**All Sites:** Unique and verified

**Ready for synthesis and cloning.** 🎯
