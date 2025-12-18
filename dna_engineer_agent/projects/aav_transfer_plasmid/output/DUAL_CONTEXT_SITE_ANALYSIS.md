# Dual-Context Restriction Site Analysis

**Date:** 2025-12-18
**Analysis:** Sites unique in BOTH RepCap and Transfer Plasmid contexts

---

## Executive Summary

Analyzed all 6 VP1 regions for restriction sites that are unique in:
1. ✅ RepCap source plasmid (AAV9-RepCap-NOCAP-WITH-RESTRICTION-SITES.gb)
2. ✅ ssAAV transfer plasmid (pGS-ssAAV-EF1A-VP1-rBG_v02.gb)
3. ✅ scAAV transfer plasmid (pGS-scAAV-EF1A-VP1-rBG_v02.gb)

**Result:** 3 of 6 regions have viable unique sites.

---

## Results by Region

### ✅ Region 4: VR4 to VR5 - VIABLE

**Site:** BsrGI (TGTACA)
- VP1 relative position: 1415
- RepCap position: 3782
- ssAAV position: 2770
- scAAV position: 1833
- **Status:** ✅ UNIQUE in all three contexts

---

### ✅ Region 5: VR5 to VR8 - VIABLE

**Sites:** BmtI/NheI (GCTAGC) - Same recognition sequence
- VP1 relative position: 1572
- RepCap position: 3939-3943
- ssAAV position: 2927-2931
- scAAV position: 1990-1994
- **Status:** ✅ UNIQUE in all three contexts

---

### ✅ Region 6: Post-VR8 - VIABLE

**Site:** BstZ17I (GTATAC)
- VP1 relative position: 1798
- RepCap position: 4167
- ssAAV position: 3155
- scAAV position: 2218
- **Status:** ✅ UNIQUE in all three contexts

---

## 🚩 Flagged Regions Requiring Manual Review

### ❌ Region 1: Rep68-stop to VP2-start (VP1 pos 0-280)

**Original Site:** SmaI (CCCGGG) @ position 141

**Problem:** All candidate enzymes fail uniqueness test

| Enzyme | RepCap Count | ssAAV Count | scAAV Count | Issue |
|--------|--------------|-------------|-------------|-------|
| **SmaI** | 1 | 5 | 5 | **ITR/backbone duplicates** |
| EcoRI | 3 | 1 | 1 | Multiple in RepCap |
| NruI | 3 | 1 | 1 | Multiple in RepCap |
| StuI | 3 | 4 | 2 | Multiple everywhere |
| XhoI | 3 | 4 | 3 | Multiple everywhere |

**Recommendation:**
- **Option A:** Skip this region, rely on Regions 4-6
- **Option B:** Use different promoter that doesn't contain these sites
- **Option C:** Accept SmaI with caveat that it appears in ITRs/backbone (requires flanking digests)

---

### ❌ Region 2: VP2-AAP Intergenic (VP1 pos 400-550)

**Original Site:** BbvCI (CCTCAGC) @ position 464

**Problem:** All candidate enzymes fail uniqueness test

| Enzyme | RepCap Count | ssAAV Count | scAAV Count | Issue |
|--------|--------------|-------------|-------------|-------|
| **BbvCI** | 1 | 2 | 1 | **EF1α promoter duplicate** |
| ApaLI | 4 | 4 | 4 | Multiple everywhere |
| StuI | 3 | 4 | 2 | Multiple everywhere |

**Recommendation:**
- **Option A:** Skip this region, rely on Regions 4-6
- **Option B:** Use BbvCI with caveat (appears in promoter, requires careful digest strategy)
- **Option C:** Use different promoter (CMV, PGK, CAG) without BbvCI site

---

### ❌ Region 3: AAP-stop to VR4 (VP1 pos 1150-1300)

**Original Site:** AgeI (ACCGGT) @ position 1219

**Problem:** Only one candidate enzyme, fails uniqueness

| Enzyme | RepCap Count | ssAAV Count | scAAV Count | Issue |
|--------|--------------|-------------|-------------|-------|
| **AgeI** | 1 | 2 | 2 | **EF1α promoter duplicate** |

**Recommendation:**
- **Option A:** Skip this region, rely on Regions 4-6
- **Option B:** Use AgeI with caveat (appears in promoter, requires flanking digests)
- **Option C:** Use different promoter without AgeI site
- **Option D:** Engineer a different site in this region (would require going back to RepCap plasmid)

---

## Overall Recommendations

### Immediate Solution: Use the 3 Viable Sites

**Working Sites:**
1. Region 4: **BsrGI** ✅
2. Region 5: **BmtI/NheI** ✅
3. Region 6: **BstZ17I** ✅

This gives you **3 unique restriction sites** across VP1 that work in all contexts:
- Sufficient for basic modular cloning
- All sites well-distributed across VP1 (positions 1415, 1572, 1798)
- Can be used for inserting peptides in VR4-VR5, VR5-VR8, and post-VR8 regions

---

### Long-term Solutions for Regions 1-3

#### Option 1: Use Alternative Promoter ⭐ **RECOMMENDED**

**Replace EF1α with a promoter that doesn't contain AgeI, BbvCI, or common sites.**

**Candidates to test:**
- **PGK promoter** (~450 bp) - No AgeI or BbvCI
- **CMV immediate-early** (~600 bp) - Check for site conflicts
- **CAG promoter** (~1,600 bp) - Hybrid CMV enhancer + chicken β-actin
- **UBC promoter** (~1,200 bp) - Ubiquitin C promoter

**Action:** Screen candidate promoters for absence of problem sites.

---

#### Option 2: Accept Limitations, Use Partial Digests

Keep current design but use strategic digest planning:
- **AgeI:** Use with flanking enzymes to isolate VP1 fragment
- **BbvCI:** Type IIS enzyme, cuts outside recognition - may still be usable
- **SmaI:** Most problematic (5 sites in ssAAV) - avoid or use very carefully

---

#### Option 3: Redesign Regions 1-3 in RepCap Source

Go back to the RepCap plasmid and engineer different sites:
- Screen EF1α promoter for absent enzymes
- Choose sites that don't appear in ITRs, promoter, or backbone
- Update regions 1-3 with new silent mutations

**Effort:** High - requires new RepCap plasmid generation

---

## Detailed Site Positions

### RepCap Plasmid (7,078 bp)
- VP1 CDS: 2365-4575

### ssAAV Transfer (6,421 bp)
- EF1α promoter: 159-1352 (1194 bp)
- VP1 CDS: 1359-3572

### scAAV Transfer (5,447 bp)
- EF1α core: 159-400 (242 bp)
- VP1 CDS: ~422-2635

---

## Conflict Details

### EF1α Promoter Contains:
- **AgeI (ACCGGT)** @ position 239 in ssAAV
- **BbvCI (CCTCAGC)** @ position 1063 in ssAAV
- Result: Regions 2 & 3 sites are not unique

### ITR/Backbone Contains:
- **SmaI (CCCGGG)** @ positions 57, 68, 3812, 3823 in ssAAV
- Result: Region 1 SmaI is not unique

---

## Decision Matrix

| Region | Current Site | Status | Keep? | Alternative Action |
|--------|--------------|--------|-------|-------------------|
| 1 | SmaI | ❌ Not unique | No | Skip or redesign |
| 2 | BbvCI | ❌ Not unique | No | Skip or redesign |
| 3 | AgeI | ❌ Not unique | No | Skip or redesign |
| 4 | BsrGI | ✅ Unique | **YES** | Use as-is |
| 5 | BmtI/NheI | ✅ Unique | **YES** | Use as-is |
| 6 | BstZ17I | ✅ Unique | **YES** | Use as-is |

---

## Next Steps

1. **Decide on strategy:**
   - Accept 3 viable sites (quickest)
   - Switch to alternative promoter (best long-term)
   - Redesign regions 1-3 (most effort)

2. **If switching promoter:**
   - Screen PGK, CMV, CAG, UBC for site conflicts
   - Update COMPONENTS in assembly script
   - Regenerate transfer plasmids as v03

3. **If keeping current design:**
   - Document limitations in lab protocols
   - Update cloning strategies to work around non-unique sites
   - Proceed with BsrGI, BmtI, BstZ17I only

---

## ADDENDUM: Original Alternatives Analysis (2025-12-18)

**Question:** Were there other enzyme options from the original six-regions analysis?

**Answer:** Yes, 20 total candidates were tested (including alternatives).

**Result:** Only **5 of 20** (25%) work in both RepCap and transfer plasmid contexts.

### Why Most Alternatives Failed

1. **Never Engineered (8 candidates):** Alternative enzymes like AvrII, EagI, FseI, NaeI, EcoNI, PshAI, BsmBI, RsrII were identified as theoretically possible but the silent mutations were never actually made in the RepCap plasmid.

2. **Promoter Conflicts (3 candidates):** AgeI and BbvCI clash with natural sites in the EF1α promoter.

3. **ITR Conflicts (1 candidate):** SmaI appears 5× in ITR/backbone regions.

### Viable Alternatives Found

| Region | Viable Enzymes | Status |
|--------|---------------|--------|
| Region 1 | None (0/3) | All fail |
| Region 2 | None (0/1) | BbvCI in promoter |
| Region 3 | None (0/4) | All fail |
| Region 4 | BsrGI (1/5) | ✅ Already in v02 |
| Region 5 | BmtI, NheI (2/2) | ✅ Already in v02 |
| Region 6 | BstZ17I, AfeI (2/3) | ✅ BstZ17I in v02, AfeI beyond 200bp |

**Conclusion:** The current v02 design with 3 sites (BsrGI, BmtI/NheI, BstZ17I) represents the maximum achievable with the existing RepCap design and EF1α promoter.

**See:** `ORIGINAL_ALTERNATIVES_ANALYSIS.md` for complete details.

---

**Analysis Complete** ✓
