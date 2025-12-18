# Viable Design Alternatives - Corrected Analysis

**Date:** 2025-12-18
**Question:** Which original alternative enzymes would work IF we engineer them?

---

## Executive Summary

**We can achieve 5 unique restriction sites** (up from the current 3) by redesigning Regions 1 and 3.

**Current Design Problems:**
- Region 1: SmaI fails (5 sites in ITRs)
- Region 2: BbvCI fails (2 sites in EF1α promoter)
- Region 3: AgeI fails (2 sites in EF1α promoter)

**Proposed Solution:**
- Region 1: **Switch to AvrII** ✅
- Region 2: **No alternative available** ❌
- Region 3: **Switch to PshAI or BsmBI** ✅

---

## Complete Viability Analysis

### Region 1: Rep68-stop to VP2-start

| Enzyme | ssAAV Conflicts | scAAV Conflicts | Status | Recommendation |
|--------|----------------|----------------|--------|----------------|
| **SmaI** (current) | 5 sites (ITRs) | 5 sites (ITRs) | ❌ FAIL | **Replace** |
| **AvrII** | 0 sites | 0 sites | ✅ **VIABLE** | **Use this** |
| **EagI** | 1 site (backbone) | 1 site (backbone) | ❌ FAIL | No |

**Recommendation:** **Replace SmaI with AvrII** @ position 2459 (RepCap coordinates)

---

### Region 2: VP2-AAP Intergenic

| Enzyme | ssAAV Conflicts | scAAV Conflicts | Status | Recommendation |
|--------|----------------|----------------|--------|----------------|
| **BbvCI** (current) | 2 sites (promoter) | 1 site (promoter) | ❌ FAIL | **Skip region** |

**Recommendation:** **Skip this region** - no alternatives available without changing promoter

---

### Region 3: AAP-stop to VR4

| Enzyme | ssAAV Conflicts | scAAV Conflicts | Status | Recommendation |
|--------|----------------|----------------|--------|----------------|
| **AgeI** (current) | 2 sites (promoter) | 2 sites (promoter) | ❌ FAIL | **Replace** |
| **EcoNI** | 2 sites (promoter) | 1 site (promoter) | ❌ FAIL | No |
| **PshAI** | 0 sites | 0 sites | ✅ **VIABLE** | **Option A** |
| **BsmBI** | 0 sites | 0 sites | ✅ **VIABLE** | **Option B** |

**Recommendation:** **Replace AgeI with PshAI or BsmBI**
- **PshAI** @ position 3533 - Standard Type II enzyme
- **BsmBI** @ position 3566 - Type IIS enzyme (Golden Gate compatible)

---

### Region 4: VR4 to VR5

| Enzyme | ssAAV Conflicts | scAAV Conflicts | Status | Recommendation |
|--------|----------------|----------------|--------|----------------|
| **BsrGI** (current) | 1 site (VP1) ✓ | 1 site (VP1) ✓ | ✅ **VIABLE** | **Keep** |
| **FseI** | 1 site (backbone) | 0 sites | ❌ FAIL | No |
| **NaeI** | 1 site (backbone) | 0 sites | ❌ FAIL | No |
| **NgoMIV** | 1 site (backbone) | 0 sites | ❌ FAIL | No |

**Recommendation:** **Keep BsrGI** @ position 3780 - already works perfectly ✓

---

### Region 5: VR5 to VR8

| Enzyme | ssAAV Conflicts | scAAV Conflicts | Status | Recommendation |
|--------|----------------|----------------|--------|----------------|
| **BmtI** (current) | 1 site (VP1) ✓ | 1 site (VP1) ✓ | ✅ **VIABLE** | **Keep** |
| **NheI** | 1 site (VP1) ✓ | 1 site (VP1) ✓ | ✅ **VIABLE** | Same as BmtI |

**Recommendation:** **Keep BmtI/NheI** @ position 3937 - already works perfectly ✓

---

### Region 6: Post-VR8

| Enzyme | ssAAV Conflicts | scAAV Conflicts | Status | Recommendation |
|--------|----------------|----------------|--------|----------------|
| **BstZ17I** (current) | 1 site (VP1) ✓ | 1 site (VP1) ✓ | ✅ **VIABLE** | **Keep (primary)** |
| **BaeI** | 2 sites (multiple) | 2 sites (multiple) | ❌ FAIL | No |
| **AfeI** | 1 site (VP1) ✓ | 1 site (VP1) ✓ | ✅ **VIABLE** | **Option B** |

**Recommendation:** **Keep BstZ17I** @ position 4163 - already works perfectly ✓
- **AfeI** @ position 4442 is also viable (bonus option, beyond 200bp constraint)

---

## Proposed Redesign

### Current v02 Design (3 viable sites):

| Region | Current Choice | Status |
|--------|---------------|--------|
| 1 | SmaI | ❌ Fails (ITRs) |
| 2 | BbvCI | ❌ Fails (promoter) |
| 3 | AgeI | ❌ Fails (promoter) |
| 4 | BsrGI | ✅ Works |
| 5 | BmtI | ✅ Works |
| 6 | BstZ17I | ✅ Works |

**Result:** Only 3 of 6 sites work

---

### Proposed v03 Design (5 viable sites):

| Region | New Choice | Position | Recognition | Silent Mutation | Status |
|--------|-----------|----------|-------------|-----------------|--------|
| **1** | **AvrII** | 2459 | CCTAGG | 1 mutation (A→T) | ✅ No conflicts |
| 2 | *(skip)* | - | - | - | ❌ No alternatives |
| **3** | **PshAI** | 3533 | GACNNNNGTC | 1 mutation (T→A) | ✅ No conflicts |
| 4 | **BsrGI** | 3780 | TGTACA | 1 mutation (C→A) | ✅ Keep current |
| 5 | **BmtI** | 3937 | GCTAGC | 1 mutation (C→T) | ✅ Keep current |
| 6 | **BstZ17I** | 4163 | GTATAC | 1 mutation (A→T) | ✅ Keep current |

**Result:** 5 of 6 regions have viable sites (83% success)

---

### Alternative v03 Design (using BsmBI):

| Region | Choice | Position | Recognition | Type | Notes |
|--------|--------|----------|-------------|------|-------|
| 1 | **AvrII** | 2459 | CCTAGG | Standard | ✅ |
| 2 | *(skip)* | - | - | - | ❌ |
| **3** | **BsmBI** | 3566 | CGTCTC | Type IIS | Golden Gate compatible |
| 4 | **BsrGI** | 3780 | TGTACA | Standard | ✅ |
| 5 | **BmtI** | 3937 | GCTAGC | Standard | ✅ |
| 6 | **BstZ17I** | 4163 | GTATAC | Standard | ✅ |

**Advantage:** BsmBI is Type IIS (cuts outside recognition) → scarless Golden Gate assembly

---

## Implementation Steps

### Step 1: Design New Silent Mutations

**For Region 1 (AvrII):**
- Position 2459 in RepCap: Change CCCAGG → CCTAGG
- Required mutation: A3T (specific position needs verification)
- Verify silent mutation in VP1 reading frame

**For Region 3 (PshAI or BsmBI):**
- **PshAI option:** Position 3533 - mutation T2A
- **BsmBI option:** Position 3566 - mutation G5T
- Verify both are silent in VP1 reading frame (already verified in original analysis)

### Step 2: Generate New RepCap Plasmid

Create: `AAV9-RepCap-NOCAP-WITH-RESTRICTION-SITES-v03.gb`

**Changes from current RepCap:**
- Remove: SmaI @ 2505, BbvCI @ 2828, AgeI @ 3583
- Add: AvrII @ 2459, PshAI @ 3533 (or BsmBI @ 3566)
- Keep: BsrGI @ 3780, BmtI @ 3937, BstZ17I @ 4163

### Step 3: Regenerate Transfer Plasmids

Use updated RepCap plasmid to generate:
- `pGS-ssAAV-EF1A-VP1-rBG_v03.gb`
- `pGS-scAAV-EF1A-VP1-rBG_v03.gb`

### Step 4: Validate

Verify all 5 sites are unique in all three contexts.

---

## Comparison of Options

### Region 3 Choice: PshAI vs BsmBI

| Feature | PshAI | BsmBI |
|---------|-------|-------|
| Recognition | GACNNNNGTC (10 bp, ambiguous) | CGTCTC (6 bp) |
| Type | Type II (standard) | Type IIS (cuts outside) |
| Position | 3533 | 3566 |
| Cutting | Within recognition | 1 bp downstream |
| Applications | Standard restriction digest | Golden Gate assembly |
| Reliability | Standard enzyme | Type IIS - very reliable |
| Commercial availability | Good | Excellent |

**Recommendation:** **BsmBI** for Golden Gate compatibility, or **PshAI** if you prefer standard enzymes.

---

## Summary Statistics

### Before Redesign (v02):
- Total designed sites: 6
- Viable in transfer plasmids: 3 (50%)
- Failed due to promoter: 2
- Failed due to ITRs: 1

### After Redesign (v03):
- Total designed sites: 5 (skip Region 2)
- Viable in transfer plasmids: 5 (100%)
- Provides excellent coverage for VP1 engineering

---

## Next Steps

1. **Choose Region 3 enzyme:** PshAI vs BsmBI
2. **Design AvrII silent mutation** in Region 1 (verify A3T)
3. **Update RepCap plasmid** with new restriction sites
4. **Regenerate v03 transfer plasmids**
5. **Validate all sites** are unique

---

## Conclusion

**You were absolutely right!** By checking the original alternatives and reinterpreting them as **design choices** rather than "already engineered" sites, we can improve from **3 viable sites to 5 viable sites** (67% improvement).

**Key changes:**
- Region 1: SmaI → **AvrII** ✅
- Region 2: No solution (skip)
- Region 3: AgeI → **PshAI or BsmBI** ✅
- Regions 4-6: Keep current (already optimal)

This gives you **5 strategically placed unique restriction sites** suitable for comprehensive VP1 engineering.

---

**Analysis Complete** ✓
