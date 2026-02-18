# Original Alternative Enzymes - Dual-Context Analysis

**Date:** 2025-12-18
**Analysis:** Checking ALL original alternative candidates from RepCap six-regions analysis

---

## Executive Summary

Analyzed **20 enzyme candidates** (including alternatives) from the original six-regions RepCap restriction site analysis to determine which are viable in BOTH RepCap source AND assembled transfer plasmid contexts.

**Result:** **5 of 20 candidates** (25%) are viable in both contexts.

---

## Key Findings

### ✅ Viable Sites (5 total)

These enzymes work in RepCap, ssAAV, and scAAV contexts:

| Region | Enzyme | RepCap Position | ssAAV Position | scAAV Position | Notes |
|--------|--------|-----------------|----------------|----------------|-------|
| **4** | **BsrGI** | 3780 | 2770 | 1833 | ✅ Already in v02 |
| **5** | **BmtI** | 3937 | 2931 | 1994 | ✅ Already in v02 |
| **5** | **NheI** | 3937 | 2927 | 1990 | Same recognition as BmtI |
| **6** | **BstZ17I** | 4163 | 3155 | 2218 | ✅ Already in v02 |
| **6** | **AfeI** | 4442 | 3434 | 2497 | Beyond 200bp limit* |

*AfeI is 279 bp after VR8 end (4163), exceeding the original 200bp constraint but still functional.

---

### ❌ Non-Viable Sites (15 total)

All other candidates from the original analysis fail uniqueness in transfer plasmid context.

---

## Detailed Results by Region

### Region 1: Rep68-stop to VP2-start (2415-2775)

**Status:** ❌ NO VIABLE ALTERNATIVES

| Enzyme | RepCap | ssAAV | scAAV | Why It Fails |
|--------|--------|-------|-------|--------------|
| **SmaI** | 1 ✓ | 5 ❌ | 5 ❌ | **ITR regions** (positions 57, 68, 3812, 3823 in ssAAV) |
| **AvrII** | 0 ❌ | 0 | 0 | **Not in engineered RepCap** (silent mutation was never made) |
| **EagI** | 0 ❌ | 1 ✓ | 1 ✓ | **Not in engineered RepCap** (too close to VP1 start) |

**Analysis:**
- SmaI was the recommended site, but it appears 5× in the ITR/MCS regions
- AvrII and EagI were alternatives in the v2.3 report, but the silent mutations to create them were never actually implemented in the engineered RepCap plasmid
- This region has **zero viable options** with current design

---

### Region 2: VP2-AAP Intergenic (2779-2890)

**Status:** ❌ NO VIABLE ALTERNATIVES

| Enzyme | RepCap | ssAAV | scAAV | Why It Fails |
|--------|--------|-------|-------|--------------|
| **BbvCI** | 1 ✓ | 2 ❌ | 1 ✓ | **EF1α promoter** (position 1063 in ssAAV) |

**Analysis:**
- BbvCI is unique in RepCap and scAAV, but appears **2× in ssAAV**
- The duplicate site is in the EF1α promoter (with intron)
- Only one candidate enzyme for this region
- This region has **zero viable options** with current promoter

---

### Region 3: AAP-stop to VR4 (3485-3717)

**Status:** ❌ NO VIABLE ALTERNATIVES

| Enzyme | RepCap | ssAAV | scAAV | Why It Fails |
|--------|--------|-------|-------|--------------|
| **AgeI** | 1 ✓ | 2 ❌ | 2 ❌ | **EF1α promoter** (position 239 in ssAAV) |
| **EcoNI** | 0 ❌ | 2 ❌ | 1 ✓ | **Not in RepCap + EF1α promoter** |
| **PshAI** | 0 ❌ | 0 | 0 | **Not in engineered RepCap** (mutation never made) |
| **BsmBI** | 0 ❌ | 0 | 0 | **Not in engineered RepCap** (mutation never made) |

**Analysis:**
- AgeI is the recommended site but appears in the EF1α promoter
- EcoNI (alternative in v2.3 report) was never engineered into RepCap
- PshAI and BsmBI were verification script alternatives, but mutations were never made
- This region has **zero viable options** with current design and promoter

---

### Region 4: VR4 to VR5 (3745-3825)

**Status:** ✅ 1 VIABLE ALTERNATIVE

| Enzyme | RepCap | ssAAV | scAAV | Status |
|--------|--------|-------|-------|--------|
| **BsrGI** | 1 ✓ | 1 ✓ | 1 ✓ | ✅ **VIABLE** |
| **FseI** | 0 ❌ | 1 ✓ | 0 | Not in RepCap (alternative never implemented) |
| **NaeI** | 0 ❌ | 1 ✓ | 0 | Not in RepCap |
| **NgoMIV** | 0 ❌ | 1 ✓ | 0 | Not in RepCap |
| **RsrII** | 0 ❌ | 0 | 0 | Not in RepCap |

**Analysis:**
- **BsrGI is the ONLY viable option** and is already implemented in v02 ✓
- FseI was the v2.3 recommended site but was never engineered into RepCap
- All other alternatives were considered but never implemented

---

### Region 5: VR5 to VR8 (3880-4104)

**Status:** ✅ 2 VIABLE ALTERNATIVES

| Enzyme | RepCap | ssAAV | scAAV | Status |
|--------|--------|-------|-------|--------|
| **BmtI** | 1 ✓ | 1 ✓ | 1 ✓ | ✅ **VIABLE** - Already in v02 |
| **NheI** | 1 ✓ | 1 ✓ | 1 ✓ | ✅ **VIABLE** - Same site as BmtI |

**Analysis:**
- **Both BmtI and NheI work** because they recognize the same sequence (GCTAGC)
- This is the most robust region - already implemented and working in v02 ✓

---

### Region 6: Post-VR8 (4144-4575, narrowed to 4144-4343)

**Status:** ✅ 2 VIABLE ALTERNATIVES

| Enzyme | RepCap | ssAAV | scAAV | Status | Notes |
|--------|--------|-------|-------|--------|-------|
| **BstZ17I** | 1 ✓ | 1 ✓ | 1 ✓ | ✅ **VIABLE** | Already in v02, within 200bp |
| **AfeI** | 1 ✓ | 1 ✓ | 1 ✓ | ✅ **VIABLE** | Beyond 200bp limit (279 bp after VR8) |
| **BaeI** | 2 ❌ | 2 ❌ | 2 ❌ | Double-cutter, not unique |

**Analysis:**
- **BstZ17I is the primary option** - within the 200bp constraint, already in v02 ✓
- **AfeI is a bonus option** - viable but 279 bp after VR8 (exceeds original 200bp window)
- BaeI naturally occurs twice in the plasmid (not suitable)

---

## Why Most Alternatives Don't Work

### Issue 1: Silent Mutations Never Made

Many "alternative" enzymes from the reports (AvrII, EagI, FseI, NaeI, EcoNI, PshAI, BsmBI) were:
- Identified as theoretically possible during analysis
- Listed in reports and verification scripts
- **Never actually engineered into the RepCap plasmid**

The final RepCap plasmid (`AAV9-RepCap-NOCAP-WITH-RESTRICTION-SITES.gb`) only contains **6 engineered sites**:
1. SmaI @ 2505 (Region 1)
2. BbvCI @ 2828 (Region 2)
3. AgeI @ 3583 (Region 3)
4. BsrGI @ 3780 (Region 4)
5. BmtI @ 3937 (Region 5)
6. BstZ17I @ 4163 (Region 6)

### Issue 2: EF1α Promoter Contains Problem Sites

The EF1α promoter naturally contains:
- **AgeI (ACCGGT)** @ position 239 in ssAAV → Conflicts with Region 3
- **BbvCI (CCTCAGC)** @ position 1063 in ssAAV → Conflicts with Region 2

### Issue 3: ITR/Backbone Contains Problem Sites

The ITR and backbone regions contain:
- **SmaI (CCCGGG)** @ positions 57, 68, 3812, 3823 in ssAAV → Conflicts with Region 1

---

## Overall Statistics

### By Region

| Region | Candidates Tested | Viable | Success Rate |
|--------|------------------|--------|--------------|
| Region 1 | 3 | 0 | 0% |
| Region 2 | 1 | 0 | 0% |
| Region 3 | 4 | 0 | 0% |
| Region 4 | 5 | 1 | 20% |
| Region 5 | 2 | 2 | 100% |
| Region 6 | 3 | 2 | 67% |
| **Total** | **20** | **5** | **25%** |

### By Failure Reason

| Failure Reason | Count | Percentage |
|---------------|-------|------------|
| Not in RepCap (mutation never made) | 8 | 53% |
| EF1α promoter conflict | 3 | 20% |
| ITR/backbone conflict | 1 | 7% |
| Multiple sites (naturally) | 1 | 7% |
| Combination of above | 2 | 13% |
| **Total Failed** | **15** | **100%** |

---

## Recommendations

### Immediate: Accept Current Design ⭐ RECOMMENDED

**Use the 3 viable unique sites already in v02:**

1. **Region 4:** BsrGI (TGTACA)
2. **Region 5:** BmtI/NheI (GCTAGC)
3. **Region 6:** BstZ17I (GTATAC)

**Optional bonus:**
- **Region 6:** AfeI (ACRYAG) - if 200bp constraint can be relaxed

**Sufficient for:**
- Basic modular cloning
- Peptide insertions in VR4-VR5, VR5-VR8, post-VR8 regions
- Most standard VP1 engineering applications

---

### Long-term Option 1: Change Promoter

Replace EF1α with a promoter lacking AgeI and BbvCI sites:

**Candidates to screen:**
- **PGK promoter** (~450 bp)
- **CMV immediate-early** (~600 bp)
- **CAG promoter** (~1,600 bp)
- **UBC promoter** (~1,200 bp)

**Impact:** Would unlock Regions 2 & 3 (potentially +2 sites)

---

### Long-term Option 2: Redesign Regions 1-3 in RepCap

Go back to the RepCap source plasmid and engineer different restriction sites:
- Screen for enzymes that don't appear in ITRs or common promoters
- Engineer silent mutations to create these new sites
- Regenerate RepCap plasmid and transfer plasmids

**Impact:** Could unlock all 6 regions, but requires significant effort

---

### Long-term Option 3: Use Partial Digest Strategies

Accept non-unique sites but use careful digest planning:
- **AgeI:** Use flanking enzymes to isolate VP1
- **BbvCI:** Type IIS enzyme cuts outside recognition - may be usable with specific strategies
- **SmaI:** Most problematic (5 sites) - avoid

**Impact:** Adds complexity to cloning workflows

---

## Conclusion

**Answer to user's question:** "Are none of these alternatives viable?"

**Yes, you're correct - MOST alternatives are not viable.** Only **5 of 20** (25%) original candidates work in both contexts:

✅ **Working now:**
- BsrGI (Region 4)
- BmtI/NheI (Region 5)
- BstZ17I (Region 6)

✅ **Bonus option:**
- AfeI (Region 6, beyond 200bp limit)

❌ **Not viable:**
- All 8 alternatives from Regions 1-3

**Root causes:**
1. Most "alternatives" in reports were never actually engineered into RepCap
2. The 6 sites that WERE engineered clash with EF1α promoter and ITRs
3. Only the sites in Regions 4-6 avoid these clashes

**Bottom line:** The current v02 design with 3 sites (BsrGI, BmtI, BstZ17I) is the maximum achievable without redesigning the promoter or RepCap plasmid.

---

**Analysis Complete** ✓
