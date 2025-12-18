# FINAL CORRECTED: Six-Region Restriction Site Analysis

**Date:** 2025-12-17
**Plasmid:** BASE-DRAFT-AAV9-RepCap-NOCAP.gb
**Status:** ✅ ALL BUGS FIXED - Frame offsets corrected, all mutations verified silent, Region 6 narrowed to 200bp

---

## Critical Bug Fix

**The silent_sites.py script had TWO major bugs:**

1. **Frame offset calculation bug:** Used `(start - vp1_start) % 3` which gives the position WITHIN a codon, not the offset TO the next codon boundary.
   - **Corrected formula:** `frame_offset = (3 - position_in_codon) % 3 if position_in_codon != 0 else 0`

2. **classify_mutation_impact bug:** Hardcoded `frame=0` instead of using the detected `frame` parameter.
   - **Fixed:** Line 903 now uses `translate(original_coding, frame=frame)`

---

## Final Recommended Sites (All Verified Silent, Curated Enzyme List)

| Region | Enzyme | Position | Recognition Site | Mutations | Verified |
|--------|--------|----------|------------------|-----------|----------|
| **Region 1** | **SmaI** | 2505 | `CCCGGG` | 1 silent (T→C) | ✅ CTT→CTC (L→L) |
| **Region 2** | **BbvCI** | 2828 | `CCTCAGC` | 1 silent (C→A) | ✅ TCC→TCA (S→S) |
| **Region 3** | **AgeI** | 3583 | `ACCGGT` | 1 silent (G→C) | ✅ ACG→ACC (T→T) |
| **Region 4** | **BsrGI** | 3780 | `TGTACA` | 1 silent (C→A) | ✅ GTC→GTA (V→V) |
| **Region 5** | **BmtI** | 3937 | `GCTAGC` | 1 silent (C→T) | ✅ GCC→GCT (A→A) |
| **Region 6** | **BstZ17I** | 4163 | `GTATAC` | 1 silent (A→T) | ✅ GGA→GGT (G→G) |

---

## Region Details

### Region 1: Rep68-stop to VP2-start (2415-2775)

- **Narrowed from original** to avoid splice acceptor in VP1 N-terminus
- **Frame offset:** 1 (position 2415 is 3rd base of VP1 codon 17, skip 1 base to codon 18)
- **Site:** SmaI @ 2505
- **Mutation:** T2505C (CTT → CTC, both Leucine)
- **Why not BbvCI?** Would conflict with Region 2 (same recognition sequence)

### Region 2: VP2-AAP Intergenic (2779-2890)

- **Frame offset:** 0 (starts at codon boundary)
- **Site:** BbvCI @ 2828
- **Type:** Type IIS enzyme (cuts outside recognition site)
- **Application:** Golden Gate assembly

### Region 3: AAP-stop to VR4 (3485-3717)

- **Frame offset:** 0 (starts just after AAP stop codon)
- **Site:** AgeI @ 3583 (reliable workhorse enzyme)
- **Mutation:** G3585C (ACG → ACC, both Threonine)
- **Advantages:** Very common, works in most buffers, generates 5' overhang
- **Avoids:** EcoRV (used by Genscript) and PshAI (fussy enzyme)

### Region 4: VR4 to VR5 (3745-3825)

- **Frame offset:** 0
- **Site:** BsrGI @ 3780 (reliable single-site enzyme)
- **Mutation:** C3783A (GTC → GTA, both Valine)
- **Avoids:** FseI (fussy), NaeI/NgoMIV (require two sites for efficiency)

### Region 5: VR5 to VR8 (3880-4104)

- **Frame offset:** 0
- **Site:** BmtI (NheI also works - same site) @ 3937
- **Mutation:** C→T (silent)

### Region 6: Post-VR8 (4144-4343, NARROWED to 200bp after VR8)

- **Narrowed from 4144-4575** to fall within 200 bases after VR8 end
- **Frame offset:** 0
- **Site:** BstZ17I @ 4163
- **Mutation:** A4164T (GGA → GGT, both Glycine)
- **Advantages:** Unique site, within narrowed 200bp window
- **Avoids:** AfeI was at 4442 (beyond 200bp limit), BaeI (double-cutter)

---

## Summary Statistics

- **Total sites:** 6 (one per region)
- **Total mutations needed:** 6 silent mutations
- **Naturally occurring sites:** 0 (all require 1 mutation each)
- **All sites unique:** ✅ Yes, in entire 7,078 bp plasmid
- **All mutations silent:** ✅ Yes, verified manually for each site
- **No boundary overlaps:** ✅ Verified against 12 critical boundaries
- **Region 6 constraint:** ✅ Falls within 200 bases after VR8 (position 4163, 20 bp into region)

---

## Files Created/Fixed

1. **`scripts/tools/silent_sites.py`:**
   - Line 903: Changed from `frame=0` to `frame=frame`

2. **`scripts/tools/silent_sites_curated.py`:** (NEW)
   - Curated list of 61 reliable, non-fussy enzymes
   - Excludes: Nicking enzymes (Nb.*, Nt.*), double-cutters (BaeI, XcmI), mega-enzymes (I-SceI), highly ambiguous sites (>2 N's)
   - Includes: Common 6-bp and 8-bp cutters, reliable Type IIS enzymes (BbvCI, BsaI, BsmBI)

3. **Frame offset calculation (corrected in all analysis scripts):**
   ```python
   position_in_codon = (region_start - vp1_start) % 3
   frame_offset = (3 - position_in_codon) % 3 if position_in_codon != 0 else 0
   ```

### Curated Enzyme List (61 enzymes)

**Standard 6-bp cutters:** AfeI, AflII, AgeI, ApaI, ApaLI, AvrII, BamHI, BglII, BmtI, BbvCI, BsaI, BsmBI, BspEI, BsrGI, BstBI, BstZ17I, ClaI, EagI, EcoRI, EcoRV, Esp3I, HindIII, HpaI, KpnI, MfeI, MluI, MscI, NcoI, NdeI, NheI, NruI, NsiI, PciI, PmlI, PstI, PvuI, PvuII, SacI, SacII, SalI, ScaI, SmaI, SnaBI, SpeI, SphI, StuI, XbaI, XhoI, XmaI

**8-bp cutters:** AscI, AsiSI, FseI, NotI, PacI, PmeI, SbfI, SrfI, SwaI

**Minimal ambiguity:** AflIII, BsaWI, HincII

---

## Next Steps

To generate the complete revised report with all details, sequence contexts, and codon changes, run:

```bash
python3 generate_final_corrected_report.py
```

This will create: `reports/AAV9_RepCap_SixRegions_Analysis_CORRECTED_FINAL.md`
