# AAV Transfer Plasmid Assembly - Prompt v2.0

**Version:** 2.0
**Date:** 2025-12-18
**Status:** Implementation Ready

---

## **TASK: Draft Two AAV Transfer Plasmids with VP1 Expression Cassette (v2.0)**

---

### **1. OBJECTIVE**

Create two draft AAV transfer plasmids for VP1/VP2/VP3 protein expression by assembling a transgene cassette between the ITRs of two backbone plasmids. **Critical update:** Use the VP1 CDS from the engineered plasmid that already contains 6 silent restriction sites.

**Input Backbones:**
1. `test_data/pGS-scAAV-ITR128-Amp-empty.gb` → Output: `pGS-scAAV-EF1A-VP1-rBG_v02.gb`
2. `test_data/pGS-ssAAV-ITR128-Amp-empty.gb` → Output: `pGS-ssAAV-EF1A-VP1-rBG_v02.gb`

---

### **2. CRITICAL CHANGE: VP1 SOURCE**

**DO NOT use native AAV9 VP1.** Instead, extract the VP1 coding sequence from:

```
Source File: test_data/AAV9-RepCap-NOCAP-WITH-RESTRICTION-SITES.gb
VP1 CDS Coordinates: 2365-4575 (2211 bp, 736 aa)
```

This VP1 contains **6 pre-engineered silent restriction sites**:

| Site # | Enzyme | Position | Recognition | Mutation | Codon Change |
|--------|--------|----------|-------------|----------|--------------|
| 1 | **SmaI** | 2505 | `CCCGGG` | T2505C | CTT→CTC (L→L) |
| 2 | **BbvCI** | 2828 | `CCTCAGC` | C2832A | TCC→TCA (S→S) |
| 3 | **AgeI** | 3583 | `ACCGGT` | G3585C | ACG→ACC (T→T) |
| 4 | **BsrGI** | 3780 | `TGTACA` | C3783A | GTC→GTA (V→V) |
| 5 | **BmtI/NheI** | 3937 | `GCTAGC` | C3939T | GCC→GCT (A→A) |
| 6 | **BstZ17I** | 4163 | `GTATAC` | A4164T | GGA→GGT (G→G) |

**When extracting VP1 for the new plasmid, these sites must be preserved in their relative positions within the CDS.**

---

### **3. ADDITIONAL UNIQUE RESTRICTION SITES**

Add **two additional unique restriction sites** at junction points:

#### **Site A: Between EF1α Promoter and VP1 CDS**
- **Location:** After EF1α intron splice acceptor, before VP1 Kozak/ATG
- **Requirements:**
  - Must be unique in the entire output plasmid
  - Must NOT be present in: EF1α promoter, VP1 CDS (including the 6 engineered sites), rBG polyA, ITRs, or backbone
  - Should NOT disrupt splice acceptor (AG) at intron 3' end
  - Should be a common, reliable enzyme (not fussy)

- **Suggested candidates** (verify uniqueness in final plasmid):
  - **NheI** (`GCTAGC`) - if not already present in backbone
  - **AflII** (`CTTAAG`) - cuts well, common
  - **NcoI** (`CCATGG`) - can overlap Kozak if positioned correctly
  - **EcoRI** (`GAATTC`) - very common, reliable
  - **PacI** (`TTAATTAA`) - 8-cutter, very likely unique

- **Annotation:**
  ```
  Feature Type: misc_feature
  Label: Upstream_Cloning_Site
  Note: "Unique [ENZYME] site for promoter swapping. [Action]: Engineered | [Risk]: Low"
  ```

#### **Site B: Between VP1 Stop Codon and 3' UTR/polyA**
- **Location:** After VP1 TAA/TAG/TGA stop codon, before 3' UTR element
- **Requirements:**
  - Must be unique in entire plasmid
  - Must NOT be in VP1 CDS, 3' UTR, polyA, ITRs, or backbone
  - Provide clean separation between ORF and 3' regulatory elements

- **Suggested candidates** (verify uniqueness):
  - **XhoI** (`CTCGAG`) - very common
  - **SalI** (`GTCGAC`) - compatible with XhoI overhangs
  - **NotI** (`GCGGCCGC`) - 8-cutter, likely unique
  - **ClaI** (`ATCGAT`) - if not already present
  - **BamHI** (`GGATCC`) - very common

- **Annotation:**
  ```
  Feature Type: misc_feature
  Label: Downstream_Cloning_Site
  Note: "Unique [ENZYME] site for 3'UTR/polyA swapping. [Action]: Engineered | [Risk]: Low"
  ```

---

### **4. 3' UTR ELEMENT**

Add a **small, common 3' UTR element** between the VP1 stop codon and the rabbit β-globin polyA signal:

#### **Option A: Minimal SV40 Late 3' UTR (Recommended)**
- **Size:** ~55-65 bp
- **Source:** SV40 late region 3' UTR (NOT the polyA signal itself)
- **Sequence example:**
  ```
  TGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACC
  ```
- **Note:** This is a minimal 3' UTR that provides mRNA stability without adding a second polyA signal

#### **Option B: Synthetic Minimal 3' UTR**
- **Size:** ~30-50 bp
- **Contains:** Spacer sequence with no cryptic splice sites or internal polyA signals
- **Purpose:** Provides spacing between stop codon and polyA for efficient termination

**Annotation:**
```
Feature Type: 3'UTR
Label: mini_3UTR
Note: "Minimal 3' UTR element for mRNA stability. [Action]: Inserted | [Risk]: Low"
```

**IMPORTANT:** The 3' UTR goes BEFORE the rBG polyA signal, not after. Order is:
```
[VP1 STOP] → [Downstream Site] → [mini 3'UTR] → [rBG polyA]
```

---

### **5. COMPREHENSIVE ANNOTATION REQUIREMENTS**

#### **5.1 ITR Annotations**

The backbone plasmids use **truncated ITR128** (128 bp) instead of full-length 145 bp ITRs.

**For each ITR, annotate:**
```
Feature Type: repeat_region
Label: ITR-L (or ITR-R)
Note: "AAV2 Inverted Terminal Repeat (truncated 128bp). [Type]: ITR128"
```

**Sub-annotations within each ITR (if identifiable):**
- `misc_feature` for RBE (Rep Binding Element) - GAGC repeats
- `misc_feature` for TRS (Terminal Resolution Site) - GGTTGA motif
- `misc_feature` for D-sequence (if present)

**For scAAV backbone:** Note if one ITR has ΔTRS mutation for self-complementary packaging.

#### **5.2 Cap Gene Annotations (VP1, VP2, VP3, AAP)**

Extract and annotate ALL overlapping ORFs from the VP1 CDS:

| Feature | Type | Coordinates (relative to VP1 start=1) | Note |
|---------|------|---------------------------------------|------|
| **VP1** | CDS | 1-2211 (full length) | Primary annotation |
| **VP2** | CDS | join(412..414, 415..2211) | ACG start (Thr), non-canonical |
| **VP3** | CDS | 607-2211 | Conventional ATG start |
| **AAP** | CDS | 527-1120 | +1 frame relative to VP1 |

**VP1 annotation:**
```
Feature Type: CDS
Label: VP1
Gene: cap
Codon_start: 1
Translation: <full 736 aa sequence>
Note: "AAV9 VP1 capsid protein. Contains 6 engineered silent restriction sites."
```

**VP2 annotation (CRITICAL - uses non-canonical ACG start):**
```
Feature Type: CDS
Label: VP2
Gene: cap
Codon_start: 1
Translation: <starting from ACG-encoded Thr>
Note: "AAV9 VP2 capsid protein. Non-canonical ACG (Thr) start codon."
```

**VP3 annotation:**
```
Feature Type: CDS
Label: VP3
Gene: cap
Codon_start: 1
Translation: <starting from ATG at position 607>
Note: "AAV9 VP3 capsid protein. Most abundant capsid component (1:1:10 ratio)."
```

**AAP annotation (CRITICAL - +1 reading frame):**
```
Feature Type: CDS
Label: AAV9_AAP
Codon_start: 1
Translation: <AAP protein sequence from +1 frame>
Note: "Assembly-Activating Protein. Translated from +1 frame of VP1. Essential for capsid assembly."
```

#### **5.3 Variable Region (VR) Annotations**

Annotate all 9 variable regions using coordinates from `AAV9-RepCap-NOCAP-WITH-RESTRICTION-SITES.gb`:

| VR | DNA Coordinates* | Amino Acids | Feature Label |
|----|-----------------|-------------|---------------|
| VR-I | 784-807 | AA 262-269 | AAV9_VR-I |
| VR-II | 979-996 | AA 327-332 | AAV9_VR-II |
| VR-III | 1144-1158 | AA 382-386 | AAV9_VR-III |
| VR-IV | 1354-1380 | AA 452-460 | AAV9_VR-IV |
| VR-V | 1462-1515 | AA 488-505 | AAV9_VR-V |
| VR-VI | 1579-1617 | AA 527-539 | AAV9_VR-VI |
| VR-VII | 1633-1674 | AA 545-558 | AAV9_VR-VII |
| VR-VIII | 1741-1779 | AA 581-593 | AAV9_VR-VIII |
| VR-IX | 2110-2142 | AA 704-714 | AAV9_VR-IX |

*Coordinates are relative to VP1 start (position 1 = first bp of VP1 ATG). Adjust to absolute plasmid coordinates.

**Annotation format:**
```
Feature Type: misc_feature
Label: AAV9_VR-IV
Note: "AAV9 Variable Region VR-IV (AA 452-460). Surface-exposed loop for capsid engineering."
```

#### **5.4 Engineered Restriction Sites (Within VP1)**

Annotate each of the 6 pre-engineered sites:
```
Feature Type: misc_feature
Label: SmaI_site (or appropriate enzyme)
Note: "Engineered SmaI restriction site (CCCGGG). Silent mutation: CTT→CTC (L→L)"
```

---

### **6. COMPLETE CASSETTE STRUCTURE**

The final transgene cassette between ITRs should be (5' → 3'):

```
[ITR-L] ─┬─ [EF1α promoter core]
         │
         ├─ [EF1α Intron A]
         │
         ├─ [UPSTREAM UNIQUE SITE] ← New (e.g., NheI)
         │
         ├─ [Kozak] ─ [VP1 ATG]
         │              │
         │              ├─ VP2 ACG start (internal)
         │              │
         │              ├─ AAP start (internal, +1 frame)
         │              │
         │              ├─ VP3 ATG start (internal)
         │              │
         │              ├─ [6 engineered restriction sites]
         │              │
         │              ├─ [9 Variable Regions]
         │              │
         │              └─ VP1/VP2/VP3/AAP stops
         │
         ├─ [DOWNSTREAM UNIQUE SITE] ← New (e.g., NotI)
         │
         ├─ [mini 3' UTR]
         │
         └─ [rBG polyA]
         │
[ITR-R] ─┘
```

---

### **7. UNIQUENESS VERIFICATION**

Before finalizing, verify that these enzymes cut **exactly once** in the entire plasmid:

**Sites that MUST be unique:**
1. Upstream Cloning Site (new)
2. Downstream Cloning Site (new)
3. SmaI (from VP1)
4. BbvCI (from VP1)
5. AgeI (from VP1)
6. BsrGI (from VP1)
7. BmtI/NheI (from VP1)
8. BstZ17I (from VP1)

**Verification script must:**
```python
from Bio.Restriction import RestrictionBatch
# Count occurrences of each site in full plasmid
# ABORT if any required site appears >1 time or 0 times
```

---

### **8. SIZE CALCULATIONS**

| Component | Approximate Size |
|-----------|------------------|
| EF1α promoter + intron | ~1,172 bp |
| Upstream site + spacer | ~10-20 bp |
| VP1 ORF (with stop) | ~2,214 bp |
| Downstream site + spacer | ~10-20 bp |
| Mini 3' UTR | ~55-65 bp |
| rBG polyA | ~127 bp |
| **TOTAL CASSETTE** | **~3,600 bp** |

⚠️ **scAAV WARNING:** This exceeds the ~2.4 kb packaging limit. Flag in report.
✅ **ssAAV OK:** Within 4.7 kb limit.

---

### **9. OUTPUT REQUIREMENTS**

**Generate:**

1. **Plasmid Files:**
   - `test_data/pGS-scAAV-EF1A-VP1-rBG_v02.gb`
   - `test_data/pGS-ssAAV-EF1A-VP1-rBG_v02.gb`

2. **Report:** `reports/EF1A-VP1-rBG_Assembly_Report_v02.md`

3. **ASCII Map:**
```
================================================================
CONSTRUCT: pGS-ssAAV-EF1A-VP1-rBG_v02.gb (XXXX bp)
================================================================
[ITR128-L]                                          [ITR128-R]
|                                                            |
v                                                            v
5'==|==[EF1A+IntronA]==[NheI]==[VP1/VP2/VP3/AAP]==[NotI]==[3'UTR]==[rBGpA]==|==3'
    |                          |___ 6 engineered sites ___|                  |
    |<----------------------- ~3600 bp cassette -------------------------->|

Internal VP1 Sites: SmaI-BbvCI-AgeI-BsrGI-BmtI-BstZ17I
Variable Regions:   VR-I through VR-IX annotated
================================================================
```

4. **Restriction Map Table** showing all unique sites and their positions

---

### **10. VERIFICATION CHECKLIST**

| Check | Requirement | Action if Failed |
|-------|-------------|------------------|
| ITR-L intact | 128 bp match to reference | ABORT |
| ITR-R intact | 128 bp match to reference | ABORT |
| VP1 ORF integrity | 2211 bp, no frameshifts | ABORT |
| VP2/VP3/AAP frames correct | Proper start codons present | ABORT |
| All 6 VP1 sites unique | Exactly 1 occurrence each | ABORT |
| Upstream site unique | Exactly 1 occurrence | ABORT |
| Downstream site unique | Exactly 1 occurrence | ABORT |
| 9 VR regions annotated | All present in output | WARN |
| AAP overlapping CDS annotated | Feature present | WARN |
| scAAV size check | >2.4 kb → flag oversized | WARN |
| ssAAV size check | Must be <4.7 kb | ABORT if exceeded |

---

## **CHANGES FROM v1.0**

1. **VP1 Source:** Now uses engineered plasmid with 6 pre-existing silent sites instead of native AAV9
2. **Junction Sites:** Added 2 new unique restriction sites for modular cloning
3. **3' UTR:** Added minimal 3' UTR element between stop and polyA
4. **ITR Type:** Explicitly specified ITR128 (truncated) instead of full-length
5. **VP2 Start:** Explicitly calls out non-canonical ACG start codon
6. **AAP Frame:** Explicitly requires +1 frame annotation
7. **VR Regions:** All 9 variable regions must be annotated
8. **Restriction Sites:** All 8 unique sites (6 internal + 2 junctions) must be verified

---

**END OF PROMPT v2.0**
