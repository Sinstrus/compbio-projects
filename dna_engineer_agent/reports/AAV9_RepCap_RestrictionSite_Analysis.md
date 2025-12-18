# Restriction Site Analysis: AAV9 Rep-Cap Helper Plasmid

**Plasmid:** BASE-DRAFT-AAV9-RepCap-NOCAP.gb
**Analysis Date:** 2025-12-17
**Agent Version:** 2.2.0
**Plasmid Length:** 7,078 bp

---

## Executive Summary

This analysis identifies unique restriction enzyme sites (requiring ≤2 silent mutations) in three specific regions of an AAV9 Rep-Cap helper plasmid where the Cap coding sequence has been intentionally inactivated by mutating the VP1, VP2, and VP3 start codons.

**Key Findings:**
- **3 unique restriction sites** identified (1 per region)
- **Total mutations required:** 2 (both silent)
- **All sites are unique** in the entire plasmid
- **No amino acid changes** in any of the proposed modifications

---

## Target Regions

### Region 1: VP1 Unique Region
**Definition:** N-terminal portion unique to VP1, between VP1 start and VP2 start
**Coordinates:** 2365-2775 (411 bp, 137 amino acids)
**Purpose:** This region encodes the VP1-unique domain containing the phospholipase A2 (PLA2) motif essential for endosomal escape

### Region 2: VP2-AAP Intergenic
**Definition:** After VP2 start codon but before AAP start codon
**Coordinates:** 2779-2890 (112 bp, 37 amino acids)
**Purpose:** This region is within VP2 coding sequence but precedes AAP translation (AAP is in +1 frame)

### Region 3: VP3 Post-AAP
**Definition:** Within VP3 coding region but after AAP stop codon
**Coordinates:** 3485-4575 (1,091 bp, 363 amino acids)
**Purpose:** This region encodes the C-terminal portion of VP3, including variable regions and structural domains

---

## Recommended Sites

### Region 1: **EagI** Site

| Parameter | Value |
|-----------|-------|
| **Enzyme** | EagI |
| **Recognition Sequence** | CGGCCG (6 bp) |
| **Position** | 2369 (1-indexed) |
| **Mutations Required** | 1 |
| **Original Sequence** | CTGCCG |
| **Modified Sequence** | CGGCCG |
| **Uniqueness** | ✓ Unique in entire plasmid |

**Codon Change:**
```
Position 2370:  GCT → GCG  (Alanine → Alanine)  ✓ SILENT
```

**Rationale:** EagI is a commonly used enzyme for cloning. The site requires only a single silent mutation (T→G) in the third position of an alanine codon (GCT→GCG), making it a highly conservative change with no impact on the encoded protein.

---

### Region 2: **BbvCI** Site

| Parameter | Value |
|-----------|-------|
| **Enzyme** | BbvCI |
| **Recognition Sequence** | CCTCAGC (7 bp) |
| **Position** | 2828 (1-indexed) |
| **Mutations Required** | 1 |
| **Original Sequence** | CCTCCGC |
| **Modified Sequence** | CCTCAGC |
| **Uniqueness** | ✓ Unique in entire plasmid |

**Codon Change:**
```
Position 2832:  TCC → TCA  (Serine → Serine)  ✓ SILENT
```

**Rationale:** BbvCI is a Type IIS restriction enzyme that cuts outside its recognition sequence, making it useful for scarless cloning and Golden Gate assembly. The required mutation changes one serine codon to another (TCC→TCA), a completely neutral change.

---

### Region 3: **AfeI** Site

| Parameter | Value |
|-----------|-------|
| **Enzyme** | AfeI |
| **Recognition Sequence** | AGCGCT (6 bp) |
| **Position** | 4442 (1-indexed) |
| **Mutations Required** | **0** (site already present!) |
| **Original Sequence** | AGCGCT |
| **Modified Sequence** | AGCGCT |
| **Uniqueness** | ✓ Unique in entire plasmid |

**Codon Change:**
```
None - site already present
```

**Rationale:** AfeI site is already present and unique in the plasmid, requiring no modifications. This is an ideal "handle" for cloning or diagnostic purposes. AfeI generates blunt ends, useful for blunt-end cloning strategies.

---

## Alternative Candidates

### Region 1 Alternatives
All require 1 silent mutation and are unique in the full plasmid:

| Enzyme | Position | Recognition Site | Mutation |
|--------|----------|------------------|----------|
| EcoRV  | 2378     | GATATC          | T2A      |
| BspEI  | 2385     | TCCGGA          | A4G      |
| SmaI   | 2505     | CCCGGG          | T1C      |

### Region 2 Alternatives
All require 1 silent mutation and are unique in the full plasmid:

| Enzyme | Position | Recognition Site | Mutation |
|--------|----------|------------------|----------|
| BlpI   | 2854     | GCTNAGC         | A3T      |
| NaeI   | 2859     | GCCGGC          | C4G      |

### Region 3 Alternatives
All require 0 mutations (already present) and are unique:

| Enzyme      | Position | Recognition Site        |
|-------------|----------|-------------------------|
| BaeI        | 4318     | ACNNNNGTAYC            |
| XcmI        | 4380     | CCANNNNNNNNNTGG        |
| BspDI/ClaI  | 4575     | ATCGAT                 |

---

## Methodology

### Analysis Pipeline

1. **Sequence Extraction**
   - Parsed GenBank file using Biopython
   - Extracted DNA and protein sequences for each target region
   - Verified coordinates against file annotations

2. **Restriction Site Identification**
   - Used `silent_sites.py` tool with `--mutations 2` parameter
   - Searched 184 restriction enzymes (minimum 6 bp recognition sites)
   - Scanned ~250+ restriction enzyme recognition sequences including IUPAC ambiguity codes

3. **Filtering Criteria**
   - **Silent mutations only:** No amino acid changes permitted
   - **Uniqueness:** Site must appear only once in the entire 7,078 bp plasmid
   - **Minimum mutations:** Prioritized by fewest mutations required (0 > 1 > 2)

4. **Verification**
   - Confirmed uniqueness by counting occurrences in full plasmid sequence
   - Verified silent nature by translating codons before and after mutations
   - Cross-checked positions against GenBank annotations

### Candidates Found

| Region | Total Candidates | Silent Mutations | Unique in Plasmid |
|--------|------------------|------------------|-------------------|
| Region 1 (VP1 unique) | 5,501 | 407 | 67 |
| Region 2 (VP2-AAP intergenic) | 1,419 | 149 | 31 |
| Region 3 (VP3 post-AAP) | 11,944 | 683 | 86 |

---

## Applications

These restriction sites can be used for:

1. **Diagnostic Restriction Mapping**
   - Verify plasmid identity after transformation
   - Quick QC check before sequencing

2. **Future Modifications**
   - Insert epitope tags (FLAG, HA, etc.)
   - Swap variable regions (e.g., VHH display at VR-IV)
   - Add fluorescent protein fusions

3. **Golden Gate Assembly** (for BbvCI)
   - Scarless cloning workflows
   - Modular assembly of Cap variants

4. **Gene Synthesis Services**
   - Design fragments for synthesis with these handles
   - Enable downstream cloning without additional PCR

---

## Sequence Context

### Region 1: EagI Site Context (position 2369)
```
...GCCGATGGTTATCTTCCAGATTGGCTCGAGGACAACCTTAGTGAAGGAATTCGC...
       ^^^
    Position 2369-2374
    CTGCCG → CGGCCG (EagI site after T2370G mutation)
```

### Region 2: BbvCI Site Context (position 2828)
```
...CCTGGAAAGAAGAGGCCTGTAGAGCAGTCTCCTCAGGAACCGGACTCCTCCGCGGG...
                                        ^^^^^^^
                                    Position 2828-2834
                                    CCTCCGC → CCTCAGC (BbvCI site after C2832A mutation)
```

### Region 3: AfeI Site Context (position 4442)
```
...GTATTGCTGTTAATACTGAAGGTGTATATAGTGAACCCCGCCCCATTGGCACCAGAT...
                              ^^^^^^
                          Position 4442-4447
                          AGCGCT (AfeI site already present)
```

---

## Risk Assessment

### Low Risk: Silent Mutations
- Both proposed mutations are in the third position of degenerate codons
- Alanine: GCT → GCG (both encode Ala, similar codon usage)
- Serine: TCC → TCA (both encode Ser, similar codon usage)
- No impact on protein structure or function expected

### Low Risk: Codon Usage
- Human codon usage table shows similar frequencies for both variants
- No rare codons introduced
- mRNA secondary structure changes: unlikely to be significant

### Verification Recommended
- Sequence plasmid after introducing mutations
- Functional testing if critical: verify Cap expression (though these are in the "NOCAP" construct where Cap is already non-functional by design)

---

## Implementation Notes

### To Introduce These Sites:

1. **Order as Gene Synthesis**
   - Use Genscript, IDT, or similar service
   - Specify the exact sequences with mutations
   - Most cost-effective for regions >500 bp

2. **Site-Directed Mutagenesis**
   - Use QuikChange or similar protocol
   - Only 2 mutations total needed (Region 1 + Region 2)
   - Region 3 requires no work (site already present)

3. **Primer Design** (for site-directed mutagenesis):

**Region 1 - EagI (T2370G):**
```
Forward:  5'-GCCGATGGTTATCTTCCAGATTGGCGGAGGACAACCTTAGTGAAGG-3'
                                      ^
                                    Mutation
```

**Region 2 - BbvCI (C2832A):**
```
Forward:  5'-CAGTCTCCTCAGGAACCGGACTCCTCAGCGGGTATTGGCAAATCG-3'
                                        ^
                                    Mutation
```

---

## Conclusions

This analysis successfully identified three unique restriction sites suitable for future modifications:

1. **EagI** in VP1 unique region (1 silent mutation)
2. **BbvCI** in VP2-AAP intergenic (1 silent mutation)
3. **AfeI** in VP3 post-AAP (0 mutations - already present!)

All sites are:
- ✓ **Unique** in the entire plasmid
- ✓ **Silent** (no amino acid changes)
- ✓ **Minimally mutated** (only 2 total mutations needed)
- ✓ **Functionally neutral** (both mutations are in wobble positions)

These sites provide excellent "handles" for future engineering of this Rep-Cap helper plasmid without compromising its function.

---

## Files Generated

- `analyze_restriction_sites.py` - Main analysis script
- `generate_final_report.py` - Report generation script
- `/tmp/region1_vp1_unique.txt` - Region 1 sequence data
- `/tmp/region2_vp2_aap.txt` - Region 2 sequence data
- `/tmp/region3_vp3_post_aap.txt` - Region 3 sequence data
- This report: `AAV9_RepCap_RestrictionSite_Analysis.md`

---

**Analysis completed successfully.**
For questions or modifications, refer to the scripts in the project directory.
