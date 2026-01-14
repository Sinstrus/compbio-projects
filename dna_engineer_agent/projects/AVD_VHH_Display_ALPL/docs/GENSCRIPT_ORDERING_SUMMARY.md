# GenScript Ordering Summary: AVD005 & AVD006
**VP1-VHH3-ALPL Synthetic Fragments**

**Date:** 2026-01-14
**Project:** AAV VHH Display Platform - Phase I
**Target:** ALPL (Alkaline Phosphatase)

---

## Summary

Two VP1-only expression constructs have been designed with anti-ALPL VHH3 nanobody inserted at variable region IV (VR-IV) using Design 2 asymmetric flexible linkers. Complete plasmid sequences and synthetic fragments for ordering are provided below.

---

## Complete Plasmids Generated

| Construct | Size | Description | GenBank File |
|-----------|------|-------------|--------------|
| **AVD005** | 6,690 bp | EF1a-VP1-VHH3-ALPL-bGH | AVD005-EF1A-VP1-VHH3-ALPL-bGH.gb |
| **AVD006** | 7,521 bp | Rep2Mut2Cap9-VP1-VHH3-ALPL | AVD006-Rep2Mut2Cap9-VP1-VHH3-ALPL.gb |

**Features preserved in GenBank files:**
- ori (replication origin)
- AmpR CDS (ampicillin resistance)
- AmpR promoter
- M13 primers
- lac operator/promoter
- All promoters (EF1a for AVD005, p5/p19/p40 for AVD006)
- Rep genes (AVD006 only)

---

## Synthetic Fragments for Ordering

### AVD005: AvrII to BsrGI Fragment

**Fragment specifications:**
- **Position:** 1676-3358 (1-indexed, GenBank coordinates)
- **Length:** 1,683 bp
- **Left boundary:** AvrII site (CCTAGG) at position 1676
- **Right boundary:** BsrGI site (TGTACA) at position 3358
- **FASTA file:** `AVD005_synthetic_fragment_FINAL.fasta`

**What's included:**
- AvrII recognition site (complete 6 bp)
- **VP1 coding sequence from AvrII to BsrGI** with:
  - VP2 knockout (ACG->ACC at codon 138, silent)
  - VP3 knockout (ATG->CTG at codon 203, non-silent)
  - N-terminal linker: (GGGGS)x4 = 20 aa (60 bp)
  - VHH3 anti-ALPL: 119 aa (357 bp, codon-optimized)
  - C-terminal linker: direct fusion (0 aa)
- BsrGI recognition site (complete 6 bp)

**Cloning strategy:**
1. Digest AVD003 parent plasmid with AvrII and BsrGI
2. Digest synthesized fragment with AvrII and BsrGI
3. Gel purify both fragments
4. Ligate and transform into competent E. coli
5. Select on ampicillin

---

### AVD006: SmaI to BsrGI Fragment

**Fragment specifications:**
- **Position:** 2519-4216 (1-indexed)
- **Length:** 1,698 bp
- **Left boundary:** SmaI site (CCCGGG) at position 2519
- **Right boundary:** BsrGI site (TGTACA) at position 4216
- **FASTA file:** `AVD006_synthetic_fragment_FINAL.fasta`

**What's included:**
- SmaI recognition site (complete 6 bp)
- **Complete VP1-VHH3 fusion region** (identical VP1-VHH sequence to AVD005)
- VP2 knockout (ACG->ACC at codon 138, silent)
- VP3 knockout (ATG->CTG at codon 203, non-silent)
- N-terminal linker: (GGGGS)x4 = 20 aa (60 bp)
- VHH3 anti-ALPL: 119 aa (357 bp, codon-optimized)
- C-terminal linker: direct fusion (0 aa)
- BsrGI recognition site (complete 6 bp)

**Why SmaI instead of HindIII?**
- SmaI (2519) is 461 bp closer to the first change (2792) than HindIII (2058)
- Saves ~$115-$161 in synthesis costs
- Both sites are unique in AVD002 and AVD006

**Cloning strategy:**
1. Digest AVD002 parent plasmid with SmaI and BsrGI
2. Digest synthesized fragment with SmaI and BsrGI
3. Gel purify both fragments
4. Ligate and transform into competent E. coli
5. Select on appropriate antibiotic

**Note:** SmaI creates blunt ends. Use high-quality enzyme and optimize ligation.

---

## Ordering Information for GenScript

### Summary Table

| Construct | Fragment | Length | Start | End | FASTA File |
|-----------|----------|--------|-------|-----|------------|
| AVD005 | AvrII to BsrGI | 1,683 bp | 1676 | 3358 | AVD005_synthetic_fragment_FINAL.fasta |
| AVD006 | SmaI to BsrGI | 1,698 bp | 2519 | 4216 | AVD006_synthetic_fragment_FINAL.fasta |

### Recommended Order

**For AVD005:**
- Order fragment: `AVD005_synthetic_fragment_FINAL.fasta` (1,683 bp)
- Boundaries: AvrII (1676) to BsrGI (3358)
- Cost estimate: ~$420-$589 for 1.7 kb synthesis

**For AVD006:**
- Order fragment: `AVD006_synthetic_fragment_FINAL.fasta` (1,698 bp)
- Boundaries: SmaI (2519) to BsrGI (4216)
- Cost estimate: ~$424-$594 for 1.7 kb synthesis

**Total estimated cost:** $844-$1,183 for both constructs

---

## Quality Control Specifications

### Sequence Verification Requirements
- **100% sequence accuracy** for:
  - VP2 knockout mutation (ACG->ACC at relative position 414)
  - VP3 knockout mutation (ATG->CTG at relative position 607)
  - VHH3 coding sequence (357 bp)
  - D2 N-terminal linker (60 bp)
  - All restriction sites (AvrII, SmaI, BsrGI)

### Codon Optimization
- **Organism:** *Homo sapiens*
- **Strategy:** High-frequency codons used throughout VHH3 sequence
- **GC content:** Should be 50-60% for optimal expression

---

## Post-Synthesis Assembly Instructions

### For AVD005:
1. Receive synthesized fragment (1,683 bp) from GenScript in provided vector
2. Digest both fragment and AVD003 parent plasmid with AvrII and BsrGI
3. Gel purify digested fragments
4. Ligate fragment into AVD003 backbone
5. Transform into competent *E. coli* (DH5a or similar)
6. Select on ampicillin plates
7. Screen colonies by:
   - Restriction digest (verify AvrII, BsrGI presence)
   - Sanger sequencing (verify VP2/VP3 knockouts, VHH insertion)
   - Full plasmid sequencing (verify no off-target mutations)

### For AVD006:
1. Receive synthesized fragment (1,698 bp) from GenScript
2. Digest AVD002 parent plasmid with SmaI and BsrGI
3. Digest synthesized fragment with SmaI and BsrGI
4. Gel purify both fragments
5. Ligate and transform (note: SmaI creates blunt ends - optimize ligation)
6. Screen as above

---

## Files Provided

### GenBank Files (Complete Plasmids)
- `AVD005-EF1A-VP1-VHH3-ALPL-bGH.gb` (6,690 bp)
- `AVD006-Rep2Mut2Cap9-VP1-VHH3-ALPL.gb` (7,521 bp)

### FASTA Files (Synthetic Fragments for Ordering)
- `AVD005_synthetic_fragment_FINAL.fasta` (1,683 bp) - AvrII to BsrGI
- `AVD006_synthetic_fragment_FINAL.fasta` (1,698 bp) - SmaI to BsrGI

### Documentation
- `SYNTHETIC_FRAGMENTS_FINAL.csv` - Ordering spreadsheet with full sequences
- `DESIGN_VERIFICATION_AVD005_AVD006.md` - Complete design verification report
- `GENSCRIPT_ORDERING_SUMMARY.md` - This file

---

## Key Design Features

### VP1-VHH3 Fusion Architecture
```
VP1-N-terminus (1-455)
    -> LINK_D2_N (GGGGS)x4 [20 aa]
    -> VHH3-anti-ALPL [119 aa]
    -> LINK_D2_C [0 aa, direct fusion]
    -> VP1-C-terminus (456-end)
```

### Mutations Introduced
| Position | Type | Original | Mutated | AA Change | Purpose |
|----------|------|----------|---------|-----------|---------|
| Codon 138 | VP2 start | ACG | ACC | Thr->Thr (silent) | VP2 knockout |
| Codon 203 | VP3 start | ATG | CTG | Met->Leu (non-silent) | VP3 knockout |

**No other mutations** - the VP1 backbone is otherwise identical to the parent plasmid.

---

## Contact & Support

For questions about:
- **Sequence design:** Refer to `DESIGN_VERIFICATION_AVD005_AVD006.md`
- **Cloning strategy:** See assembly instructions above
- **Synthesis issues:** Contact GenScript technical support
- **Project timeline:** Mice experiments end of March 2026

---

**STATUS: READY FOR SYNTHESIS ORDER**

All sequences verified, fragments defined, and documentation complete.
