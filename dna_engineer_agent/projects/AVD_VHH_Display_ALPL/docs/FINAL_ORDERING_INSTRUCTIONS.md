# FINAL Ordering Instructions: AVD005 & AVD006
**Synthetic Fragments for GenScript Order**

**Date:** 2026-01-14
**Status:** READY FOR GENSCRIPT ORDER

---

## Final Synthetic Fragments for GenScript Order

| Construct | Fragment Size | Boundaries | Cost Estimate | FASTA File |
|-----------|---------------|------------|---------------|------------|
| **AVD005** | **1,683 bp** | AvrII (1676) -> BsrGI (3358) | **$420-$589** | AVD005_synthetic_fragment_FINAL.fasta |
| **AVD006** | **1,698 bp** | SmaI (2519) -> BsrGI (4216) | **$424-$594** | AVD006_synthetic_fragment_FINAL.fasta |
| **TOTAL** | **3,381 bp** | - | **$844-$1,183** | Both fragments |

---

## AVD005 Fragment Details

**Fragment:** AvrII to BsrGI (1,683 bp)

**Coordinates:** 1676-3358 (1-indexed, as in GenBank)

**What's included:**
- AvrII recognition site (CCTAGG) - complete 6 bp
- **VP1 coding region from AvrII to BsrGI** including:
  - VP2 knockout: ACG->ACC at codon 138 (silent)
  - VP3 knockout: ATG->CTG at codon 203 (non-silent)
  - VHH3 insertion at VR-IV (aa 456) with D2 linkers
  - N-linker: (GGGGS)x4 = 60 bp
  - VHH3: codon-optimized, 357 bp (119 aa)
  - C-linker: direct fusion (0 bp)
- BsrGI recognition site (TGTACA) - complete 6 bp

**Why AvrII to BsrGI?**
- These are the closest unique restriction sites FLANKING the modified region
- Captures all nucleotide changes between AVD003 and AVD005
- Smaller fragment = cheaper synthesis

**Cloning strategy:**
```
1. Digest AVD003 parent with AvrII and BsrGI
2. Digest synthetic fragment with AvrII and BsrGI
3. Gel purify
4. Ligate synthetic fragment into AVD003 backbone
5. Transform and select on ampicillin
```

---

## AVD006 Fragment Details

**Fragment:** SmaI to BsrGI (1,698 bp)

**Coordinates:** 2519-4216 (1-indexed)

**What's included:**
- SmaI recognition site (CCCGGG) - complete 6 bp
- **Complete VP1-VHH3 fusion region** with:
  - VP2 knockout: ACG->ACC at codon 138 (silent)
  - VP3 knockout: ATG->CTG at codon 203 (non-silent)
  - VHH3 insertion at VR-IV (aa 456) with D2 linkers
- BsrGI recognition site (TGTACA) - complete 6 bp

**Why SmaI to BsrGI?**
- SmaI is the closest unique site UPSTREAM of the first change (at position 2792)
- BsrGI is the closest unique site DOWNSTREAM of the last change (at position 4163)
- **461 bp smaller than HindIII boundary** - saves ~$115-$161

**Cloning strategy:**
```
1. Digest AVD002 parent with SmaI and BsrGI
2. Digest synthetic fragment with SmaI and BsrGI
3. Gel purify
4. Ligate synthetic fragment into AVD002 backbone
5. Transform and select on appropriate antibiotic
```

**Note:** SmaI creates blunt ends. Ensure high-quality enzyme and optimize ligation conditions.

---

## File Organization

### FINAL Files (Use These!)
```
AVD005-EF1A-VP1-VHH3-ALPL-bGH.gb              # Complete plasmid (6,690 bp)
AVD005_synthetic_fragment_FINAL.fasta         # Order this (1,683 bp)

AVD006-Rep2Mut2Cap9-VP1-VHH3-ALPL.gb          # Complete plasmid (7,521 bp)
AVD006_synthetic_fragment_FINAL.fasta         # Order this (2,159 bp)

SYNTHETIC_FRAGMENTS_FINAL.csv                 # Ordering spreadsheet with sequences
```

---

## GenScript Order Form

**Copy this information into GenScript order form:**

### Fragment 1: AVD005
- **Construct name:** AVD005-VP1-VHH3-ALPL-AvrII-BsrGI
- **Sequence length:** 1,683 bp
- **Sequence file:** Upload `AVD005_synthetic_fragment_FINAL.fasta`
- **Cloning sites:** AvrII (5') and BsrGI (3')
- **Expression system:** Mammalian (HEK293)
- **Optimization:** Already codon-optimized for *Homo sapiens*
- **Vector:** Provide in standard cloning vector (pUC57 or similar)
- **Quantity:** 4-10 ug purified plasmid DNA
- **QC:** Full Sanger sequencing (both strands)

### Fragment 2: AVD006
- **Construct name:** AVD006-VP1-VHH3-ALPL-SmaI-BsrGI
- **Sequence length:** 1,698 bp
- **Sequence file:** Upload `AVD006_synthetic_fragment_FINAL.fasta`
- **Cloning sites:** SmaI (5') and BsrGI (3')
- **Expression system:** Mammalian (HEK293)
- **Optimization:** Already codon-optimized for *Homo sapiens*
- **Vector:** Provide in standard cloning vector (pUC57 or similar)
- **Quantity:** 4-10 ug purified plasmid DNA
- **QC:** Full Sanger sequencing (both strands)

---

## Quality Control Checklist

Upon receipt of synthetic fragments, verify:

### Critical Mutations
- [ ] VP2 knockout: ACG->ACC at codon 138 (position 414 relative to VP1 start)
- [ ] VP3 knockout: ATG->CTG at codon 203 (position 607 relative to VP1 start)
- [ ] VP1 start: ATG preserved at position 0

### VHH3 Insertion (verify all of these)
- [ ] N-linker present: (GGGGS)x4 = 60 bp at correct position
- [ ] VHH3 sequence: 357 bp coding for 119 aa anti-ALPL nanobody
- [ ] Insertion at aa 456 of VP1 (bp 1368)
- [ ] No stop codons in VHH3 or linker regions

### Restriction Sites
- [ ] AvrII present at 5' end of AVD005 fragment: CCTAGG
- [ ] SmaI present at 5' end of AVD006 fragment: CCCGGG
- [ ] BsrGI present at 3' end of both fragments: TGTACA

### Full Plasmid Assembly
- [ ] Digest check: AvrII/BsrGI (AVD005) or SmaI/BsrGI (AVD006) yields expected fragment size
- [ ] Colony PCR: Correct insert size
- [ ] Sanger sequencing: 100% match to designed sequence
- [ ] Midi/maxi prep: >1 ug/uL concentration
- [ ] Sequence entire plasmid to confirm no off-target mutations

---

## Timeline

**Week 1:** Submit order to GenScript with FASTA files
**Week 2-4:** Synthesis (GenScript turnaround: 10-15 business days)
**Week 5:** Receive fragments, perform QC
**Week 6:** Clone into final vectors (AVD003/AVD002 backbones)
**Week 7:** Sequence verification and midi prep
**Week 8:** Ready for transfection experiments

**Target:** Mice by end of March 2026

---

## Cost Breakdown

| Item | Cost |
|------|------|
| AVD005 synthesis (1,683 bp) | $420-$589 |
| AVD006 synthesis (1,698 bp) | $424-$594 |
| Sanger sequencing (4 reactions @ $15 each) | $60 |
| Cloning supplies (enzymes, primers) | ~$200 |
| **TOTAL** | **$1,104-$1,443** |

*Synthesis cost estimates based on $0.25-0.35 per bp*

**Savings vs previous design (HindIII boundary):** $117-$162 (461 bp smaller)

---

## Sequence Comparison Summary

**Changes from parent plasmid to final construct:**

| Mutation | Codon | Change | AA | Type | Purpose |
|----------|-------|--------|-----|------|---------|
| VP2 knockout | 138 | ACG->ACC | T->T | SILENT | Prevent VP2 expression |
| VP3 knockout | 203 | ATG->CTG | M->L | NON-SILENT | Prevent VP3 expression |

**Insertion at aa 456 of VP1:**
- 417 bp insertion = (GGGGS)x4 linker (60 bp) + VHH3 (357 bp)

**No other mutations** - the VP1 backbone is otherwise identical to the parent plasmid.

---

## Support

For questions about:
- **Sequence design:** See `DESIGN_VERIFICATION_AVD005_AVD006.md`
- **Fragment boundaries:** See this document
- **Full plasmids:** See GenBank files (.gb)
- **GenScript order:** See FASTA files and this document

---

**Last updated:** 2026-01-14
**Status:** READY FOR ORDER
**Action required:** Upload FASTA files to GenScript order form
