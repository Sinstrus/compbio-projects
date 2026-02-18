# Design Verification Report: AVD005-AVD011

**Date:** 2026-01-22
**Revision:** v1.1 (Fixed BUG-006: N-terminal insertion positions corrected)

**Base plasmid:** AVD002-Rep2Mut2Cap9-6R-wt.dna (7,104 bp)

---

## ⚠️ Important Note: BUG-006 Fix

This report reflects **CORRECTED** designs for AVD008, AVD009, and AVD010. The initial designs incorrectly placed N-terminal VHH3 insertions after a single amino acid (M) instead of after the M-A dipeptide.

**What was fixed:**
- **AVD008/009:** VP1 N-terminal insertion moved from bp 2381 → **bp 2384** (after M-A)
- **AVD010:** VP2 N-terminal insertion moved from bp 2792 → **bp 2795** (after M-A)

**Corrected sequences:**
- AVD008/009: M-A-[VHH3]-[GGGGS5]-**A**-D-G-Y-L... ✓ (previously missing first A)
- AVD010: M-A-[VHH3]-[GGGGS5]-**P**-G-K-K-R... ✓ (previously missing A before P)

**For details, see:**
- LESSONS_LEARNED.md: BUG-006 entry
- Verification script: `scripts/verify_flanking_sequences.py`

---

## Summary of Constructs

| ID | Size (bp) | VP1 | VP2 | VP3 | Description |
|----|-----------|-----|-----|-----|-------------|
| AVD005 | 7,104 | KO (AAG) | Native (ACG) | Native (ATG) | VP1 knockout helper |
| AVD006 | 7,551 | VHH3-VR4 | Native (ACG) | Native (ATG) | VHH3 at VR4, VP2/3 functional |
| AVD007 | 7,551 | VHH3-VR4 | KO (ACC) | KO (CTG) | VHH3 at VR4, VP2/3 KO |
| AVD008 | 7,536 | VHH3-Nterm | Native (ACG) | Native (ATG) | VHH3 at VP1 N-term, VP2/3 functional |
| AVD009 | 7,536 | VHH3-Nterm | KO (ACC) | KO (CTG) | VHH3 at VP1 N-term, VP2/3 KO |
| AVD010 | 7,536 | Native (ATG) | VHH3-Nterm | Native (ATG) | VHH3 at VP2 N-term, VP1/3 functional |
| AVD011 | 7,104 | Native (ATG) | KO (ACC) | Native (ATG) | VP2 knockout helper |

---

## AVD005: VP1 knockout helper

**Size:** 7,104 bp

**Modifications:**
- VP1 knockout: ATG→AAG at bp 2379-2381

**Verification:**
- [x] VP1 = AAG (knockout)
- [x] VP2 = ACG (native)
- [x] VP3 = ATG (native)
- [x] Total size: 7,104 bp

---

## AVD006: VHH3 at VR4, VP2/3 functional

**Size:** 7,551 bp

**Modifications:**
- VHH3 insert at VR4: 447 bp (GGGGS5-VHH3-GGGGS1) at bp 3743

**Verification:**
- [x] VP1 = ATG (native, has VHH3 insert)
- [x] VP2 = ACG (native)
- [x] VP3 = ATG (native)
- [x] VHH3 insert at bp 3743: 447 bp
- [x] Total size: 7,551 bp

---

## AVD007: VHH3 at VR4, VP2/3 KO

**Size:** 7,551 bp

**Modifications:**
- VP2 knockout: ACG→ACC at bp 2790-2792 (silent)
- VP3 knockout: ATG→CTG at bp 2985-2987 (non-silent)
- VHH3 insert at VR4: 447 bp (GGGGS5-VHH3-GGGGS1) at bp 3743

**Verification:**
- [x] VP1 = ATG (native, has VHH3 insert)
- [x] VP2 = ACC (knockout)
- [x] VP3 = CTG (knockout)
- [x] VHH3 insert at bp 3743: 447 bp
- [x] Total size: 7,551 bp

---

## AVD008: VHH3 at VP1 N-term, VP2/3 functional

**Size:** 7,536 bp

**Modifications:**
- VHH3 insert at VP1 N-term: 432 bp (VHH3-GGGGS5) at bp 2384
- Sequence: M-A-[VHH3]-[GGGGS5]-A-D-G-Y-L...

**Verification:**
- [x] VP1 N-term VHH3 insert at bp 2384 (after M-A): 432 bp
- [x] VP2 = ACG (native)
- [x] VP3 = ATG (native)
- [x] Total size: 7,536 bp

---

## AVD009: VHH3 at VP1 N-term, VP2/3 KO

**Size:** 7,536 bp

**Modifications:**
- VP2 knockout: ACG→ACC at bp 2790-2792 (silent)
- VP3 knockout: ATG→CTG at bp 2985-2987 (non-silent)
- VHH3 insert at VP1 N-term: 432 bp (VHH3-GGGGS5) at bp 2384
- Sequence: M-A-[VHH3]-[GGGGS5]-A-D-G-Y-L...

**Verification:**
- [x] VP1 N-term VHH3 insert at bp 2384 (after M-A): 432 bp
- [x] VP2 = ACC (knockout)
- [x] VP3 = CTG (knockout)
- [x] Total size: 7,536 bp

---

## AVD010: VHH3 at VP2 N-term, VP1/3 functional

**Size:** 7,536 bp

**Modifications:**
- VHH3 insert at VP2 N-term: 432 bp (VHH3-GGGGS5) at bp 2795
- Sequence: M-A-[VHH3]-[GGGGS5]-P-G-K-K-R... (ACG→M at translation)

**Verification:**
- [x] VP1 = ATG (native)
- [x] VP2 N-term VHH3 insert at bp 2795 (after M-A): 432 bp
- [x] VP3 = ATG (native)
- [x] Total size: 7,536 bp

---

## AVD011: VP2 knockout helper

**Size:** 7,104 bp

**Modifications:**
- VP2 knockout: ACG→ACC at bp 2790-2792 (silent)

**Verification:**
- [x] VP1 = ATG (native)
- [x] VP2 = ACC (knockout)
- [x] VP3 = ATG (native)
- [x] Total size: 7,104 bp

---

## VHH3 Insert Sequences

**VHH3 protein sequence** (119 aa):
```
EVQLVESGGGLVQPGGSLRLSCAASGFTFSTADMGWFRQAPGKGRELVAAVSGSGFSTYSDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAKATISLYYAMDVWGQGTTVTVSS
```

**VHH3 DNA sequence** (357 bp):
```
GAGGTGCAACTGGTTGAAAGCGGCGGAGGACTTGTTCAACCCGGCGGCAGCCTTAGGCTTTCTTGCGCTGCCAGCGGCTTCACCTTTAGCACCGCCGACATGGGCTGGTTTAGGCAAGCTCCCGGAAAAGGCAGGGAACTTGTTGCCGCTGTGAGCGGCAGCGGCTTCAGCACCTACTCTGATAGCGTTGAGGGCAGGTTCACCATCAGCAGGGACAACGCCAAGAGGATGGTGTACCTGCAGATGAACAGCTTGAGGGCCGAGGACACCGCCGTGTACTACTGCGCCAAGGCCACAATTAGCCTGTACTACGCCATGGATGTGTGGGGACAGGGCACCACCGTGACCGTGAGCAGC
```

**GGGGS5 linker** (25 aa, 75 bp):
```
GGTGGAGGCGGATCTGGAGGCGGTGGTTCAGGCGGTGGAGGAAGTGGTGGCGGAGGTTCTGGTGGAGGCGGTTCT
```

**GGGGS1 linker** (5 aa, 15 bp):
```
GGTGGAGGCGGTTCT
```

## Next Steps

1. Sequence verify all plasmids (Sanger sequencing)
2. Production screening:
   - AVD005 + AVD007: VP1-VHH3 only display
   - AVD005 + AVD009: VP1-VHH3 N-term only display
   - Compare display efficiency across constructs
3. Functional assays:
   - Western blot (VP expression, VHH3 detection)
   - ELISA (anti-ALPL binding)
   - Transduction efficiency
   - ALPL-targeted cell binding

---
*Report generated: 2026-01-22 13:40:08*
