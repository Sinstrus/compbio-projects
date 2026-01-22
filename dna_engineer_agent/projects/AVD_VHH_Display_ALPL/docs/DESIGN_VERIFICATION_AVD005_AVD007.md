# Design Verification Report: AVD005 and AVD007

**Date:** 2026-01-20

**Base plasmid:** AVD002-Rep2Mut2Cap9-6R-wt.dna (7,104 bp)

## AVD005: VP1 Knockout Helper Plasmid

### Design Summary
- **Size:** 7104 bp
- **Modification:** VP1 start codon knockout (ATG→AAG at bp 2379-2381)
- **Purpose:** Trans-complementation helper for AVD007

### Verification Checklist
- [x] VP1 start codon = AAG (knockout)
- [x] VP2 start codon = ACG (preserved)
- [x] VP3 start codon = ATG (preserved)
- [x] Total size = 7,104 bp

## AVD007: VHH3 Display with D2 Linkers

### Design Summary
- **Size:** 7551 bp
- **Modifications:**
  - VP2 knockout: ACG→ACC at bp 2790-2792 (silent)
  - VP3 knockout: ATG→CTG at bp 2985-2987 (non-silent)
  - VHH3 insertion: 447 bp at VR-IV region (bp 3743, after SKTINGSG)
- **Insert composition:**
  - N-terminal linker: 5×GGGGS (25 aa, 75 bp)
  - VHH3: 119 aa, 357 bp
  - C-terminal linker: 1×GGGGS (5 aa, 15 bp)
  - **Total insert:** 149 aa, 447 bp

### Verification Checklist
- [x] VP2 knockout = ACC
- [x] VP3 knockout = CTG
- [x] Insert length = 447 bp
- [x] No stop codons in insert
- [x] Insert is multiple of 3
- [x] Total size = 7,551 bp

### Insert Sequence Analysis

- **GC content:** 61.3%
- **Translation:** GGGGSGGGGSGGGGSGGGGSGGGGSEVQLVESGGGLVQPGGSLRLSCAASGFTFSTADMGWFRQAPGKGRELVAAVSGSGFSTYSDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAKATISLYYAMDVWGQGTTVTVSSGGGGS

### Synthetic Fragment for Ordering

```
GGTGGAGGCGGATCTGGAGGCGGTGGTTCAGGCGGTGGAGGAAGTGGTGGCGGAGGTTCTGGTGGAGGCGGTTCTGAGGTGCAACTGGTTGAAAGCGGCGGAGGACTTGTTCAACCCGGCGGCAGCCTTAGGCTTTCTTGCGCTGCCAGCGGCTTCACCTTTAGCACCGCCGACATGGGCTGGTTTAGGCAAGCTCCCGGAAAAGGCAGGGAACTTGTTGCCGCTGTGAGCGGCAGCGGCTTCAGCACCTACTCTGATAGCGTTGAGGGCAGGTTCACCATCAGCAGGGACAACGCCAAGAGGATGGTGTACCTGCAGATGAACAGCTTGAGGGCCGAGGACACCGCCGTGTACTACTGCGCCAAGGCCACAATTAGCCTGTACTACGCCATGGATGTGTGGGGACAGGGCACCACCGTGACCGTGAGCAGCGGTGGAGGCGGTTCT
```

## Next Steps

1. Order synthetic fragment from GenScript (~$200-400)
2. Clone into AVD002 using appropriate restriction sites
3. Sequence verify both plasmids (Sanger sequencing)
4. Co-transfect AVD005 (helper) + AVD007 (VHH display)
5. Functional assays: Western blot, ELISA, transduction

---
*Report generated: 2026-01-20 16:56:30*
