# Checkpoints 1–7: Pre-Synthesis Verification

> Part of: [Agent Instructions](../README.md)

**Execute these checkpoints BEFORE finalizing any synthesis design.**

---

## Checkpoint 1: ITR Integrity (AAV-specific)
- [ ] 5' ITR sequence matches reference (100% identity or known functional variants)
- [ ] 3' ITR sequence matches reference
- [ ] ITR length correct for type (128bp, 145bp, etc.)
- [ ] No mutations within 10bp of ITR boundaries
- [ ] Secondary structure predictions match expected stem-loops

## Checkpoint 2: Regulatory Element Validation
- [ ] Promoter matches published sequence (≥95% identity)
- [ ] Kozak consensus optimal (GCCACCATGG) or acceptable (maintain -3 A/G and +4 G)
- [ ] PolyA signal present and correct (AATAAA hexamer preserved)
- [ ] PolyA signal 50-200bp downstream of stop codon
- [ ] WPRE (if present) matches reference sequence

## Checkpoint 3: Coding Sequence Integrity
- [ ] CDS in-frame from start (ATG) to stop (TAA/TAG/TGA)
- [ ] No internal stop codons
- [ ] Start codon in good Kozak context
- [ ] Protein BLAST confirms correct protein identity
- [ ] No frameshifts or indels vs. reference

## Checkpoint 4: Packaging Size Limit (AAV-specific)
- [ ] Distance between ITRs calculated
- [ ] scAAV: ≤ ~2.2 kb (must form dimer to package)
- [ ] ssAAV: ≤ ~4.7 kb (single-strand packaging limit)
- [ ] If borderline, flag for experimental validation

## Checkpoint 5: Cryptic Feature Scan
- [ ] No cryptic splice sites with MaxEntScan score > 6.0
- [ ] No upstream ATGs in strong Kozak context (unless intentional uORF)
- [ ] No cryptic polyA signals (AATAAA) before intended polyA
- [ ] No obvious secondary structures that might impede transcription/translation

## Checkpoint 6: Restriction Site Scan (Basic)
- [ ] No sites that would interfere with standard cloning
- [ ] Annotate all restriction sites for reference
- [ ] Flag rare restriction sites (8bp cutters) for potential use

## Checkpoint 7: Homology to Mammalian Genome
- [ ] BLAST CDS against human/mouse/rat genomes
- [ ] Flag high-identity regions (>18bp perfect match) that might cause:
  - Off-target integration
  - Homology-based recombination
  - Immune recognition (if endogenous sequence)
