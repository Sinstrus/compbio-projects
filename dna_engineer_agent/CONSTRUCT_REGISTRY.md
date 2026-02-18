# Construct Registry

Authoritative manifest of all AVD-series DNA constructs.
All files are stored in `constructs/` (flat directory, no subdirectories).

**Last updated:** 2026-02-18
**Total constructs:** 101 active + 2 archived

---

## Parent Plasmids (AVD001-004)

| ID | Filename | Description |
|----|----------|-------------|
| AVD001 | AVD001-Rep2Mut2Cap9-6R-NoCap.dna | Rep2Mut2-Cap9 backbone with six unique restriction sites engineered in; capsid ORFs removed. Negative control — no VP expression. |
| AVD002 | AVD002-Rep2Mut2Cap9-6R-wt.dna | Wild-type Cap9 on the 6R backbone. Primary build template for all VHH display and splicing constructs (7,104 bp). |
| AVD002 | AVD002-Rep2Mut2Cap9-6R-wt-5fold.gb | AVD002 variant with 5-fold axis residues annotated. GenBank format for downstream build scripts. |
| AVD002 | AVD002-Rep2Mut2Cap9-6R-wt-5fold-hoffmann.gb | AVD002 variant using Hoffmann numbering convention for 5-fold axis residues. |
| AVD003 | AVD003-EF1A-VPall-bGH.dna | Standalone VP1/VP2/VP3 expression cassette under EF1a promoter with bGH polyA. Trans-complementation helper for mosaic capsid experiments. |
| AVD004 | AVD004-ssAAV-EF1A-VPall-bGH.dna | ITR-flanked single-stranded AAV version of AVD003. Experimental packaging variant for trans-complementation. |

## Anti-ALPL VHH3 Display (AVD005-011)

| ID | Filename | Description |
|----|----------|-------------|
| AVD005 | AVD005-Rep2Mut2Cap9-VP1ko.gb | Identical to AVD002 except VP1 start ATG>AAG. No VHH insert. VP1 knockout control confirming VP2/VP3-only capsids assemble. |
| AVD006 | AVD006-Rep2Mut2Cap9-VP1-VHH3-VR4.gb | VHH3 anti-ALPL at VR-IV (pos 456) with D2 asymmetric linkers (5xGGGGS N-term, 1xGGGGS C-term). VP1 ATG>AAG; VP2/VP3 intact. Primary VR-IV display construct. |
| AVD007 | AVD007-Rep2Mut2Cap9-VP1-VHH3-D2.gb | Same VR-IV insertion as AVD006 plus VP2 knockout (ACG>ACC) and VP3 knockout (ATG>CTG). VP1-only expression for mosaic capsid production with helper. |
| AVD008 | AVD008-Rep2Mut2Cap9-VP1n-VHH3.gb | VHH3 fused at VP1 N-terminus after MA dipeptide with 5xGGGGS C-terminal linker. VP2/VP3 intact. Tests N-terminal display topology. |
| AVD009 | AVD009-Rep2Mut2Cap9-VP1n-VHH3-D2.gb | Same VP1 N-terminal fusion as AVD008, plus VP2/VP3 knockouts. VP1-only mosaic variant. |
| AVD010 | AVD010-Rep2Mut2Cap9-VP2n-VHH3.gb | VHH3 fused at VP2 N-terminus after MA dipeptide. Tests VP2-mediated display. VP1/VP3 intact. |
| AVD011 | AVD011-Rep2Mut2Cap9-VP2ko.gb | Identical to AVD002 except VP2 start ACG>ACC. No VHH insert. VP2 knockout control confirming VP1+VP3 capsids assemble. |

## Anti-TfR1 VHH Display (AVD012-026)

Three clones (TfR-A-SS, TfR-A-PP, TfR-B) x five insertion sites (VR4, D2, VP1n, VP1n-D2, VP2n). Same architecture as AVD006-010 but with anti-TfR1 nanobodies.

| ID | Filename | Description |
|----|----------|-------------|
| AVD012 | AVD012-Rep2Mut2Cap9-VP1-TfR-A-SS-VR4.gb | Anti-TfR1 clone A (short spacer) at VR-IV. VP1 ATG>AAG. Same insertion architecture as AVD006. |
| AVD013 | AVD013-Rep2Mut2Cap9-VP1-TfR-A-SS-D2.gb | Clone A-SS at VR-IV with VP1/VP2/VP3 all knocked out. Mosaic variant like AVD007. |
| AVD014 | AVD014-Rep2Mut2Cap9-VP1n-TfR-A-SS.gb | Clone A-SS at VP1 N-terminus. Like AVD008 with TfR-A-SS instead of VHH3. |
| AVD015 | AVD015-Rep2Mut2Cap9-VP1n-TfR-A-SS-D2.gb | Clone A-SS at VP1 N-term, VP2/VP3 KO. Like AVD009. |
| AVD016 | AVD016-Rep2Mut2Cap9-VP2n-TfR-A-SS.gb | Clone A-SS at VP2 N-terminus. Like AVD010. |
| AVD017 | AVD017-Rep2Mut2Cap9-VP1-TfR-A-PP-VR4.gb | Anti-TfR1 clone A (proline peptide variant) at VR-IV. VP1 ATG>AAG. |
| AVD018 | AVD018-Rep2Mut2Cap9-VP1-TfR-A-PP-D2.gb | Clone A-PP at VR-IV with all VP knockouts. |
| AVD019 | AVD019-Rep2Mut2Cap9-VP1n-TfR-A-PP.gb | Clone A-PP at VP1 N-terminus. |
| AVD020 | AVD020-Rep2Mut2Cap9-VP1n-TfR-A-PP-D2.gb | Clone A-PP at VP1 N-term, VP2/VP3 KO. |
| AVD021 | AVD021-Rep2Mut2Cap9-VP2n-TfR-A-PP.gb | Clone A-PP at VP2 N-terminus. |
| AVD022 | AVD022-Rep2Mut2Cap9-VP1-TfR-B-VR4.gb | Anti-TfR1 clone B at VR-IV. VP1 ATG>AAG. |
| AVD023 | AVD023-Rep2Mut2Cap9-VP1-TfR-B-D2.gb | Clone B at VR-IV with all VP knockouts. |
| AVD024 | AVD024-Rep2Mut2Cap9-VP1n-TfR-B.gb | Clone B at VP1 N-terminus. |
| AVD025 | AVD025-Rep2Mut2Cap9-VP1n-TfR-B-D2.gb | Clone B at VP1 N-term, VP2/VP3 KO. |
| AVD026 | AVD026-Rep2Mut2Cap9-VP2n-TfR-B.gb | Clone B at VP2 N-terminus. |

## Splicing Reporters (AVD027-058)

32 dual-reporter plasmids on TRX155 backbone testing splice donor strength across 7 tiers. Each has a unique 11-mer 5' splice donor. Intron retained = NanoLuc active; intron spliced = Firefly Luc active.

| ID | Filename | Description |
|----|----------|-------------|
| AVD027 | AVD027-SplicingReporter-ACGGTAAAGCA-Tier1.gb | Tier 1 (strongest H-bond). Splice donor ACGGTAAAGCA. Dual Nluc/Fluc frame-shift reporter with VHH3 intron body. |
| AVD028 | AVD028-SplicingReporter-GTGGTAAGGTG-Tier1.gb | Tier 1. Splice donor GTGGTAAGGTG. |
| AVD029 | AVD029-SplicingReporter-ACGGTAAATCA-Tier1.gb | Tier 1. Splice donor ACGGTAAATCA. |
| AVD030 | AVD030-SplicingReporter-ACGGTAAACCG-Tier1.gb | Tier 1. Splice donor ACGGTAAACCG. |
| AVD031 | AVD031-SplicingReporter-ACGGTAAAGCC-Tier2.gb | Tier 2. Splice donor ACGGTAAAGCC. |
| AVD032 | AVD032-SplicingReporter-GTGGTAAGGTA-Tier2.gb | Tier 2. Splice donor GTGGTAAGGTA. |
| AVD033 | AVD033-SplicingReporter-AAGGTAAATCC-Tier2.gb | Tier 2. Splice donor AAGGTAAATCC. |
| AVD034 | AVD034-SplicingReporter-CCGGTAAACAA-Tier2.gb | Tier 2. Splice donor CCGGTAAACAA. |
| AVD035 | AVD035-SplicingReporter-ACGGTAAAGGC-Tier3.gb | Tier 3. Splice donor ACGGTAAAGGC. |
| AVD036 | AVD036-SplicingReporter-TTGGTAAGGTG-Tier3.gb | Tier 3. Splice donor TTGGTAAGGTG. |
| AVD037 | AVD037-SplicingReporter-ACGGTAAATGG-Tier3.gb | Tier 3. Splice donor ACGGTAAATGG. |
| AVD038 | AVD038-SplicingReporter-ACGGTAACCGG-Tier3.gb | Tier 3. Splice donor ACGGTAACCGG. |
| AVD039 | AVD039-SplicingReporter-ACGGTAAGAGG-Tier3.gb | Tier 3. Splice donor ACGGTAAGAGG. |
| AVD040 | AVD040-SplicingReporter-ACGGTAAAGAA-Tier4.gb | Tier 4. Splice donor ACGGTAAAGAA. |
| AVD041 | AVD041-SplicingReporter-GTGGTAAGGAT-Tier4.gb | Tier 4. Splice donor GTGGTAAGGAT. |
| AVD042 | AVD042-SplicingReporter-AGGGTAATTCC-Tier4.gb | Tier 4. Splice donor AGGGTAATTCC. |
| AVD043 | AVD043-SplicingReporter-AGGGTAATCCT-Tier4.gb | Tier 4. Splice donor AGGGTAATCCT. |
| AVD044 | AVD044-SplicingReporter-AGGGTAATCTC-Tier4.gb | Tier 4. Splice donor AGGGTAATCTC. |
| AVD045 | AVD045-SplicingReporter-CAGGTAAAGCC-Tier5.gb | Tier 5. Splice donor CAGGTAAAGCC. |
| AVD046 | AVD046-SplicingReporter-GTGGTAAGGGG-Tier5.gb | Tier 5. Splice donor GTGGTAAGGGG. |
| AVD047 | AVD047-SplicingReporter-CAGGTAAATAC-Tier5.gb | Tier 5. Splice donor CAGGTAAATAC. |
| AVD048 | AVD048-SplicingReporter-CAGGTAAACAT-Tier5.gb | Tier 5. Splice donor CAGGTAAACAT. |
| AVD049 | AVD049-SplicingReporter-CAGGTAAACGC-Tier5.gb | Tier 5. Splice donor CAGGTAAACGC. |
| AVD050 | AVD050-SplicingReporter-AAGGTAAAGGG-Tier6.gb | Tier 6. Splice donor AAGGTAAAGGG. |
| AVD051 | AVD051-SplicingReporter-TTGGTAATGTT-Tier6.gb | Tier 6. Splice donor TTGGTAATGTT. |
| AVD052 | AVD052-SplicingReporter-AGGGTAATTGG-Tier6.gb | Tier 6. Splice donor AGGGTAATTGG. |
| AVD053 | AVD053-SplicingReporter-GGGGTAATCGG-Tier6.gb | Tier 6. Splice donor GGGGTAATCGG. |
| AVD054 | AVD054-SplicingReporter-CAGGTAATTCA-Tier6.gb | Tier 6. Splice donor CAGGTAATTCA. |
| AVD055 | AVD055-SplicingReporter-CAGGTAATGAG-Tier7.gb | Tier 7 (weakest). Splice donor CAGGTAATGAG. |
| AVD056 | AVD056-SplicingReporter-GGGGTAAGGGA-Tier7.gb | Tier 7. Splice donor GGGGTAAGGGA. |
| AVD057 | AVD057-SplicingReporter-CAGGTAATTGG-Tier7.gb | Tier 7. Splice donor CAGGTAATTGG. |
| AVD058 | AVD058-SplicingReporter-CAGGTAATCAG-Tier7.gb | Tier 7 (weakest). Splice donor CAGGTAATCAG. |

## Intron Retention Cassettes (AVD059-099)

41 constructs on AVD002 backbone inserting gene-derived splice donor/acceptor pairs at VR-IV, flanking GT/AG-depleted VHH3. Spliced form = near-WT capsid (6aa scar). Retained form = full VHH3 display. Selected for predicted log2FC >= 1.0 in HEK293.

| ID | Filename | Description |
|----|----------|-------------|
| AVD059 | AVD059-Rep2Mut2Cap9-VR4-IR-NABP2.gb | Intron retention cassette using NABP2 splice donor/acceptor at VR-IV. |
| AVD060 | AVD060-Rep2Mut2Cap9-VR4-IR-NCSTN.gb | Intron retention cassette using NCSTN splice donor/acceptor at VR-IV. |
| AVD061 | AVD061-Rep2Mut2Cap9-VR4-IR-QARS1.gb | Intron retention cassette using QARS1 splice donor/acceptor at VR-IV. |
| AVD062 | AVD062-Rep2Mut2Cap9-VR4-IR-GPS1.gb | Intron retention cassette using GPS1 splice donor/acceptor at VR-IV. |
| AVD063 | AVD063-Rep2Mut2Cap9-VR4-IR-PSMD12.gb | Intron retention cassette using PSMD12 splice donor/acceptor at VR-IV. |
| AVD064 | AVD064-Rep2Mut2Cap9-VR4-IR-PLOD3.gb | Intron retention cassette using PLOD3 splice donor/acceptor at VR-IV. |
| AVD065 | AVD065-Rep2Mut2Cap9-VR4-IR-ATG9A.gb | Intron retention cassette using ATG9A splice donor/acceptor at VR-IV. |
| AVD066 | AVD066-Rep2Mut2Cap9-VR4-IR-RPS2.gb | Intron retention cassette using RPS2 splice donor/acceptor at VR-IV. |
| AVD067 | AVD067-Rep2Mut2Cap9-VR4-IR-YIPF3.gb | Intron retention cassette using YIPF3 splice donor/acceptor at VR-IV. |
| AVD068 | AVD068-Rep2Mut2Cap9-VR4-IR-SLC6A9.gb | Intron retention cassette using SLC6A9 splice donor/acceptor at VR-IV. |
| AVD069 | AVD069-Rep2Mut2Cap9-VR4-IR-FTSJ1.gb | Intron retention cassette using FTSJ1 splice donor/acceptor at VR-IV. |
| AVD070 | AVD070-Rep2Mut2Cap9-VR4-IR-PSMD4.gb | Intron retention cassette using PSMD4 splice donor/acceptor at VR-IV. |
| AVD071 | AVD071-Rep2Mut2Cap9-VR4-IR-DRAP1.gb | Intron retention cassette using DRAP1 splice donor/acceptor at VR-IV. |
| AVD072 | AVD072-Rep2Mut2Cap9-VR4-IR-ELOA.gb | Intron retention cassette using ELOA splice donor/acceptor at VR-IV. |
| AVD073 | AVD073-Rep2Mut2Cap9-VR4-IR-SPTLC1.gb | Intron retention cassette using SPTLC1 splice donor/acceptor at VR-IV. |
| AVD074 | AVD074-Rep2Mut2Cap9-VR4-IR-RELA.gb | Intron retention cassette using RELA splice donor/acceptor at VR-IV. |
| AVD075 | AVD075-Rep2Mut2Cap9-VR4-IR-SLC25A11.gb | Intron retention cassette using SLC25A11 splice donor/acceptor at VR-IV. |
| AVD076 | AVD076-Rep2Mut2Cap9-VR4-IR-SCAMP3.gb | Intron retention cassette using SCAMP3 splice donor/acceptor at VR-IV. |
| AVD077 | AVD077-Rep2Mut2Cap9-VR4-IR-TGFB1I1.gb | Intron retention cassette using TGFB1I1 splice donor/acceptor at VR-IV. |
| AVD078 | AVD078-Rep2Mut2Cap9-VR4-IR-DDX3X.gb | Intron retention cassette using DDX3X splice donor/acceptor at VR-IV. |
| AVD079 | AVD079-Rep2Mut2Cap9-VR4-IR-UQCRC1.gb | Intron retention cassette using UQCRC1 splice donor/acceptor at VR-IV. |
| AVD080 | AVD080-Rep2Mut2Cap9-VR4-IR-EIF3I.gb | Intron retention cassette using EIF3I splice donor/acceptor at VR-IV. |
| AVD081 | AVD081-Rep2Mut2Cap9-VR4-IR-TIMM17B.gb | Intron retention cassette using TIMM17B splice donor/acceptor at VR-IV. |
| AVD082 | AVD082-Rep2Mut2Cap9-VR4-IR-ANKRD28.gb | Intron retention cassette using ANKRD28 splice donor/acceptor at VR-IV. |
| AVD083 | AVD083-Rep2Mut2Cap9-VR4-IR-DDX23.gb | Intron retention cassette using DDX23 splice donor/acceptor at VR-IV. |
| AVD084 | AVD084-Rep2Mut2Cap9-VR4-IR-PRRC2A.gb | Intron retention cassette using PRRC2A splice donor/acceptor at VR-IV. |
| AVD085 | AVD085-Rep2Mut2Cap9-VR4-IR-ATP6V0B.gb | Intron retention cassette using ATP6V0B splice donor/acceptor at VR-IV. |
| AVD086 | AVD086-Rep2Mut2Cap9-VR4-IR-MTA1.gb | Intron retention cassette using MTA1 splice donor/acceptor at VR-IV. |
| AVD087 | AVD087-Rep2Mut2Cap9-VR4-IR-PIH1D1.gb | Intron retention cassette using PIH1D1 splice donor/acceptor at VR-IV. |
| AVD088 | AVD088-Rep2Mut2Cap9-VR4-IR-PABPC4.gb | Intron retention cassette using PABPC4 splice donor/acceptor at VR-IV. |
| AVD089 | AVD089-Rep2Mut2Cap9-VR4-IR-FLII.gb | Intron retention cassette using FLII splice donor/acceptor at VR-IV. |
| AVD090 | AVD090-Rep2Mut2Cap9-VR4-IR-OTUB1.gb | Intron retention cassette using OTUB1 splice donor/acceptor at VR-IV. |
| AVD091 | AVD091-Rep2Mut2Cap9-VR4-IR-ACADVL.gb | Intron retention cassette using ACADVL splice donor/acceptor at VR-IV. |
| AVD092 | AVD092-Rep2Mut2Cap9-VR4-IR-CTSA.gb | Intron retention cassette using CTSA splice donor/acceptor at VR-IV. |
| AVD093 | AVD093-Rep2Mut2Cap9-VR4-IR-RPL13.gb | Intron retention cassette using RPL13 splice donor/acceptor at VR-IV. |
| AVD094 | AVD094-Rep2Mut2Cap9-VR4-IR-MT2A.gb | Intron retention cassette using MT2A splice donor/acceptor at VR-IV. |
| AVD095 | AVD095-Rep2Mut2Cap9-VR4-IR-NUCB1.gb | Intron retention cassette using NUCB1 splice donor/acceptor at VR-IV. |
| AVD096 | AVD096-Rep2Mut2Cap9-VR4-IR-3B.gb | Intron retention cassette using 3B splice donor/acceptor at VR-IV. |
| AVD097 | AVD097-Rep2Mut2Cap9-VR4-IR-RPS6KB2.gb | Intron retention cassette using RPS6KB2 splice donor/acceptor at VR-IV. |
| AVD098 | AVD098-Rep2Mut2Cap9-VR4-IR-POLD2.gb | Intron retention cassette using POLD2 splice donor/acceptor at VR-IV. |
| AVD099 | AVD099-Rep2Mut2Cap9-VR4-IR-KRI1.gb | Intron retention cassette using KRI1 splice donor/acceptor at VR-IV. |

## Archived (Superseded)

| ID | Filename | Description |
|----|----------|-------------|
| AVD005 (v1) | archived/AVD005-EF1A-VP1-VHH3-ALPL-bGH.gb | Original VP1-VHH3-ALPL design on EF1a backbone. Superseded by AVD005 (VP1ko on Rep2Mut2 backbone). |
| AVD006 (v1) | archived/AVD006-Rep2Mut2Cap9-VP1-VHH3-ALPL.gb | Original VHH3 display design before linker architecture revision. Superseded by current AVD006 with D2 asymmetric linkers. |
