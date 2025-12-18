# REVISED: Comprehensive Restriction Site Analysis - 6 Regions

**Plasmid:** BASE-DRAFT-AAV9-RepCap-NOCAP.gb
**Analysis Date:** 2025-12-17
**Version:** 2.3.0 (Region 1 REVISED)
**Plasmid Length:** 7,078 bp

---

## Revision Notes

**Region 1 has been revised** to avoid the splice acceptor in the VP1 unique N-terminus:

- **Previous Region 1:** 2365-2775 (VP1 start to VP2 start) - EagI @ 2369
- **Revised Region 1:** 2415-2775 (Rep68 stop to VP2 start) - AvrII @ 2459

The revised Region 1 is narrowed by 50 bp to start after the Rep68 stop codon (position 2414),
avoiding potential splice acceptor sites in the VP1 unique N-terminal region.

---

## Executive Summary

This analysis identifies unique restriction enzyme sites (requiring ≤2 silent mutations) in six
strategically chosen regions of an AAV9 Rep-Cap helper plasmid. These sites enable future
modifications while avoiding overlap with critical boundaries (start codons, variable regions).

**Key Findings:**
- **6 unique restriction sites** identified (1 per region)
- **Total mutations required:** 5 silent mutations + 1 site already present
- **All sites verified** to not overlap critical boundaries
- **All sites are unique** in the entire 7,078 bp plasmid
- **No amino acid changes** in any proposed modification

---

## Summary Table

| Region | Enzyme | Position | Recognition Site | Mutations | Status |
|--------|--------|----------|------------------|-----------|--------|
| Region 1 | **AvrII** | 2459 | `CCTAGG` | 1 silent | ✓ |
| Region 2 | **BbvCI** | 2828 | `CCTCAGC` | 1 silent | ✓ |
| Region 3 | **EcoNI** | 3487 | `CCTCAGTAAGG` | 1 silent | ✓ |
| Region 4 | **FseI** | 3759 | `GGCCGGCC` | 1 silent | ✓ |
| Region 5 | **BmtI** | 3937 | `GCTAGC` | 1 silent | ✓ |
| Region 6 | **BaeI** | 4318 | `ACACCTGTACC` | 0 (exists) | ✓ |

---

## Detailed Analysis by Region

### Region 1: Rep68-stop to VP2-start

**Region Boundaries:** 2415-2775 (361 bp)
**Description:** After Rep68 stop, before VP2 start

**⚠️ REVISED:** This region was narrowed from the original analysis to start after Rep68 stop codon,
avoiding the splice acceptor in the VP1 unique N-terminus.

#### Recommended Site

| Parameter | Value |
|-----------|-------|
| **Enzyme** | AvrII |
| **Recognition Sequence** | CCTAGG |
| **Position** | 2459 (1-indexed) |
| **Mutations Required** | 1 |
| **Original Sequence** | `CCCAGG` |
| **Modified Sequence** | `CCTAGG` |
| **Uniqueness** | ✓ Unique in entire plasmid |
| **Boundary Check** | ✓ No overlaps |

#### Codon Changes

```
Position 2461: CCT → CCC (P → P) ✓ SILENT
```

**Mutation Details:**
- A3T

#### Sequence Context

```
Position 2424-2500
CGAGTGGTGGGCTTTGAAACCTGGAGCCCCTCAACCCAAGGCAAATCAACAACATCAAGACAACGCTCGAGGTCTT
                                   ^^^^^^
AvrII site @ 2459
```

---

### Region 2: VP2-AAP Intergenic

**Region Boundaries:** 2779-2890 (112 bp)
**Description:** After VP2 start, before AAP start

#### Recommended Site

| Parameter | Value |
|-----------|-------|
| **Enzyme** | BbvCI |
| **Recognition Sequence** | CCTCAGC |
| **Position** | 2828 (1-indexed) |
| **Mutations Required** | 1 |
| **Original Sequence** | `CCTCCGC` |
| **Modified Sequence** | `CCTCAGC` |
| **Uniqueness** | ✓ Unique in entire plasmid |
| **Boundary Check** | ✓ No overlaps |

#### Codon Changes

```
Position 2832: TCC → TCA (S → S) ✓ SILENT
```

**Mutation Details:**
- C5A

#### Sequence Context

```
Position 2793-2870
GAGGCCTGTAGAGCAGTCTCCTCAGGAACCGGACTCCTCCGCGGGTATTGGCAAATCGGGTGCACAGCCCGCTAAAA
                                   ^^^^^^^
BbvCI site @ 2828
```

---

### Region 3: AAP-stop to VR4

**Region Boundaries:** 3485-3717 (233 bp)
**Description:** Between AAP stop and VR4 start

#### Recommended Site

| Parameter | Value |
|-----------|-------|
| **Enzyme** | EcoNI |
| **Recognition Sequence** | CCTCAGTAAGG |
| **Position** | 3487 (1-indexed) |
| **Mutations Required** | 1 |
| **Original Sequence** | `CCTCAGTACGG` |
| **Modified Sequence** | `CCTCAGTAAGG` |
| **Uniqueness** | ✓ Unique in entire plasmid |
| **Boundary Check** | ✓ No overlaps |

#### Codon Changes

```
Position 3495-3497: CGG → AGG (R → R) ✓ SILENT
```

**Mutation Details:**
- C9A

#### Sequence Context

```
Position 3452-3533
GCCTCCCGCCGTTCCCAGCGGACGTTTTCATGATTCCTCAGTACGGGTATCTGACGCTTAATGATGGAAGCCAGGCCGTGG
                                   ^^^^^^^^^^^
EcoNI site @ 3487
```

---

### Region 4: VR4 to VR5

**Region Boundaries:** 3745-3825 (81 bp)
**Description:** Between VR4 end and VR5 start

#### Recommended Site

| Parameter | Value |
|-----------|-------|
| **Enzyme** | FseI |
| **Recognition Sequence** | GGCCGGCC |
| **Position** | 3759 (1-indexed) |
| **Mutations Required** | 1 |
| **Original Sequence** | `GGCCGGAC` |
| **Modified Sequence** | `GGCCGGCC` |
| **Uniqueness** | ✓ Unique in entire plasmid |
| **Boundary Check** | ✓ No overlaps |

#### Codon Changes

```
Position 3765: GGA → GGC (G → G) ✓ SILENT
```

**Mutation Details:**
- A7C

#### Sequence Context

```
Position 3724-3802
TCTGGACAGAATCAACAAACGCTAAAATTCAGTGTGGCCGGACCCAGCAACATGGCTGTCCAGGGAAGAAACTACATA
                                   ^^^^^^^^
FseI site @ 3759
```

---

### Region 5: VR5 to VR8

**Region Boundaries:** 3880-4104 (225 bp)
**Description:** Between VR5 end and VR8 start

#### Recommended Site

| Parameter | Value |
|-----------|-------|
| **Enzyme** | BmtI |
| **Recognition Sequence** | GCTAGC |
| **Position** | 3937 (1-indexed) |
| **Mutations Required** | 1 |
| **Original Sequence** | `GCCAGC` |
| **Modified Sequence** | `GCTAGC` |
| **Uniqueness** | ✓ Unique in entire plasmid |
| **Boundary Check** | ✓ No overlaps |

#### Codon Changes

```
Position 3939: GCC → GCT (A → A) ✓ SILENT
```

**Mutation Details:**
- C3T

#### Sequence Context

```
Position 3902-3978
GACGTAATAGCTTGATGAATCCTGGACCTGCTATGGCCAGCCACAAAGAAGGAGAGGACCGTTTCTTTCCTTTGTC
                                   ^^^^^^
BmtI site @ 3937
```

---

### Region 6: Post-VR8

**Region Boundaries:** 4144-4575 (432 bp)
**Description:** After VR8 end

#### Recommended Site

| Parameter | Value |
|-----------|-------|
| **Enzyme** | BaeI |
| **Recognition Sequence** | ACACCTGTACC |
| **Position** | 4318 (1-indexed) |
| **Mutations Required** | 0 |
| **Original Sequence** | `ACACCTGTACC` |
| **Modified Sequence** | `ACACCTGTACC` |
| **Uniqueness** | ✓ Unique in entire plasmid |
| **Boundary Check** | ✓ No overlaps |

**No mutations required** - site already present in the plasmid.

#### Sequence Context

```
Position 4283-4364
TGAAGCACCCGCCTCCTCAGATCCTCATCAAAAACACACCTGTACCTGCGGATCCTCCAACGGCCTTCAACAAGGACAAGC
                                   ^^^^^^^^^^^
BaeI site @ 4318
```

---

## Applications

These restriction sites provide excellent handles for:

1. **Variable Region Engineering**
   - Insert VHH nanobodies at VR4/VR5 using flanking sites
   - Swap entire variable region cassettes
   - Region 3 (AAP-stop to VR4) enables insertions just before VR4
   - Region 4 (VR4 to VR5) enables dual insertions at VR4 and VR5

2. **Epitope Tagging**
   - Add FLAG, HA, or His tags in VP1 unique region
   - Insert fluorescent proteins for trafficking studies

3. **Golden Gate Assembly** (BbvCI in Region 2)
   - Scarless cloning workflows
   - Modular assembly of Cap variants

4. **Diagnostic Restriction Mapping**
   - Quick QC after cloning
   - Verify plasmid identity

---

## Conclusions

This revised analysis successfully identified **6 unique restriction sites** across strategic regions
of the AAV9 Rep-Cap plasmid:

✓ **All sites unique** in entire plasmid
✓ **All mutations silent** (no amino acid changes)
✓ **No boundary overlaps** (verified for 11 critical boundaries)
✓ **Minimal mutations** (5 total across 6 sites, 1 site exists naturally)
✓ **Strategically positioned** for VR engineering
✓ **Region 1 revised** to avoid splice acceptor in VP1 N-terminus

These sites provide excellent molecular handles for:
- VHH display at variable regions
- Epitope tagging in VP1 unique domain
- Modular assembly via Golden Gate
- Diagnostic restriction mapping

All modifications are low-risk and preserve the functional architecture of the Cap gene.

---

**Analysis completed successfully.**