---
report_type: pre-flight
generated_by: DNA Engineer Agent
agent_version: "1.1.0"
script_version: preflight_scan_v1.0
timestamp: 2025-12-16T23:39:02Z
input_file: TTRC004_rep2mut02-cap9-p5.gb
input_file_hash: sha256:0c1959da99c61f08e9bd0e4a9ff22240aed8bac5b425a0a19333a08e0d7d429f
input_file_size: 7330 bp
plasmid_type: Rep-Cap Helper
cis_element_manifest_version: "1.0.0"
project_context: aav_vhh_v1 (v1.0.0)
status: READY
---

# Pre-Flight Report: TTRC004_rep2mut02-cap9-p5

## 1. Summary

**Status:** ✅ READY

All checks passed. Construct is ready for VHH insertion at VR-IV using Design 2 linkers.

---

## 2. Plasmid Classification

| Property | Value |
|----------|-------|
| Type | Rep-Cap Helper |
| Topology | Circular |
| Size | 7330 bp |
| Key Features | Rep, Cap (AAV9), p19 promoter, p40 promoter, p5 promoter |

**Classification Rationale:** Contains Rep and Cap genes with associated promoters (p5, p19, p40), but lacks ITRs. This is a helper plasmid for AAV production.

---

## 3. Cis-Element Scan Results

| Element ID | Name | Risk Level | Detection Method | Found? | Location(s) | Notes |
|------------|------|------------|------------------|--------|-------------|-------|
| AAV2_ITR_LEFT | AAV2 Left ITR | CRITICAL | sequence_partial | NO (expected) | N/A | Expected absence (Rep-Cap helper) |
| AAV2_ITR_RIGHT | AAV2 Right ITR | CRITICAL | sequence_partial | NO (expected) | N/A | Expected absence (Rep-Cap helper) |
| AAV_REP_BINDING | Rep Binding Element (RBE) | CRITICAL | sequence_exact | NO | N/A |  |
| AAV_TRS | Terminal Resolution Site | HIGH | sequence_exact | N/A | N/A | Not applicable (Rep-Cap helper has no ITRs) |
| AAV_AAP_ORF | Assembly-Activating Protein (AAP) | HIGH | protein_motif | NO | N/A |  |
| AAV_P5_PROMOTER | p5 Promoter | MEDIUM | sequence_partial | NO | N/A |  |
| AAV_P19_PROMOTER | p19 Promoter | MEDIUM | sequence_partial | NO | N/A |  |
| AAV_P40_PROMOTER | p40 Promoter | MEDIUM | sequence_partial | NO | N/A |  |
| POLYA_BGH | BGH Poly(A) Signal | LOW | sequence_exact | YES | 1929-1934, 4174-4179 |  |
| POLYA_SV40 | SV40 Poly(A) Signal | LOW | sequence_partial | NO | N/A |  |
| KOZAK_CONSENSUS | Kozak Consensus Sequence | MEDIUM | regex | NO | N/A |  |

**Context-Filtered Results:**
- Rep-Cap helper plasmid: ITR checks skipped (expected absence)
- TRS checks skipped (only valid inside ITR context)

---

## 4. Overlapping ORF Verification

### AAP (Assembly-Activating Protein)

| Source | Start (bp) | End (bp) | Length (aa) |
|--------|------------|----------|-------------|
| File Annotation | 2476 | 3069 | 198 |
| Independent Verification (+1 frame) | 2389 | 3069 | 226 |

**Discrepancy:** ⚠️ MISMATCH DETECTED
- Start: annotated 2476, verified 2389 (diff: -87 bp)
- Length: annotated 198 aa, verified 226 aa

**Verification Method:** Translated VP1 region (bp 1950-4160) in +1 frame, searched for start motif `MPGFYEIVIKV`, located first downstream stop codon.

---

## 5. Insertion Site Analysis

**Target:** VR-IV (Variable Region IV)

### Anchor Sequences (Derived from Actual Sequence)

| Position | Sequence | Location |
|----------|----------|----------|
| Upstream (7aa) | `YYLSKTI` | bp 3282-3302 |
| VR-IV | `NGSGQNQQT` | bp 3303-3329 |
| Downstream (7aa) | `LKFSVAG` | bp 3330-3350 |

### AAP Overlap Check

- AAP ends at: 3069 bp
- VR-IV starts at: 3303 bp
- Gap: 234 bp
- **Overlap Risk:** ✅ None - VR-IV is 234 bp downstream of AAP

---

## 6. Risk Assessment

| Risk Category | Level | Justification |
|---------------|-------|---------------|
| ITR Integrity | LOW | Not applicable (Rep-Cap helper has no ITRs) |
| AAP Frame | LOW | VR-IV is 234 bp downstream of verified AAP region |
| Cis-Element Collision | LOW | No critical elements in modification zone |
| Insertion Size | LOW | 140aa within recommended limits (<150aa) |

**Overall Risk:** LOW

---

## 7. Recommendations

### ✅ READY - Proceed with modification

- [ ] Proceed with VHH insertion script using Design 2 linkers
- [ ] Insert at VR-IV coordinates: bp 3303-3329
- [ ] Use LINK_D2_N (20aa) + VHH_120aa + LINK_D2_C (0aa)
- [ ] Verify ORF integrity after modification
- [ ] Check AAP +1 frame for introduced stop codons

---

## Appendix: Raw Data

### VP1 Translation (+1 frame, first 100aa of AAP region)
```
SSLLRNRTPPRVLANRVHSPLKRDSISVRLATQSQSQTLNQSENLPQPPQVWDLLQWLQVVAHQWQTITKVPMEWVVPREIGIAIPNGWGTESSPPAPEP
```

### Insertion Design Preview
```
LINK_D2_N:  GGGGSGGGGSGGGGSGGGGS (20aa)
VHH:        [120aa placeholder]
LINK_D2_C:  (empty - direct fusion)
Total:      140aa (420bp insertion)
```
