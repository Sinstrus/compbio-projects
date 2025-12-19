---
report_type: pre-flight
generated_by: DNA Engineer Agent
agent_version: "1.1.0"
script_version: preflight_scan_v1.0
timestamp: YYYY-MM-DDTHH:MM:SSZ
input_file: <filename>.gb
input_file_hash: sha256:<hash>
input_file_size: <N> bp
plasmid_type: <Transfer | Rep-Cap Helper | Expression | Other>
cis_element_manifest_version: "1.0.0"
project_context: <project_name> (v<version>) | null
status: <READY | BLOCKED | WARNINGS>
---

# Pre-Flight Report: <CONSTRUCT_NAME>

## 1. Summary

**Status:** ✅ READY | ❌ BLOCKED | ⚠️ WARNINGS

<One-line summary of critical findings>

---

## 2. Plasmid Classification

| Property | Value |
|----------|-------|
| Type | <Transfer / Rep-Cap Helper / Expression / Other> |
| Topology | <Circular / Linear> |
| Size | <N> bp |
| Key Features | <Rep, Cap, ITRs, Transgene, etc.> |

**Classification Rationale:** <Why this plasmid was classified this way>

---

## 3. Cis-Element Scan Results

| Element ID | Name | Risk Level | Detection Method | Found? | Location(s) | Notes |
|------------|------|------------|------------------|--------|-------------|-------|
| AAV2_ITR_LEFT | Left ITR | CRITICAL | sequence_partial | ✅/❌ | N/A or bp range | |
| AAV_AAP_ORF | AAP | HIGH | protein_motif | ✅/❌ | bp range | |
| ... | ... | ... | ... | ... | ... | ... |

**Context-Filtered Results:**
- <N> spurious TRS hits outside ITR context → IGNORED
- <Other context-based filtering decisions>

---

## 4. Overlapping ORF Verification

### AAP (Assembly-Activating Protein)

| Source | Start (bp) | End (bp) | Length (aa) |
|--------|------------|----------|-------------|
| File Annotation | <N> | <N> | <N> |
| Independent Verification (+1 frame) | <N> | <N> | <N> |

**Discrepancy:** ✅ None | ⚠️ <Description of difference>

**Verification Method:** Translated VP1 region in +1 frame, searched for start motif `MPGFYEIVIKVPSD`, located first downstream stop codon.

---

## 5. Insertion Site Analysis

**Target:** <VR-IV / VR-VIII / Other>

### Anchor Sequences (Derived from Actual Sequence)

| Position | Sequence | VP1 Residue |
|----------|----------|-------------|
| Upstream | <AA sequence> | <N>-<N> |
| Target Region | <AA sequence> | <N>-<N> |
| Downstream | <AA sequence> | <N>-<N> |

### AAP Overlap Check

- AAP ends at: <N> bp
- Insertion site starts at: <N> bp
- Gap: <N> bp
- **Overlap Risk:** ✅ None | ⚠️ <Description>

---

## 6. Risk Assessment

| Risk Category | Level | Justification |
|---------------|-------|---------------|
| ITR Integrity | ✅ Low | Not applicable (Rep-Cap helper) / ITRs distant from modification |
| AAP Frame | ✅ Low / ⚠️ Medium / ❌ High | <Reason> |
| Cis-Element Collision | ✅ Low | No critical elements in modification zone |
| Insertion Size | ✅ Low / ⚠️ Medium | <N> aa within / exceeds recommended limits |

**Overall Risk:** ✅ LOW | ⚠️ MEDIUM | ❌ HIGH

---

## 7. Recommendations

### If READY:
- [ ] Proceed with modification script
- [ ] Verify insertion coordinates match plan
- [ ] Run bio-sanity checks after modification

### If BLOCKED:
- [ ] <Specific issue that must be resolved>
- [ ] <Action required before proceeding>

### If WARNINGS:
- [ ] <Issue requiring user acknowledgment>
- [ ] <Recommended verification step>

---

## Appendix: Raw Data

<Optional: Include translated sequences, coordinate calculations, or other supporting data that informed the analysis>
