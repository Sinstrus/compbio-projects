# DNA Engineer Agent — Lessons Learned

This document captures critical bugs, design issues, and best practices discovered during agent-assisted DNA engineering projects.

---

## BUG-005: Parent Sequence Mismatch in Multi-Construct Builds

**Severity:** CRITICAL
**Discovered:** 2026-01-14
**Project:** AVD005/AVD006 VP1-VHH fusion design

### The Problem

When building two related constructs (AVD005 from AVD003, AVD006 from AVD002), the build script extracted VP1 from AVD003 and used it for BOTH constructs:

```python
# WRONG: Uses AVD003's VP1 for both constructs
vp1_original = str(avd003.seq[1520:3731])
vp1_vhh_fusion = build_vp1_vhh_fusion(vp1_original)
avd005 = create_avd005(avd003, vp1_vhh_fusion)
avd006 = create_avd006(avd002, vp1_vhh_fusion)  # WRONG! Should use AVD002's VP1
```

### The Consequence

AVD003 and AVD002 had 7 pre-existing silent mutations in their VP1 sequences (from prior 6R site engineering). When AVD003's VP1 was inserted into AVD006, those 7 differences appeared as "unexplained" mutations.

Comparing AVD002 → AVD006 showed **9 mutations** instead of the expected **2**:
- 2 intentional (VP2/VP3 knockouts)
- 7 unintentional (inherited from AVD003's VP1)

### The Fix

Each construct must use the VP1 sequence from its OWN parent plasmid:

```python
# CORRECT: Each construct uses VP1 from its own parent
vp1_from_avd003 = str(avd003.seq[1520:3731])
vp1_from_avd002 = str(avd002.seq[2378:4589])

vp1_vhh_for_avd005 = build_vp1_vhh_fusion(vp1_from_avd003)
vp1_vhh_for_avd006 = build_vp1_vhh_fusion(vp1_from_avd002)

avd005 = create_avd005(avd003, vp1_vhh_for_avd005)
avd006 = create_avd006(avd002, vp1_vhh_for_avd006)
```

### The Lesson

**ALWAYS start from the actual parent sequence, not a reference sequence from a different plasmid.**

Different plasmids of the same "type" may have:
- Silent mutations from codon optimization
- Restriction site modifications (like 6R engineering)
- Vendor-specific sequence variants
- Historical sequence drift

### Prevention Checklist

- [ ] When building construct B from parent A, extract the modified region from A's actual sequence
- [ ] Never assume two "similar" plasmids have identical sequences
- [ ] After building, perform a full diff between parent and child to verify ONLY intentional changes exist
- [ ] If the diff shows unexpected mutations, investigate the source before proceeding

---

## DESIGN-005: Synthetic Fragment Boundary Selection

**Type:** Best Practice
**Discovered:** 2026-01-14
**Project:** AVD005/AVD006 synthetic fragment ordering

### The Problem

Initial fragment boundaries used restriction sites that were INSIDE the modified VP1-VHH region:
- BsrGI at position 3353 = 1,833 bp INTO VP1-VHH
- BmtI at position 3510 = 1,990 bp INTO VP1-VHH
- BstZ17I at position 3736 = 2,216 bp INTO VP1-VHH

These sites were from the 6R engineering panel and move WITH the VP1 sequence. Using them as cloning boundaries would create incomplete fragments.

### The Solution

Find the closest unique restriction sites that FLANK (are OUTSIDE) the modified region:

1. **Identify the changed region:**
   - Compare parent vs child sequence
   - Find the first and last nucleotide that differs
   - This defines the "modified region"

2. **Find flanking sites:**
   - Search UPSTREAM of first change for unique restriction sites
   - Search DOWNSTREAM of last change for unique restriction sites
   - Choose the closest unique sites

3. **For AVD005:**
   - Modified region: VP2 knockout (bp 413), VP3 knockout (bp 606), VHH insertion (bp 1368-1785)
   - Upstream flank: AvrII (CCTAGG) at position 1676
   - Downstream flank: BsrGI (TGTACA) at position 3358
   - Fragment: 1,683 bp

4. **For AVD006:**
   - Same modified region within VP1
   - Upstream flank: HindIII (AAGCTT) at position 2058
   - Downstream flank: BsrGI (TGTACA) at position 4216
   - Fragment: 2,159 bp

### The Lesson

**Synthetic fragment boundaries must use sites that FLANK the modified region, not sites within it.**

### Prevention Checklist

- [ ] Before selecting cloning sites, map ALL restriction sites in the construct
- [ ] Identify which sites are INSIDE the modified region vs OUTSIDE
- [ ] Choose the closest UNIQUE sites that are OUTSIDE the modified region
- [ ] Verify the fragment includes ALL nucleotide changes
- [ ] Smaller fragments are cheaper — don't extend unnecessarily

---

## DESIGN-006: VP1/VP2/VP3 Knockout Strategy

**Type:** Design Pattern
**Discovered:** 2026-01-14
**Project:** AVD005/AVD006 VP1-only expression

### Background

AAV capsid genes encode three proteins (VP1, VP2, VP3) from the same mRNA via alternative start codons:
- VP1: Full length, starts at first ATG
- VP2: Shorter, starts at internal ACG (non-canonical start)
- VP3: Shortest, starts at internal ATG

For VP1-only expression (e.g., VHH display), VP2 and VP3 must be knocked out while preserving VP1.

### The Strategy

**VP2 Knockout (Silent):**
- VP2 starts at ACG (encodes Thr) at codon 138 of VP1
- Mutation: ACG → ACC (both encode Thr)
- Effect: Removes VP2 start codon context, still silent in VP1

**VP3 Knockout (Non-Silent):**
- VP3 starts at ATG (encodes Met) at codon 203 of VP1
- Mutation: ATG → CTG (Met → Leu)
- Effect: Eliminates VP3 start codon, changes VP1 amino acid

### Important Notes

1. **VP2 knockout IS silent** — the mutation ACG→ACC doesn't change the amino acid in VP1
2. **VP3 knockout is NOT silent** — Met→Leu changes VP1 at position 203
3. **Position is relative to VP1 start** — codon numbering starts at VP1 ATG
4. **AAP reading frame is preserved** — mutations are in VP1 frame, not +1 frame

### Prevention Checklist

- [ ] Verify VP2 knockout mutation is silent (codon 138: ACG→ACC = Thr→Thr)
- [ ] Acknowledge VP3 knockout is non-silent (codon 203: ATG→CTG = Met→Leu)
- [ ] Check that AAP reading frame (+1 relative to VP1) is not disrupted
- [ ] Confirm VP1 start codon (position 1) is preserved

---

## DESIGN-007: VHH Insertion at Variable Regions with Asymmetric Linkers

**Type:** Design Pattern
**Discovered:** 2026-01-14
**Project:** AVD005/AVD006 anti-ALPL VHH display
**Reference:** See `docs/DESIGN_PATTERN_Asymmetric_Linkers.md` for full analysis

### Background

AAV capsid has 9 variable regions (VR-I through VR-IX) that tolerate insertions:
- **VR-IV (aa 452-460):** Located at 3-fold spike apex - **OPTIMAL for VHH display**
- VR-VIII (aa 586-591): Also used but disrupts HSPG binding in AAV2

### Design Parameters Used

**Insertion Site:**
- Location: VR-IV at amino acid 456 of VP1
- Position: bp 1368 (456 × 3) from VP1 start
- Rationale: Highest solvent exposure, best receptor engagement

**Linker Design (D2 = Asymmetric Flexible) - FROM BIOGEN PATENT:**
- **N-terminal linker:** (GGGGS)×4 = 20 amino acids = 60 bp
- **C-terminal linker:** (GGGGS)×1 = 5 amino acids = 15 bp (SHORT, not direct fusion)
- **Total linker contribution:** 75 bp
- **Asymmetry rationale:**
  - Long N-terminal = rotational freedom ("flagpole" dynamics)
  - Short C-terminal = structural anchor (nucleates folding)
  - Empirically validated: >10-fold enhancement in CNS transduction vs symmetric designs

**VHH Insert:**
- Target: ALPL (Alkaline Phosphatase)
- VHH clone: VHH3 anti-ALPL
- Size: 118 amino acids = 354 bp (dnachisel-optimized for synthesis)

**Total Insertion:**
- 432 bp = 60 bp (N-linker) + 354 bp (VHH) + 15 bp (C-linker)
- 143 amino acids = 20 aa (N-linker) + 118 aa (VHH) + 5 aa (C-linker)

### Why Asymmetric? (Critical Design Insight)

The asymmetric linker is **intentional and evidence-based**, NOT a mistake:

1. **Polymer physics:** Long N-terminal linker acts as entropic spring, allows VHH to "scan" for receptors
2. **Structural stability:** C-terminal anchor prevents VHH from collapsing onto capsid surface
3. **Empirical validation:** Biogen patent data shows >10× improvement over symmetric linkers
4. **Best for high-affinity VHHs:** Design D2 is optimal when KD < 10 nM

**Common misconception:** "Both sides should have linkers for symmetry" OR "Asymmetric means direct fusion"
**Reality:** Biogen D2 specifies ASYMMETRY = long N-term (20 aa) + SHORT C-term (5 aa), NOT zero C-term. The 4:1 ratio (not 4:0) provides optimal performance.

### The Lesson

Document all insertion parameters AND the biophysical rationale:
1. Insertion site (VR region, amino acid position, bp position)
2. Linker design (type, sequence, length, **asymmetry reasoning**)
3. Insert sequence (target, size, codon optimization)
4. Total size impact
5. **Biophysical mechanism** (why this design works)

---

## Checkpoint 10: Parent-Child Sequence Verification (NEW)

**Purpose:** After building any construct, verify that ONLY intentional changes exist by comparing to the direct parent plasmid.

### Procedure

1. **Load parent and child sequences:**
   ```
   parent = parent_plasmid.seq
   child = child_plasmid.seq
   ```

2. **Compare backbone regions (outside modified area):**
   ```
   - Before modified region: should be IDENTICAL
   - After modified region: should be IDENTICAL
   ```

3. **Compare modified region:**
   ```
   - List ALL nucleotide differences
   - For each difference, classify as:
     - INTENTIONAL: Matches design specification (e.g., VP2/VP3 knockout)
     - UNINTENTIONAL: Not in design specification (ERROR!)
   ```

4. **Acceptance Criteria:**
   ```
   ✅ PASS: All differences are INTENTIONAL
   ❌ FAIL: Any UNINTENTIONAL differences found
   ```

5. **If FAIL:**
   ```
   - Investigate source of unintentional changes
   - Check if wrong parent sequence was used
   - Check for copy-paste errors in sequence assembly
   - Rebuild from correct parent sequence
   ```

### Example Output

```
=== CHECKPOINT 10: PARENT-CHILD VERIFICATION ===

Parent: AVD002 (7,104 bp)
Child: AVD006 (7,521 bp)
Expected difference: +417 bp (VHH insertion)

BACKBONE COMPARISON:
Before VP1 (0-2378): IDENTICAL ✅
After VP1 (5006-end): IDENTICAL ✅

MODIFIED REGION COMPARISON:
Found 2 point mutations + 1 insertion:

| Position | Parent | Child | Type | Status |
|----------|--------|-------|------|--------|
| 2792 | G | C | VP2 knockout (ACG→ACC) | ✅ INTENTIONAL |
| 2985 | A | C | VP3 knockout (ATG→CTG) | ✅ INTENTIONAL |
| 3747 | - | +417bp | VHH insertion | ✅ INTENTIONAL |

OVERALL: ✅ PASS — Only intentional changes detected
```

---

## DESIGN-008: Codon Optimization for Synthesis Compatibility

**Type:** Best Practice
**Discovered:** 2026-01-14
**Project:** AVD005/AVD006 synthetic fragment ordering

### The Problem

Initial codon optimization used simple lookup tables preferring "high-frequency" human codons:

```python
CODON_USAGE_HUMAN = {
    'G': 'GGC',  # 75% GC!
    'S': 'AGC',  # 67% GC!
    ...
}
```

This created synthesis problems:
- (GGGGS)×4 linker: **93.3% GC** (using GGC-GGC-GGC-GGC-AGC)
- VHH sequence: **70.3% GC**
- **718 direct repeats** >= 9bp from repetitive linker pattern

GenScript flagged these as synthesis issues (high GC + inverted repeats).

### The Solution

Use dnachisel to optimize with explicit constraints:

```python
from dnachisel import *

problem = DnaOptimizationProblem(
    sequence=original_sequence,
    constraints=[
        EnforceTranslation(),                 # Preserve protein sequence
        EnforceGCContent(mini=0.35, maxi=0.65, window=50),  # GC in range
        AvoidHairpins(stem_size=9, hairpin_window=200),     # No hairpins
        AvoidPattern(enzyme_pattern("BsaI")),  # Avoid cloning sites
    ],
    objectives=[
        CodonOptimize(species="h_sapiens"),   # Prefer common codons
    ]
)
problem.resolve_constraints()
problem.optimize()
```

For repetitive linkers like (GGGGS)×4, use varied codons manually:

```python
# OLD: Repetitive high-GC (93% GC)
"GGCGGCGGCGGCAGCGGCGGCGGCGGCAGCGGCGGCGGCGGCAGCGGCGGCGGCGGCAGC"

# NEW: Varied codons (67% GC)
"GGTGGAGGCGGATCTGGAGGCGGTGGTTCAGGCGGTGGAGGAAGTGGTGGCGGAGGTTCT"
```

### Results

| Metric | Original | Optimized | Change |
|--------|----------|-----------|--------|
| Linker GC | 93.3% | 66.7% | -26.6% |
| VHH GC | 70.3% | 59.9% | -10.4% |
| Insert GC | 73.2% | 60.9% | -12.3% |
| Direct repeats | 718 | 8 | -99% |

### The Lesson

**Simple codon optimization can create synthesis problems. Use dnachisel with explicit GC and hairpin constraints, especially for repetitive sequences like glycine-serine linkers.**

### Prevention Checklist

- [ ] Before ordering synthesis, check GC content in 50bp windows (target: 35-65%)
- [ ] Check for direct repeats >= 9bp
- [ ] Check for inverted repeats (potential hairpins) >= 9bp
- [ ] For repetitive sequences (GGGGS linkers), use varied codon synonyms
- [ ] Use dnachisel with EnforceGCContent and AvoidHairpins constraints
- [ ] Run sequence through synthesis vendor's complexity checker before ordering

---

## Summary of Lessons

| ID | Type | Summary | Prevention |
|----|------|---------|------------|
| BUG-005 | Critical Bug | Parent sequence mismatch in multi-construct builds | Always use VP1 from each construct's own parent |
| DESIGN-005 | Best Practice | Synthetic fragment boundary selection | Use flanking sites, not internal sites |
| DESIGN-006 | Design Pattern | VP2/VP3 knockout strategy | ACG→ACC (silent), ATG→CTG (non-silent) |
| DESIGN-007 | Design Pattern | VHH insertion at variable regions | Document site, linkers, insert clearly |
| DESIGN-008 | Best Practice | Codon optimization for synthesis | Use dnachisel with GC and hairpin constraints |
| Checkpoint 10 | Verification | Parent-child sequence comparison | Verify ONLY intentional changes exist |

---

## Version History

- **2026-01-14:** Added BUG-005, DESIGN-005, DESIGN-006, DESIGN-007, DESIGN-008, Checkpoint 10 from AVD005/AVD006 project
- **2026-01-14:** Added DESIGN-008 (dnachisel optimization for synthesis compatibility)
