# Checkpoint 10: Parent-Child Sequence Verification

> Part of: [Agent Instructions](../README.md)

**Purpose:** After building any construct, verify that ONLY intentional changes exist by comparing to the direct parent plasmid. Prevents BUG-005 (parent sequence mismatch).

**CRITICAL RULE:** When building multiple related constructs, ALWAYS use the sequence from each construct's OWN parent plasmid. Never assume two "similar" plasmids have identical sequences.

**This checkpoint is CRITICAL after:**
- Building a new construct from a parent plasmid
- Inserting VHH/peptide fusions into capsid proteins
- Applying VP2/VP3 knockouts
- Any multi-construct build project

---

## Procedure

1. **Load Parent and Child Sequences**

2. **Identify the Modified Region**
   - Start/end of modified region (e.g., VP1 start/end + insertion)
   - Expected size difference (e.g., +417 bp for VHH insertion)

3. **Compare Backbone Regions (Outside Modified Area)**
   - Both "before" regions should be IDENTICAL
   - Both "after" regions should be IDENTICAL

4. **Compare Modified Region**
   - For each difference, classify as INTENTIONAL or UNINTENTIONAL

5. **Acceptance Criteria**
   - PASS: Backbone regions identical, all differences intentional, size matches
   - FAIL: Any backbone differences, any unintentional differences, size mismatch

6. **Output Example:**
   ```
   === CHECKPOINT 10: PARENT-CHILD VERIFICATION ===

   Parent: AVD002-Rep2Mut2Cap9-6R-wt.dna (7,104 bp)
   Child: AVD006-Rep2Mut2Cap9-VP1-VHH3-ALPL.gb (7,521 bp)
   Expected difference: +417 bp (VHH insertion)
   Actual difference: +417 bp

   BACKBONE COMPARISON:
   Before modified region: IDENTICAL
   After modified region: IDENTICAL

   MODIFIED REGION COMPARISON:
   | Position | Parent | Child | Classification | Status |
   |----------|--------|-------|----------------|--------|
   | 2792     | G      | C     | VP2 knockout   | INTENTIONAL |
   | 2985     | A      | C     | VP3 knockout   | INTENTIONAL |
   | 3747     | -      | +417bp| VHH insertion  | INTENTIONAL |

   OVERALL: PASS — Only intentional changes detected
   ```

7. **If FAIL:**
   - DO NOT proceed with synthesis
   - Investigate: wrong parent used? Sequences from different plasmids mixed?
   - Common cause: Using VP1 from plasmid A when building from plasmid B
   - Rebuild using the correct parent sequence
   - Re-run Checkpoint 10

**References:** BUG-005 in Lessons_learned.md
