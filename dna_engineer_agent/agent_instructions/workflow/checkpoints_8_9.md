# Checkpoints 8–9: Silent Mutations & Cloning Site Uniqueness

> Part of: [Agent Instructions](../README.md)

---

## Checkpoint 8: Silent Mutation Verification

**Purpose:** Verify that all nucleotide mutations introduced to remove restriction sites or optimize codon usage are truly SILENT (synonymous) and do not change the amino acid sequence.

**This checkpoint is CRITICAL after:**
- Removing internal HindIII/XbaI sites from CDS
- Removing Type IIS sites (BsaI, BsmBI, BbsI) for Golden Gate compatibility
- Any codon optimization operations

### Procedure

1. **Identify All Mutations**
   - Compare original CDS to final optimized CDS
   - List all nucleotide positions that differ

2. **For Each Mutation:**
   ```
   a. Determine codon position:
      - Calculate frame offset: (mutation_position - CDS_start) % 3
      - Frame 0 = first position of codon
      - Frame 1 = second position (middle)
      - Frame 2 = third position (wobble)

   b. Extract affected codon:
      - Original codon (3 bases including mutation site)
      - Mutated codon (3 bases with new nucleotide)

   c. Translate both codons:
      - original_AA = translate(original_codon)
      - mutated_AA = translate(mutated_codon)

   d. Verify silence:
      - PASS: original_AA == mutated_AA
      - FAIL: original_AA != mutated_AA (non-synonymous!)
      - CRITICAL: mutation creates STOP codon
   ```

3. **Special Cases:**
   - **Start codon (ATG):** Should NEVER be mutated (Met has only one codon)
   - **Stop codons (TAA/TAG/TGA):** Can be mutated between stop codons (all silent), but NEVER to sense codon
   - **Wobble position (frame 2):** Usually safe, but verify translation
   - **Second position (frame 1):** ALWAYS changes amino acid — mutation here is NEVER silent

4. **Acceptance Criteria:**
   - PASS if ALL mutations are silent (100%)
   - REVIEW if >95% silent (check the non-silent ones)
   - FAIL if <95% silent or ANY stop codon created

5. **Output Example:**
   ```
   === CHECKPOINT 8: SILENT MUTATION VERIFICATION ===

   Total mutations: 15
   Silent mutations: 15 (100%)
   Non-silent mutations: 0
   Stop codons created: 0

   MUTATIONS ANALYZED:
   Position | Original Codon | New Codon | Original AA | New AA | Status
   ---------|---------------|-----------|-------------|--------|--------
   342      | GAA           | GAG       | Glu         | Glu    | SILENT
   567      | CTT           | CTC       | Leu         | Leu    | SILENT

   OVERALL: PASS — All mutations are silent
   ```

6. **If Failures Found:**
   - DO NOT PROCEED with synthesis
   - Redesign the mutations to be truly silent
   - Common fixes:
     * For position 0 or 1: try mutating position 2 (wobble) instead
     * For Leu, Ser, Arg (6 codons): many silent alternatives
     * For Met, Trp (1 codon): cannot remove restriction site at these positions
   - Re-run Checkpoint 8 after fixes

**Test Coverage:** `scripts/tools/tests/test_silent_classification.py`

---

## Checkpoint 9: Cloning Site Uniqueness Verification

**Purpose:** Verify that restriction sites used for cloning (HindIII, XbaI, EcoRV) appear EXACTLY ONCE in the final assembled construct (backbone + insert).

**This checkpoint prevents:**
- **DESIGN-001:** Internal cloning sites that prevent successful restriction digest/ligation
- **BUG-003:** Incorrect uniqueness counting (must check both strands for palindromic sites)

### Procedure

1. **Load Project Specification** — Backbone sequence, cloning sites, insert sequence from Phase 0

2. **Assemble Full Construct**
   ```
   full_construct = backbone_5prime + insert + backbone_3prime
   ```

3. **Count Cloning Sites in Full Construct**
   ```
   For each cloning enzyme (e.g., HindIII, XbaI):

   a. Get recognition sequence:
      HindIII: AAGCTT | XbaI: TCTAGA | EcoRV: GATATC

   b. Count on BOTH strands:
      For PALINDROMIC sites: count = occurrences in forward strand
      For NON-PALINDROMIC sites: count = forward + reverse complement

   c. Record all positions
   ```

4. **Apply Uniqueness Criteria**
   - UNIQUE: Exactly 1 occurrence in full construct
   - ABSENT: 0 occurrences (cloning won't work)
   - NOT UNIQUE: 2+ occurrences (digest creates multiple fragments)

5. **Context Analysis (if NOT UNIQUE)**
   - Identify source of extra sites (backbone, insert CDS, regulatory elements)
   - NEVER remove the backbone cloning site
   - Remove internal sites via silent mutations
   - Check exclusion zones before mutating
   - Re-run Checkpoints 8 and 9 after fixes

6. **Output Example:**
   ```
   === CHECKPOINT 9: CLONING SITE UNIQUENESS ===

   Backbone: pGS-ssAAV-ITR128-Amp-empty (5010 bp)
   Insert: EF1a-transgene-bGHpA (1523 bp)
   Cloning method: HindIII / XbaI

   Enzyme   | Recognition | Count | Positions | Status
   ---------|-------------|-------|-----------|--------
   HindIII  | AAGCTT      | 1     | [185]     | UNIQUE
   XbaI     | TCTAGA      | 1     | [2276]    | UNIQUE

   OVERALL: PASS — All cloning sites are unique
   ```

**Note:** Checkpoints 8 and 9 often run in a loop:
1. Run CP9 to find non-unique sites
2. Design silent mutations to remove internal sites
3. Run CP8 to verify mutations are silent
4. Re-run CP9 to confirm uniqueness
5. Repeat until both pass

**Test Coverage:** `scripts/tools/tests/test_uniqueness_counting.py`
