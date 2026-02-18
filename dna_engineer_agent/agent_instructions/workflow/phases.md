# Workflow Phases (0–5)

> Part of: [Agent Instructions](../README.md)

---

## Phase 0: Project Setup and System Identification

Before any analysis or design, establish the project parameters.

### Step 0.1: Identify the Biological System

**Questions to answer:**
1. What biological GOAL does this sequence serve?
2. What SYSTEM does it belong to?
3. What is the END USE (verification, synthesis design, optimization)?

**How to determine this:**
- Ask the user if not clear
- BLAST the whole sequence against NCBI to find similar constructs
- Look for diagnostic features (ITRs → AAV, LTRs → Lentivirus, etc.)

**Load the appropriate knowledge base:**
```
knowledge_base/systems/{system_id}.json
```

Example systems:
- `recombinant_aav_production` — For AAV transfer plasmids, Rep-Cap helpers
- (Future: `lentivirus_production`, `crispr_delivery`, etc.)

### Step 0.2: Select Synthesis Workflow

**If this is a SYNTHESIS design project**, determine:

1. **Backbone Selection**
   - Load backbone catalog: `backbones/genscript/BACKBONE_CATALOG.json`
   - Choose backbone based on:
     - System requirements (AAV → pGS-AAV, FLASH → pUC57)
     - ITR type (ITR128 vs ITR145, scAAV vs ssAAV)
     - Selection marker (Amp vs Kan)
     - Downstream workflow (Golden Gate → use FLASH "Key IIS Free")

2. **Cloning Method**
   - Load workflows: `knowledge_base/synthesis/genscript_workflows.json`
   - Determine cloning strategy:
     - **HindIII/XbaI** (for pGS-AAV backbones)
     - **EcoRV blunt** (for FLASH backbones)
   - Document cloning sites that MUST be unique in final construct

3. **Exclusion Zones**
   - Load: `knowledge_base/exclusion_zones.json`
   - Identify regions where mutations are forbidden or restricted:
     - ITR boundaries (±10bp)
     - Kozak consensus
     - PolyA signals
     - Promoter cores
     - Splice sites

4. **Enzyme Constraints**
   - Load: `knowledge_base/enzyme_metadata.json`
   - If Golden Gate downstream: identify Type IIS sites to avoid (BsaI, BsmBI, BbsI)
   - Document fussy enzymes to avoid
   - **CRITICAL**: Note that ALL pGS-AAV backbones contain BbsI site just downstream of 3' ITR

**Output:** Project specification document:
```
=== PROJECT SPECIFICATION ===
System: recombinant_aav_production / transfer_plasmid
End Use: Synthesis design for GenScript

BACKBONE:
- Name: pGS-ssAAV-ITR128-Amp-empty
- Length: 5010 bp
- ITRs: 5' ITR128 (27-154), 3' ITR128 (2312-2439)
- Selection: Ampicillin

CLONING METHOD:
- Strategy: HindIII (AAGCTT) / XbaI (TCTAGA)
- 5' cloning site: HindIII at position 185
- 3' cloning site: XbaI at position 2276

CONSTRAINTS:
- NO internal HindIII sites (required for cloning)
- NO internal XbaI sites (required for cloning)
- Backbone contains BbsI at position 2299 (incompatible with BbsI-based Golden Gate)
- Preserve ITR sequences (27-154, 2312-2439)

EXCLUSION ZONES:
- ITR boundaries: ±10bp from positions 154 and 2312
- Promoter core (if applicable)
- Kozak consensus around ATG
- PolyA signal region

CHECKPOINTS TO EXECUTE:
[All standard checkpoints 1-7, plus new:]
[ ] Checkpoint 8: Silent Mutation Verification
[ ] Checkpoint 9: Cloning Site Uniqueness Verification
```

---

## Phase 1: Derive Requirements from the Knowledge Base

Once you know the system, the knowledge base tells you what molecular components are required.

**For each system, extract:**

1. **Required Proteins** — What proteins must be produced? Function? Reference sequences?
2. **Required Cis-Elements** — Promoters, regulatory elements, packaging/replication signals?
3. **Construct Type Classification** — What type? What should/shouldn't it contain?

**Output:** A dynamically generated manifest of things to find and verify.

**CRITICAL: Create an explicit checklist:**

```
=== VERIFICATION CHECKLIST ===
Derived from: recombinant_aav_production / transfer_plasmid

PROTEINS TO VERIFY:
[ ] Transgene CDS — User-specified protein

CIS-ELEMENTS TO VERIFY:
[ ] 5' ITR (verify length, serotype)
[ ] 3' ITR (verify length, serotype)
[ ] Promoter (verify identity, position)
[ ] Kozak sequence (verify consensus)
[ ] PolyA signal (verify sequence, position)
[ ] WPRE (if present)

MUST NOT CONTAIN:
[ ] Confirm NO Rep genes
[ ] Confirm NO Cap genes
[ ] Confirm NO packaging signals beyond ITRs

STRUCTURAL RULES:
[ ] Insert size within AAV packaging limit (scAAV: ~2.2kb, ssAAV: ~4.7kb total between ITRs)
[ ] Kozak consensus optimal or acceptable
[ ] No cryptic splice sites in transgene
[ ] PolyA signal 50-200bp after CDS stop codon
```

**This checklist drives all subsequent phases. Every box must be checked before generating the final report.**

---

## Phase 1.5: Checklist Execution Loop (MANDATORY)

This loop ensures every derived requirement is verified. **Do not skip items because the user didn't mention them.**

```
FOR EACH item in VERIFICATION CHECKLIST:
    1. Fetch reference sequence (from UniProt/NCBI/knowledge base)
    2. Find corresponding region in input sequence
    3. Align reference to input
    4. Apply acceptance criteria:
       - Identity ≥ 85% → ✅ VERIFIED
       - Identity 70-84% → ⚠️ PARTIAL
       - Identity < 70% → ❌ FAILED
       - Reference unavailable → ⬜ UNVERIFIED (document why)
    5. Record: coordinates, identity %, status
    6. Mark checkbox complete: [✅], [⚠️], [❌], or [⬜]

CONTINUE until ALL checkboxes are marked.
ONLY THEN proceed to checkpoint phases.
```

---

## Phase 2: BLAST for Gross Anatomy

Before looking for individual elements, understand the overall structure.

**Step 2.1: Whole-sequence BLAST**
- BLAST the entire input sequence against NCBI nr/nt
- What known construct is this most similar to? Overall identity?

**Step 2.2: ORF Identification**
- Find all open reading frames > 100 amino acids
- BLAST each ORF against NCBI protein database
- What known proteins does each ORF match?

**Output:** A map of "this region matches X protein at Y% identity"

---

## Phase 3: Verify Each Required Element by Homology

For each element in your derived manifest:

1. **Fetch Reference Sequence** — From knowledge_base/references/ or external databases
2. **Align Reference to Input** — Record: coordinates, identity %, coverage
3. **Apply Acceptance Criteria** — ≥85% VERIFIED, 70-84% PARTIAL, <70% FAILED
4. **Compare to File Annotations** — Flag discrepancies > 30bp

**Critical Rule:** The alignment result is the truth, not the file annotation.

---

## Phase 4: Validate Structural Rules

After verifying individual elements, check how they work together.

**For AAV transfer plasmids:**
1. Insert Size ≤ packaging limit
2. Element Order: 5'ITR → promoter → Kozak → CDS → polyA → 3'ITR
3. Appropriate spacing between elements
4. CDS in-frame with start codon

**For Rep-Cap helpers:**
1. Promoter Architecture: p5, p19, p40 in correct positions
2. Splicing Pattern: Rep68/40 splice sites intact
3. Frame Relationships: AAP in +1 frame relative to VP1
4. Nesting: VP1/VP2/VP3 share C-terminus

---

## Phase 5: Generate Report

**Report Structure:**
```
=== DNA SEQUENCE VERIFICATION REPORT ===
Generated: [timestamp]
Input: [filename or sequence ID]
System: [biological system identified]
Status: [PASS / PARTIAL / FAIL]

--- SUMMARY ---
[2-3 sentence summary]

--- VERIFIED ELEMENTS ---
--- STRUCTURAL VALIDATION ---
--- CHECKPOINTS ---
--- WARNINGS/FLAGS ---
--- RECOMMENDATIONS ---
```

**Status Determination:**
- ✅ **PASS**: All elements ≥85% identity, all checkpoints passed
- ⚠️ **PARTIAL**: Some elements 70-84%, or minor checkpoint warnings
- ❌ **FAIL**: Any element <70%, or critical checkpoint failure
