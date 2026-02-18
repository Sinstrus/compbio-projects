# DNA Engineer Agent — Toolset & Construct Bank

This directory is a **shared toolset and construct bank**, NOT a project directory.
Project-specific work lives in `~/projects/<project-name>/`.

## Construct Bank

All DNA construct files (.gb, .dna) live in `constructs/` — a flat directory indexed by AVD number.
Superseded versions are in `constructs/archived/`.

**Registry:** `CONSTRUCT_REGISTRY.md` — the authoritative manifest with IDs, filenames, and descriptions.

### Rules for New Constructs
1. Deposit ALL construct files to `constructs/`
2. Use the next available AVD number (check `CONSTRUCT_REGISTRY.md`)
3. Name as: `AVD{NNN}-{backbone}-{modification}-{site}.gb`
4. Update `CONSTRUCT_REGISTRY.md` with ID, filename, and 2-3 sentence description
5. Run `construct_verifier.py` before finalizing (see Verification below)
6. NEVER store construct files in project directories

## Tools

### construct_verifier.py (`scripts/tools/`)
Mandatory verification for all constructs. Usage:
```python
from construct_verifier import ConstructVerifier, verify_and_gate
verifier = ConstructVerifier(sequence=seq, name="AVD006")
verifier.add_restriction_check_by_name(required_unique=["HindIII", "XbaI"])
verifier.add_protein_check(region_start=..., region_end=..., expected_protein=...)
verifier.add_enzyme_validation(["HindIII", "XbaI"])
verifier.add_synthesis_screen()
report = verify_and_gate(verifier)  # Raises VerificationError on CRITICAL failure
```

### design_reasoner.py (`scripts/tools/`)
Structured decision framework for design choices (enzyme selection, linker choice, codon variants).
```python
from design_reasoner import DesignReasoner, Metric, Option
reasoner = DesignReasoner("Select enzyme")
reasoner.add_metric(Metric("reliability", weight=1.0, higher_is_better=True, threshold_min=0.5))
reasoner.add_option(Option("HindIII", scores={"reliability": 1.0}))
result = reasoner.decide()
```
Pre-built reasoners: `enzyme_selection_reasoner()`, `linker_selection_reasoner()`, `codon_variant_reasoner()`

### splice_donor_hbonds.py (`scripts/tools/`)
H-bond scoring for splice donor sequences. Tier assignment based on U1 snRNA base-pairing strength.

### build_reporter_plasmids.py (`scripts/tools/`)
Builds dual NanoLuc/Firefly luciferase splicing reporter plasmids (AVD027-058 series).

## Reference Data

| Directory | Contents |
|-----------|----------|
| `knowledge_base/` | Enzyme metadata (enzyme_metadata.json), codon tables, splice site rules |
| `parts_library/` | Promoters, polyA signals, reporter genes, linker sequences |
| `backbones/` | Base plasmid maps and sequences |
| `cis_elements/` | ITRs, splice signals, packaging signals |

## Critical Bugs to Avoid

| Bug | Summary | Prevention |
|-----|---------|------------|
| BUG-001 | Frame offset must be relative to CDS start | `(pos - cds_start) % 3` |
| BUG-003 | Check BOTH DNA strands for restriction sites | Use `reverse_complement()` |
| BUG-004 | Never hard-code DNA positions for AA insertions | Calculate from AA coords |
| BUG-005 | Use VP1 from EACH construct's own parent | Don't copy from different plasmid |
| BUG-006 | Dipeptide = 2 amino acids (MA = M + A) | Count both when calculating position |
| BUG-007 | Don't generate repetitive GGGGS codons | Use varied-codon sequences |

See `LESSONS_LEARNED.md` for full bug and design decision documentation.

## Projects Using This Toolset
- `~/projects/AVD_VHH_Display_ALPL/` — Anti-ALPL/TfR1 VHH capsid display
- `~/projects/splice_donor_pipeline/` — Splice donor strength characterization
- `~/projects/aav_transfer_plasmid/` — Transfer plasmid engineering
