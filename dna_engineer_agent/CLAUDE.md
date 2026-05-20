# DNA Engineer Agent — Toolset & Construct Bank

This is the **unified DNA engineering workspace** — toolset, construct bank, build scripts, and design docs.
Build scripts from former satellite projects (AVD_VHH_Display_ALPL, splice_donor_pipeline) are now in `scripts/`.

## Construct Bank

All DNA construct files (.gb, .dna) live in `constructs/` — a flat directory indexed by AVD number.
Superseded versions are in `constructs/archived/`.

**Registry:** `CONSTRUCT_REGISTRY.md` — the authoritative manifest with IDs, filenames, and descriptions.

### Rules for New Constructs
1. Deposit ALL construct files to `constructs/`
2. Use the next available AVD number (check `CONSTRUCT_REGISTRY.md`)
3. **AVD numbers MUST be unique to a single DNA sequence. Never assign an existing AVD number to a different sequence. Skipping numbers is fine; duplication is not. Multiple file formats or annotation variants of the same sequence share one AVD number.**
4. Name as: `AVD{NNN}-{backbone}-{modification}-{site}.gb`
5. Update `CONSTRUCT_REGISTRY.md` with ID, filename, and 2-3 sentence description
6. Run `construct_verifier.py` before finalizing (see Verification below)
7. NEVER store construct files in project directories

## Concept Tracker

Pre-pipeline concepts are tracked in `CONCEPT_TRACKER.md` — a status board for ideas that haven't entered the formal design → order → plasmid pipeline yet.

**When to consult:** At the start of any design or order-related task. Check whether a relevant concept is already queued before starting from scratch.

**When to update:**
- Concept moves to ORDERED → add the ORD barcode to the entry
- Concept is DELIVERED → add the AVD + BAT barcode to the entry
- Concept is SCRAPPED or BACK-BURNERED → update status in place

**DQ numbering:** Simple sequential integers (DQ001, DQ002…), assigned by the human or Claude when adding an entry. Check the last DQ number in the file — no automated tool required.

**Entry length discipline:** Each entry is ≤15 lines — concept (2–3 sentences), status, blockers, cross-refs. Deep design reasoning, codon tables, MaxEntScan scores, and mechanistic analysis belong in `DESIGN_RATIONALE.md`, not here. Add a link to the relevant DESIGN_RATIONALE.md once it exists.

**Pre-order design reasoning:** If a concept needs substantive design work before an ORD number exists, create `orders/staging/DQ###-slug.md` and link to it from the CONCEPT_TRACKER entry. When the order is placed, move that file to `orders/ORD-YYYYMMDD-NNN-slug/DESIGN_RATIONALE.md`.

**Never delete SCRAPPED or DELIVERED entries.** Update status in place to preserve the audit trail.

---

## Two distinct construction workflows

DNA constructs are built by **one of two paths**, tracked by separate registries. Pick the right path before starting work.

| Path | Tracker | When to use | Typical turnaround |
|---|---|---|---|
| **Vendor synthesis** (Genscript FLASH / Premium) | `orders/ORDER_REGISTRY.md` (`ORD-YYYYMMDD-NNN`) | A novel synthetic fragment is needed; vendor builds the plasmid | FLASH: ~3 wk; Premium: 4–6 wk; ITR-containing: 6–8 wk |
| **In-house cloning** (RE digest / Gibson / Golden Gate / SDM / PCR-overhang) | `cloning/CLONING_REGISTRY.md` (`CLN-YYYYMMDD-NNN`) | All DNA pieces already exist in-house; speed matters; ITR plasmids that would otherwise be 6–8 wk vendor-blocked | 5–10 days bench |

If a build needs **both** a synthesized fragment **and** an in-house assembly step, run the ORD first and chain a CLN to assemble the final. The CLN plan references the upstream ORD as a "Donor" input.

See `cloning/CLAUDE.md` for in-house cloning rules. Below is the vendor synthesis workflow.

---

## Orders

Synthesis orders are tracked in `orders/` with date-based barcoding.
Registry: `orders/ORDER_REGISTRY.md`

### Genscript Synthesis Workflows
New plasmids are almost always built by modifying an existing Genscript-made plasmid — cut out a region with restriction enzymes and ligate in a newly synthesized fragment. Two workflows:

- **FLASH Gene Synthesis** (fast, cheap): Requires (1) starting plasmid stored in VectorArk or Vector Onboarding, (2) synthetic fragment rated "simple" by Genscript complexity screening, (3) no ITR-containing plasmids (ITRs are ineligible for VectorArk storage).
- **Premium Gene Synthesis** (slower, handles anything): Can start from any existing Genscript-made plasmid, handles "complex" sequences and ITR-containing plasmids.
- **DESIGN-017 (applies to both FLASH and Premium):** Only the backbone is digested — internal enzyme sites in the synthetic fragment are irrelevant. The fragment is synthesized de novo and never exposed to restriction enzymes. Always use the closest unique enzyme pair flanking the modification site.

### VectorArk Inventory
Registry: `orders/VECTORARK_REGISTRY.md` — tracks constructs stored in Genscript VectorArk (28 total: 25 TRX + 3 AVD). VectorArk is required for FLASH orders; ITR plasmids are ineligible.

### Rules for New Orders
1. Use barcode format `ORD-YYYYMMDD-NNN-descriptive-slug` — NNN = sequence number for that date (starting at 001), then up to 5 lowercase hyphenated keywords describing the panel (e.g. `ORD-20260417-001-wave2-VR4-IR-panel`)
2. Create a subdirectory `orders/ORD-YYYYMMDD-NNN-descriptive-slug/` for order files (CSVs, notes)
3. Update `orders/ORDER_REGISTRY.md` with Order ID, date, vendor, method, backbone, constructs, status
4. Update `CONSTRUCT_REGISTRY.md` Order History section with cross-reference
5. Every order MUST include a `Narrative` column (2-3 sentences: what hypothesis, what panels, what the order tests)
6. Every amendment MUST be logged in the Changelog section of ORDER_REGISTRY.md with date and reason
7. NEVER scatter order files in project directories or constructs/
8. Store vendor quote zips in the order subdirectory (`orders/ORD-YYYYMMDD-NNN-descriptive-slug/`). The full `Order_*_Files.zip` (with Quote PDF + Order XLSX) is REQUIRED, not just the FASTA submission file (DESIGN-028).
9. Before confirming any order: (a) generate `quote_metadata.json` via `python scripts/tools/parse_genscript_quote.py orders/<ORD-folder>` (or `inventory_tracker.py ingest-quote <ORD>`), then (b) run `python scripts/tools/audit_order.py <vendor_quote.zip>` — Phase 0 BLOCKS if quote zip / metadata missing (DESIGN-028); Phase 1 verifies vendor sequences against submitted CSV and registries.
10. When placing a FLASH order, verify the target vector is in `orders/VECTORARK_REGISTRY.md`
11. For Premium orders modifying an existing plasmid, record the source plasmid in the order CSV
12. After submitting a new VectorArk storage order, update `orders/VECTORARK_REGISTRY.md`
13. **Every order MUST have a `DESIGN_RATIONALE.md` in its folder.** For simple single-construct amendments, a one-paragraph stub is sufficient. For non-trivial panels, cover: scientific hypothesis, panel structure, key design decisions, wave-1 comparators, and predicted outcomes. This is the authoritative design context — all deep reasoning belongs here, not in `CONCEPT_TRACKER.md`.
14. **Every Narrative MUST end with prep-scale summary** in the form `<N>×<scale> at <conc>` (e.g. "6×1mg at 1mg/mL TE", "44×0.2mg at 1mg/mL TE"). Heterogeneous scales: list breakdown ("4×1mg + 2×0.1mg at 1mg/mL TE"). Mirrors `quote_metadata.json` for at-a-glance review.

### Order CSV column format (mandatory)
The CSV **must** be generated by script (never manual transcription) and follow this exact column order:
1. **`Name`** — `.gb` construct filename minus the `.gb` extension (e.g. `AVD210-Rep2Mut2Cap9-VR4-VHH27-D2v4-mosaic`). Look up the actual filename in `constructs/`.
2. **`Sequence`** — synthetic fragment sequence read directly from the fragment `.txt` file (plain uppercase, no FASTA header, no line breaks).
3. Remaining metadata in order: `Backbone`, `VectorArk_ID`, `Enzyme_5`, `Enzyme_3`, `Fragment_Size_bp`, `Final_Construct_Size_bp`, `Partner_Plasmid`, `Production_Strategy`, `Narrative`

Use `csv.DictWriter` with `quoting=csv.QUOTE_ALL` and `encoding='utf-8-sig'` (Excel-safe).
Also produce a companion `ORD-YYYYMMDD-NNN_sequences_for_submission.txt` FASTA from the same script.

## Tools

### construct_verifier.py (`scripts/tools/`)
Mandatory verification for all constructs. Usage:
```python
from construct_verifier import ConstructVerifier, verify_and_gate
verifier = ConstructVerifier(sequence=seq, name="AVD006")
verifier.add_avd_uniqueness_check()
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

### asset_registry.py (`scripts/tools/`)
Unified registry utility for FIG###, PPP###, and AVD### barcoded asset banks. Parses markdown-table registries, assigns next available numbers, and audits bank directories for consistency.
```python
from asset_registry import parse_registry, get_next_number, validate_bank

# Get next available FIG number
next_fig = get_next_number("../figure_bank/FIGURE_REGISTRY.md", "FIG")  # → "FIG007"

# Get next available PPP number
next_ppp = get_next_number("../voyager_slide_decks/PRESENTATION_REGISTRY.md", "PPP")

# Audit figure bank for orphans, missing files, duplicates
issues = validate_bank("../figure_bank/FIGURE_REGISTRY.md", "../figure_bank/", "FIG")

# Parse registry into structured dict
entries = parse_registry("CONSTRUCT_REGISTRY.md", "AVD")  # {1: {barcode: "AVD001", ...}, ...}
```

### select_flash_backbone.py (`scripts/tools/`)
Selects the optimal VectorArk backbone (AVD001 or AVD002) and enzyme pair for
FLASH synthesis of AAV capsid display constructs. Ensures VP start codon positions
outside the enzyme pair match the backbone (DESIGN-018).
```python
from select_flash_backbone import select_flash_backbone, detect_vp_nucleotides, explain_selection

# Select backbone + enzyme pair for a construct
result = select_flash_backbone(construct_seq, insertion_site="VP2n")
print(result["backbone"])   # "AVD001" or "AVD002"
print(result["enzyme_5"])   # e.g. "SmaI"
print(result["enzyme_3"])   # e.g. "BbvCI"
print(explain_selection(result))

# Detect VP start codon nucleotides
vp = detect_vp_nucleotides(construct_seq)  # {'VP1': 'A', 'VP2': 'G', 'VP3': 'C'}
```
Key functions: `select_flash_backbone()`, `detect_vp_nucleotides()`, `extract_fragment()`, `explain_selection()`

### optimize_intron_body_codons.py (`scripts/tools/`)
3-stage intron body codon optimizer (GGGGS5 + VHH3 + GGGGS1, 447bp).
v4 GGGGS design uses `KMER_SIZE=9, MAX_KMER_REPEATS=1, MAXENT_THRESHOLD=0.0`.
Key functions: `design_ggggs_maxent_guided()`, `_collect_kmers()`, `_has_duplicate_kmers()`.
**Note:** Old `*_8mers` functions were renamed to `*_kmers` and parameterized on `KMER_SIZE` in v4.

### parse_genscript_quote.py (`scripts/tools/`)
Ingest a Genscript `Order_*_Files.zip` into structured `quote_metadata.json` for an
ORD folder. Captures per-construct ordered prep scale (mg), concentration,
supercoil/endotoxin spec, unit price, and template lot ID — the inputs needed to
reason about plasmid inventory before delivery and to backfill the ledger after
delivery.

Primary parsing path: `openpyxl` on the `Order-*.xlsx` "Plasmid unified" sheet
(cleanest, machine-structured). Fallback: regex on `GenScript_Order_*.txt` "Item >N"
blocks. Multi-zip ORD folders (e.g. wave2 panel split across two quotes) are merged
into one `quote_metadata.json`.

```bash
# Auto-locate the zip in an ORD folder
python scripts/tools/parse_genscript_quote.py orders/ORD-20260429-001-vhh38-cdr3mut-neg-ctrl

# Override path / re-generate
python scripts/tools/parse_genscript_quote.py <ORD-folder> --zip <path> --force

# Self-test against ORD-20260429-001 known-good zip
python scripts/tools/parse_genscript_quote.py --selftest

# Convenience wrapper from inventory_tracker
python scripts/tools/inventory_tracker.py ingest-quote ORD-20260429-001
```

Importable: `parse_quote_zip()`, `find_quote_zips()`, `extract_avd_id()`,
`parse_quantity_mg()`, `write_quote_metadata()`. Required by `audit_order.py`
Phase 0 (DESIGN-028 BLOCK).

### genscript_coa_pull.py (`scripts/tools/`)
Batch download and parse Genscript CoA PDFs to extract concentration/QC data for inventory backfill.
Four subcommands form a sequential pipeline:
```bash
# Step 1 — get JS snippet to run in browser console (while logged into Genscript orders page)
python scripts/tools/genscript_coa_pull.py print-js

# Step 2 — paste JS output into browser, save JSON to file, copy document.cookie string

# Step 3 — download all CoA PDFs
python scripts/tools/genscript_coa_pull.py download orders/coa_urls.json "COOKIE_STRING" \
    --out orders/coa_pdfs/

# Step 4 — parse PDFs for amount/concentration/A260-280
python scripts/tools/genscript_coa_pull.py parse orders/coa_pdfs/

# Step 5 — generate inventory_tracker backfill commands
python scripts/tools/genscript_coa_pull.py backfill orders/coa_pdfs/ \
    --ord ORD-YYYYMMDD-NNN --date YYYY-MM-DD
```
Key limitation: the JS snippet uses DOM scraping (Strategy 1 = anchor tags, Strategy 2 = onclick/data attributes). If CoA links are loaded dynamically via XHR, the snippet will warn and you'll need to capture the API call from the Network tab instead.

### audit_order.py (`scripts/tools/`)
Pre-submission audit for synthesis orders against Genscript quote zips. Three-phase:
- **Phase 0 (BLOCKING, DESIGN-028):** verifies `Order_*_Files.zip` is present in the ORD folder.
  If `quote_metadata.json` is missing, **auto-generates it** from the zip(s) before continuing.
  Verifies every submitted CSV AVD has an `ordered_mg`. Reports prep-scale summary and total cost.
  Exits 1 only if the zip itself is missing or coverage is incomplete (override with
  `--allow-missing-quote` for legitimate exceptions).
- **Phase 1 (programmatic):** vendor insert vs. submitted sequence, registry existence,
  synthesis screen, flanking sites, .gb spot-check.
- **Phase 2 (optional LLM):** Claude Code review of narrative accuracy + scientific coherence.

```bash
python scripts/tools/audit_order.py /path/to/Order_*.zip          # full audit
python scripts/tools/audit_order.py /path/to/Order_*.zip --no-llm  # programmatic only
python scripts/tools/audit_order.py /path/to/Order_*.zip --allow-missing-quote  # legitimate exception only
```

### plasmidsaurus_checker.py (`scripts/tools/`)
Aligns Plasmidsaurus long-read sequencing results (zip format) against reference
constructs (.gb/.dna) and produces structured PASS/WARN/FAIL/CRITICAL_FAIL reports
per sample. Handles circular sequence alignment (local-anchor + global strategy),
per-base quality annotation (coverage, homopolymer flags, Dam/Dcm methylation sites),
and ITR-aware variant interpretation. A structural deletion ≥20 bp in an ITR feature
is always CRITICAL_FAIL — never an artifact.

Known Nanopore pitfalls encoded: (1) homopolymers ≥5 nt → LIKELY_ARTIFACT for INS/DEL;
(2) Dam (GATC) and Dcm (CCAGG/CCTGG) sites → INTERPRET_WITH_CAUTION within ±2 bp;
(3) predicted_methylation_site from TSV → INTERPRET_WITH_CAUTION; (4) coverage <10
reads → INTERPRET_WITH_CAUTION; (5) ITR large deletions → always CRITICAL_FAIL.

```bash
# Check all samples in a zip against one reference
python scripts/tools/plasmidsaurus_checker.py check results.zip --ref constructs/AVD548.gb

# Batch: auto-match each sample by name prefix (e.g., AVD548-1 → constructs/AVD548*.gb)
python scripts/tools/plasmidsaurus_checker.py batch results.zip

# Override auto-match for specific samples
python scripts/tools/plasmidsaurus_checker.py batch results.zip \
    --ref TP688-146202=constructs/TP688-CBA-hGBint-EGFP-T2A-AkaLuc.dna

# JSON output
python scripts/tools/plasmidsaurus_checker.py check results.zip --ref constructs/AVD548.gb --json
```
Per-sample per-base TSV (coverage, homopolymer, methylation) is indexed by FASTA position.
Requires: `biopython`, `pandas`. Handles .gb/.gbk (BioPython genbank) and .dna (BioPython snapgene).

### gel_render.py (`scripts/tools/`)
SDS-PAGE / Western blot SVG figure generator. Renders predicted gel band patterns
from lane architecture data. Supports Coomassie, silver, and Western stain types
with color inversion. Includes MW calculation from amino acid sequences and
construct file CDS extraction via BioPython (handles both .gb and .dna SnapGene files).
Key functions: `render_gel_svg()`, `quick_gel()`, `extract_proteins_from_file()`,
`calc_mw_from_sequence()`, `check_comigration()`, `format_lane_table()`, `suggest_gel_pct()`.
MW ladder presets: Chameleon 700, Precision Plus, PageRuler.

### build_reporter_plasmids.py (`scripts/tools/`)
Builds dual NanoLuc/Firefly luciferase splicing reporter plasmids (AVD027-058 series).

### audit_slides.py (`scripts/tools/`)
**Mandatory post-build visual audit** for PPTX slides. Uses Claude Vision (Haiku by
default) to detect text overflow, collisions, figure sizing issues, and spacing
problems on rendered slides.

**FIG SVGs are NOT audited here.** Figure-level layout is enforced at render time
by the SCF### scaffold system in `figure_bank/scaffolds/` — slot bounds raise
exceptions before a figure can be saved. See `figure_bank/CLAUDE.md`.

**MANDATORY for slide decks:** Never run `node build.js` directly — use
`bash ../../build_and_audit.sh` which builds then audits automatically. Fix all
errors before presenting. Max 3 fix-rebuild cycles before escalating to the user.

**Caching:** Slide audits cache results by content hash. Cache files
(`*.audit-cache.json`) live in the deck directory. Unchanged passing slides are
skipped automatically. Use `--no-cache` to force re-audit.

```bash
# Normal workflow (from inside a deck directory):
bash ../../build_and_audit.sh

# Audit a PPTX directly (LibreOffice → PDF → per-slide PNGs → Claude Vision)
python scripts/tools/audit_slides.py path/to/PPPxxx-file.pptx

# Audit specific slides only
python scripts/tools/audit_slides.py path/to/PPPxxx-file.pptx --slides 3,5,7

# Use Sonnet for edge cases (default is Haiku)
python scripts/tools/audit_slides.py path/to/PPPxxx-file.pptx --model claude-sonnet-4-6

# Skip fix suggestions (faster)
python scripts/tools/audit_slides.py path/to/PPPxxx-file.pptx --no-fix-suggestions
```
Requires: `anthropic`, `pymupdf`, `python-pptx`, `Pillow`, `CairoSVG` (all in `.venv`).
API key auto-loaded from `~/projects/second-brain/.env`.
Key functions: `audit_pptx()`, `SlideRenderer`, `VisionAuditor`.

### design_primers.py (`scripts/tools/`)
RT-PCR primer design and QC tool. Uses **primer3-py v2.3.0** (native Python bindings)
for thermodynamic design and **MFEprimer v3.3.1** (`~/bin/mfeprimer`) for specificity
checking against a template sequence.

Two subcommands:
- `design`: design new primer pairs for an amplicon flanking a variable region
- `qc`: run thermodynamic QC on an existing primer pair

```bash
# Design primers flanking a variable region (e.g., VR4 AgeI=3596, BsrGI=3794)
python scripts/tools/design_primers.py design <construct.gb> \
    --left-boundary 3596 --right-boundary 3794 \
    --amplicon-min 200 --amplicon-max 400 \
    --tm-min 57 --tm-max 65 --num-return 5

# QC an existing primer pair (thermodynamics only)
python scripts/tools/design_primers.py qc ATATTTCCCGTCGCAAATGC TCCATTGAGAGCCCAAGAAG

# QC with MFEprimer specificity check
python scripts/tools/design_primers.py qc ATATTTCCCGTCGCAAATGC TCCATTGAGAGCCCAAGAAG \
    --construct constructs/AVD002-Rep2Mut2Cap9-6R-wt-5fold.gb
```

Key importable functions: `design_primers(sequence, left_boundary, right_boundary, ...)`,
`score_existing_pair(fwd, rev)`, `check_specificity(pairs, template_fasta, ...)`.

Thresholds: hairpin/homodimer/heterodimer Tm < 47°C = PASS; < 40°C = WARN (PASS).
MFEprimer: PASS = exactly 1 predicted amplicon (no off-targets in template).
Dependencies: `pip install primer3-py biopython` (already in `.venv`); `~/bin/mfeprimer` binary.

### inventory_tracker.py (`scripts/tools/`)
Plasmid inventory tracking with batch-level ledger. Tracks physical DNA tubes from
Genscript deliveries through experimental consumption. Batch barcode: `BAT-YYYYMMDD-NNN`.
Ledger: `orders/INVENTORY_LEDGER.md` (9 columns: Batch, Construct, Order, Delivered,
**Ordered (mg)**, Conc, Starting (ug), Remaining (ug), Status — `Ordered (mg)` sourced
from `quote_metadata.json` per ORD folder; `--` = unknown/legacy).
```bash
# Check stock levels
python scripts/tools/inventory_tracker.py status [--low] [AVD###]

# Ingest a Genscript quote zip (writes ORD-*/quote_metadata.json)
python scripts/tools/inventory_tracker.py ingest-quote <ORD> [--zip <path>] [--force]

# Auto-route prep-only restock zips dropped in orders/staging/ — scaffolds new
# ORD folder + registry row + DESIGN_RATIONALE.md stub when no existing ORD matches.
# (For sequence-synthesis orders, the ORD folder is created up-front by the build
# script and the quote zip goes directly there — bypass staging.)
python scripts/tools/inventory_tracker.py ingest-staging

# Register batches from a delivery (defaults to ordered_mg from quote_metadata.json)
python scripts/tools/inventory_tracker.py backfill <ORD> <delivery-date> [--ug <ug>] [--conc <mg/mL>]

# Re-populate Ordered (mg) on every batch from all quote_metadata.json files (idempotent)
python scripts/tools/inventory_tracker.py refresh-ordered-mg

# Parse experiment consumption (dry run)
python scripts/tools/inventory_tracker.py parse-experiment <EXP>

# Debit inventory from an experiment
python scripts/tools/inventory_tracker.py debit-experiment <EXP> [--date YYYY-MM-DD]

# Validate ledger integrity
python scripts/tools/inventory_tracker.py validate
```
Key importable functions: `parse_ledger()`, `write_ledger()`, `receive_batch()`,
`consume_from_batch()`, `consume_construct()` (FIFO), `get_stock()`, `get_low_stock()`,
`parse_experiment_consumption()`, `load_quote_metadata()`, `validate_ledger()`.

## Reference Data

| Directory | Contents |
|-----------|----------|
| `knowledge_base/` | Enzyme metadata (enzyme_metadata.json), codon tables, splice site rules, Kozak efficiency table (kozak_wang2014.csv) |
| `parts_library/` | Promoters, polyA signals, reporter genes, linker sequences |
| `backbones/` | Base plasmid maps and sequences |
| `cis_elements/` | ITRs, splice signals, packaging signals |
| `orders/` | Synthesis order tracking (ORD-YYYYMMDD-NNN barcoding) |

### Kozak Efficiency Reference (`knowledge_base/kozak_wang2014.csv`)
Wang et al. 2014 exhaustive Kozak dataset: 65,536 rows covering all 4^8 possible 8-mer contexts
around ATG. Columns: `sequence RNA` (11-mer RNA, `nnnnnnAUGnn` — 6 upstream + ATG + 2 downstream),
`sequence DNA` (same in DNA form), `efficiency` (relative translation strength, higher = stronger),
`lower.bound`, `upper.bound` (95% CI). GT/AG columns indicate position (1-based int) of a GT or AG
dinucleotide in the 11-mer context, or `#VALUE!` if absent.

Note: `sequence RNA` column has a UTF-8 BOM — use `encoding='utf-8-sig'` or pandas `utf-8-sig`.
Use pandas with `na_values="#VALUE!"` and `encoding='utf-8-sig'` to load cleanly:

```python
import pandas as pd
kozak_db = pd.read_csv(
    "knowledge_base/kozak_wang2014.csv",
    na_values="#VALUE!",
    encoding="utf-8-sig"
)
# Query by DNA sequence context (11-mer, positions -6 to +5 around ATG):
# Kozak consensus has GCCACC upstream + ATG + G at +4; search by upstream context:
hits = kozak_db[kozak_db["sequence DNA"].str.startswith("GCCACCATG")]
# Or exact 11-mer lookup:
row = kozak_db[kozak_db["sequence DNA"] == "GCCACCATGCA"]
```

## Critical Bugs to Avoid

| Bug | Summary | Prevention |
|-----|---------|------------|
| BUG-001 | Frame offset must be relative to CDS start | `(pos - cds_start) % 3` |
| BUG-003 | Check BOTH DNA strands for restriction sites | Use `reverse_complement()` |
| BUG-004 | Never hard-code DNA positions for AA insertions | Calculate from AA coords |
| BUG-005 | Use VP1 from EACH construct's own parent | Don't copy from different plasmid |
| BUG-006 | Dipeptide = 2 amino acids (MA = M + A) | Count both when calculating position |
| BUG-007 | Don't generate repetitive GGGGS codons | Use varied-codon sequences with k-mer uniqueness |
| BUG-008 | Never reuse AVD numbers for different sequences | Use `add_avd_uniqueness_check()` |
| BUG-011 | **AVD100-131 are DEAD. Never use, never use as design starting point, never copy sequence from them.** Old repetitive GGGGS + β-globin acceptor v1 encodes LINLFILLFSHPQ hydrophobic stretch → retained VP3-VHH crashes to pellet (confirmed EXP26000390). Use AVD250-281 (FP v3) for all VR4 IR VHH3 designs. | Check AVD number before copying any FP VR4 sequence |
| BUG-010 | GT/AG-null controls MUST use a Tier 6 or Tier 7 parent. A Tier 1 parent (~0-5% splicing) is already retention-dominant — knocking out GT/AG changes nothing measurable. The entire point of a GT/AG-null is to abolish splicing in a construct that *would otherwise splice efficiently*. **Before designing any GT/AG-null: (1) confirm parent is Tier 6/7, (2) confirm parent and null differ at exactly 2 positions (the GT and the AG), (3) confirm same construct series — same GGGGS version, same acceptor design, same HBB pad status.** | Check tier in tiered_panel_32_v7_summary.txt before picking parent |
| BUG-012 | **AATAAA is the poly(A) signal hexamer (CPSF recognition), NOT a spliceosomal branch point.** Branch point consensus = YTRAY (position 1 must be pyrimidine Y = C or T). AATAAA has A at position 1 (purine) — fails YTRAY. Mislabeling caused v1+v2 FP constructs (AVD134-165, AVD216-249) to be built with a biologically invalid "branch point." v3 fix: CTAAC (C=Y,T=T,A=R,A=branch-A,C=Y — YTRAY-compliant). All FP designs must use AVD250-281 (v3) or later. | Before annotating any splice acceptor branch point, verify it satisfies YTRAY — pos 1 must be C or T. Never assert AATAAA as a branch point. |
| BUG-013 | **GenBank KEYWORDS fields are assumed to be unreliable carryover from template files.** KEYWORDS are copied forward when cloning plasmids in Benchling and frequently describe a prior construct's design type, not the current one. Using KEYWORDS to determine construct identity caused selection of the wrong VCAP-195 lot (VOY-1637, miniZfr regulatable) instead of the correct AlwaysON lot (VOY-1631) in EXP26000387 run 3 — a wasted experiment. | **Never use KEYWORDS to determine construct identity, design type, or transgene name.** Always verify against authoritative production records: Construct-Information.csv, ORDER_REGISTRY, Benchling Vector Production entry, or equivalent. The transgene name in the production CSV is ground truth. |

See `LESSONS_LEARNED.md` for full bug and design decision documentation.

## Key Design Principles

| Principle | Summary | Reference |
|-----------|---------|-----------|
| DESIGN-010 | Do-no-harm safety checks should use tolerable thresholds, not absolute zero-regression rules | `LESSONS_LEARNED.md` |
| DESIGN-012 | When k-mer uniqueness conflicts with other constraints (splice sites, GC%), relax k by +1. Unique (k+1)-mers is strictly less restrictive — more codon choices, still eliminates repeat patterns vendors flag. v4 GGGGS: k=8→9 achieved all MaxEntScan < 0.0 while keeping Genscript FLASH "simple". | `LESSONS_LEARNED.md` |
| DESIGN-017 | Synthesized fragment internal restriction sites are irrelevant (both FLASH and Premium) — only the backbone is digested. Always use the closest unique enzyme pair flanking the modification. VP1n=KpnI/SmaI, VP2n=SmaI/BbvCI, VR4=AgeI/BsrGI. | `LESSONS_LEARNED.md` |
| DESIGN-018 | VP start codons outside the enzyme pair must match the backbone. Use `select_flash_backbone()` to pick optimal backbone + pair — never hardcode backbone assignments. | `LESSONS_LEARNED.md` |
| DESIGN-020 | Enzyme uniqueness: **prefer globally unique** (1 site in backbone AND ≤1 in final construct). Genscript rejects enzymes with multiple sites even in discarded fragment. Fall back to "conditionally unique" (1 outside excision) only when necessary. Upstream enzyme overlap into replacement region is OK (serves as junction). | `LESSONS_LEARNED.md` |
| DESIGN-026 | **GT/AG cryptic splice filters apply ONLY to intronic (non-splice-intended) sequence.** The exonic pad immediately downstream of the splice acceptor AG is exonic — an AG there cannot function as a competing acceptor (no room for branch point/PPT in 1-2 nt). **Never apply `"AG" not in exonic_pad` asserts.** GT is still forbidden in the exonic pad (cryptic donor = exon skipping risk). | `LESSONS_LEARNED.md` |
| DESIGN-027 | For any isolated ORF (no overlap with another ORF), always use **TAA** as the stop codon. TGA has near-cognate suppression risk in mammalian cells; TAG has intermediate read-through. TAA is first choice unless a specific constraint overrides (overlapping ORF, cloning junction). | `LESSONS_LEARNED.md` |
| DESIGN-028 | **Quote zip ingestion is mandatory before submission.** Every order folder MUST contain `Order_*_Files.zip` AND `quote_metadata.json` (generated by `parse_genscript_quote.py`). The xlsx "Plasmid unified" sheet is the authoritative source for per-construct prep scale (mg), concentration, supercoil/endotoxin spec, and unit price. `audit_order.py` Phase 0 BLOCKS submission (exit 1) without these. The submitted FASTA alone is insufficient. | `parse_genscript_quote.py`, `audit_order.py` |

## Additional Directories

| Directory | Contents |
|-----------|----------|
| `scripts/` | Build scripts for all construct campaigns (build_avd*.py, build_first_principles*.py, etc.) |
| `archive/` | Superseded scripts, reference PDFs/notes, historical synthetic fragments |
