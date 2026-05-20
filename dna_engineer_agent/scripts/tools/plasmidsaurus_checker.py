#!/usr/bin/env python3
"""
plasmidsaurus_checker.py — Align Plasmidsaurus long-read sequencing results against
reference constructs and produce structured PASS/WARN/FAIL/CRITICAL_FAIL reports.

Usage:
    # Check all samples in a zip against one reference
    python scripts/tools/plasmidsaurus_checker.py check results.zip --ref constructs/AVD548.gb

    # Batch mode: auto-match each sample to its reference by name prefix in constructs/
    python scripts/tools/plasmidsaurus_checker.py batch results.zip

    # Override auto-match for specific samples
    python scripts/tools/plasmidsaurus_checker.py batch results.zip \\
        --ref TP688-146202=constructs/TP688-CBA-hGBint-EGFP-T2A-AkaLuc.dna

    # Emit JSON for further analysis
    python scripts/tools/plasmidsaurus_checker.py check results.zip \\
        --ref constructs/AVD548.gb --json

KNOWN NANOPORE / PLASMIDSAURUS SEQUENCING PITFALLS:

  1. Homopolymers (>=5 nt): Nanopore systematically miscounts run length. DEL/INS variants
     in homopolymer context are LIKELY_ARTIFACT. The TSV `homopolymer` column (integer
     run-length) flags these positions. Never fail a sample on homopolymer indels alone.

  2. Dam/Dcm methylation sites: Dam = GATC (6mA on adenine); Dcm = CCWGG (CCAGG/CCTGG,
     5mC on internal cytosine). Both cause systematic Nanopore base-call shifts. Variants
     within ±2 bp of any Dam or Dcm site are flagged as INTERPRET_WITH_CAUTION.

  3. ITR palindromes: AAV ITR hairpins (high GC, palindromic secondary structure) cause
     (a) reduced read coverage and (b) alignment ambiguity between flip/flop orientations.
     CRITICAL: A deletion >=20 bp within an ITR region is CRITICAL_FAIL — this is a
     structural defect (e.g., c+c' arm ~42 bp deletion), not a sequencing artifact. Do
     not rationalize ITR deletions as noise.

  4. Low coverage (<10 reads): base calls are unreliable regardless of other flags.
     Variants at coverage < 10 are INTERPRET_WITH_CAUTION.

  5. 6mA / 5mC predicted methylation sites (TSV `predicted_methylation_site` column):
     methylation shifts basecalling distributions; flagged for human review.

  6. Short insertions in low-complexity sequence: typically homopolymer overcalls.
     Trust high-coverage non-homopolymer indels.
"""

import argparse
import json
import re
import sys
import tempfile
import zipfile
from dataclasses import asdict, dataclass, field
from pathlib import Path

try:
    import pandas as pd
    _HAS_PANDAS = True
except ImportError:
    _HAS_PANDAS = False

try:
    from Bio import Align, SeqIO
    from Bio.Seq import Seq
    _HAS_BIOPYTHON = True
except ImportError:
    _HAS_BIOPYTHON = False

import datetime
import html as html_module

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
AGENT_DIR = SCRIPT_DIR.parent.parent
CONSTRUCTS_DIR = AGENT_DIR / "constructs"
SEQUENCING_DIR = AGENT_DIR / "sequencing"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
LOW_COVERAGE_THRESHOLD = 10
MIN_LARGE_DELETION = 20        # bp — structural deletion threshold
CIRCULAR_ROTATIONS = 4         # rotation points to try for circular alignment
HOMOPOLYMER_MIN = 5            # minimum run length to flag
METHYLATION_RADIUS = 2         # ±bp around Dam/Dcm site to flag

DAM_PATTERN = "GATC"
DCM_PATTERNS = ["CCAGG", "CCTGG"]

VERDICT_ORDER = ["PASS", "WARN", "FAIL", "CRITICAL_FAIL"]


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------
@dataclass
class Variant:
    ref_pos: int           # 1-based in reference (-1 for insertions without ref pos)
    query_pos: int         # 1-based in FASTA consensus (-1 for deletions)
    ref_base: str
    query_base: str
    change_type: str       # SUB | DEL | INS
    coverage: int          # reads at query_pos (0 for DEL)
    mismatch_rate: float
    in_homopolymer: bool
    homopolymer_len: int
    predicted_methylation: str    # "" | "6mA" | "5mC"
    at_dam_dcm: bool
    feature_name: str
    reliability: str       # RELIABLE | INTERPRET_WITH_CAUTION | LIKELY_ARTIFACT


@dataclass
class DeletionSpan:
    ref_start: int
    ref_end: int
    size: int
    feature_name: str


@dataclass
class RegionSummary:
    feature_name: str
    ref_start: int
    ref_end: int
    n_variants: int
    n_reliable_variants: int
    coverage_min: int
    coverage_mean: float
    deletions: list = field(default_factory=list)
    verdict: str = "INTACT"
    notes: str = ""


@dataclass
class SampleResult:
    sample_name: str
    reference_file: str
    query_length: int
    ref_length: int
    circular_offset: int
    identity_pct: float
    quality: dict
    variants: list = field(default_factory=list)
    large_deletions: list = field(default_factory=list)
    region_summaries: list = field(default_factory=list)
    verdict: str = "PASS"
    verdict_reason: str = ""


# ---------------------------------------------------------------------------
# Zip parsing
# ---------------------------------------------------------------------------
_TYPE_MAP = {
    "fasta-files": "fasta",
    "per-base-data": "per_base",
    "summary-files": "summary",
    "genbank-files": "gbk",
}


def parse_zip(zip_path: Path) -> dict:
    """
    Extract Plasmidsaurus results zip and return {sample_name: {type_key: Path}}.
    File naming: {RUN_ID}_{N}_{SAMPLE_NAME}.{ext} inside {RUN_ID}_{type-key}/ directories.
    """
    tmpdir = Path(tempfile.mkdtemp(prefix="psaurus_"))
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(tmpdir)

    samples = {}
    for fpath in sorted(tmpdir.rglob("*")):
        if not fpath.is_file():
            continue
        parent = fpath.parent.name  # e.g. "22VPNM_fasta-files"
        file_type = None
        for suffix, key in _TYPE_MAP.items():
            if parent.endswith("_" + suffix):
                file_type = key
                break
        if file_type is None:
            continue
        parts = fpath.stem.split("_", 2)   # ["22VPNM", "1", "AVD548-1"]
        if len(parts) < 3:
            continue
        sample_name = parts[2]
        samples.setdefault(sample_name, {})[file_type] = fpath

    return samples


# ---------------------------------------------------------------------------
# Reference loading
# ---------------------------------------------------------------------------
def load_reference(ref_path: Path) -> tuple:
    """
    Load .gb/.gbk or .dna reference file via BioPython.
    Returns (seq_str_upper, features) where features = [(label, start_1based, end_1based, strand)].
    """
    ext = ref_path.suffix.lower()
    fmt = "genbank" if ext in (".gb", ".gbk", ".genbank") else "snapgene"
    record = SeqIO.read(str(ref_path), fmt)
    seq = str(record.seq).upper()

    features = []
    for feat in record.features:
        if feat.type == "source":
            continue
        label = ""
        for qk in ("label", "gene", "product", "note"):
            if qk in feat.qualifiers:
                label = feat.qualifiers[qk][0]
                break
        if not label:
            label = feat.type
        s = int(feat.location.start)  # 0-based
        e = int(feat.location.end)    # 0-based exclusive
        features.append((label, s + 1, e, feat.location.strand))  # 1-based inclusive

    return seq, features


# ---------------------------------------------------------------------------
# Summary parsing
# ---------------------------------------------------------------------------
def parse_summary(txt_path: Path) -> dict:
    text = txt_path.read_text(encoding="utf-8")
    q = {"monomer_pct_moles": None, "dimer_pct_moles": None, "ecoli_pct": None}
    m = re.search(r"moles\s+([\d.]+)\s+([\d.]+)", text)
    if m:
        q["monomer_pct_moles"] = float(m.group(1))
        q["dimer_pct_moles"] = float(m.group(2))
    m = re.search(r"E\.\s*coli\s+genomic\s+contamination:\s*([\d.]+)%", text)
    if m:
        q["ecoli_pct"] = float(m.group(1))
    return q


# ---------------------------------------------------------------------------
# Per-base TSV
# ---------------------------------------------------------------------------
def parse_per_base(tsv_path: Path):
    """Parse per-base TSV; return DataFrame indexed by pos (1-based int)."""
    df = pd.read_csv(str(tsv_path), sep="\t", dtype=str)
    for col in ["pos", "reads_all", "matches", "mismatches", "deletions",
                "insertions", "A", "C", "T", "G"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)
    if "homopolymer" in df.columns:
        df["homopolymer"] = pd.to_numeric(df["homopolymer"], errors="coerce").fillna(0).astype(int)
    if "low_conf" in df.columns:
        df["low_conf"] = df["low_conf"].str.strip().str.lower() == "true"
    if "predicted_methylation_site" in df.columns:
        df["predicted_methylation_site"] = df["predicted_methylation_site"].fillna("")
    df.set_index("pos", inplace=True)
    return df


def _tsv_row(df, qpos: int) -> dict:
    """Safely fetch per-base stats for a query position."""
    default = {"reads_all": 0, "homopolymer": 0, "mismatches": 0,
               "predicted_methylation_site": ""}
    if df is None or qpos not in df.index:
        return default
    row = df.loc[qpos]
    return {
        "reads_all": int(row.get("reads_all", 0)),
        "homopolymer": int(row.get("homopolymer", 0)),
        "mismatches": int(row.get("mismatches", 0)),
        "predicted_methylation_site": str(row.get("predicted_methylation_site", "")).strip(),
    }


# ---------------------------------------------------------------------------
# Methylation site detection
# ---------------------------------------------------------------------------
def find_methylation_sites(ref_seq: str) -> set:
    """
    Return set of 1-based ref positions within ±METHYLATION_RADIUS bp of any
    Dam (GATC) or Dcm (CCAGG, CCTGG) site on either strand.
    """
    ref_len = len(ref_seq)
    flagged = set()
    rev_comp = str(Seq(ref_seq).reverse_complement())

    for pattern in [DAM_PATTERN] + DCM_PATTERNS:
        for is_rev, strand_seq in enumerate([ref_seq, rev_comp]):
            start = 0
            while True:
                pos = strand_seq.find(pattern, start)
                if pos == -1:
                    break
                fwd_pos = (ref_len - pos - len(pattern)) if is_rev else pos  # 0-based
                for off in range(-METHYLATION_RADIUS, len(pattern) + METHYLATION_RADIUS):
                    actual = fwd_pos + off + 1  # 1-based
                    if 1 <= actual <= ref_len:
                        flagged.add(actual)
                start = pos + 1

    return flagged


# ---------------------------------------------------------------------------
# Circular alignment
# ---------------------------------------------------------------------------
def _make_aligner():
    a = Align.PairwiseAligner()
    a.mode = "global"
    a.match_score = 2
    a.mismatch_score = -1
    a.open_gap_score = -8
    a.extend_gap_score = -0.3
    return a


def _reconstruct_gapped(alignment, seq_a: str, seq_b: str) -> tuple:
    """
    Reconstruct gapped sequences from a PairwiseAlignment.
    alignment.aligned[0] = block coordinates in seq_a (first arg to aligner.align)
    alignment.aligned[1] = block coordinates in seq_b (second arg)
    Returns (gapped_a, gapped_b) where gaps are '-'.
    """
    a_coords = alignment.aligned[0]
    b_coords = alignment.aligned[1]
    ga, gb = [], []
    a_prev = b_prev = 0

    for (a_s, a_e), (b_s, b_e) in zip(a_coords, b_coords):
        if b_s > b_prev:           # gap in a (b has content a skipped)
            ga.append("-" * (b_s - b_prev))
            gb.append(seq_b[b_prev:b_s])
        if a_s > a_prev:           # gap in b (a has content b skipped)
            ga.append(seq_a[a_prev:a_s])
            gb.append("-" * (a_s - a_prev))
        ga.append(seq_a[a_s:a_e])
        gb.append(seq_b[b_s:b_e])
        a_prev, b_prev = a_e, b_e

    if b_prev < len(seq_b):
        ga.append("-" * (len(seq_b) - b_prev))
        gb.append(seq_b[b_prev:])
    if a_prev < len(seq_a):
        ga.append(seq_a[a_prev:])
        gb.append("-" * (len(seq_a) - a_prev))

    return "".join(ga), "".join(gb)


def align_circular(query: str, ref: str) -> tuple:
    """
    Align query to a circular reference.

    Step 1: Local alignment of the first 300 bp of query against doubled reference
            to pinpoint where the query starts in the reference (rotation offset).
            Using the first 300 bp avoids any large deletion in the anchor region.
    Step 2: Rotate reference to that offset and perform global alignment.

    This two-step approach eliminates the circular-junction artefact that occurs
    when the reference wrap-point falls in the middle of a gap region.

    Returns (aligned_query, aligned_ref, offset, identity_pct).
    aligned_ref uses the ROTATED reference; offset converts rotated positions back
    to unrotated reference positions via: actual = (rotated - 1 + offset) % len(ref) + 1.
    """
    ref_len = len(ref)
    doubled = ref + ref
    aligner = _make_aligner()

    # --- Step 1: find rotation offset via local anchor alignment ---
    anchor = query[:min(300, len(query))]
    aligner.mode = "local"
    local_alns = aligner.align(anchor, doubled)
    best_local = next(iter(local_alns))
    # aligned[1] = blocks in doubled (second argument to aligner.align)
    doubled_start = best_local.aligned[1][0][0]  # 0-based start in doubled ref
    offset = doubled_start % ref_len

    # --- Step 2: global alignment against rotated reference ---
    rotated_ref = ref[offset:] + ref[:offset]
    aligner.mode = "global"
    alignments = aligner.align(query, rotated_ref)
    best = next(iter(alignments))

    aligned_q, aligned_r = _reconstruct_gapped(best, query, rotated_ref)

    matches = sum(1 for a, b in zip(aligned_q, aligned_r) if a == b and a != "-")
    total = sum(1 for a, b in zip(aligned_q, aligned_r) if not (a == "-" and b == "-"))
    identity = (matches / total * 100) if total > 0 else 0.0

    return aligned_q, aligned_r, offset, identity


# ---------------------------------------------------------------------------
# Position mapping
# ---------------------------------------------------------------------------
def build_position_map(aligned_q: str, aligned_r: str, offset: int, ref_len: int) -> tuple:
    """
    Build bidirectional maps between FASTA (query) positions and actual reference positions.
    offset: the circular rotation applied to ref before alignment.
    Returns (q_to_r, r_to_q):
        q_to_r[query_pos_1based] = actual_ref_pos_1based  (None if insertion)
        r_to_q[actual_ref_pos_1based] = query_pos_1based  (None if deletion)
    """
    q_to_r, r_to_q = {}, {}
    q_pos = rotated_r_pos = 0

    for q_char, r_char in zip(aligned_q, aligned_r):
        if q_char != "-":
            q_pos += 1
        if r_char != "-":
            rotated_r_pos += 1

        actual_r = (rotated_r_pos - 1 + offset) % ref_len + 1

        if q_char != "-" and r_char != "-":
            q_to_r[q_pos] = actual_r
            r_to_q[actual_r] = q_pos
        elif q_char != "-":   # insertion in query
            q_to_r[q_pos] = None
        else:                  # deletion in query (r_char != '-')
            r_to_q[actual_r] = None

    return q_to_r, r_to_q


# ---------------------------------------------------------------------------
# Large deletion detection
# ---------------------------------------------------------------------------
def find_large_deletions(aligned_q: str, aligned_r: str, offset: int, ref_len: int,
                          ref_features: list) -> list:
    """
    Find runs of >=MIN_LARGE_DELETION consecutive gaps in aligned_q (deletions from query).
    Returns list of DeletionSpan with actual (unrotated) reference coordinates.
    """
    rotated_r_pos = 0
    col_ref = []   # actual ref pos for each column; -1 for insertion columns
    for r_char in aligned_r:
        if r_char != "-":
            rotated_r_pos += 1
            col_ref.append((rotated_r_pos - 1 + offset) % ref_len + 1)
        else:
            col_ref.append(-1)

    result = []
    i = 0
    n = len(aligned_q)
    while i < n:
        if aligned_q[i] == "-" and col_ref[i] != -1:
            j = i
            while j < n and aligned_q[j] == "-" and col_ref[j] != -1:
                j += 1
            del_refs = [col_ref[k] for k in range(i, j)]
            size = len(del_refs)
            if size >= MIN_LARGE_DELETION:
                feat = _feature_at(del_refs[0], ref_features) or _feature_at(del_refs[-1], ref_features)
                result.append(DeletionSpan(del_refs[0], del_refs[-1], size, feat or ""))
            i = j
        else:
            i += 1

    return result


def _feature_at(ref_pos: int, features: list) -> str:
    """Return name of first feature whose range contains ref_pos (1-based)."""
    for label, start, end, _ in features:
        if start <= ref_pos <= end:
            return label
    return ""


# ---------------------------------------------------------------------------
# Variant calling
# ---------------------------------------------------------------------------
def _large_del_set(large_deletions: list) -> set:
    """Build set of ref positions covered by large deletions (for exclusion)."""
    covered = set()
    for span in large_deletions:
        for p in range(span.ref_start, span.ref_end + 1):
            covered.add(p)
    return covered


def call_variants(aligned_q: str, aligned_r: str, per_base_df,
                  offset: int, ref_len: int, ref_features: list,
                  methylation_sites: set, large_deletions: list) -> list:
    """
    Walk aligned strings and call SUB / small DEL / INS variants.
    Positions covered by large_deletions are excluded (reported separately).
    Returns list of Variant.
    """
    large_del_ref_pos = _large_del_set(large_deletions)
    variants = []
    q_pos = rotated_r_pos = 0

    for q_char, r_char in zip(aligned_q, aligned_r):
        if q_char != "-":
            q_pos += 1
        if r_char != "-":
            rotated_r_pos += 1

        if q_char == r_char:
            continue  # exact match

        actual_r = (rotated_r_pos - 1 + offset) % ref_len + 1 if rotated_r_pos > 0 else -1

        if q_char == "-":
            change = "DEL"
            use_qpos = -1
        elif r_char == "-":
            change = "INS"
            use_qpos = q_pos
        else:
            change = "SUB"
            use_qpos = q_pos

        # Skip positions covered by a large deletion
        if actual_r in large_del_ref_pos:
            continue

        # Per-base TSV lookup
        row = _tsv_row(per_base_df, use_qpos) if use_qpos > 0 else _tsv_row(per_base_df, q_pos)
        coverage = row["reads_all"]
        hplen = row["homopolymer"]
        pred_meth = row["predicted_methylation_site"]
        mismatch_rate = (row["mismatches"] / coverage) if coverage > 0 else 0.0

        in_homopolymer = hplen >= HOMOPOLYMER_MIN
        at_dam_dcm = actual_r in methylation_sites

        feat = _feature_at(actual_r, ref_features) if actual_r > 0 else ""

        if in_homopolymer and change in ("DEL", "INS"):
            reliability = "LIKELY_ARTIFACT"
        elif coverage < LOW_COVERAGE_THRESHOLD:
            reliability = "INTERPRET_WITH_CAUTION"
        elif at_dam_dcm or pred_meth:
            reliability = "INTERPRET_WITH_CAUTION"
        else:
            reliability = "RELIABLE"

        variants.append(Variant(
            ref_pos=actual_r,
            query_pos=use_qpos,
            ref_base=r_char,
            query_base=q_char,
            change_type=change,
            coverage=coverage,
            mismatch_rate=round(mismatch_rate, 3),
            in_homopolymer=in_homopolymer,
            homopolymer_len=hplen,
            predicted_methylation=pred_meth,
            at_dam_dcm=at_dam_dcm,
            feature_name=feat,
            reliability=reliability,
        ))

    return variants


# ---------------------------------------------------------------------------
# Region summaries
# ---------------------------------------------------------------------------
def summarize_regions(ref_features: list, variants: list, large_deletions: list,
                      per_base_df, r_to_q: dict, ref_len: int) -> list:
    """
    Summarize alignment quality per annotated feature region.
    Returns list of RegionSummary.
    """
    seen = set()
    results = []

    for label, feat_start, feat_end, _ in ref_features:
        if label in seen:
            continue
        seen.add(label)

        feat_variants = [v for v in variants if feat_start <= v.ref_pos <= feat_end]
        n_reliable = sum(1 for v in feat_variants if v.reliability == "RELIABLE")

        feat_dels = [d for d in large_deletions
                     if d.ref_start <= feat_end and d.ref_end >= feat_start]

        # Coverage stats from TSV via r_to_q mapping
        coverages = []
        for rp in range(feat_start, feat_end + 1):
            qp = r_to_q.get(rp)
            if qp is not None and per_base_df is not None and qp in per_base_df.index:
                coverages.append(int(per_base_df.loc[qp, "reads_all"]))

        cov_min = min(coverages) if coverages else 0
        cov_mean = round(sum(coverages) / len(coverages), 1) if coverages else 0.0

        # Verdict
        if feat_dels:
            verdict = "CRITICAL_FAIL"
            notes = (f"Structural deletion ({feat_dels[0].size} bp, "
                     f"ref {feat_dels[0].ref_start}–{feat_dels[0].ref_end})")
        elif n_reliable > 3:
            verdict = "FAIL"
            notes = f"{n_reliable} reliable variants"
        elif n_reliable > 0:
            verdict = "WARN"
            notes = f"{n_reliable} reliable variant(s)"
        elif feat_variants and all(v.reliability != "RELIABLE" for v in feat_variants):
            verdict = "WARN"
            notes = f"{len(feat_variants)} ambiguous/artifact variant(s)"
        elif not coverages:
            verdict = "WARN"
            notes = "No coverage data for this region"
        else:
            verdict = "INTACT"
            notes = ""

        results.append(RegionSummary(
            feature_name=label,
            ref_start=feat_start,
            ref_end=feat_end,
            n_variants=len(feat_variants),
            n_reliable_variants=n_reliable,
            coverage_min=cov_min,
            coverage_mean=cov_mean,
            deletions=feat_dels,
            verdict=verdict,
            notes=notes,
        ))

    return results


# ---------------------------------------------------------------------------
# Overall verdict
# ---------------------------------------------------------------------------
def _upgrade_verdict(current: str, candidate: str) -> str:
    i_cur = VERDICT_ORDER.index(current) if current in VERDICT_ORDER else 0
    i_cand = VERDICT_ORDER.index(candidate) if candidate in VERDICT_ORDER else 0
    return VERDICT_ORDER[max(i_cur, i_cand)]


def determine_verdict(region_summaries: list, identity_pct: float) -> tuple:
    verdict = "PASS"
    reason = f"Identity {identity_pct:.2f}%, all regions intact"
    first_warn_reason = None

    for rs in region_summaries:
        verdict = _upgrade_verdict(verdict, rs.verdict)
        if rs.verdict == "CRITICAL_FAIL":
            reason = f"Structural deletion in '{rs.feature_name}': {rs.notes}"
            return verdict, reason
        if rs.verdict in ("FAIL", "WARN") and first_warn_reason is None:
            first_warn_reason = f"{rs.verdict} in '{rs.feature_name}': {rs.notes}"

    if verdict == "PASS" and identity_pct < 97.0:
        verdict = "FAIL"
        reason = f"Identity {identity_pct:.2f}% < 97% threshold"
    elif identity_pct < 99.0:
        verdict = _upgrade_verdict(verdict, "WARN")
        reason = f"Identity {identity_pct:.2f}% < 99%"
    elif first_warn_reason:
        reason = first_warn_reason

    return verdict, reason


# ---------------------------------------------------------------------------
# Stdout summary (brief)
# ---------------------------------------------------------------------------
_VERDICT_ICON = {
    "PASS":          "✓",
    "WARN":          "~",
    "FAIL":          "✗",
    "CRITICAL_FAIL": "!!",
}


def print_verdict_summary(results: list) -> None:
    """Print a brief per-sample verdict table to stdout."""
    print()
    print(f"  {'Sample':<25} {'Identity':>9}  {'Verdict'}")
    print(f"  {'-'*25} {'-'*9}  {'-'*30}")
    for r in results:
        icon = _VERDICT_ICON.get(r.verdict, "?")
        print(f"  {r.sample_name:<25} {r.identity_pct:>8.2f}%  {icon} {r.verdict}  —  {r.verdict_reason}")
    print()


def _is_itr(feature_name: str) -> bool:
    n = feature_name.upper()
    return "ITR" in n or "INVERTED TERMINAL" in n


# ---------------------------------------------------------------------------
# Output path helpers
# ---------------------------------------------------------------------------
def _run_id_from_zip(zip_path: Path) -> str:
    """Extract Plasmidsaurus run ID from zip filename: '22VPNM_results.zip' → '22VPNM'."""
    stem = zip_path.stem  # e.g. "22VPNM_results"
    return stem.split("_")[0]


def _run_slug(run_id: str, results: list, single_ref: Path = None) -> str:
    """
    Build a short descriptive slug for this run's folder and filename.
    check mode:  22VPNM-8x-AVD548
    batch mode:  FQL4T3-TP688-TP737-TP738  (up to 3 unique prefixes, then +N)
    """
    n = len(results)
    if single_ref is not None:
        ref_stem = single_ref.stem.split("-")[0]  # "AVD548" from "AVD548-CBA-..."
        return f"{run_id}-{n}x-{ref_stem}"
    else:
        prefixes = []
        seen = set()
        for r in results:
            p = r.sample_name.split("-")[0]
            if p not in seen:
                prefixes.append(p)
                seen.add(p)
        if len(prefixes) <= 3:
            slug = "-".join(prefixes)
        else:
            slug = "-".join(prefixes[:3]) + f"+{len(prefixes)-3}"
        return f"{run_id}-{slug}"


def _report_filename(run_id: str, results: list, single_ref: Path = None) -> str:
    """Build the HTML report filename (slug + .html)."""
    return _run_slug(run_id, results, single_ref) + ".html"


def _outdir(zip_path: Path, results: list = None, single_ref: Path = None,
            out_override: str = None) -> Path:
    """Return (and create) the output directory for this run."""
    if out_override:
        d = Path(out_override)
    else:
        run_id = _run_id_from_zip(zip_path)
        if results:
            slug = _run_slug(run_id, results, single_ref)
            d = SEQUENCING_DIR / slug
        else:
            d = SEQUENCING_DIR / run_id
    d.mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------
_VERDICT_COLOR = {
    "PASS":          ("#166534", "#dcfce7"),   # dark green text, light green bg
    "WARN":          ("#92400e", "#fef3c7"),   # amber
    "FAIL":          ("#991b1b", "#fee2e2"),   # red
    "CRITICAL_FAIL": ("#ffffff", "#7f1d1d"),   # white text, dark red bg
}

_REGION_COLOR = {
    "INTACT":        ("#166534", "#f0fdf4"),
    "WARN":          ("#92400e", "#fffbeb"),
    "FAIL":          ("#991b1b", "#fff1f2"),
    "CRITICAL_FAIL": ("#ffffff", "#7f1d1d"),
}

_RELIABILITY_COLOR = {
    "RELIABLE":                ("#166534", "#f0fdf4"),
    "INTERPRET_WITH_CAUTION":  ("#92400e", "#fffbeb"),
    "LIKELY_ARTIFACT":         ("#6b7280", "#f9fafb"),
}

_CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: ui-sans-serif, system-ui, -apple-system, sans-serif;
  font-size: 14px; color: #111827; background: #f9fafb;
  padding: 24px;
}
h1 { font-size: 22px; font-weight: 700; color: #1e3a5f; margin-bottom: 4px; }
h2 { font-size: 15px; font-weight: 600; color: #374151; margin: 20px 0 8px; }
.meta { color: #6b7280; font-size: 12px; margin-bottom: 24px; }
.meta span { margin-right: 16px; }

/* Summary table */
.summary-table { width: 100%; border-collapse: collapse; margin-bottom: 28px;
  background: white; border-radius: 8px; overflow: hidden;
  box-shadow: 0 1px 3px rgba(0,0,0,.08); }
.summary-table th { background: #1e3a5f; color: white; text-align: left;
  padding: 8px 12px; font-size: 12px; font-weight: 600; text-transform: uppercase;
  letter-spacing: .05em; }
.summary-table td { padding: 8px 12px; border-bottom: 1px solid #e5e7eb;
  vertical-align: middle; }
.summary-table tr:last-child td { border-bottom: none; }
.summary-table tr:hover td { background: #f0f9ff; }

/* Verdict badge */
.badge { display: inline-block; padding: 2px 8px; border-radius: 9999px;
  font-size: 11px; font-weight: 700; letter-spacing: .03em; white-space: nowrap; }

/* Sample cards */
details { background: white; border-radius: 8px; margin-bottom: 12px;
  box-shadow: 0 1px 3px rgba(0,0,0,.08); overflow: hidden; }
details[open] { box-shadow: 0 2px 8px rgba(0,0,0,.12); }
summary { display: flex; align-items: center; gap: 10px; padding: 12px 16px;
  cursor: pointer; user-select: none; font-weight: 600; font-size: 14px;
  list-style: none; }
summary::-webkit-details-marker { display: none; }
summary::before { content: "▶"; font-size: 10px; color: #9ca3af;
  transition: transform .15s; min-width: 10px; }
details[open] summary::before { transform: rotate(90deg); }
summary:hover { background: #f8fafc; }
.sample-body { padding: 16px; border-top: 1px solid #e5e7eb; }

/* Alert box for structural deletions */
.deletion-alert { background: #fef2f2; border: 1px solid #fecaca;
  border-left: 4px solid #dc2626; border-radius: 6px;
  padding: 12px 14px; margin-bottom: 14px; }
.deletion-alert h3 { font-size: 13px; color: #991b1b; margin-bottom: 6px; font-weight: 700; }
.deletion-alert .del-row { font-family: ui-monospace, monospace; font-size: 12px;
  color: #7f1d1d; margin: 3px 0; }
.deletion-alert .itr-note { font-size: 11px; color: #b91c1c; font-style: italic;
  margin-top: 6px; }

/* No-issues box */
.no-issues { background: #f0fdf4; border: 1px solid #bbf7d0;
  border-left: 4px solid #16a34a; border-radius: 6px;
  padding: 10px 14px; margin-bottom: 14px; font-size: 12px; color: #166534; }

/* Meta grid */
.meta-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(180px,1fr));
  gap: 10px; margin-bottom: 14px; }
.meta-card { background: #f8fafc; border: 1px solid #e5e7eb; border-radius: 6px;
  padding: 8px 12px; }
.meta-card label { font-size: 10px; font-weight: 600; color: #6b7280;
  text-transform: uppercase; letter-spacing: .05em; display: block; }
.meta-card .val { font-size: 13px; font-weight: 600; color: #111827; margin-top: 2px; }

/* Data tables */
.data-table { width: 100%; border-collapse: collapse; font-size: 12px;
  margin-bottom: 14px; }
.data-table th { background: #f3f4f6; color: #374151; font-weight: 600;
  padding: 6px 10px; text-align: left; border-bottom: 1px solid #e5e7eb; }
.data-table td { padding: 5px 10px; border-bottom: 1px solid #f3f4f6;
  font-family: ui-monospace, monospace; }
.data-table tr:hover td { background: #fafafa; }
.data-table tr:last-child td { border-bottom: none; }
.section-label { font-size: 11px; font-weight: 700; color: #6b7280;
  text-transform: uppercase; letter-spacing: .06em; margin: 14px 0 6px; }
"""


def _badge(verdict: str) -> str:
    tc, bg = _VERDICT_COLOR.get(verdict, ("#111827", "#f3f4f6"))
    label = verdict.replace("_", " ")
    return (f'<span class="badge" style="color:{tc};background:{bg}">'
            f'{label}</span>')


def _region_cell(verdict: str) -> str:
    tc, bg = _REGION_COLOR.get(verdict, ("#111827", "#f3f4f6"))
    return (f'<span style="color:{tc};background:{bg};padding:1px 6px;'
            f'border-radius:4px;font-size:11px;font-weight:600">{verdict}</span>')


def _reliability_cell(r: str) -> str:
    tc, bg = _RELIABILITY_COLOR.get(r, ("#111827", "#f3f4f6"))
    return (f'<span style="color:{tc};background:{bg};padding:1px 6px;'
            f'border-radius:4px;font-size:11px">{r}</span>')


def _h(s) -> str:
    """HTML-escape a value."""
    return html_module.escape(str(s))


def _sample_html(r: SampleResult, open_default: bool = False) -> str:
    q = r.quality
    mono = f"{q.get('monomer_pct_moles', '?')}%" if q.get("monomer_pct_moles") is not None else "—"
    ecoli = f"{q.get('ecoli_pct', '?')}%" if q.get("ecoli_pct") is not None else "—"
    ref_name = Path(r.reference_file).name
    open_attr = " open" if open_default else ""

    # Verdict color for summary bar
    tc, bg = _VERDICT_COLOR.get(r.verdict, ("#111827", "#f3f4f6"))
    border_color = bg if r.verdict == "CRITICAL_FAIL" else tc

    parts = [
        f'<details{open_attr} style="border-left:4px solid {border_color}">',
        f'<summary>{_badge(r.verdict)}'
        f'<span style="flex:1">{_h(r.sample_name)}</span>'
        f'<span style="font-size:12px;color:#6b7280;font-weight:400">{_h(r.verdict_reason)}</span>'
        f'</summary>',
        '<div class="sample-body">',
    ]

    # Meta grid
    parts.append('<div class="meta-grid">')
    for label, val in [
        ("Reference", ref_name),
        ("Query length", f"{r.query_length:,} bp"),
        ("Ref length", f"{r.ref_length:,} bp"),
        ("Identity", f"{r.identity_pct:.2f}%"),
        ("Circ. offset", f"{r.circular_offset:,} bp"),
        ("Monomer purity", mono),
        ("E. coli contam", ecoli),
    ]:
        parts.append(
            f'<div class="meta-card"><label>{_h(label)}</label>'
            f'<div class="val">{_h(val)}</div></div>'
        )
    parts.append('</div>')

    # Structural deletions
    if r.large_deletions:
        parts.append('<div class="deletion-alert">')
        parts.append(f'<h3>⚠ Structural Deletion{"s" if len(r.large_deletions)>1 else ""}'
                     f' (&ge;{MIN_LARGE_DELETION} bp)</h3>')
        for span in r.large_deletions:
            feat_str = f" [{span.feature_name}]" if span.feature_name else ""
            parts.append(
                f'<div class="del-row">DEL  ref {span.ref_start}–{span.ref_end}'
                f'  {span.size} bp{_h(feat_str)}</div>'
            )
            if span.feature_name and _is_itr(span.feature_name):
                parts.append(
                    '<div class="itr-note">ITR structural deletion — '
                    'not a sequencing artifact; confirmed real</div>'
                )
        parts.append('</div>')
    else:
        parts.append(
            f'<div class="no-issues">No structural deletions (&ge;{MIN_LARGE_DELETION} bp)</div>'
        )

    # Variants
    non_art = [v for v in r.variants if v.reliability != "LIKELY_ARTIFACT"]
    art = [v for v in r.variants if v.reliability == "LIKELY_ARTIFACT"]
    if r.variants:
        parts.append(
            f'<div class="section-label">Variants '
            f'({len(r.variants)} total · {len(non_art)} non-artifact · '
            f'{len(art)} LIKELY_ARTIFACT)</div>'
        )
        parts.append(
            '<table class="data-table"><thead><tr>'
            '<th>ref pos</th><th>type</th><th>ref</th><th>query</th>'
            '<th>cov</th><th>homopol</th><th>methyl</th><th>feature</th><th>reliability</th>'
            '</tr></thead><tbody>'
        )
        for v in sorted(r.variants, key=lambda x: x.ref_pos):
            hp = str(v.homopolymer_len) if v.in_homopolymer else "—"
            meth = v.predicted_methylation or ("Dam/Dcm" if v.at_dam_dcm else "—")
            parts.append(
                f'<tr><td>{v.ref_pos if v.ref_pos > 0 else "—"}</td>'
                f'<td><b>{_h(v.change_type)}</b></td>'
                f'<td>{_h(v.ref_base)}</td><td>{_h(v.query_base)}</td>'
                f'<td>{v.coverage}</td><td>{hp}</td><td>{_h(meth)}</td>'
                f'<td>{_h(v.feature_name)}</td>'
                f'<td>{_reliability_cell(v.reliability)}</td></tr>'
            )
        parts.append('</tbody></table>')
    else:
        parts.append('<div class="no-issues">No single-base variants</div>')

    # Region summaries
    non_intact = [rs for rs in r.region_summaries if rs.verdict != "INTACT"]
    intact = [rs for rs in r.region_summaries if rs.verdict == "INTACT"]

    parts.append(
        f'<div class="section-label">Region Summaries '
        f'({len(r.region_summaries)} features)</div>'
    )
    parts.append(
        '<table class="data-table"><thead><tr>'
        '<th>Feature</th><th>Verdict</th><th>cov min</th><th>cov mean</th>'
        '<th>Variants</th><th>Notes</th>'
        '</tr></thead><tbody>'
    )
    # Non-intact first, then intact (collapsed visually by ordering)
    for rs in non_intact + intact:
        parts.append(
            f'<tr><td>{_h(rs.feature_name)}</td>'
            f'<td>{_region_cell(rs.verdict)}</td>'
            f'<td>{rs.coverage_min}</td><td>{rs.coverage_mean:.1f}</td>'
            f'<td>{rs.n_variants}</td><td>{_h(rs.notes)}</td></tr>'
        )
    parts.append('</tbody></table>')

    parts.append('</div></details>')
    return "\n".join(parts)


def generate_html_report(results: list, run_id: str, zip_path: Path) -> str:
    """Generate a self-contained HTML report for a Plasmidsaurus run."""
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    n_cf = sum(1 for r in results if r.verdict == "CRITICAL_FAIL")
    n_fail = sum(1 for r in results if r.verdict == "FAIL")
    n_warn = sum(1 for r in results if r.verdict == "WARN")
    n_pass = sum(1 for r in results if r.verdict == "PASS")

    # Summary counts banner
    counts = []
    if n_cf:
        counts.append(f'<span style="color:#7f1d1d;font-weight:700">'
                      f'{n_cf} CRITICAL_FAIL</span>')
    if n_fail:
        counts.append(f'<span style="color:#991b1b;font-weight:700">{n_fail} FAIL</span>')
    if n_warn:
        counts.append(f'<span style="color:#92400e;font-weight:700">{n_warn} WARN</span>')
    if n_pass:
        counts.append(f'<span style="color:#166534;font-weight:700">{n_pass} PASS</span>')
    counts_html = " &nbsp;·&nbsp; ".join(counts) if counts else "—"

    summary_rows = []
    for r in results:
        ref_name = Path(r.reference_file).name
        summary_rows.append(
            f'<tr><td>{_h(r.sample_name)}</td>'
            f'<td style="font-size:11px;color:#6b7280">{_h(ref_name)}</td>'
            f'<td style="text-align:right">{r.identity_pct:.2f}%</td>'
            f'<td>{_badge(r.verdict)}</td>'
            f'<td style="color:#6b7280;font-size:11px">{_h(r.verdict_reason)}</td></tr>'
        )

    # Auto-open CRITICAL_FAIL / FAIL samples; collapse WARN / PASS
    sample_htmls = []
    for r in results:
        auto_open = r.verdict in ("CRITICAL_FAIL", "FAIL")
        sample_htmls.append(_sample_html(r, open_default=auto_open))

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Plasmidsaurus Report — {_h(run_id)}</title>
<style>{_CSS}</style>
</head>
<body>

<h1>Plasmidsaurus Sequencing Report</h1>
<div class="meta">
  <span><b>Run ID:</b> {_h(run_id)}</span>
  <span><b>Zip:</b> {_h(zip_path.name)}</span>
  <span><b>Samples:</b> {len(results)}</span>
  <span><b>Generated:</b> {now}</span>
  <span>{counts_html}</span>
</div>

<h2>Summary</h2>
<table class="summary-table">
  <thead><tr>
    <th>Sample</th><th>Reference</th><th style="text-align:right">Identity</th>
    <th>Verdict</th><th>Reason</th>
  </tr></thead>
  <tbody>{"".join(summary_rows)}</tbody>
</table>

<h2>Sample Details</h2>
{"".join(sample_htmls)}

</body>
</html>"""


# ---------------------------------------------------------------------------
# Auto-match reference
# ---------------------------------------------------------------------------
def auto_match_reference(sample_name: str, overrides: dict,
                          constructs_dir: Path = CONSTRUCTS_DIR) -> Path:
    """
    Match sample to a reference .gb/.dna by prefix.
    'AVD548-1' -> prefix 'AVD548' -> constructs/AVD548*.gb
    overrides: {sample_name: Path} from --ref flags
    """
    if sample_name in overrides:
        p = Path(overrides[sample_name])
        if not p.exists():
            raise FileNotFoundError(f"Override reference not found: {p}")
        return p

    prefix = sample_name.split("-")[0]
    matches = []
    for ext in (".gb", ".gbk", ".dna"):
        matches.extend(constructs_dir.glob(f"{prefix}*{ext}"))
    matches = sorted(set(matches))

    if not matches:
        raise ValueError(
            f"No reference for '{sample_name}' (prefix '{prefix}') in {constructs_dir}. "
            f"Use --ref {sample_name}=PATH to specify."
        )
    if len(matches) > 1:
        names = [m.name for m in matches]
        raise ValueError(
            f"Ambiguous reference for '{sample_name}' (prefix '{prefix}'): {names}. "
            f"Use --ref {sample_name}=PATH to specify."
        )
    return matches[0]


# ---------------------------------------------------------------------------
# Core pipeline: process one sample
# ---------------------------------------------------------------------------
def process_sample(sample_name: str, files: dict, ref_path: Path) -> SampleResult:
    """Full pipeline for one sample: align, annotate, summarize, verdict."""
    # Load FASTA consensus
    fasta_path = files.get("fasta")
    if not fasta_path:
        raise FileNotFoundError(f"No FASTA file for sample '{sample_name}'")
    fasta_record = SeqIO.read(str(fasta_path), "fasta")
    query_seq = str(fasta_record.seq).upper()

    # Load reference
    ref_seq, ref_features = load_reference(ref_path)

    # Parse quality summary (optional)
    quality = {}
    if "summary" in files:
        try:
            quality = parse_summary(files["summary"])
        except Exception:
            pass

    # Parse per-base TSV (optional)
    per_base_df = None
    if _HAS_PANDAS and "per_base" in files:
        try:
            per_base_df = parse_per_base(files["per_base"])
        except Exception:
            pass

    # Methylation sites in reference
    methylation_sites = find_methylation_sites(ref_seq)

    # Circular alignment
    aligned_q, aligned_r, offset, identity = align_circular(query_seq, ref_seq)

    # Position maps
    q_to_r, r_to_q = build_position_map(aligned_q, aligned_r, offset, len(ref_seq))

    # Large deletions
    large_dels = find_large_deletions(aligned_q, aligned_r, offset, len(ref_seq), ref_features)

    # Single-base variants (excluding large deletion positions)
    variants = call_variants(aligned_q, aligned_r, per_base_df,
                             offset, len(ref_seq), ref_features,
                             methylation_sites, large_dels)

    # Region summaries
    region_summaries = summarize_regions(ref_features, variants, large_dels,
                                         per_base_df, r_to_q, len(ref_seq))

    # Overall verdict
    verdict, reason = determine_verdict(region_summaries, identity)

    return SampleResult(
        sample_name=sample_name,
        reference_file=str(ref_path),
        query_length=len(query_seq),
        ref_length=len(ref_seq),
        circular_offset=offset,
        identity_pct=round(identity, 4),
        quality=quality,
        variants=variants,
        large_deletions=large_dels,
        region_summaries=region_summaries,
        verdict=verdict,
        verdict_reason=reason,
    )


# ---------------------------------------------------------------------------
# CLI commands
# ---------------------------------------------------------------------------
def cmd_check(args) -> int:
    """check subcommand: all samples in zip vs. single reference."""
    zip_path = Path(args.zip)
    ref_path = Path(args.ref)

    if not zip_path.exists():
        print(f"ERROR: zip not found: {zip_path}", file=sys.stderr)
        return 1
    if not ref_path.exists():
        print(f"ERROR: reference not found: {ref_path}", file=sys.stderr)
        return 1

    print(f"Loading zip: {zip_path.name}")
    samples = parse_zip(zip_path)
    if not samples:
        print("ERROR: no samples found in zip", file=sys.stderr)
        return 1

    print(f"Found {len(samples)} sample(s): {', '.join(sorted(samples))}")
    print()

    results = []
    exit_code = 0
    for sample_name in sorted(samples):
        print(f"  Processing {sample_name}...")
        try:
            result = process_sample(sample_name, samples[sample_name], ref_path)
        except Exception as e:
            print(f"  ERROR: {e}", file=sys.stderr)
            exit_code = 1
            continue
        results.append(result)
        if result.verdict in ("FAIL", "CRITICAL_FAIL"):
            exit_code = 1

    print()
    if args.json:
        print(json.dumps([_result_to_dict(r) for r in results], indent=2))
    else:
        print_verdict_summary(results)
        run_id = _run_id_from_zip(zip_path)
        outdir = _outdir(zip_path, results, single_ref=ref_path, out_override=getattr(args, "out", None))
        fname = _report_filename(run_id, results, single_ref=ref_path)
        html_path = outdir / fname
        html_path.write_text(generate_html_report(results, run_id, zip_path), encoding="utf-8")
        print(f"Report saved → {html_path}")

    return exit_code


def cmd_batch(args) -> int:
    """batch subcommand: auto-match each sample to its reference by prefix."""
    zip_path = Path(args.zip)
    if not zip_path.exists():
        print(f"ERROR: zip not found: {zip_path}", file=sys.stderr)
        return 1

    # Parse --ref overrides: "SAMPLENAME=path" pairs
    overrides = {}
    for raw in (args.ref or []):
        if "=" in raw:
            sname, path = raw.split("=", 1)
            overrides[sname.strip()] = path.strip()
        else:
            print(f"ERROR: --ref override must be SAMPLENAME=PATH, got: {raw}", file=sys.stderr)
            return 1

    constructs_dir = Path(args.constructs) if args.constructs else CONSTRUCTS_DIR

    print(f"Loading zip: {zip_path.name}")
    samples = parse_zip(zip_path)
    if not samples:
        print("ERROR: no samples found in zip", file=sys.stderr)
        return 1

    print(f"Found {len(samples)} sample(s): {', '.join(sorted(samples))}")
    print()

    results = []
    exit_code = 0
    for sample_name in sorted(samples):
        try:
            ref_path = auto_match_reference(sample_name, overrides, constructs_dir)
        except (ValueError, FileNotFoundError) as e:
            print(f"  SKIP {sample_name}: {e}", file=sys.stderr)
            exit_code = 1
            continue

        print(f"  Processing {sample_name} → {ref_path.name}...")
        try:
            result = process_sample(sample_name, samples[sample_name], ref_path)
        except Exception as e:
            print(f"  ERROR {sample_name}: {e}", file=sys.stderr)
            exit_code = 1
            continue
        results.append(result)
        if result.verdict in ("FAIL", "CRITICAL_FAIL"):
            exit_code = 1

    print()
    if args.json:
        print(json.dumps([_result_to_dict(r) for r in results], indent=2))
    else:
        print_verdict_summary(results)
        run_id = _run_id_from_zip(zip_path)
        outdir = _outdir(zip_path, results, out_override=getattr(args, "out", None))
        fname = _report_filename(run_id, results)
        html_path = outdir / fname
        html_path.write_text(generate_html_report(results, run_id, zip_path), encoding="utf-8")
        print(f"Report saved → {html_path}")

    return exit_code


def _result_to_dict(result: SampleResult) -> dict:
    """Convert SampleResult to JSON-serializable dict."""
    d = asdict(result)
    # Convert DeletionSpan/RegionSummary nested dataclasses (asdict handles recursion)
    return d


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main():
    if not _HAS_BIOPYTHON:
        print("ERROR: biopython is required. Install with: pip install biopython", file=sys.stderr)
        sys.exit(1)
    if not _HAS_PANDAS:
        print("WARNING: pandas not available — per-base quality data will be skipped.", file=sys.stderr)

    parser = argparse.ArgumentParser(
        description="Align Plasmidsaurus long-read results against reference constructs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # check subcommand
    p_check = sub.add_parser("check", help="Check all samples in zip vs. one reference")
    p_check.add_argument("zip", help="Plasmidsaurus results zip")
    p_check.add_argument("--ref", required=True, help="Reference construct (.gb or .dna)")
    p_check.add_argument("--json", action="store_true", help="Output JSON instead of text")
    p_check.add_argument("--out", metavar="DIR", help="Output directory for HTML report (default: sequencing/{RUN_ID}/)")

    # batch subcommand
    p_batch = sub.add_parser("batch",
        help="Auto-match each sample to reference by name prefix in constructs/")
    p_batch.add_argument("zip", help="Plasmidsaurus results zip")
    p_batch.add_argument("--ref", action="append", metavar="SAMPLE=PATH",
                         help="Override auto-match for specific sample (repeatable)")
    p_batch.add_argument("--constructs", help=f"Constructs directory (default: {CONSTRUCTS_DIR})")
    p_batch.add_argument("--json", action="store_true", help="Output JSON instead of text")
    p_batch.add_argument("--out", metavar="DIR", help="Output directory for HTML report (default: sequencing/{RUN_ID}/)")

    args = parser.parse_args()

    if args.command == "check":
        sys.exit(cmd_check(args))
    elif args.command == "batch":
        sys.exit(cmd_batch(args))


if __name__ == "__main__":
    main()
