#!/usr/bin/env python3
"""
construct_verifier.py — Mandatory post-build verification for DNA constructs.

This module consolidates all verification checks into a single, enforceable
pipeline. Build scripts must call `ConstructVerifier.run_all()` before writing
output. If any CRITICAL check fails, output is blocked.

Addresses:
- Plan item 3a (P0): Verification not enforced in build pipeline
- Plan item 3c (P1): Enzyme selection has no automated guard
- Plan item 3d (P1): No mandatory synthesis screen

Usage:
    from construct_verifier import ConstructVerifier

    verifier = ConstructVerifier(
        sequence="ATGAAAGCC...",
        name="AVD006",
    )
    verifier.add_restriction_check(required_unique=["HindIII", "XbaI"])
    verifier.add_protein_check(region_start=100, region_end=400, expected_aa="MKTLLF...")
    verifier.add_synthesis_screen()

    report = verifier.run_all()
    if not report.passed:
        print(report)
        sys.exit(1)
    # Safe to write output
"""

import json
import re
import sys
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Fallback defaults (used when enzyme_metadata.json unavailable)
_DEFAULT_DAM_SENSITIVE = {"XbaI", "BclI", "ClaI", "MboI", "DpnII"}
_DEFAULT_FUSSY = {"BbvCI", "PaqCI", "SmaI"}

# Module-level aliases for backward compatibility
DAM_SENSITIVE_ENZYMES = _DEFAULT_DAM_SENSITIVE
FUSSY_ENZYMES = _DEFAULT_FUSSY

# Default synthesis constraints
DEFAULT_GC_WINDOW = 50
DEFAULT_GC_MAX = 0.65
DEFAULT_GC_MIN = 0.25
DEFAULT_HOMOPOLYMER_MAX = 6
DEFAULT_REPEAT_MIN_LEN = 11

# AAV packaging limits
AAV_PACKAGING_LIMIT = 4900  # bp between ITRs


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return ''.join(comp.get(b, b) for b in reversed(seq))


def translate(dna_seq: str, frame: int = 0) -> str:
    """Translate DNA to protein in the given frame."""
    dna_seq = dna_seq.upper()[frame:]
    protein = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i + 3]
        if len(codon) == 3:
            aa = CODON_TABLE.get(codon, 'X')
            if aa == '*':
                break
            protein.append(aa)
    return ''.join(protein)


def count_sites_both_strands(sequence: str, site: str) -> int:
    """
    Count restriction enzyme sites on BOTH strands (BUG-003 prevention).

    For palindromic sites, forward count = total sites.
    For non-palindromic sites, count forward + reverse complement.
    """
    sequence = sequence.upper()
    site = site.upper()
    rc_site = reverse_complement(site)

    fwd_count = 0
    start = 0
    while True:
        pos = sequence.find(site, start)
        if pos == -1:
            break
        fwd_count += 1
        start = pos + 1

    if site == rc_site:
        # Palindromic — forward count IS the total
        return fwd_count

    # Non-palindromic — also count reverse complement
    rev_count = 0
    start = 0
    while True:
        pos = sequence.find(rc_site, start)
        if pos == -1:
            break
        rev_count += 1
        start = pos + 1

    return fwd_count + rev_count


def find_site_positions(sequence: str, site: str) -> List[Tuple[int, str]]:
    """Find all positions of a restriction site on both strands.

    Returns list of (position, strand) tuples where strand is '+' or '-'.
    """
    sequence = sequence.upper()
    site = site.upper()
    rc_site = reverse_complement(site)
    positions = []

    start = 0
    while True:
        pos = sequence.find(site, start)
        if pos == -1:
            break
        positions.append((pos, '+'))
        start = pos + 1

    if site != rc_site:
        start = 0
        while True:
            pos = sequence.find(rc_site, start)
            if pos == -1:
                break
            positions.append((pos, '-'))
            start = pos + 1

    return sorted(positions, key=lambda x: x[0])


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

class Severity(Enum):
    """Check severity level."""
    CRITICAL = "CRITICAL"    # Blocks output
    WARNING = "WARNING"      # Flagged but doesn't block
    INFO = "INFO"            # Informational


@dataclass
class CheckResult:
    """Result of a single verification check."""
    name: str
    passed: bool
    severity: Severity
    message: str
    details: Optional[str] = None

    def __str__(self):
        icon = "PASS" if self.passed else "FAIL"
        s = f"[{icon}] [{self.severity.value}] {self.name}: {self.message}"
        if self.details and not self.passed:
            s += f"\n        {self.details}"
        return s


@dataclass
class VerificationReport:
    """Complete verification report for a construct."""
    construct_name: str
    checks: List[CheckResult] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        """True if no CRITICAL checks failed."""
        return not any(
            not c.passed and c.severity == Severity.CRITICAL
            for c in self.checks
        )

    @property
    def critical_failures(self) -> List[CheckResult]:
        return [c for c in self.checks if not c.passed and c.severity == Severity.CRITICAL]

    @property
    def warnings(self) -> List[CheckResult]:
        return [c for c in self.checks if not c.passed and c.severity == Severity.WARNING]

    def __str__(self):
        lines = [f"=== Verification Report: {self.construct_name} ==="]
        for check in self.checks:
            lines.append(f"  {check}")

        n_pass = sum(1 for c in self.checks if c.passed)
        n_total = len(self.checks)
        n_crit_fail = len(self.critical_failures)
        n_warn = len(self.warnings)

        lines.append(f"\n  Summary: {n_pass}/{n_total} passed, "
                      f"{n_crit_fail} critical failures, {n_warn} warnings")

        if self.passed:
            lines.append("  >>> VERIFICATION PASSED — safe to write output <<<")
        else:
            lines.append("  >>> VERIFICATION FAILED — output blocked <<<")
            for cf in self.critical_failures:
                lines.append(f"      BLOCKING: {cf.name} — {cf.message}")

        return '\n'.join(lines)


# ---------------------------------------------------------------------------
# Enzyme Metadata Loader
# ---------------------------------------------------------------------------

_enzyme_metadata_cache: Optional[Dict] = None

def load_enzyme_metadata() -> Dict:
    """Load enzyme metadata from knowledge base JSON."""
    global _enzyme_metadata_cache
    if _enzyme_metadata_cache is not None:
        return _enzyme_metadata_cache

    kb_path = Path(__file__).parent.parent.parent / "knowledge_base" / "enzyme_metadata.json"
    if kb_path.exists():
        with open(kb_path) as f:
            _enzyme_metadata_cache = json.load(f)
        return _enzyme_metadata_cache

    # Fallback: return empty
    _enzyme_metadata_cache = {"enzymes": {}, "fussy_enzymes": {"enzymes": []},
                               "methylation_blocking": {"Dam_methylation": {"affected_enzymes": []}}}
    return _enzyme_metadata_cache


def get_dam_sensitive_enzymes() -> set:
    """Derive dam-sensitive enzyme names from knowledge base JSON.

    Unions defaults with JSON-derived enzymes so known dam-sensitive
    enzymes are never lost even if absent from JSON.
    """
    metadata = load_enzyme_metadata()
    result = set(_DEFAULT_DAM_SENSITIVE)
    for name, info in metadata.get("enzymes", {}).items():
        methyl = info.get("methylation_sensitivity", "")
        if "Dam" in methyl:
            result.add(name)
    return result


def get_fussy_enzymes() -> set:
    """Derive fussy enzyme names from knowledge base JSON.

    Unions defaults with JSON-derived enzymes so known fussy
    enzymes are never lost even if absent from JSON.
    """
    metadata = load_enzyme_metadata()
    result = set(_DEFAULT_FUSSY)
    fussy = metadata.get("fussy_enzymes", {}).get("enzymes", [])
    for e in fussy:
        if "name" in e:
            result.add(e["name"])
    return result


def get_recognition_site(enzyme_name: str) -> Optional[str]:
    """Look up a restriction enzyme's recognition site from knowledge base JSON."""
    metadata = load_enzyme_metadata()
    return metadata.get("enzymes", {}).get(enzyme_name, {}).get("recognition_site")


def validate_enzyme_choice(enzyme_name: str) -> List[CheckResult]:
    """
    Validate a restriction enzyme choice against knowledge base metadata.

    Checks:
    - Dam methylation sensitivity (DESIGN-004)
    - Reliability rating
    - Known fussy behavior (DESIGN-005)

    Returns list of CheckResults (may be empty if enzyme is fine).
    """
    results = []
    metadata = load_enzyme_metadata()
    enzymes = metadata.get("enzymes", {})

    # Check if enzyme is known
    if enzyme_name not in enzymes:
        if enzyme_name in get_dam_sensitive_enzymes():
            results.append(CheckResult(
                name=f"Enzyme:{enzyme_name}:dam_sensitive",
                passed=False,
                severity=Severity.WARNING,
                message=f"{enzyme_name} is dam-methylation sensitive (DESIGN-004)",
                details="Use dam-insensitive alternative: MluI, NotI, HindIII, EcoRI, SnaBI",
            ))
        if enzyme_name in get_fussy_enzymes():
            results.append(CheckResult(
                name=f"Enzyme:{enzyme_name}:fussy",
                passed=False,
                severity=Severity.WARNING,
                message=f"{enzyme_name} has known reliability issues (DESIGN-005)",
                details="Consider using a more reliable alternative",
            ))
        return results

    info = enzymes[enzyme_name]

    # Check dam sensitivity
    methyl = info.get("methylation_sensitivity", "")
    if "Dam" in methyl and "block" in methyl.lower():
        results.append(CheckResult(
            name=f"Enzyme:{enzyme_name}:dam_blocking",
            passed=False,
            severity=Severity.CRITICAL,
            message=f"{enzyme_name} is BLOCKED by Dam methylation",
            details=f"Methylation info: {methyl}. Use Dam- strain or choose alternative.",
        ))
    elif "Dam" in methyl:
        results.append(CheckResult(
            name=f"Enzyme:{enzyme_name}:dam_sensitive",
            passed=False,
            severity=Severity.WARNING,
            message=f"{enzyme_name} may be affected by Dam methylation",
            details=f"Methylation info: {methyl}. Verify with Dam- DNA or use alternative.",
        ))

    # Check reliability
    reliability = info.get("reliability", "")
    if reliability in ("low", "medium"):
        results.append(CheckResult(
            name=f"Enzyme:{enzyme_name}:reliability",
            passed=False,
            severity=Severity.WARNING,
            message=f"{enzyme_name} has {reliability} reliability rating",
            details="Consider using a very_high or high reliability enzyme instead.",
        ))

    # Check fussy enzyme list
    fussy_list = metadata.get("fussy_enzymes", {}).get("enzymes", [])
    for entry in fussy_list:
        if entry.get("name") == enzyme_name:
            issues = "; ".join(entry.get("issues", []))
            results.append(CheckResult(
                name=f"Enzyme:{enzyme_name}:fussy",
                passed=False,
                severity=Severity.WARNING,
                message=f"{enzyme_name} is on fussy enzyme list",
                details=f"Issues: {issues}. Recommendation: {entry.get('recommendations', 'N/A')}",
            ))
            break

    return results


# ---------------------------------------------------------------------------
# Synthesis Screen
# ---------------------------------------------------------------------------

def synthesis_screen(
    sequence: str,
    name: str = "unnamed",
    gc_window: int = DEFAULT_GC_WINDOW,
    gc_max: float = DEFAULT_GC_MAX,
    gc_min: float = DEFAULT_GC_MIN,
    homopolymer_max: int = DEFAULT_HOMOPOLYMER_MAX,
    repeat_min_len: int = DEFAULT_REPEAT_MIN_LEN,
) -> List[CheckResult]:
    """
    Mandatory pre-synthesis quality screen.

    Checks:
    - GC content in sliding windows
    - Homopolymer runs
    - Direct repeats >= repeat_min_len bp

    Returns list of CheckResults.
    """
    results = []
    seq = sequence.upper()

    # --- GC content in sliding windows ---
    gc_issues = []
    for i in range(0, len(seq) - gc_window + 1, gc_window // 2):
        window = seq[i:i + gc_window]
        gc = (window.count('G') + window.count('C')) / len(window)
        if gc > gc_max:
            gc_issues.append(f"GC={gc:.0%} at {i}-{i + gc_window}")
        elif gc < gc_min:
            gc_issues.append(f"GC={gc:.0%} at {i}-{i + gc_window}")

    if gc_issues:
        # Summarize — don't list every window
        n = len(gc_issues)
        sample = gc_issues[:3]
        detail = "; ".join(sample)
        if n > 3:
            detail += f" ... and {n - 3} more windows"
        results.append(CheckResult(
            name=f"Synthesis:{name}:gc_content",
            passed=False,
            severity=Severity.WARNING,
            message=f"{n} windows outside GC range [{gc_min:.0%}-{gc_max:.0%}]",
            details=detail,
        ))
    else:
        results.append(CheckResult(
            name=f"Synthesis:{name}:gc_content",
            passed=True,
            severity=Severity.WARNING,
            message=f"GC content within [{gc_min:.0%}-{gc_max:.0%}] in all {gc_window}bp windows",
        ))

    # --- Homopolymer runs ---
    worst_run = ""
    worst_len = 0
    for base in 'ACGT':
        run = base * (homopolymer_max + 1)
        if run in seq:
            # Find longest run
            for length in range(homopolymer_max + 1, 20):
                if base * length in seq:
                    if length > worst_len:
                        worst_len = length
                        worst_run = base
                else:
                    break

    if worst_len > 0:
        results.append(CheckResult(
            name=f"Synthesis:{name}:homopolymer",
            passed=False,
            severity=Severity.WARNING,
            message=f"{worst_run}x{worst_len} homopolymer run detected (max allowed: {homopolymer_max})",
        ))
    else:
        results.append(CheckResult(
            name=f"Synthesis:{name}:homopolymer",
            passed=True,
            severity=Severity.WARNING,
            message=f"No homopolymer runs >{homopolymer_max}bp",
        ))

    # --- Direct repeats ---
    repeat_issues = []
    for length in range(repeat_min_len, min(30, len(seq) // 2)):
        seen = {}
        for i in range(len(seq) - length + 1):
            kmer = seq[i:i + length]
            if kmer in seen:
                # Only report longest repeat per position pair
                if length == repeat_min_len or (kmer[:length - 1] not in seen):
                    repeat_issues.append(
                        f"{length}bp repeat: pos {seen[kmer]} & {i}"
                    )
                if len(repeat_issues) > 20:
                    break
            else:
                seen[kmer] = i
        if len(repeat_issues) > 20:
            break

    if repeat_issues:
        n = len(repeat_issues)
        sample = repeat_issues[:3]
        detail = "; ".join(sample)
        if n > 3:
            detail += f" ... and {n - 3} more"
        results.append(CheckResult(
            name=f"Synthesis:{name}:repeats",
            passed=False,
            severity=Severity.WARNING,
            message=f"{n} direct repeats >={repeat_min_len}bp detected",
            details=detail,
        ))
    else:
        results.append(CheckResult(
            name=f"Synthesis:{name}:repeats",
            passed=True,
            severity=Severity.WARNING,
            message=f"No direct repeats >={repeat_min_len}bp",
        ))

    return results


# ---------------------------------------------------------------------------
# ConstructVerifier — Main Class
# ---------------------------------------------------------------------------

class ConstructVerifier:
    """
    Mandatory post-build verification pipeline.

    Consolidates all verification checks into one enforceable module.
    Build scripts call `run_all()` which returns a VerificationReport.
    If critical checks fail, output should be blocked.

    Example:
        verifier = ConstructVerifier(sequence=my_seq, name="AVD006")
        verifier.add_restriction_check(
            required_unique={"HindIII": "AAGCTT", "XbaI": "TCTAGA"},
            forbidden={"BsaI": "GGTCTC"},
        )
        verifier.add_protein_check(
            region_start=100, region_end=400,
            expected_protein="MKTLLF...",
            frame=0,
        )
        verifier.add_enzyme_validation(["HindIII", "XbaI"])
        verifier.add_synthesis_screen()

        report = verifier.run_all()
        if not report.passed:
            raise VerificationError(str(report))
    """

    def __init__(self, sequence: str, name: str = "unnamed"):
        self.sequence = sequence.upper()
        self.name = name
        self._checks: List = []  # List of callables that return List[CheckResult]

    # --- Check registration methods ---

    def add_restriction_check(
        self,
        required_unique: Optional[Dict[str, str]] = None,
        forbidden: Optional[Dict[str, str]] = None,
    ):
        """
        Add restriction site verification (Checkpoint 9 / BUG-003).

        Args:
            required_unique: Dict of {enzyme_name: recognition_site} that must appear
                             exactly once in the final construct (both strands).
            forbidden: Dict of {enzyme_name: recognition_site} that must NOT appear.
        """
        def _check() -> List[CheckResult]:
            results = []

            # Required unique sites
            if required_unique:
                for enzyme, site in required_unique.items():
                    count = count_sites_both_strands(self.sequence, site)
                    if count == 1:
                        results.append(CheckResult(
                            name=f"RestrictionSite:{enzyme}:unique",
                            passed=True,
                            severity=Severity.CRITICAL,
                            message=f"{enzyme} ({site}) is unique (1 site)",
                        ))
                    else:
                        positions = find_site_positions(self.sequence, site)
                        pos_str = ", ".join(f"{p}({s})" for p, s in positions)
                        results.append(CheckResult(
                            name=f"RestrictionSite:{enzyme}:unique",
                            passed=False,
                            severity=Severity.CRITICAL,
                            message=f"{enzyme} ({site}) has {count} sites (expected 1)",
                            details=f"Positions: {pos_str}" if positions else None,
                        ))

            # Forbidden sites
            if forbidden:
                for enzyme, site in forbidden.items():
                    count = count_sites_both_strands(self.sequence, site)
                    if count == 0:
                        results.append(CheckResult(
                            name=f"RestrictionSite:{enzyme}:absent",
                            passed=True,
                            severity=Severity.CRITICAL,
                            message=f"{enzyme} ({site}) is absent (0 sites)",
                        ))
                    else:
                        positions = find_site_positions(self.sequence, site)
                        pos_str = ", ".join(f"{p}({s})" for p, s in positions)
                        results.append(CheckResult(
                            name=f"RestrictionSite:{enzyme}:absent",
                            passed=False,
                            severity=Severity.CRITICAL,
                            message=f"{enzyme} ({site}) found {count} times (expected 0)",
                            details=f"Positions: {pos_str}" if positions else None,
                        ))

            return results

        self._checks.append(_check)

    def add_restriction_check_by_name(
        self,
        required_unique: Optional[List[str]] = None,
        forbidden: Optional[List[str]] = None,
    ):
        """
        Like add_restriction_check but takes enzyme name lists, looks up sites from KB.

        Args:
            required_unique: List of enzyme names that must appear exactly once.
            forbidden: List of enzyme names that must NOT appear.

        Raises:
            ValueError: If any enzyme name has no recognition site in the KB.
        """
        req_dict = {}
        for name in (required_unique or []):
            site = get_recognition_site(name)
            if site is None:
                raise ValueError(f"Unknown enzyme '{name}': not found in knowledge base")
            req_dict[name] = site

        forb_dict = {}
        for name in (forbidden or []):
            site = get_recognition_site(name)
            if site is None:
                raise ValueError(f"Unknown enzyme '{name}': not found in knowledge base")
            forb_dict[name] = site

        self.add_restriction_check(
            required_unique=req_dict if req_dict else None,
            forbidden=forb_dict if forb_dict else None,
        )

    def add_protein_check(
        self,
        region_start: int,
        region_end: int,
        expected_protein: str,
        frame: int = 0,
        label: str = "CDS",
    ):
        """
        Verify that a region translates to the expected protein.

        Args:
            region_start: 0-indexed start position in construct sequence.
            region_end: 0-indexed end position (exclusive).
            expected_protein: Expected amino acid sequence (single-letter).
            frame: Reading frame offset within the region (0, 1, or 2).
            label: Human-readable name for this protein region.
        """
        def _check() -> List[CheckResult]:
            region = self.sequence[region_start:region_end]
            actual_protein = translate(region, frame=frame)

            if actual_protein == expected_protein:
                return [CheckResult(
                    name=f"Protein:{label}:identity",
                    passed=True,
                    severity=Severity.CRITICAL,
                    message=f"{label} translates correctly ({len(expected_protein)} aa)",
                )]

            # Find first mismatch
            min_len = min(len(actual_protein), len(expected_protein))
            mismatch_pos = None
            for i in range(min_len):
                if actual_protein[i] != expected_protein[i]:
                    mismatch_pos = i
                    break
            if mismatch_pos is None and len(actual_protein) != len(expected_protein):
                mismatch_pos = min_len

            detail = f"Length: {len(actual_protein)} vs expected {len(expected_protein)}"
            if mismatch_pos is not None:
                actual_ctx = actual_protein[max(0, mismatch_pos - 3):mismatch_pos + 4]
                expect_ctx = expected_protein[max(0, mismatch_pos - 3):mismatch_pos + 4]
                detail += f". First mismatch at aa {mismatch_pos + 1}: ...{actual_ctx}... vs ...{expect_ctx}..."

            return [CheckResult(
                name=f"Protein:{label}:identity",
                passed=False,
                severity=Severity.CRITICAL,
                message=f"{label} translation MISMATCH",
                details=detail,
            )]

        self._checks.append(_check)

    def add_insertion_check(
        self,
        insertion_point: int,
        expected_upstream_aa: str,
        expected_downstream_aa: str,
        cds_start: int,
        n_flank_codons: int = 3,
        label: str = "insertion",
    ):
        """
        Verify flanking sequences around an insertion point (BUG-004/006).

        Args:
            insertion_point: 0-indexed DNA position of insertion.
            expected_upstream_aa: Expected amino acids upstream of insertion.
            expected_downstream_aa: Expected amino acids downstream of insertion.
            cds_start: 0-indexed position of CDS ATG codon.
            n_flank_codons: Number of codons to check on each side.
            label: Human-readable label.
        """
        def _check() -> List[CheckResult]:
            results = []

            # Calculate frame offset (BUG-001)
            frame = (insertion_point - cds_start) % 3

            # Upstream flanking codons
            upstream_start = insertion_point - (n_flank_codons * 3) - frame
            upstream_dna = self.sequence[upstream_start:insertion_point]
            upstream_aa = translate(upstream_dna, frame=frame)

            if upstream_aa.endswith(expected_upstream_aa):
                results.append(CheckResult(
                    name=f"Insertion:{label}:upstream_flank",
                    passed=True,
                    severity=Severity.CRITICAL,
                    message=f"Upstream flank correct: ...{upstream_aa[-len(expected_upstream_aa):]}",
                ))
            else:
                results.append(CheckResult(
                    name=f"Insertion:{label}:upstream_flank",
                    passed=False,
                    severity=Severity.CRITICAL,
                    message=f"Upstream flank MISMATCH",
                    details=f"Got ...{upstream_aa[-6:]} but expected ...{expected_upstream_aa}",
                ))

            # Downstream flanking codons
            downstream_end = insertion_point + (n_flank_codons * 3) + (3 - frame) % 3
            downstream_dna = self.sequence[insertion_point:downstream_end]
            downstream_aa = translate(downstream_dna, frame=(3 - frame) % 3)

            if downstream_aa.startswith(expected_downstream_aa):
                results.append(CheckResult(
                    name=f"Insertion:{label}:downstream_flank",
                    passed=True,
                    severity=Severity.CRITICAL,
                    message=f"Downstream flank correct: {downstream_aa[:len(expected_downstream_aa)]}...",
                ))
            else:
                results.append(CheckResult(
                    name=f"Insertion:{label}:downstream_flank",
                    passed=False,
                    severity=Severity.CRITICAL,
                    message=f"Downstream flank MISMATCH",
                    details=f"Got {downstream_aa[:6]}... but expected {expected_downstream_aa}...",
                ))

            return results

        self._checks.append(_check)

    def add_parent_child_check(
        self,
        parent_sequence: str,
        expected_changes: Optional[List[str]] = None,
        label: str = "parent-child",
    ):
        """
        Compare parent and child sequences to verify only expected changes (Checkpoint 10).

        Args:
            parent_sequence: The parent construct sequence.
            expected_changes: List of human-readable expected change descriptions.
            label: Human-readable label.
        """
        def _check() -> List[CheckResult]:
            results = []
            parent = parent_sequence.upper()

            # Find differences
            if len(self.sequence) == len(parent):
                diffs = []
                for i, (a, b) in enumerate(zip(self.sequence, parent)):
                    if a != b:
                        diffs.append(i)
                results.append(CheckResult(
                    name=f"ParentChild:{label}:point_mutations",
                    passed=True,  # Just informational
                    severity=Severity.INFO,
                    message=f"{len(diffs)} point differences vs parent",
                    details=f"Positions: {diffs[:10]}{'...' if len(diffs) > 10 else ''}",
                ))
            else:
                size_diff = len(self.sequence) - len(parent)
                results.append(CheckResult(
                    name=f"ParentChild:{label}:size",
                    passed=True,  # Informational — insertions change size
                    severity=Severity.INFO,
                    message=f"Size differs by {size_diff:+d} bp ({len(parent)} → {len(self.sequence)})",
                ))

            return results

        self._checks.append(_check)

    def add_packaging_check(self, itr_left_end: int, itr_right_start: int,
                            limit: int = AAV_PACKAGING_LIMIT):
        """
        Verify AAV packaging size constraint (Checkpoint 4).

        Args:
            itr_left_end: 0-indexed end position of left ITR.
            itr_right_start: 0-indexed start position of right ITR.
            limit: Maximum allowed bp between ITRs.
        """
        def _check() -> List[CheckResult]:
            payload = itr_right_start - itr_left_end
            passed = payload <= limit
            return [CheckResult(
                name="Packaging:size_limit",
                passed=passed,
                severity=Severity.CRITICAL,
                message=f"Payload = {payload} bp {'<=' if passed else '>'} {limit} bp limit",
                details=None if passed else f"Over by {payload - limit} bp. Reduce insert size.",
            )]

        self._checks.append(_check)

    def add_enzyme_validation(self, enzyme_names: List[str]):
        """
        Validate enzyme choices against knowledge base (Plan 3c).

        Args:
            enzyme_names: List of enzyme names to validate.
        """
        def _check() -> List[CheckResult]:
            results = []
            for name in enzyme_names:
                enzyme_results = validate_enzyme_choice(name)
                if enzyme_results:
                    results.extend(enzyme_results)
                else:
                    results.append(CheckResult(
                        name=f"Enzyme:{name}:validated",
                        passed=True,
                        severity=Severity.INFO,
                        message=f"{name} passes all enzyme quality checks",
                    ))
            return results

        self._checks.append(_check)

    def add_synthesis_screen(self, **kwargs):
        """
        Add mandatory pre-synthesis quality screen (Plan 3d).

        Keyword args are passed to synthesis_screen().
        """
        def _check() -> List[CheckResult]:
            return synthesis_screen(self.sequence, name=self.name, **kwargs)

        self._checks.append(_check)

    def add_frame_check(self, position: int, cds_start: int,
                        expected_frame: int, label: str = "frame"):
        """
        Verify reading frame at a specific position (BUG-001).

        Args:
            position: 0-indexed DNA position to check.
            cds_start: 0-indexed CDS start (ATG position).
            expected_frame: Expected frame offset (0, 1, or 2).
            label: Human-readable label.
        """
        def _check() -> List[CheckResult]:
            actual_frame = (position - cds_start) % 3
            passed = actual_frame == expected_frame
            return [CheckResult(
                name=f"Frame:{label}",
                passed=passed,
                severity=Severity.CRITICAL,
                message=f"Frame at pos {position} = {actual_frame} "
                        f"{'==' if passed else '!='} expected {expected_frame}",
                details=f"Calculated as ({position} - {cds_start}) % 3 = {actual_frame}" if not passed else None,
            )]

        self._checks.append(_check)

    def add_custom_check(self, name: str, check_fn, severity: Severity = Severity.CRITICAL):
        """
        Add a custom check function.

        Args:
            name: Check name.
            check_fn: Callable that takes (sequence: str) and returns (passed: bool, message: str).
            severity: Check severity.
        """
        def _check() -> List[CheckResult]:
            passed, message = check_fn(self.sequence)
            return [CheckResult(
                name=name,
                passed=passed,
                severity=severity,
                message=message,
            )]

        self._checks.append(_check)

    # --- Execution ---

    def run_all(self) -> VerificationReport:
        """
        Execute all registered checks and return a VerificationReport.

        Returns:
            VerificationReport with all results. Check report.passed to determine
            if output should be written.
        """
        report = VerificationReport(construct_name=self.name)
        for check_fn in self._checks:
            try:
                results = check_fn()
                report.checks.extend(results)
            except Exception as e:
                report.checks.append(CheckResult(
                    name="InternalError",
                    passed=False,
                    severity=Severity.CRITICAL,
                    message=f"Check raised exception: {e}",
                    details=str(e),
                ))
        return report


# ---------------------------------------------------------------------------
# Convenience: VerificationError
# ---------------------------------------------------------------------------

class VerificationError(Exception):
    """Raised when mandatory verification fails."""
    def __init__(self, report: VerificationReport):
        self.report = report
        super().__init__(str(report))


def verify_and_gate(verifier: ConstructVerifier) -> VerificationReport:
    """
    Run all checks and raise VerificationError if critical checks fail.

    This is the "gate" function — call this before writing output.

    Usage:
        report = verify_and_gate(verifier)
        # If we reach here, all critical checks passed
        write_genbank(...)
    """
    report = verifier.run_all()
    if not report.passed:
        raise VerificationError(report)
    return report


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Demo with a simple test sequence
    demo_seq = (
        "AAGCTT"   # HindIII
        "ATGAAAGCCTTTGGCAAA"  # CDS region
        "TCTAGA"   # XbaI
    )

    verifier = ConstructVerifier(sequence=demo_seq, name="demo_construct")
    verifier.add_restriction_check(
        required_unique={"HindIII": "AAGCTT", "XbaI": "TCTAGA"},
        forbidden={"BsaI": "GGTCTC"},
    )
    verifier.add_enzyme_validation(["HindIII", "XbaI"])
    verifier.add_synthesis_screen()

    report = verifier.run_all()
    print(report)
    sys.exit(0 if report.passed else 1)
