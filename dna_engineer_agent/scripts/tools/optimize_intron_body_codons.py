#!/usr/bin/env python3
"""
Full Intron Body Codon Optimization with DNAchisel + MaxEntScan

Three-stage optimization of the AVD059-099 intron body (447 bp = GGGGS5 + VHH3 + GGGGS1):

Stage 0: Design GT/AG-free GGGGS5 and GGGGS1 linkers with varied codons.
Stage 1: DNAchisel global optimization of VHH3 to minimize GT/AG while preserving protein.
Stage 2: MaxEntScan context refinement of remaining GT donors and AG acceptors.

All changes are synonymous — the full body protein sequence is preserved exactly.
All 149 codons are mutable in Stage 2 (GGGGS + VHH3).

Usage:
    python scripts/tools/optimize_intron_body_codons.py
"""

import sys
from itertools import product
from pathlib import Path
from typing import Dict, List, Tuple

# Add tools directory to path for imports
TOOLS_DIR = Path(__file__).parent
sys.path.insert(0, str(TOOLS_DIR))

from build_splice_cassettes import (
    AA_TO_CODONS,
    CODON_TABLE,
    score_donor_maxentscan,
    score_acceptor_maxentscan,
    classify_maxent_score,
    translate,
)
from splice_donor_hbonds import calculate_hbonds
from u1_duplex_dg import calculate_u1_duplex_dg

# ============================================================================
# Constants — AVD006 canonical input sequences (before optimization)
# ============================================================================

GGGGS5_ORIGINAL = "GGTGGAGGCGGATCTGGAGGCGGTGGTTCAGGCGGTGGAGGAAGTGGTGGCGGAGGTTCTGGTGGAGGCGGTTCT"

VHH3_ORIGINAL = (
    "GAGGTGCAACTGGTTGAAAGCGGCGGAGGACTTGTTCAACCCGGCGGCAGCCTTAGGCTTTCTTGCGCTGCCAGC"
    "GGCTTCACCTTTAGCACCGCCGACATGGGCTGGTTTAGGCAAGCTCCCGGAAAAGGCAGGGAACTTGTTGCCGCT"
    "GTGAGCGGCAGCGGCTTCAGCACCTACTCTGATAGCGTTGAGGGCAGGTTCACCATCAGCAGGGACAACGCCAAG"
    "AGGATGGTGTACCTGCAGATGAACAGCTTGAGGGCCGAGGACACCGCCGTGTACTACTGCGCCAAGGCCACAATT"
    "AGCCTGTACTACGCCATGGATGTGTGGGGACAGGGCACCACCGTGACCGTGAGCAGC"
)

GGGGS1_ORIGINAL = "GGTGGAGGCGGTTCT"

VHH3_PROTEIN = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFSTADMGWFRQAPGKGR"
    "ELVAAVSGSGFSTYSDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAKATISLYY"
    "AMDVWGQGTTVTVSS"
)

# Human codon usage frequencies (per thousand) for tiebreaking
# Source: Kazusa codon usage database, Homo sapiens
HUMAN_CODON_FREQ = {
    "TTT": 17.6, "TTC": 20.3, "TTA": 7.7, "TTG": 12.9,
    "TCT": 15.2, "TCC": 17.7, "TCA": 12.2, "TCG": 4.4,
    "TAT": 12.2, "TAC": 15.3, "TAA": 1.0, "TAG": 0.8,
    "TGT": 10.6, "TGC": 12.6, "TGA": 1.6, "TGG": 13.2,
    "CTT": 13.2, "CTC": 19.6, "CTA": 7.2, "CTG": 39.6,
    "CCT": 17.5, "CCC": 19.8, "CCA": 16.9, "CCG": 6.9,
    "CAT": 10.9, "CAC": 15.1, "CAA": 12.3, "CAG": 34.2,
    "CGT": 4.5, "CGC": 10.4, "CGA": 6.2, "CGG": 11.4,
    "ATT": 16.0, "ATC": 20.8, "ATA": 7.5, "ATG": 22.0,
    "ACT": 13.1, "ACC": 18.9, "ACA": 15.1, "ACG": 6.1,
    "AAT": 17.0, "AAC": 19.1, "AAA": 24.4, "AAG": 31.9,
    "AGT": 12.1, "AGC": 19.5, "AGA": 12.2, "AGG": 12.0,
    "GTT": 11.0, "GTC": 14.5, "GTA": 7.1, "GTG": 28.1,
    "GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4,
    "GAT": 21.8, "GAC": 25.1, "GAA": 29.0, "GAG": 39.6,
    "GGT": 10.8, "GGC": 22.2, "GGA": 16.5, "GGG": 16.5,
}


# Restriction sites that must not be created (Golden Gate compatibility)
FORBIDDEN_RESTRICTION_SITES = {
    "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC",
    "BbsI": "GAAGAC",
}


# ============================================================================
# Utility functions
# ============================================================================

def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def count_restriction_sites(seq: str) -> int:
    """Count total forbidden restriction sites on both strands."""
    total = 0
    for name, site in FORBIDDEN_RESTRICTION_SITES.items():
        rc = reverse_complement(site)
        total += seq.count(site)
        if site != rc:
            total += seq.count(rc)
    return total


def has_new_restriction_sites(original_seq: str, modified_seq: str) -> bool:
    """Check if modification introduced new restriction sites."""
    return count_restriction_sites(modified_seq) > count_restriction_sites(original_seq)


def split_codons(dna: str) -> List[str]:
    """Split DNA into codons (triplets)."""
    return [dna[i:i+3] for i in range(0, len(dna), 3)]


def join_codons(codons: List[str]) -> str:
    """Join codons into a DNA string."""
    return "".join(codons)


def find_dinuc_positions(seq: str, dinuc: str) -> List[int]:
    """Find all 0-indexed positions of a dinucleotide."""
    positions = []
    for i in range(len(seq) - 1):
        if seq[i:i+2] == dinuc:
            positions.append(i)
    return positions


def has_gt_or_ag_internal(codon: str) -> bool:
    """Check if a single codon contains GT or AG internally."""
    return "GT" in codon or "AG" in codon


def has_forbidden_dinuc_at_boundary(codons: List[str], idx: int) -> bool:
    """
    Check if the codon at idx creates GT or AG at its boundaries with neighbors.

    Checks:
    - Internal GT/AG within codon[idx]
    - Junction with left neighbor: codons[idx-1][2] + codons[idx][0]
    - Junction with right neighbor: codons[idx][2] + codons[idx+1][0]
    """
    codon = codons[idx]

    # Internal
    if has_gt_or_ag_internal(codon):
        return True

    # Left boundary
    if idx > 0:
        junction_left = codons[idx - 1][2] + codon[0]
        if junction_left in ("GT", "AG"):
            return True

    # Right boundary
    if idx < len(codons) - 1:
        junction_right = codon[2] + codons[idx + 1][0]
        if junction_right in ("GT", "AG"):
            return True

    return False


def count_dinuc(seq: str, dinuc: str) -> int:
    """Count occurrences of a dinucleotide in a sequence."""
    return len(find_dinuc_positions(seq, dinuc))


# ============================================================================
# Scoring functions (extended for AG acceptors)
# ============================================================================

def score_gt_at_pos(seq: str, gt_pos: int) -> Dict:
    """
    Score a GT dinucleotide at position gt_pos in seq using MaxEntScan and U1 dG.

    Returns dict with maxent_score and u1_dg (or None if window out of bounds).
    """
    result = {"maxent_score": None, "u1_dg": None}

    # 9-mer for MaxEntScan: pos-3 to pos+6
    s9 = gt_pos - 3
    e9 = gt_pos + 6
    if 0 <= s9 and e9 <= len(seq):
        niner = seq[s9:e9]
        result["maxent_score"] = score_donor_maxentscan(niner)

    # 11-mer for U1 dG: pos-3 to pos+8
    s11 = gt_pos - 3
    e11 = gt_pos + 8
    if 0 <= s11 and e11 <= len(seq):
        eleven_mer = seq[s11:e11]
        dg_result = calculate_u1_duplex_dg(eleven_mer)
        result["u1_dg"] = dg_result.delta_g

    return result


def score_ag_at_pos(seq: str, ag_pos: int) -> Dict:
    """
    Score an AG dinucleotide at position ag_pos using MaxEntScan acceptor scoring.

    23-mer window: seq[ag_pos-18 : ag_pos+5] where AG is at indices 18-19.
    Valid when ag_pos >= 18 and ag_pos + 5 <= len(seq).
    """
    result = {"acceptor_score": None}
    start = ag_pos - 18
    end = ag_pos + 5
    if 0 <= start and end <= len(seq):
        mer23 = seq[start:end]
        result["acceptor_score"] = score_acceptor_maxentscan(mer23)
    return result


def score_all_gt_ag_in_region(seq: str, start: int, end: int) -> List[Dict]:
    """
    Score all GT and AG dinucleotides in seq[start:end].

    GT sites get MaxEntScan donor scores and U1 dG.
    AG sites get MaxEntScan acceptor scores.
    """
    results = []
    for i in range(max(0, start), min(end, len(seq) - 1)):
        dinuc = seq[i:i+2]
        if dinuc == "GT":
            scores = score_gt_at_pos(seq, i)
            results.append({
                "pos": i,
                "dinuc": "GT",
                "maxent_score": scores["maxent_score"],
                "u1_dg": scores["u1_dg"],
            })
        elif dinuc == "AG":
            ag_scores = score_ag_at_pos(seq, i)
            results.append({
                "pos": i,
                "dinuc": "AG",
                "maxent_score": ag_scores["acceptor_score"],
                "u1_dg": None,
            })
    return results


def do_no_harm_check(
    original_seq: str,
    modified_seq: str,
    change_center: int,
    scan_radius: int = 25,
    safe_ceiling: float = 0.0,
) -> bool:
    """
    Safety check: ensure no GT or AG site in the region around a change scores
    higher (stronger splice site) than it did before the change.

    Also rejects if:
    - New GT or AG dinucleotides appear that weren't there before
    - New forbidden restriction sites (BsaI, BsmBI, BbsI) are introduced

    The safe_ceiling parameter (default 0.0) relaxes the regression check:
    if a site's MaxEntScan score increases but remains below safe_ceiling,
    the change is allowed. Sites below 0.0 are "Not Detected" by the
    spliceosome and are functionally irrelevant.

    Returns True if the change is safe, False if it violates "do no harm".
    """
    # Check for new restriction sites (full sequence check)
    if has_new_restriction_sites(original_seq, modified_seq):
        return False

    start = max(0, change_center - scan_radius)
    end = min(len(original_seq), change_center + scan_radius)

    orig_sites = score_all_gt_ag_in_region(original_seq, start, end)
    mod_sites = score_all_gt_ag_in_region(modified_seq, start, end)

    # Build lookup for original scores
    orig_by_pos = {}
    for site in orig_sites:
        orig_by_pos[(site["pos"], site["dinuc"])] = site

    # Check each site in modified sequence
    for site in mod_sites:
        key = (site["pos"], site["dinuc"])
        if key not in orig_by_pos:
            # New GT or AG appeared — reject
            return False
        orig = orig_by_pos[key]
        # For both GT and AG: check MaxEntScan score didn't increase
        if (site["maxent_score"] is not None and orig["maxent_score"] is not None
                and site["maxent_score"] > orig["maxent_score"] + 0.01):
            # Allow regression if the site stays below the safe ceiling (Not Detected)
            if site["maxent_score"] >= safe_ceiling:
                return False
        # For GT sites only: also check U1 dG (more negative = stronger)
        if site["dinuc"] == "GT":
            if (site["u1_dg"] is not None and orig["u1_dg"] is not None
                    and site["u1_dg"] < orig["u1_dg"] - 0.01):
                return False

    return True


# ============================================================================
# Stage 0: Design GT/AG-free GGGGS linkers
# ============================================================================

def design_ggggs_zero_gtag() -> Tuple[str, str]:
    """
    Design GGGGS5 (75 bp) and GGGGS1 (15 bp) with zero GT, zero AG,
    and varied codons to avoid synthesis issues (BUG-007).

    Codon rules for zero GT/AG in each GGGGS unit (G-G-G-G-S):
      G1-G3 (Gly before Gly): {GGC, GGG}
        - GGT has internal GT
        - GGA creates junction AG with next G-codon
      G4 (Gly before Ser):    {GGC, GGA}
        - GGT has internal GT
        - GGG creates junction GT with T-starting Ser codons
      S (Ser, mid-chain):     {TCT, TCC, TCG}
        - TCA creates junction AG with next G-codon
        - AGT has internal GT + AG
        - AGC has internal AG

    48 valid encodings per unit. 5 different units chosen for codon diversity.

    Returns (ggggs5_dna, ggggs1_dna).
    """
    # 5 varied GGGGS units for GGGGS5 (each unique for codon diversity)
    units = [
        #           G1     G2     G3     G4     S
        ["GGC", "GGG", "GGC", "GGA", "TCT"],  # Unit 1
        ["GGG", "GGC", "GGG", "GGC", "TCC"],  # Unit 2
        ["GGC", "GGC", "GGG", "GGA", "TCG"],  # Unit 3
        ["GGG", "GGG", "GGC", "GGC", "TCT"],  # Unit 4
        ["GGC", "GGG", "GGG", "GGA", "TCC"],  # Unit 5
    ]

    # GGGGS1 (different from all 5 units above)
    ggggs1_unit = ["GGG", "GGC", "GGC", "GGA", "TCG"]

    # Build DNA strings
    ggggs5 = "".join(codon for unit in units for codon in unit)
    ggggs1 = "".join(ggggs1_unit)

    # --- Verification ---
    assert len(ggggs5) == 75, f"GGGGS5 length {len(ggggs5)} != 75"
    assert len(ggggs1) == 15, f"GGGGS1 length {len(ggggs1)} != 15"

    assert count_dinuc(ggggs5, "GT") == 0, f"GGGGS5 has {count_dinuc(ggggs5, 'GT')} GT!"
    assert count_dinuc(ggggs5, "AG") == 0, f"GGGGS5 has {count_dinuc(ggggs5, 'AG')} AG!"
    assert count_dinuc(ggggs1, "GT") == 0, f"GGGGS1 has {count_dinuc(ggggs1, 'GT')} GT!"
    assert count_dinuc(ggggs1, "AG") == 0, f"GGGGS1 has {count_dinuc(ggggs1, 'AG')} AG!"

    assert translate(ggggs5) == "GGGGS" * 5, f"GGGGS5 protein mismatch!"
    assert translate(ggggs1) == "GGGGS", f"GGGGS1 protein mismatch!"

    # Verify no unit-to-unit junctions create GT/AG
    for i in range(len(units) - 1):
        last_nt = units[i][-1][-1]   # Last nt of S codon
        first_nt = units[i+1][0][0]  # First nt of next G1 codon
        junction = last_nt + first_nt
        assert junction not in ("GT", "AG"), \
            f"Junction GT/AG between units {i+1} and {i+2}: {junction}"

    # Verify no 15-mer repeats (BUG-007)
    unit_strs = [join_codons(u) for u in units]
    assert len(set(unit_strs)) == 5, "GGGGS5 has repeated units!"
    assert join_codons(ggggs1_unit) not in unit_strs, "GGGGS1 matches a GGGGS5 unit!"

    # Verify no restriction sites
    assert count_restriction_sites(ggggs5) == 0, "GGGGS5 has restriction sites!"
    assert count_restriction_sites(ggggs1) == 0, "GGGGS1 has restriction sites!"

    return ggggs5, ggggs1


# ============================================================================
# Stage 1: DNAchisel optimization of VHH3
# ============================================================================

def optimize_vhh3_dnachisel(vhh3_dna: str) -> Tuple[str, Dict]:
    """
    Use DNAchisel to minimize GT/AG in VHH3 while preserving the protein.

    Constraints:
    - EnforceTranslation (protein must be identical)
    - No BsaI/BsmBI/BbsI restriction sites (both strands)

    Objectives (weighted):
    - Minimize GT dinucleotides (boost=5.0)
    - Minimize AG dinucleotides (boost=5.0)
    - Match human codon usage (boost=0.5, tiebreaker)

    Returns (optimized_vhh3_dna, stats_dict).
    """
    from dnachisel import (
        DnaOptimizationProblem, AvoidPattern, EnforceTranslation,
        CodonOptimize,
    )

    problem = DnaOptimizationProblem(
        sequence=vhh3_dna,
        constraints=[
            EnforceTranslation(),
            AvoidPattern("GGTCTC"),    # BsaI
            AvoidPattern("GAGACC"),    # BsaI rc
            AvoidPattern("CGTCTC"),    # BsmBI
            AvoidPattern("GAGACG"),    # BsmBI rc
            AvoidPattern("GAAGAC"),    # BbsI
            AvoidPattern("GTCTTC"),    # BbsI rc
        ],
        objectives=[
            AvoidPattern("GT", boost=5.0),
            AvoidPattern("AG", boost=5.0),
            CodonOptimize(species='h_sapiens', method='match_codon_usage', boost=0.5),
        ],
    )

    # Seed for reproducibility (DNAchisel uses numpy random internally)
    import numpy as np
    np.random.seed(42)

    problem.resolve_constraints()
    problem.optimize()

    optimized = problem.sequence

    # Verify protein preserved
    assert translate(optimized) == translate(vhh3_dna), "DNAchisel broke VHH3 protein!"

    stats = {
        "before_gt": count_dinuc(vhh3_dna, "GT"),
        "before_ag": count_dinuc(vhh3_dna, "AG"),
        "after_gt": count_dinuc(optimized, "GT"),
        "after_ag": count_dinuc(optimized, "AG"),
    }

    return optimized, stats


def fix_boundary_gtag(
    ggggs5: str, vhh3_codons: List[str], ggggs1: str,
) -> Tuple[List[str], List[Dict]]:
    """
    Fix any GT/AG at GGGGS↔VHH3 boundaries after DNAchisel optimization.

    Checks:
    - GGGGS5 last nt + VHH3 first nt (must not be GT or AG)
    - VHH3 last nt + GGGGS1 first nt (must not be GT or AG)

    Returns (fixed_vhh3_codons, changes).
    """
    vhh3_codons = list(vhh3_codons)
    changes = []

    # Check GGGGS5 → VHH3 junction
    last_nt_ggggs5 = ggggs5[-1]
    junction_left = last_nt_ggggs5 + vhh3_codons[0][0]
    if junction_left in ("GT", "AG"):
        aa = CODON_TABLE[vhh3_codons[0]]
        old = vhh3_codons[0]
        # Pick best synonymous codon that fixes the junction
        candidates = sorted(
            AA_TO_CODONS[aa],
            key=lambda c: HUMAN_CODON_FREQ.get(c, 0),
            reverse=True,
        )
        for alt in candidates:
            new_junction = last_nt_ggggs5 + alt[0]
            if new_junction not in ("GT", "AG") and not has_gt_or_ag_internal(alt):
                vhh3_codons[0] = alt
                if alt != old:
                    changes.append({
                        "codon_idx": 0, "aa": aa,
                        "old_codon": old, "new_codon": alt,
                        "phase": "1b",
                        "reason": f"Fix GGGGS5->VHH3 junction {junction_left}",
                    })
                break

    # Check VHH3 → GGGGS1 junction
    first_nt_ggggs1 = ggggs1[0]
    junction_right = vhh3_codons[-1][2] + first_nt_ggggs1
    if junction_right in ("GT", "AG"):
        aa = CODON_TABLE[vhh3_codons[-1]]
        old = vhh3_codons[-1]
        candidates = sorted(
            AA_TO_CODONS[aa],
            key=lambda c: HUMAN_CODON_FREQ.get(c, 0),
            reverse=True,
        )
        for alt in candidates:
            new_junction = alt[2] + first_nt_ggggs1
            if new_junction not in ("GT", "AG") and not has_gt_or_ag_internal(alt):
                vhh3_codons[-1] = alt
                if alt != old:
                    changes.append({
                        "codon_idx": len(vhh3_codons) - 1, "aa": aa,
                        "old_codon": old, "new_codon": alt,
                        "phase": "1b",
                        "reason": f"Fix VHH3->GGGGS1 junction {junction_right}",
                    })
                break

    return vhh3_codons, changes


# ============================================================================
# Stage 2a: Optimize context at remaining GT donor sites (MaxEntScan)
# ============================================================================

def stage2a_optimize_gt_context(
    codons: List[str],
    mutable_range: Tuple[int, int],
) -> Tuple[List[str], List[Dict]]:
    """
    For each remaining GT site, enumerate synonymous codon combinations
    in the flanking window and pick the combination that minimizes MaxEntScan
    score while passing the do-no-harm safety check.

    All 149 codons are mutable (GGGGS + VHH3).

    Returns (modified codons, list of changes made).
    """
    codons = list(codons)
    changes = []
    mutable_start, mutable_end = mutable_range

    seq = join_codons(codons)
    gt_positions = find_dinuc_positions(seq, "GT")

    for gt_pos in gt_positions:
        seq = join_codons(codons)  # refresh after possible changes

        # Identify which codons contribute to the 9-mer window (pos-3 to pos+5)
        window_start = max(0, gt_pos - 3)
        window_end = min(len(seq), gt_pos + 6)

        # Find codon indices that overlap this window
        codon_indices = set()
        for bp in range(window_start, window_end):
            codon_indices.add(bp // 3)
        codon_indices = sorted(codon_indices)

        # Filter to mutable codons only, and cap at 3 closest to GT
        mutable_in_window = [i for i in codon_indices if mutable_start <= i < mutable_end]
        gt_codon = gt_pos // 3
        mutable_in_window.sort(key=lambda i: abs(i - gt_codon))
        mutable_in_window = mutable_in_window[:3]

        if not mutable_in_window:
            continue

        # Get current scores
        current_scores = score_gt_at_pos(seq, gt_pos)
        current_maxent = current_scores["maxent_score"]
        current_u1_dg = current_scores["u1_dg"]

        if current_maxent is None:
            continue

        # Build list of synonym options for each mutable codon
        codon_options = []
        for ci in mutable_in_window:
            aa = CODON_TABLE.get(codons[ci], "?")
            if aa == "?" or aa == "*":
                codon_options.append([codons[ci]])
            else:
                codon_options.append(AA_TO_CODONS.get(aa, [codons[ci]]))

        best_combo = None
        best_maxent = current_maxent
        best_u1_dg = current_u1_dg if current_u1_dg is not None else 0.0
        best_freq = -1.0

        for combo in product(*codon_options):
            # Apply combination
            test_codons = list(codons)
            for ci, alt in zip(mutable_in_window, combo):
                test_codons[ci] = alt

            test_seq = join_codons(test_codons)

            # Do no harm check (rejects new GT/AG, score regression, restriction sites)
            change_center = mutable_in_window[len(mutable_in_window) // 2] * 3 + 1
            if not do_no_harm_check(seq, test_seq, change_center):
                continue

            # Score this combination
            if gt_pos < len(test_seq) - 1 and test_seq[gt_pos:gt_pos+2] == "GT":
                new_scores = score_gt_at_pos(test_seq, gt_pos)
                new_maxent = new_scores["maxent_score"] if new_scores["maxent_score"] is not None else current_maxent
                new_u1_dg = new_scores["u1_dg"] if new_scores["u1_dg"] is not None else 0.0
            else:
                # GT was eliminated — best possible outcome
                new_maxent = -999.0
                new_u1_dg = 0.0

            # Codon usage frequency (tiebreaker)
            freq = sum(HUMAN_CODON_FREQ.get(alt, 0) for alt in combo)

            # Compare: minimize MaxEntScan (primary), maximize U1 dG (secondary),
            # maximize codon frequency (tiebreaker)
            if (new_maxent < best_maxent - 0.01 or
                (abs(new_maxent - best_maxent) < 0.01 and new_u1_dg > best_u1_dg + 0.01) or
                (abs(new_maxent - best_maxent) < 0.01 and abs(new_u1_dg - best_u1_dg) < 0.01
                 and freq > best_freq)):
                best_combo = combo
                best_maxent = new_maxent
                best_u1_dg = new_u1_dg
                best_freq = freq

        # Apply best combination if it improves on current
        if best_combo is not None and best_maxent < current_maxent - 0.01:
            for ci, alt in zip(mutable_in_window, best_combo):
                if alt != codons[ci]:
                    old_codon = codons[ci]
                    codons[ci] = alt
                    changes.append({
                        "codon_idx": ci,
                        "aa": CODON_TABLE.get(old_codon, "?"),
                        "old_codon": old_codon,
                        "new_codon": alt,
                        "phase": "2a",
                        "reason": f"Optimize GT context at pos {gt_pos} "
                                  f"(MaxEnt {current_maxent:.2f} -> {best_maxent:.2f})",
                    })

    return codons, changes


# ============================================================================
# Stage 2b: Optimize context at remaining AG acceptor sites (MaxEntScan)
# ============================================================================

def stage2b_optimize_ag_context(
    codons: List[str],
    mutable_range: Tuple[int, int],
) -> Tuple[List[str], List[Dict]]:
    """
    For each remaining AG site, enumerate synonymous codon combinations
    in the flanking window and pick the combination that minimizes MaxEntScan
    acceptor score while passing the do-no-harm safety check.

    Returns (modified codons, list of changes made).
    """
    codons = list(codons)
    changes = []
    mutable_start, mutable_end = mutable_range

    seq = join_codons(codons)
    ag_positions = find_dinuc_positions(seq, "AG")

    for ag_pos in ag_positions:
        seq = join_codons(codons)  # refresh after possible changes

        # Get current acceptor score
        current_scores = score_ag_at_pos(seq, ag_pos)
        current_acceptor = current_scores["acceptor_score"]

        if current_acceptor is None:
            continue

        # Identify codons near the AG for optimization
        # Use a window of ±3 bp around the AG dinucleotide
        window_start = max(0, ag_pos - 3)
        window_end = min(len(seq), ag_pos + 5)

        codon_indices = set()
        for bp in range(window_start, window_end):
            codon_indices.add(bp // 3)
        codon_indices = sorted(codon_indices)

        # Filter to mutable codons, cap at 3 closest
        mutable_in_window = [i for i in codon_indices if mutable_start <= i < mutable_end]
        ag_codon = ag_pos // 3
        mutable_in_window.sort(key=lambda i: abs(i - ag_codon))
        mutable_in_window = mutable_in_window[:3]

        if not mutable_in_window:
            continue

        # Build synonym options
        codon_options = []
        for ci in mutable_in_window:
            aa = CODON_TABLE.get(codons[ci], "?")
            if aa in ("?", "*"):
                codon_options.append([codons[ci]])
            else:
                codon_options.append(AA_TO_CODONS.get(aa, [codons[ci]]))

        best_combo = None
        best_acceptor = current_acceptor
        best_freq = -1.0

        for combo in product(*codon_options):
            test_codons = list(codons)
            for ci, alt in zip(mutable_in_window, combo):
                test_codons[ci] = alt

            test_seq = join_codons(test_codons)

            # Do no harm check
            change_center = mutable_in_window[len(mutable_in_window) // 2] * 3 + 1
            if not do_no_harm_check(seq, test_seq, change_center):
                continue

            # Score the AG site
            if ag_pos < len(test_seq) - 1 and test_seq[ag_pos:ag_pos+2] == "AG":
                new_scores = score_ag_at_pos(test_seq, ag_pos)
                new_acceptor = new_scores["acceptor_score"] if new_scores["acceptor_score"] is not None else current_acceptor
            else:
                # AG was eliminated — best possible outcome
                new_acceptor = -999.0

            # Codon frequency tiebreaker
            freq = sum(HUMAN_CODON_FREQ.get(alt, 0) for alt in combo)

            if (new_acceptor < best_acceptor - 0.01 or
                (abs(new_acceptor - best_acceptor) < 0.01 and freq > best_freq)):
                best_combo = combo
                best_acceptor = new_acceptor
                best_freq = freq

        # Apply best combination if it improves on current
        if best_combo is not None and best_acceptor < current_acceptor - 0.01:
            for ci, alt in zip(mutable_in_window, best_combo):
                if alt != codons[ci]:
                    old_codon = codons[ci]
                    codons[ci] = alt
                    changes.append({
                        "codon_idx": ci,
                        "aa": CODON_TABLE.get(old_codon, "?"),
                        "old_codon": old_codon,
                        "new_codon": alt,
                        "phase": "2b",
                        "reason": f"Optimize AG acceptor at pos {ag_pos} "
                                  f"(AccScore {current_acceptor:.2f} -> {best_acceptor:.2f})",
                    })

    return codons, changes


# ============================================================================
# Reporting
# ============================================================================

def print_gt_ag_summary(label: str, seq: str):
    """Print GT/AG counts for a sequence."""
    gt = count_dinuc(seq, "GT")
    ag = count_dinuc(seq, "AG")
    print(f"  {label}: GT={gt}, AG={ag} (total forbidden={gt + ag})")


def print_per_site_scores(seq: str, label: str):
    """Print MaxEntScan and U1 dG for every GT site, and acceptor scores for AG sites."""
    gt_positions = find_dinuc_positions(seq, "GT")
    ag_positions = find_dinuc_positions(seq, "AG")

    # GT donor scores
    print(f"\n  {label} GT donor scores ({len(gt_positions)} sites):")
    print(f"  {'Pos':>4}  {'9-mer':>11}  {'MaxEnt':>7}  {'Class':>12}  {'U1 dG':>7}  {'Fn':>3}")
    print(f"  {'-'*55}")

    for pos in gt_positions:
        s9 = pos - 3
        e9 = pos + 6
        s11 = pos - 3
        e11 = pos + 8

        if 0 <= s9 and e9 <= len(seq):
            niner = seq[s9:e9]
            maxent = score_donor_maxentscan(niner)
            maxent_cls = classify_maxent_score(maxent)
        else:
            niner = "N/A"
            maxent = -999.0
            maxent_cls = "N/A"

        if 0 <= s11 and e11 <= len(seq):
            eleven = seq[s11:e11]
            dg = calculate_u1_duplex_dg(eleven)
            u1_dg = dg.delta_g
            fn = "Y" if dg.is_functional else "n"
        else:
            u1_dg = 0.0
            fn = "?"

        print(f"  {pos:>4}  {niner:>11}  {maxent:>7.2f}  {maxent_cls:>12}  {u1_dg:>7.2f}  {fn:>3}")

    # AG acceptor scores
    print(f"\n  {label} AG acceptor scores ({len(ag_positions)} sites):")
    print(f"  {'Pos':>4}  {'AccScore':>8}  {'Type':>6}  {'Region':>8}")
    print(f"  {'-'*35}")

    for pos in ag_positions:
        # Determine location type
        codon_idx = pos // 3
        pos_in_codon = pos % 3
        if pos_in_codon == 2 and pos + 1 < len(seq):
            loc = "junct"
        else:
            loc = "inter"

        # Determine region
        if pos < 75:
            region = "GGGGS5"
        elif pos < 432:
            region = "VHH3"
        else:
            region = "GGGGS1"

        # Score
        ag_scores = score_ag_at_pos(seq, pos)
        acc_score = ag_scores["acceptor_score"]

        if acc_score is not None:
            print(f"  {pos:>4}  {acc_score:>8.2f}  {loc:>6}  {region:>8}")
        else:
            print(f"  {pos:>4}  {'N/A':>8}  {loc:>6}  {region:>8}")


def print_codon_diff(before_codons: List[str], after_codons: List[str]):
    """Print a table of codon changes between original and final sequences."""
    changed = []
    for i, (old, new) in enumerate(zip(before_codons, after_codons)):
        if old != new:
            changed.append((i, old, new))

    if not changed:
        print("\n  No codon changes made.")
        return

    print(f"\n  Codon diff ({len(changed)} codons changed):")
    print(f"  {'Idx':>4}  {'AA':>3}  {'Before':>6}  {'After':>6}  {'Region':>8}")
    print(f"  {'-'*40}")

    for idx, old, new in changed:
        aa = CODON_TABLE.get(old, "?")
        if idx < 25:
            region = "GGGGS5"
        elif idx < 144:
            region = "VHH3"
        else:
            region = "GGGGS1"
        print(f"  {idx:>4}  {aa:>3}  {old:>6}  {new:>6}  {region:>8}")


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 80)
    print("Full Intron Body Codon Optimization (DNAchisel + MaxEntScan)")
    print("=" * 80)

    # Verify starting protein
    vhh3_protein_check = translate(VHH3_ORIGINAL)
    assert vhh3_protein_check == VHH3_PROTEIN, "VHH3 protein mismatch before optimization!"

    original_body = GGGGS5_ORIGINAL + VHH3_ORIGINAL + GGGGS1_ORIGINAL
    print(f"\nInput: {len(original_body)} bp = GGGGS5({len(GGGGS5_ORIGINAL)}) + "
          f"VHH3({len(VHH3_ORIGINAL)}) + GGGGS1({len(GGGGS1_ORIGINAL)})")
    print_gt_ag_summary("Original body", original_body)

    # ========== STAGE 0: GGGGS Design ==========
    print("\n" + "=" * 80)
    print("STAGE 0: Design GT/AG-free GGGGS linkers")
    print("=" * 80)

    ggggs5_new, ggggs1_new = design_ggggs_zero_gtag()

    print(f"\n  GGGGS5 new (75 bp): {ggggs5_new}")
    print(f"  GGGGS1 new (15 bp): {ggggs1_new}")
    print_gt_ag_summary("GGGGS5 original", GGGGS5_ORIGINAL)
    print_gt_ag_summary("GGGGS5 new", ggggs5_new)
    print_gt_ag_summary("GGGGS1 original", GGGGS1_ORIGINAL)
    print_gt_ag_summary("GGGGS1 new", ggggs1_new)

    # Show protein verification
    print(f"\n  GGGGS5 protein: {translate(ggggs5_new)}")
    print(f"  GGGGS1 protein: {translate(ggggs1_new)}")

    # ========== STAGE 1: DNAchisel VHH3 ==========
    print("\n" + "=" * 80)
    print("STAGE 1: DNAchisel optimization of VHH3")
    print("=" * 80)

    vhh3_opt, dnachisel_stats = optimize_vhh3_dnachisel(VHH3_ORIGINAL)

    print(f"\n  DNAchisel results (VHH3 only, {len(VHH3_ORIGINAL)} bp):")
    print(f"    GT: {dnachisel_stats['before_gt']} -> {dnachisel_stats['after_gt']}")
    print(f"    AG: {dnachisel_stats['before_ag']} -> {dnachisel_stats['after_ag']}")

    # Fix boundaries
    vhh3_codons = split_codons(vhh3_opt)
    vhh3_codons, boundary_changes = fix_boundary_gtag(ggggs5_new, vhh3_codons, ggggs1_new)
    vhh3_opt = join_codons(vhh3_codons)

    if boundary_changes:
        print(f"\n  Boundary fixes: {len(boundary_changes)}")
        for c in boundary_changes:
            print(f"    VHH3 codon {c['codon_idx']} ({c['aa']}): "
                  f"{c['old_codon']} -> {c['new_codon']}  [{c['reason']}]")
    else:
        print(f"\n  No boundary fixes needed.")

    # Verify VHH3 protein
    assert translate(vhh3_opt) == VHH3_PROTEIN, "Stage 1 broke VHH3 protein!"
    print(f"  VHH3 protein: PASS")

    # Build full body for Stage 2
    ggggs5_codons = split_codons(ggggs5_new)
    vhh3_codons = split_codons(vhh3_opt)
    ggggs1_codons = split_codons(ggggs1_new)
    all_codons = ggggs5_codons + vhh3_codons + ggggs1_codons
    assert len(all_codons) == 149, f"Expected 149 codons, got {len(all_codons)}"

    stage1_seq = join_codons(all_codons)
    print_gt_ag_summary("After Stage 0+1 (full body)", stage1_seq)

    # ========== STAGE 2a: GT MaxEntScan refinement ==========
    print("\n" + "=" * 80)
    print("STAGE 2a: Optimize context at remaining GT donor sites")
    print("=" * 80)

    MUTABLE_START = 0
    MUTABLE_END = 149

    all_codons, p2a_changes = stage2a_optimize_gt_context(
        all_codons, (MUTABLE_START, MUTABLE_END),
    )

    p2a_seq = join_codons(all_codons)
    print(f"\n  Stage 2a changes: {len(p2a_changes)}")
    for c in p2a_changes:
        print(f"    Codon {c['codon_idx']} ({c['aa']}): "
              f"{c['old_codon']} -> {c['new_codon']}  [{c['reason']}]")
    print_gt_ag_summary("After Stage 2a", p2a_seq)

    # Verify protein preserved
    vhh3_after_2a = join_codons(all_codons[25:144])
    assert translate(vhh3_after_2a) == VHH3_PROTEIN, "Stage 2a broke VHH3 protein!"
    assert translate(join_codons(all_codons[:25])) == "GGGGS" * 5, "Stage 2a broke GGGGS5 protein!"
    assert translate(join_codons(all_codons[144:])) == "GGGGS", "Stage 2a broke GGGGS1 protein!"

    # ========== STAGE 2b: AG MaxEntScan refinement ==========
    print("\n" + "=" * 80)
    print("STAGE 2b: Optimize context at remaining AG acceptor sites")
    print("=" * 80)

    all_codons, p2b_changes = stage2b_optimize_ag_context(
        all_codons, (MUTABLE_START, MUTABLE_END),
    )

    p2b_seq = join_codons(all_codons)
    print(f"\n  Stage 2b changes: {len(p2b_changes)}")
    for c in p2b_changes:
        print(f"    Codon {c['codon_idx']} ({c['aa']}): "
              f"{c['old_codon']} -> {c['new_codon']}  [{c['reason']}]")
    print_gt_ag_summary("After Stage 2b", p2b_seq)

    # Verify protein preserved
    vhh3_after_2b = join_codons(all_codons[25:144])
    assert translate(vhh3_after_2b) == VHH3_PROTEIN, "Stage 2b broke VHH3 protein!"
    assert translate(join_codons(all_codons[:25])) == "GGGGS" * 5, "Stage 2b broke GGGGS5 protein!"
    assert translate(join_codons(all_codons[144:])) == "GGGGS", "Stage 2b broke GGGGS1 protein!"

    # ========== FINAL REPORT ==========
    print("\n" + "=" * 80)
    print("FINAL REPORT")
    print("=" * 80)

    final_seq = join_codons(all_codons)
    final_vhh3 = join_codons(all_codons[25:144])

    # Before / After counts
    orig_gt = count_dinuc(original_body, "GT")
    orig_ag = count_dinuc(original_body, "AG")
    final_gt = count_dinuc(final_seq, "GT")
    final_ag = count_dinuc(final_seq, "AG")

    print(f"\n  GT count: {orig_gt} -> {final_gt} (removed {orig_gt - final_gt})")
    print(f"  AG count: {orig_ag} -> {final_ag} (removed {orig_ag - final_ag})")
    print(f"  Total forbidden: {orig_gt + orig_ag} -> {final_gt + final_ag}")

    # Break down by region
    print(f"\n  By region:")
    for name, start, end in [("GGGGS5", 0, 75), ("VHH3", 75, 432), ("GGGGS1", 432, 447)]:
        orig_region = original_body[start:end]
        final_region = final_seq[start:end]
        print(f"    {name:>6}: GT {count_dinuc(orig_region, 'GT'):>2} -> {count_dinuc(final_region, 'GT'):>2}, "
              f"AG {count_dinuc(orig_region, 'AG'):>2} -> {count_dinuc(final_region, 'AG'):>2}")

    # Per-site scores
    print_per_site_scores(original_body, "BEFORE")
    print_per_site_scores(final_seq, "AFTER")

    # Codon diff table
    original_codons = split_codons(original_body)
    print_codon_diff(original_codons, all_codons)

    # Amino acid preservation verification
    print(f"\n  Amino acid preservation:")
    print(f"    Original VHH3 protein: {VHH3_PROTEIN[:40]}...")
    final_protein = translate(final_vhh3)
    print(f"    Final VHH3 protein:    {final_protein[:40]}...")
    protein_match = final_protein == VHH3_PROTEIN
    print(f"    VHH3 match: {'PASS' if protein_match else 'FAIL'}")
    assert protein_match, "CRITICAL: VHH3 protein sequence changed!"

    # Full body protein
    ggggs5_protein = translate(join_codons(all_codons[:25]))
    ggggs1_protein = translate(join_codons(all_codons[144:]))
    full_protein_ok = (ggggs5_protein == "GGGGS" * 5 and
                       ggggs1_protein == "GGGGS" and
                       protein_match)
    print(f"    GGGGS5 protein: {'PASS' if ggggs5_protein == 'GGGGS' * 5 else 'FAIL'}")
    print(f"    GGGGS1 protein: {'PASS' if ggggs1_protein == 'GGGGS' else 'FAIL'}")

    # Restriction sites
    n_restriction = count_restriction_sites(final_seq)
    print(f"\n  Restriction sites (BsaI/BsmBI/BbsI): {n_restriction} {'PASS' if n_restriction == 0 else 'FAIL'}")

    # No new GT/AG check (compare Stage 1 → final)
    stage1_gt_set = set(find_dinuc_positions(stage1_seq, "GT"))
    stage1_ag_set = set(find_dinuc_positions(stage1_seq, "AG"))
    final_gt_set = set(find_dinuc_positions(final_seq, "GT"))
    final_ag_set = set(find_dinuc_positions(final_seq, "AG"))

    new_gt = final_gt_set - stage1_gt_set
    new_ag = final_ag_set - stage1_ag_set
    print(f"\n  New GT (vs Stage 1): {len(new_gt)} {'PASS' if len(new_gt) == 0 else 'FAIL'}")
    if new_gt:
        for pos in sorted(new_gt):
            print(f"    NEW GT at position {pos}: {final_seq[max(0,pos-3):pos+5]}")
    print(f"  New AG (vs Stage 1): {len(new_ag)} {'PASS' if len(new_ag) == 0 else 'FAIL'}")
    if new_ag:
        for pos in sorted(new_ag):
            print(f"    NEW AG at position {pos}: {final_seq[max(0,pos-3):pos+5]}")

    # Score regression check
    print(f"\n  Score regression check (no GT site made worse vs Stage 1):")
    gt_all_pass = True
    for pos in sorted(final_gt_set & stage1_gt_set):
        stage1_scores = score_gt_at_pos(stage1_seq, pos)
        final_scores = score_gt_at_pos(final_seq, pos)
        if (stage1_scores["maxent_score"] is not None and final_scores["maxent_score"] is not None):
            if final_scores["maxent_score"] > stage1_scores["maxent_score"] + 0.01:
                print(f"    FAIL: GT at {pos} MaxEnt {stage1_scores['maxent_score']:.2f} "
                      f"-> {final_scores['maxent_score']:.2f} (WORSE)")
                gt_all_pass = False
    if gt_all_pass:
        print(f"    All remaining GT sites: MaxEnt scores <= Stage 1 — PASS")

    # Output the final sequences
    final_ggggs5 = join_codons(all_codons[:25])
    final_ggggs1 = join_codons(all_codons[144:])

    print(f"\n{'='*80}")
    print(f"OPTIMIZED SEQUENCES")
    print(f"{'='*80}")

    print(f"\nGGGGS5_OPT = \"{final_ggggs5}\"")
    print(f"\nVHH3_OPT = (")
    for i in range(0, len(final_vhh3), 60):
        chunk = final_vhh3[i:i+60]
        print(f'    "{chunk}"')
    print(f")")
    print(f"\nGGGGS1_OPT = \"{final_ggggs1}\"")

    print(f"\n{'='*80}")
    print(f"OPTIMIZED FULL INTRON BODY ({len(final_seq)} bp)")
    print(f"{'='*80}")
    print(f"\nINTRON_BODY = GGGGS5_OPT + VHH3_OPT + GGGGS1_OPT")

    # Summary verdict
    print(f"\n{'='*80}")
    print(f"VERDICT")
    print(f"{'='*80}")
    checks = [
        ("VHH3 protein preserved", protein_match),
        ("GGGGS5 protein correct", ggggs5_protein == "GGGGS" * 5),
        ("GGGGS1 protein correct", ggggs1_protein == "GGGGS"),
        ("GGGGS5 zero GT/AG", count_dinuc(join_codons(all_codons[:25]), "GT") == 0 and
                               count_dinuc(join_codons(all_codons[:25]), "AG") == 0),
        ("GGGGS1 zero GT/AG", count_dinuc(join_codons(all_codons[144:]), "GT") == 0 and
                               count_dinuc(join_codons(all_codons[144:]), "AG") == 0),
        ("GT reduced", final_gt < orig_gt),
        ("AG reduced", final_ag < orig_ag),
        ("No new GT (vs Stage 1)", len(new_gt) == 0),
        ("No new AG (vs Stage 1)", len(new_ag) == 0),
        ("No score regression", gt_all_pass),
        ("No restriction sites", n_restriction == 0),
    ]
    all_ok = True
    for name, passed in checks:
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_ok = False
        print(f"  [{status}] {name}")

    if all_ok:
        print(f"\n  ALL CHECKS PASSED — optimized sequence is ready for use.")
    else:
        print(f"\n  SOME CHECKS FAILED — review issues above.")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
