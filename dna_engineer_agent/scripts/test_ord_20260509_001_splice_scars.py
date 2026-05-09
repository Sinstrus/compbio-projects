#!/usr/bin/env python3
"""
Splice scar test for ORD-20260509-001 (AVD428-448, AVD465-485).

For every construct in the order:
  1. Simulate splicing: remove intron body [ins_pt+3 : ins_pt+504] from each construct
  2. Assert scar ≤ 5 aa at insertion point (expected: exactly 3 aa — QES or KES)
  3. Assert plasmid sequence upstream of scar is identical to AVD002
  4. Assert plasmid sequence downstream of scar is identical to AVD002
  5. Report exact scar amino acid sequence per construct

Intron structure (510 bp total insert):
  [0:3]     exonic context (CAG/AAG) — stays in spliced mRNA
  [3:504]   intron body (GT...AG, 501 bp) — spliced OUT
  [504:510] HBB exon 3 pad (GAATCC) — stays in spliced mRNA
  => scar = exon ctx (3 bp) + exon pad (6 bp) = 9 bp = 3 aa
"""

import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

ROOT_DIR = Path(__file__).parent.parent
CONSTRUCTS_DIR = ROOT_DIR / "constructs"

VP1_START    = 2378
EXON_CTX     = 3     # bp retained upstream of intron
INTRON_BODY  = 501   # bp spliced out
EXON_PAD     = 6     # bp retained downstream of intron
INTRON_SIZE  = EXON_CTX + INTRON_BODY + EXON_PAD   # 510
SCAR_DNA     = EXON_CTX + EXON_PAD                 # 9 bp → 3 aa
MAX_SCAR_AA  = 5

SITES = [
    # (aa_pos, region)
    (490, "VR5"),
    (491, "VR5"),
    (492, "VR5"),
    (493, "VR5"),
    (494, "VR5"),
    (495, "VR5"),
    (528, "VR6pre"),
    (529, "VR6pre"),
    (583, "VR8"),
    (585, "VR8"),
    (586, "VR8"),
    (587, "VR8"),
    (588, "VR8"),
    (589, "VR8"),
    (591, "VR8"),
    (592, "VR8"),
    (593, "VR8"),
    (594, "VR8post"),
    (595, "VR8post"),
    (661, "Ctloop"),
    (706, "VR9"),
]

assert len(SITES) == 21

MODULES = [
    # (avd_start, donor_11mer, label, expected_scar_aa)
    (428, "CAGGTAATTGG", "Tier7 baseline",    "QES"),
    (465, "AAGGTAAGGCA", "Module1 (AVD388)",  "KES"),
]


def load_seq(avd_num: int) -> str:
    matches = list(CONSTRUCTS_DIR.glob(f"AVD{avd_num:03d}-*.gb"))
    if not matches:
        raise FileNotFoundError(f"AVD{avd_num:03d} not found in {CONSTRUCTS_DIR}")
    return str(SeqIO.read(str(matches[0]), "genbank").seq).upper()


def simulate_splice(seq: str, ins_pt: int) -> str:
    """Remove the 501 bp intron body at ins_pt+3..ins_pt+504, keep exon ctx and pad."""
    return seq[: ins_pt + EXON_CTX] + seq[ins_pt + EXON_CTX + INTRON_BODY :]


def get_scar_dna(seq: str, ins_pt: int) -> str:
    return seq[ins_pt : ins_pt + EXON_CTX] + seq[ins_pt + EXON_CTX + INTRON_BODY : ins_pt + INTRON_SIZE]


def main():
    print("Loading AVD002 reference...")
    avd002 = load_seq(2)
    avd002_len = len(avd002)

    failures = []
    results = []

    print(f"\n{'AVD':<6} {'Region':<8} {'aa':<5} {'Module':<22} {'Scar DNA':<12} {'Scar AA':<8} "
          f"{'Scar≤5?':<8} {'Upstream':<10} {'Downstream'}")
    print("-" * 100)

    for avd_start, donor, label, expected_scar in MODULES:
        for idx, (aa, region) in enumerate(SITES):
            avd_num = avd_start + idx
            ins_pt  = VP1_START + aa * 3

            seq = load_seq(avd_num)

            # Verify construct is the expected size (AVD002 + intron)
            assert len(seq) == avd002_len + INTRON_SIZE, (
                f"AVD{avd_num}: unexpected length {len(seq)} "
                f"(expected {avd002_len + INTRON_SIZE})"
            )

            spliced   = simulate_splice(seq, ins_pt)
            scar_dna  = get_scar_dna(seq, ins_pt)
            scar_aa   = str(Seq(scar_dna).translate())

            # Core assertions
            scar_ok       = len(scar_aa) <= MAX_SCAR_AA
            upstream_ok   = spliced[:ins_pt] == avd002[:ins_pt]
            downstream_ok = spliced[ins_pt + SCAR_DNA :] == avd002[ins_pt:]

            all_ok = scar_ok and upstream_ok and downstream_ok

            tag = "PASS" if all_ok else "FAIL"
            if not all_ok:
                failures.append({
                    "avd": avd_num, "aa": aa, "region": region, "label": label,
                    "scar_ok": scar_ok, "upstream_ok": upstream_ok, "downstream_ok": downstream_ok,
                    "scar_aa": scar_aa, "expected": expected_scar,
                })

            results.append((avd_num, region, aa, label, scar_dna, scar_aa,
                             scar_ok, upstream_ok, downstream_ok, expected_scar))

            up_str   = "OK" if upstream_ok   else "FAIL"
            down_str = "OK" if downstream_ok else "FAIL"
            scar_str = "OK" if scar_ok       else f"FAIL({len(scar_aa)}aa)"
            match    = "✓" if scar_aa == expected_scar else f"✗ (expected {expected_scar})"

            print(f"AVD{avd_num:<3} {region:<8} aa{aa:<4} {label:<22} {scar_dna:<12} "
                  f"{scar_aa} {match:<8} {scar_str:<8} {up_str:<10} {down_str}")

    # Summary
    n_total = len(results)
    n_pass  = n_total - len(failures)
    print(f"\n{'='*100}")
    print(f"Results: {n_pass}/{n_total} PASS")

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  AVD{f['avd']} aa{f['aa']} {f['region']} — "
                  f"scar_ok={f['scar_ok']}, upstream_ok={f['upstream_ok']}, "
                  f"downstream_ok={f['downstream_ok']}, scar={f['scar_aa']}")
        sys.exit(1)
    else:
        print("\nAll 42 constructs: intron splices cleanly, ≤5 aa scar, plasmid identity verified.")

    # Per-site scar table
    print(f"\n{'─'*60}")
    print("Scar summary (one row per site, both modules):")
    print(f"{'aa':<5} {'Region':<8}  {'Baseline scar (AVD428+)':<28}  {'Module1 scar (AVD465+)'}")
    print(f"{'─'*5} {'─'*8}  {'─'*28}  {'─'*24}")
    for idx, (aa, region) in enumerate(SITES):
        r_base = results[idx]
        r_mod1 = results[len(SITES) + idx]
        print(f"aa{aa:<4} {region:<8}  {r_base[5]} ({r_base[4]})               "
              f"{r_mod1[5]} ({r_mod1[4]})")


if __name__ == "__main__":
    main()
