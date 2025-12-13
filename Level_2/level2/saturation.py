# level2/saturation.py
import typer
import csv
import time
import sys
from pathlib import Path
# We use standard print for pipeline safety
from rich import print as rprint

from . import regions, mutate
from .model_backends import openspliceai

app = typer.Typer(add_completion=False)

DNA = ["A", "C", "G", "T"]

@app.command()
def main(
    input_x: str = typer.Option(..., "--input-x"),
    model_path: Path = typer.Option(..., "--model-path"),
    out_file: Path = typer.Option("saturation_mutation_results.tsv", "--out"),
    quartet_out: Path = typer.Option("quartet_positions.tsv", "--quartet-out"),
    baseline_out: Path = typer.Option(None, "--baseline-out"),
    
    gain_delta_thresh: float = typer.Option(0.10, "--gain-delta"),
    
    left_flank: str = typer.Option("", "--left-flank"),
    right_flank: str = typer.Option("", "--right-flank"),
    donor_thresh: float = typer.Option(0.5, "--donor-thresh"),
    acceptor_thresh: float = typer.Option(0.5, "--acceptor-thresh"),
    
    batch_size: int = typer.Option(1000, "--batch-size"),
):
    t_start = time.time()
    
    # High-level progress message
    print(f"  > Starting Saturation for model: {model_path.name}", flush=True)
    
    try:
        X_clean, _ = regions.parse_brackets(input_x)
    except Exception as e:
        print(f"Error parsing input: {e}", file=sys.stderr)
        sys.exit(1)
    
    # 1. Baseline
    print("    [1/3] Scoring baseline...", flush=True)
    try:
        res0 = openspliceai.score_batch([X_clean], left_flank, right_flank, model_path)[0]
    except Exception as e:
        print(f"Error running OpenSpliceAI: {e}", file=sys.stderr)
        sys.exit(1)

    quartet = mutate.pick_quartet(res0.donors, res0.acceptors, donor_thresh, acceptor_thresh, 0.05, 0.01, 3)
    
    if not quartet:
        print(f"Error: No valid quartet found in SeqX", file=sys.stderr)
        sys.exit(1)
        
    q_pos = {
        "dA": quartet.exonA_donor.pos, "aP": quartet.psex_acceptor.pos,
        "dP": quartet.psex_donor.pos, "aB": quartet.exonB_acceptor.pos
    }

    # Baseline Maps
    base_donors = {i: 0.0 for i in range(len(X_clean))}
    base_acceptors = {i: 0.0 for i in range(len(X_clean))}
    for s in res0.donors: base_donors[s.pos] = s.prob
    for s in res0.acceptors: base_acceptors[s.pos] = s.prob

    if baseline_out:
        with open(baseline_out, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["Type", "Position", "Probability"])
            for pos in sorted(base_donors.keys()): writer.writerow(["Donor", pos, f"{base_donors[pos]:.4f}"])
            for pos in sorted(base_acceptors.keys()): writer.writerow(["Acceptor", pos, f"{base_acceptors[pos]:.4f}"])

    with open(quartet_out, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Site", "Position"])
        writer.writerow(["dA", quartet.exonA_donor.pos]); writer.writerow(["aP", quartet.psex_acceptor.pos])
        writer.writerow(["dP", quartet.psex_donor.pos]); writer.writerow(["aB", quartet.exonB_acceptor.pos])

    # 2. Generate Variants
    print("    [2/3] Generating variants...", flush=True)
    variants = [] 
    variants.append({"pos": -1, "mut_str": "Original", "seq": X_clean, "new_gtr": False, "new_acc": False, "type": "Original"})

    for i in range(len(X_clean)):
        base = X_clean[i]
        for b in DNA:
            if b == base: continue
            mut_str = f"sub:{i} {base}>{b}"
            seq = mutate.apply_mutation_step(X_clean, mut_str)
            ng, na = mutate.check_new_motifs(X_clean, mut_str)
            variants.append({"pos": i, "mut_str": mut_str, "seq": seq, "new_gtr": ng, "new_acc": na, "type": "sub"})
        
        mut_str = f"del:{i} {base}"
        seq = mutate.apply_mutation_step(X_clean, mut_str)
        ng, na = mutate.check_new_motifs(X_clean, mut_str)
        variants.append({"pos": i, "mut_str": mut_str, "seq": seq, "new_gtr": ng, "new_acc": na, "type": "del"})
        
        for b in DNA:
            mut_str = f"ins:{i} {b}"
            seq = mutate.apply_mutation_step(X_clean, mut_str)
            ng, na = mutate.check_new_motifs(X_clean, mut_str)
            variants.append({"pos": i, "mut_str": mut_str, "seq": seq, "new_gtr": ng, "new_acc": na, "type": "ins"})

    # 3. Batch Scoring
    print(f"    [3/3] Scoring {len(variants)} variants...", flush=True)
    all_seqs = [v["seq"] for v in variants]
    
    with open(out_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        header = ["Position", "Mutation", "Type", "Sequence", "New_Acceptor_Motif", "New_GTR_Motif", "p_dA", "p_aP", "p_dP", "p_aB", "New_Site_Info", "Max_New_Prob"]
        writer.writerow(header)
        
        for i in range(0, len(variants), batch_size):
            chunk_end = min(i + batch_size, len(variants))
            chunk_seqs = all_seqs[i:chunk_end]
            
            # THIS is the message you want to see:
            print(f"      - Processing batch {i}-{chunk_end}...", flush=True)
            
            try:
                chunk_results = openspliceai.score_batch(chunk_seqs, left_flank, right_flank, model_path)
            except Exception as e:
                print(f"[red]Error scoring batch {i}: {e}[/red]", file=sys.stderr)
                sys.exit(1)
            
            for j, res in enumerate(chunk_results):
                idx = i + j
                v = variants[idx]
                
                def map_to_orig(mut_pos):
                    m_type, m_idx = v["type"], v["pos"]
                    if m_type == "ins": return m_idx if mut_pos == m_idx else (mut_pos if mut_pos < m_idx else mut_pos - 1)
                    elif m_type == "del": return mut_pos if mut_pos < m_idx else mut_pos + 1
                    return mut_pos

                def is_quartet(phys_pos, kind):
                    targets = [q_pos['dA'], q_pos['dP']] if kind == "Donor" else [q_pos['aP'], q_pos['aB']]
                    for t in targets:
                        if abs(phys_pos - t) <= 2: return True
                    return False

                new_sites = []; max_p = 0.0
                
                for s in res.donors:
                    orig_pos = map_to_orig(s.pos)
                    if is_quartet(orig_pos, "Donor"): continue
                    base_p = base_donors.get(orig_pos, 0.0) if orig_pos != -999 else 0.0
                    if (s.prob - base_p) > gain_delta_thresh:
                        new_sites.append(f"D@{s.pos}({s.prob:.2f})")
                        max_p = max(max_p, s.prob)

                for s in res.acceptors:
                    orig_pos = map_to_orig(s.pos)
                    if is_quartet(orig_pos, "Acceptor"): continue
                    base_p = base_acceptors.get(orig_pos, 0.0) if orig_pos != -999 else 0.0
                    if (s.prob - base_p) > gain_delta_thresh:
                        new_sites.append(f"A@{s.pos}({s.prob:.2f})")
                        max_p = max(max_p, s.prob)

                new_site_str = ";".join(new_sites) if new_sites else "None"

                def get_site_p(site_list, target_pos):
                    mapped = target_pos
                    if v["type"] == "ins" and v["pos"] <= target_pos: mapped += 1
                    elif v["type"] == "del" and v["pos"] < target_pos: mapped -= 1
                    matches = [s.prob for s in site_list if abs(s.pos - mapped) <= 5]
                    return max(matches) if matches else 0.0

                pdA = get_site_p(res.donors, q_pos['dA'])
                paP = get_site_p(res.acceptors, q_pos['aP'])
                pdP = get_site_p(res.donors, q_pos['dP'])
                paB = get_site_p(res.acceptors, q_pos['aB'])
                
                writer.writerow([
                    v["pos"], v["mut_str"], v["type"], v["seq"],
                    v["new_acc"], v["new_gtr"], 
                    f"{pdA:.4f}", f"{paP:.4f}", f"{pdP:.4f}", f"{paB:.4f}",
                    new_site_str, f"{max_p:.4f}"
                ])

    print(f"    Done! Results in {out_file}", flush=True)

if __name__ == "__main__":
    app()