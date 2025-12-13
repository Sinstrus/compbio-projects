# level2/evolution.py
import typer
import csv
import time
from pathlib import Path
from rich.console import Console
from rich import print as rprint

# Import our modules
from . import regions, mutate
from .model_backends import openspliceai

app = typer.Typer(add_completion=False)
console = Console()

# --- Helper Functions (GTR / Acceptor Motifs) ---
def is_gtr(seq: str, pos: int) -> bool:
    if pos < 0 or pos + 3 > len(seq): return False
    return seq[pos] == "G" and seq[pos+1] == "T" and seq[pos+2] in "AG"

def is_acceptor_motif(seq: str, pos: int) -> bool:
    if pos < 0 or pos + 7 > len(seq): return False
    for i in [0, 1, 2, 4]:
        if seq[pos+i] not in "CT": return False
    if seq[pos+5] != "A": return False
    if seq[pos+6] != "G": return False
    return True

def check_new_motifs(parent_seq: str, mutations_str: str) -> bool:
    if not mutations_str: return False
    indices = list(range(len(parent_seq)))
    temp_indices = indices[:]
    steps = mutations_str.split(";")
    for step in steps:
        if not step.strip(): continue
        parts = step.strip().split()
        op_part = parts[0]
        op, pos_s = op_part.split(":")
        pos = int(pos_s)
        if op == "ins": temp_indices.insert(pos, -1)
        elif op == "del":
            if pos < len(temp_indices): temp_indices.pop(pos)

    child_seq = mutate.apply_mutation_step(parent_seq, steps[0])
    for s in steps[1:]:
        child_seq = mutate.apply_mutation_step(child_seq, s)

    # Scan GTR
    for i in range(len(child_seq) - 2):
        if is_gtr(child_seq, i):
            origin = temp_indices[i : i+3]
            if -1 in origin: return True 
            if not (origin[1] == origin[0] + 1 and origin[2] == origin[1] + 1): return True
            parent_idx = origin[0]
            if not is_gtr(parent_seq, parent_idx): return True

    # Scan Acceptor
    for i in range(len(child_seq) - 6):
        if is_acceptor_motif(child_seq, i):
            origin = temp_indices[i : i+7]
            if -1 in origin: return True
            is_continuous = True
            for k in range(6):
                if origin[k+1] != origin[k] + 1:
                    is_continuous = False; break
            if not is_continuous: return True
            parent_idx = origin[0]
            if not is_acceptor_motif(parent_seq, parent_idx): return True
    return False

@app.command()
def main(
    input_x: str = typer.Option(..., "--input-x"),
    regions_str: str = typer.Option("", "--regions"),
    model_path: Path = typer.Option(..., "--model-path"),
    out_file: Path = typer.Option("evolution_history.tsv", "--out"),
    all_variants_out: Path = typer.Option("evolution_all_variants.tsv", "--all-variants-out", help="Detailed log of every scored variant"),
    fasta_out: Path = typer.Option(None, "--fasta-out", help="Optional path to save winner sequences as FASTA"),
    
    left_flank: str = "",
    right_flank: str = "",
    
    rounds: int = typer.Option(3, "--rounds"),
    N: int = typer.Option(50, "--n"), 
    M: int = typer.Option(3, "--m"),
    seed: int = 42,
    weights: str = "0.8,0.1,0.1", 
    
    donor_thresh: float = 0.5,
    acceptor_thresh: float = 0.5,
):
    """
    Evolutionary Walk with Full Logging, Round 0 tracking, and FASTA output.
    """
    t_start_main = time.time()
    OFF_TARGET_THRESH = 0.008

    # --- Setup ---
    X_clean, intervals_brackets = regions.parse_brackets(input_x)
    intervals_regions = []
    if regions_str:
        try: intervals_regions = regions.apply_region_patterns(X_clean, regions_str)
        except ValueError: pass
    all_intervals = intervals_brackets + intervals_regions
    editable = [(iv.start, iv.end) for iv in all_intervals]
    editable = sorted(list(set(editable)))
    w_parts = [float(x) for x in weights.split(",")]

    if not editable:
        rprint("[bold red]Error:[/bold red] No editable regions found.")
        raise typer.Exit(1)

    rprint(f"[bold green]Initializing Evolution on X (len={len(X_clean)})...[/bold green]")
    
    # Initial Quartet
    res0 = openspliceai.score_batch([X_clean], left_flank, right_flank, model_path)[0]
    quartet = mutate.pick_quartet(res0.donors, res0.acceptors, donor_thresh, acceptor_thresh, 0.05, 0.01, 3)
    if not quartet:
        rprint("[bold red]No valid quartet found in SeqX![/bold red]")
        raise typer.Exit(1)
    
    rprint(f"Target Quartet:\n  dA: {quartet.exonA_donor.pos}\n  aP: {quartet.psex_acceptor.pos}\n  dP: {quartet.psex_donor.pos}\n  aB: {quartet.exonB_acceptor.pos}")

    current_seq = X_clean
    current_name = "SeqX"
    history = [] 
    all_variants_log = []

    # --- Log Round 0 (Initial State) ---
    # We add the starting point to the log so it appears in plots
    all_variants_log.append({
        "Round": 0, "Variant_Idx": 0, "Status": "Initial", "Is_Winner": True,
        "p_dA": quartet.exonA_donor.prob,
        "p_aP": quartet.psex_acceptor.prob,
        "p_dP": quartet.psex_donor.prob,
        "p_aB": quartet.exonB_acceptor.prob,
        "d_aP": 0.0, "d_dP": 0.0, "Metric": 0.0, "Mutations": ""
    })

    def get_quartet_probs(res, q):
        def get_p(sites, pos):
            matches = [s.prob for s in sites if abs(s.pos - pos) <= 5]
            return max(matches) if matches else 0.0
        return {
            "dA": get_p(res.donors, q.exonA_donor.pos),
            "aP": get_p(res.acceptors, q.psex_acceptor.pos),
            "dP": get_p(res.donors, q.psex_donor.pos),
            "aB": get_p(res.acceptors, q.exonB_acceptor.pos),
        }
    
    def is_quartet_site(pos, kind):
        targets = []
        if kind == "Donor": targets = [quartet.exonA_donor.pos, quartet.psex_donor.pos]
        else: targets = [quartet.psex_acceptor.pos, quartet.exonB_acceptor.pos]
        for t in targets:
            if abs(pos - t) <= 5: return True
        return False

    # --- Evolutionary Loop ---
    for r in range(1, rounds + 1):
        rprint(f"\n[bold blue]=== ROUND {r} (Parent: {current_name}) ===[/bold blue]")
        
        attempts = 0
        max_attempts = 2
        winner = None
        round_candidates_log = []
        
        res_parent = openspliceai.score_batch([current_seq], left_flank, right_flank, model_path)[0]
        probs_parent = get_quartet_probs(res_parent, quartet)
        parent_donors_map = {s.pos: s.prob for s in res_parent.donors}
        parent_acceptors_map = {s.pos: s.prob for s in res_parent.acceptors}
        
        while attempts < max_attempts:
            attempts += 1
            current_seed = seed + (r * 100) + attempts 
            round_candidates_log = [] 

            # Generate
            rprint(f"  > Generating {N} candidates (Attempt {attempts})...")
            valid_walks = []
            gen_seed = current_seed
            
            # Track Generation Stats
            total_generated = 0
            motif_rejected = 0
            
            while len(valid_walks) < N:
                chunk = mutate.self_avoiding_walks(current_seq, editable, 20, M, tuple(w_parts), gen_seed)
                gen_seed += 1
                total_generated += len(chunk)
                for w in chunk:
                    if check_new_motifs(current_seq, w.mutations):
                        motif_rejected += 1
                    else:
                        valid_walks.append(w)
                    if len(valid_walks) >= N: break
            
            # Print Generation Stats
            rprint(f"    [dim]Gen Stats: {total_generated} generated, {motif_rejected} rejected by Motif-Filter ({(motif_rejected/total_generated)*100:.1f}%).[/dim]")
            
            # Score
            batch_seqs = [w.seq for w in valid_walks]
            batch_results = openspliceai.score_batch(batch_seqs, left_flank, right_flank, model_path)
            
            candidates = []
            ai_rejected_count = 0
            
            for idx, (w, res_child) in enumerate(zip(valid_walks, batch_results)):
                probs_child = get_quartet_probs(res_child, quartet)
                status = "Passed"
                
                # 1. AI Filter
                has_off_target = False
                for s in res_child.donors:
                    if s.prob > OFF_TARGET_THRESH:
                        if parent_donors_map.get(s.pos, 0.0) < OFF_TARGET_THRESH and not is_quartet_site(s.pos, "Donor"):
                            has_off_target = True; break
                if not has_off_target:
                    for s in res_child.acceptors:
                        if s.prob > OFF_TARGET_THRESH:
                            if parent_acceptors_map.get(s.pos, 0.0) < OFF_TARGET_THRESH and not is_quartet_site(s.pos, "Acceptor"):
                                has_off_target = True; break
                
                if has_off_target:
                    status = "AI_Reject"
                    ai_rejected_count += 1
                
                # 2. Phenotype Filter (Only check if not already rejected)
                d_aP = probs_child["aP"] - probs_parent["aP"]
                d_dP = probs_child["dP"] - probs_parent["dP"]
                metric = abs(d_aP - d_dP)

                if status == "Passed":
                    if not (d_aP > 0.01 and d_dP < -0.01 and probs_child["aP"] > 0.01 and probs_child["dP"] > 0.01):
                        status = "Pheno_Reject"
                
                variant_entry = {
                    "Round": r, "Variant_Idx": idx, "Status": status, "Is_Winner": False,
                    "p_dA": probs_child["dA"], "p_aP": probs_child["aP"], "p_dP": probs_child["dP"], "p_aB": probs_child["aB"],
                    "d_aP": d_aP, "d_dP": d_dP, "Metric": metric, "Mutations": w.mutations
                }
                round_candidates_log.append(variant_entry)

                if status == "Passed":
                    candidates.append({**variant_entry, "walk": w, "probs": probs_child, "deltas": {"aP": d_aP, "dP": d_dP}})

            # Print Evaluation Stats
            rprint(f"    [blue]Stats:[/blue] {ai_rejected_count} rejected by AI-Filter, {len(candidates)} passed Phenotype.")

            if candidates:
                candidates.sort(key=lambda x: x["Metric"], reverse=True)
                winner = candidates[0]
                for v in round_candidates_log:
                    if v["Variant_Idx"] == winner["Variant_Idx"] and v["Status"] == "Passed":
                         if v["Mutations"] == winner["Mutations"]:
                            v["Is_Winner"] = True; v["Status"] = "Winner"; break
                break
            else:
                rprint(f"    [yellow]No candidates passed filters.[/yellow]")
        
        all_variants_log.extend(round_candidates_log)

        if not winner:
            rprint(f"[bold red]Evolution stopped early at Round {r}.[/bold red]")
            break
            
        child_name = f"Seq_R{r}"
        w_data = winner
        rprint(f"  [green]Winner found![/green] {child_name} (Metric: {w_data['Metric']:.4f})")
        
        history.append({
            "Round": r, "Parent": current_name, "Child": child_name,
            "Mutations": w_data["walk"].mutations, "Seq": w_data["walk"].seq,
            "pParent_dA": probs_parent["dA"], "pParent_aP": probs_parent["aP"], "pParent_dP": probs_parent["dP"], "pParent_aB": probs_parent["aB"],
            "pChild_dA": w_data["probs"]["dA"], "pChild_aP": w_data["probs"]["aP"], "pChild_dP": w_data["probs"]["dP"], "pChild_aB": w_data["probs"]["aB"],
            "delta_aP": w_data["deltas"]["aP"], "delta_dP": w_data["deltas"]["dP"]
        })
        current_seq = w_data["walk"].seq
        current_name = child_name

    # Write Outputs
    if history:
        rprint(f"Writing history to {out_file}...")
        with open(out_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=[
                "Round", "Parent", "Child", "Mutations", "Seq",
                "pParent_dA", "pParent_aP", "pParent_dP", "pParent_aB",
                "pChild_dA", "pChild_aP", "pChild_dP", "pChild_aB",
                "delta_aP", "delta_dP"
            ], delimiter="\t")
            writer.writeheader()
            writer.writerows(history)

        if fasta_out:
            rprint(f"Writing winners to {fasta_out}...")
            with open(fasta_out, "w") as f:
                f.write(f">SeqX_Initial\n{X_clean}\n")
                for entry in history:
                    safe_muts = entry["Mutations"].replace(" ", "")
                    header = f"{entry['Child']}_Round{entry['Round']}_{safe_muts}"
                    f.write(f">{header}\n{entry['Seq']}\n")
            
    if all_variants_log:
        rprint(f"Writing detailed variants log to {all_variants_out}...")
        with open(all_variants_out, "w", newline="") as f:
            fnames = ["Round", "Variant_Idx", "Status", "Is_Winner", "p_dA", "p_aP", "p_dP", "p_aB", "d_aP", "d_dP", "Metric", "Mutations"]
            writer = csv.DictWriter(f, fieldnames=fnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(all_variants_log)

    rprint(f"[bold green]Done! (Total time: {time.time() - t_start_main:.2f}s)[/bold green]")

if __name__ == "__main__":
    app()