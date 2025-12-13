# level2/filter_and_generate.py
import typer
import pandas as pd
import random
import csv
import time
from pathlib import Path
from rich.console import Console
from rich import print as rprint

from . import regions, mutate
from .model_backends import openspliceai

app = typer.Typer(add_completion=False)
console = Console()

@app.command()
def main(
    input_x: str = typer.Option(..., "--input-x"),
    tsv_file: Path = typer.Argument(..., help="Path to saturation_mutation_results.tsv"),
    model_path: Path = typer.Option(..., "--model-path"),
    out_file: Path = typer.Option("combinatorial_results.tsv", "--out"),
    quartet_file: Path = typer.Option("quartet_positions.tsv", "--quartet", help="Path to quartet coordinates"),
    
    # Filters
    p_dA_min: float = typer.Option(0.90, "--p_dA-min"),
    p_aB_min: float = typer.Option(0.90, "--p_aB-min"),
    p_aP_min: float = typer.Option(0.10, "--p_aP-min"),
    p_dP_min: float = typer.Option(0.01, "--p_dP-min"),
    d_aP_min: float = typer.Option(0.05, "--d_aP-min"),
    d_dP_max: float = typer.Option(0.05, "--d_dP-max"),
    max_new_site_prob: float = typer.Option(0.05, "--max-new-site-prob", help="Reject variant if AI sees ANY new site > this prob"),
    
    # Generation
    G: int = typer.Option(50, "--g"),
    pick_n: int = typer.Option(3, "--pick-n"),
    seed: int = 42,
    include_original: bool = typer.Option(False, "--include-original"),
    
    left_flank: str = "",
    right_flank: str = "",
):
    """
    Filter saturation hits (using AI & Motif checks) and generate combinatorial variants.
    """
    t_start = time.time()
    rng = random.Random(seed)
    
    X_clean, _ = regions.parse_brackets(input_x)
    df = pd.read_csv(tsv_file, sep="\t")
    
    q_pos = {}
    if quartet_file.exists():
        q_df = pd.read_csv(quartet_file, sep="\t")
        for _, row in q_df.iterrows():
            q_pos[row["Site"]] = int(row["Position"])
    else:
        rprint("[red]Error: Quartet positions file needed.[/red]")
        raise typer.Exit(1)

    baseline_row = df[df["Position"] == -1]
    if baseline_row.empty:
        rprint("[red]Error: Baseline row not found.[/red]")
        raise typer.Exit(1)
    
    p_aP_orig = baseline_row.iloc[0]["p_aP"]
    p_dP_orig = baseline_row.iloc[0]["p_dP"]

    # Calculate deltas
    df["d_aP"] = df["p_aP"] - p_aP_orig
    df["d_dP"] = df["p_dP"] - p_dP_orig
    
    if "Max_New_Prob" not in df.columns:
        df["Max_New_Prob"] = 0.0 
    else:
        df["Max_New_Prob"] = df["Max_New_Prob"].fillna(0.0)

    # Apply Filters
    candidates = df[
        (df["p_dA"] > p_dA_min) &
        (df["p_aB"] > p_aB_min) &
        (df["p_aP"] > p_aP_min) &
        (df["p_dP"] > p_dP_min) &
        (df["d_aP"] > d_aP_min) &
        (df["d_dP"] < d_dP_max) &
        (df["New_Acceptor_Motif"] == False) &
        (df["New_GTR_Motif"] == False) &
        (df["Max_New_Prob"] < max_new_site_prob) &
        (df["Position"] >= 0)
    ]
    
    valid_muts = candidates["Mutation"].tolist()
    if include_original: valid_muts.append("Original")
        
    rprint(f"[bold]Found {len(valid_muts)} mutations matching criteria.[/bold]")
    if len(valid_muts) < pick_n:
        rprint(f"[red]Not enough candidates ({len(valid_muts)}) to pick {pick_n}.[/red]")
        raise typer.Exit(1)

    generated_data = [] 
    for i in range(G):
        combo = rng.sample(valid_muts, pick_n)
        parsed_ops = []
        for m in combo:
            if m == "Original": continue
            parts = m.split()
            op_raw = parts[0]
            op, pos_s = op_raw.split(":")
            pos = int(pos_s)
            parsed_ops.append({"op": op, "pos": pos, "full": m})
            
        parsed_ops.sort(key=lambda x: x["pos"], reverse=True)
        
        seq = X_clean
        ordered_muts_str = []
        for item in parsed_ops:
            seq = mutate.apply_mutation_step(seq, item["full"])
            ordered_muts_str.append(item["full"])
            
        final_mut_str = ";".join(ordered_muts_str[::-1])
        if not final_mut_str: final_mut_str = "Original"
        
        generated_data.append({"seq": seq, "mutations": final_mut_str})

    rprint(f"[blue]Scoring {len(generated_data)} combinatorial variants...[/blue]")
    batch_seqs = [d["seq"] for d in generated_data]
    results = openspliceai.score_batch(batch_seqs, left_flank, right_flank, model_path)
    
    with open(out_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        header = ["Variant_Idx", "Mutations", "Seq_Len", "Sequence", "p_dA", "p_aP", "p_dP", "p_aB", "d_aP", "d_dP"]
        writer.writerow(header)
        
        for i, (d, res) in enumerate(zip(generated_data, results)):
            def get_p(site_list, orig_pos):
                matches = [s.prob for s in site_list if abs(s.pos - orig_pos) <= 10]
                return max(matches) if matches else 0.0

            pdA = get_p(res.donors, q_pos['dA'])
            paP = get_p(res.acceptors, q_pos['aP'])
            pdP = get_p(res.donors, q_pos['dP'])
            paB = get_p(res.acceptors, q_pos['aB'])
            
            delta_aP = paP - p_aP_orig
            delta_dP = pdP - p_dP_orig
            
            writer.writerow([
                i, d["mutations"], len(d["seq"]), d["seq"],
                f"{pdA:.4f}", f"{paP:.4f}", f"{pdP:.4f}", f"{paB:.4f}",
                f"{delta_aP:.4f}", f"{delta_dP:.4f}"
            ])

    rprint(f"[green]Done! Results written to {out_file}[/green]")

if __name__ == "__main__":
    app()