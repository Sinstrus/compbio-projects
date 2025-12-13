# level2/step_analysis.py
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

@app.command()
def main(
    input_x: str = typer.Option(..., "--input-x"),
    regions_str: str = typer.Option("", "--regions"),
    model_path: Path = typer.Option(..., "--model-path"),
    out_file: Path = typer.Option("FilterPass_Seqs_step_analysis.tsv", "--out"),
    
    # Flanks
    left_flank: str = "",
    right_flank: str = "",
    
    # Config
    N: int = 50,
    M: int = 3,
    seed: int = 42,
    weights: str = "0.8,0.1,0.1", 
    
    # Selection thresholds
    donor_thresh: float = 0.5,
    acceptor_thresh: float = 0.5,
):
    """
    Decompose walks that trigger the 'Psex Donor Drop / Psex Acceptor Stable' phenotype.
    """
    t_start = time.time()
    
    # 1. Parse Inputs
    X_clean, intervals_brackets = regions.parse_brackets(input_x)
    intervals_regions = []
    if regions_str:
        try:
            intervals_regions = regions.apply_region_patterns(X_clean, regions_str)
        except ValueError:
            pass
            
    all_intervals = intervals_brackets + intervals_regions
    editable = [(iv.start, iv.end) for iv in all_intervals]
    editable = sorted(list(set(editable)))
    w_parts = [float(x) for x in weights.split(",")]

    # 2. Generate Random Walks
    rprint(f"[bold]1. Generating {N} random walks...[/bold]")
    walks = mutate.self_avoiding_walks(X_clean, editable, N, M, tuple(w_parts), seed)

    # 3. Batch Score (Pass 1)
    all_seqs = [X_clean] + [w.seq for w in walks]
    rprint(f"[bold blue]2. Scoring initial batch ({len(all_seqs)} seqs)...[/bold blue]")
    results_pass1 = openspliceai.score_batch(all_seqs, left_flank, right_flank, model_path)
    
    res0 = results_pass1[0]
    walk_results = results_pass1[1:]
    
    # 4. Find Quartet
    quartet = mutate.pick_quartet(
        res0.donors, res0.acceptors, 
        d_thr=donor_thresh, a_thr=acceptor_thresh,
        step=0.05, floor=0.01, min_sep=3
    )
    if not quartet:
        rprint("[bold red]No valid quartet found![/bold red]")
        raise typer.Exit(code=1)

    # 5. Filter and Decompose
    # Criteria: d_psex_acceptor > 0.01 AND d_psex_donor < -0.01
    
    pX_dP = quartet.psex_donor.prob
    pX_aP = quartet.psex_acceptor.prob
    
    sequences_to_score = [] # List of (seq, metadata_dict)
    
    rprint("[bold]3. Filtering for rare phenotype...[/bold]")
    pass_count = 0
    
    def get_prob(sites, pos):
        # Simple window +/- 5
        matches = [s.prob for s in sites if abs(s.pos - pos) <= 5]
        return max(matches) if matches else 0.0

    for i, (w, resY) in enumerate(zip(walks, walk_results)):
        pY_dP = get_prob(resY.donors, quartet.psex_donor.pos)
        pY_aP = get_prob(resY.acceptors, quartet.psex_acceptor.pos)
        
        d_dP = pY_dP - pX_dP
        d_aP = pY_aP - pX_aP
        
        if d_aP > 0.01 and d_dP < -0.01:
            pass_count += 1
            # Analyze steps
            muts = w.mutations.split(";")
            for step_idx, m in enumerate(muts):
                # Apply THIS single mutation to original X
                seq_step = mutate.apply_mutation_step(X_clean, m)
                
                meta = {
                    "parent_walk_idx": i,
                    "full_mutations": w.mutations,
                    "step_number": step_idx + 1,
                    "step_mutation": m
                }
                sequences_to_score.append((seq_step, meta))

    if pass_count == 0:
        rprint("[yellow]No sequences passed the phenotype filter.[/yellow]")
        return

    rprint(f"[green]Found {pass_count} sequences passing filter![/green]")
    rprint(f"[bold blue]4. Scoring {len(sequences_to_score)} decomposition steps...[/bold blue]")
    
    # 6. Batch Score (Pass 2 - The Decomposed Steps)
    seqs_only = [x[0] for x in sequences_to_score]
    results_pass2 = openspliceai.score_batch(seqs_only, left_flank, right_flank, model_path)
    
    # 7. Write Output
    with open(out_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        header = [
            "parent_walk_idx", "step_number", "step_mutation",
            "full_walk_mutations", "seq_len", "Y_seq",
            "pX_psex_donor", "pY_psex_donor", "d_psex_donor",
            "pX_psex_acceptor", "pY_psex_acceptor", "d_psex_acceptor"
        ]
        writer.writerow(header)
        
        for (seq, meta), resY in zip(sequences_to_score, results_pass2):
            pY_dP = get_prob(resY.donors, quartet.psex_donor.pos)
            pY_aP = get_prob(resY.acceptors, quartet.psex_acceptor.pos)
            
            writer.writerow([
                meta["parent_walk_idx"],
                meta["step_number"],
                meta["step_mutation"],
                meta["full_mutations"],
                len(seq),
                seq,
                f"{pX_dP:.4f}", f"{pY_dP:.4f}", f"{pY_dP - pX_dP:.4f}",
                f"{pX_aP:.4f}", f"{pY_aP:.4f}", f"{pY_aP - pX_aP:.4f}",
            ])
            
    rprint(f"[bold green]Done! Detailed step analysis written to {out_file}[/bold green]")

if __name__ == "__main__":
    app()