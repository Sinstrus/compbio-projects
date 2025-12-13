# level2/cli.py
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
    out_file: Path = typer.Option("results.tsv", "--out"),
    fasta_out: Path = typer.Option(None, "--fasta-out", help="Optional path to save sorted FASTA sequences"),
    
    # Flanks
    left_flank: str = "",
    right_flank: str = "",
    
    # Config
    N: int = 50,
    M: int = 3,
    seed: int = 42,
    weights: str = "0.8,0.1,0.1", # sub,del,ins
    
    # Selection thresholds
    donor_thresh: float = 0.5,
    acceptor_thresh: float = 0.5,
):
    """
    Level 2: Mutate editable regions and track SpliceAI quartet probabilities.
    """
    t_start_main = time.time()
    
    # 1. Parse Inputs
    X_clean, intervals_brackets = regions.parse_brackets(input_x)
    
    intervals_regions = []
    if regions_str:
        try:
            intervals_regions = regions.apply_region_patterns(X_clean, regions_str)
        except ValueError as e:
            rprint(f"[bold red]Error parsing regions:[/bold red] {e}")
            raise typer.Exit(code=1)

    all_intervals = intervals_brackets + intervals_regions
    editable = [(iv.start, iv.end) for iv in all_intervals]
    editable = sorted(list(set(editable)))

    if not editable:
        rprint("[bold yellow]Warning:[/bold yellow] No editable regions found! Mutations will not occur.")
    
    w_parts = [float(x) for x in weights.split(",")]
    
    # 2. Generate Walks (Pre-calculation)
    rprint(f"[bold]Generating {N} walks (M={M} steps)...[/bold]")
    walks = mutate.self_avoiding_walks(X_clean, editable, N, M, tuple(w_parts), seed)
    
    # 3. Single Pass Scoring
    all_seqs = [X_clean] + [w.seq for w in walks]
    
    rprint(f"[bold blue]Sending bundled batch of {len(all_seqs)} sequences to OpenSpliceAI...[/bold blue]")
    t_batch_start = time.time()
    
    all_results = openspliceai.score_batch(all_seqs, left_flank, right_flank, model_path)
    
    batch_duration = time.time() - t_batch_start
    rprint(f"  > Single-pass scoring finished in {batch_duration:.2f}s")

    # Separate results
    res0 = all_results[0]      # The baseline X
    batch_results = all_results[1:] # The mutated Ys

    # 4. Identify Quartet (using res0)
    quartet = mutate.pick_quartet(
        res0.donors, res0.acceptors, 
        d_thr=donor_thresh, a_thr=acceptor_thresh,
        step=0.05, floor=0.01, min_sep=3
    )
    
    if not quartet:
        rprint("[bold red]No valid quartet found in baseline sequence![/bold red]")
        raise typer.Exit(code=1)
        
    rprint(f"Quartet found:\n  dA: {quartet.exonA_donor.pos}\n  aP: {quartet.psex_acceptor.pos}\n  dP: {quartet.psex_donor.pos}\n  aB: {quartet.exonB_acceptor.pos}")

    # 5. Write Output & Collect FASTA Records
    rprint(f"Writing results to {out_file}...")
    
    fasta_records = [] # To store (header, sequence, sort_metric)

    with open(out_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        
        header = [
            "walk_idx", "mutations", "seq_len", "Y_seq",
            "pX_exonA_donor", "pX_psex_acceptor", "pX_psex_donor", "pX_exonB_acceptor",
            "pY_exonA_donor", "pY_psex_acceptor", "pY_psex_donor", "pY_exonB_acceptor",
            "d_exonA_donor", "d_psex_acceptor", "d_psex_donor", "d_exonB_acceptor",
        ]
        writer.writerow(header)
        
        pX_dA = quartet.exonA_donor.prob
        pX_aP = quartet.psex_acceptor.prob
        pX_dP = quartet.psex_donor.prob
        pX_aB = quartet.exonB_acceptor.prob
        
        for i, (w, resY) in enumerate(zip(walks, batch_results)):
            
            def get_prob(sites, original_pos):
                matches = [s.prob for s in sites if abs(s.pos - original_pos) <= 5]
                return max(matches) if matches else 0.0

            pY_dA = get_prob(resY.donors, quartet.exonA_donor.pos)
            pY_aP = get_prob(resY.acceptors, quartet.psex_acceptor.pos)
            pY_dP = get_prob(resY.donors, quartet.psex_donor.pos)
            pY_aB = get_prob(resY.acceptors, quartet.exonB_acceptor.pos)
            
            # Delta for sorting
            d_psex_donor = pY_dP - pX_dP

            # Collect for FASTA
            fasta_header = f"{i}_{w.mutations}"
            fasta_records.append({
                "header": fasta_header,
                "seq": w.seq,
                "sort_val": d_psex_donor
            })

            writer.writerow([
                i, w.mutations, len(w.seq), w.seq,
                f"{pX_dA:.4f}", f"{pX_aP:.4f}", f"{pX_dP:.4f}", f"{pX_aB:.4f}",
                f"{pY_dA:.4f}", f"{pY_aP:.4f}", f"{pY_dP:.4f}", f"{pY_aB:.4f}",
                f"{pY_dA - pX_dA:.4f}",
                f"{pY_aP - pX_aP:.4f}",
                f"{d_psex_donor:.4f}", # d_psex_donor
                f"{pY_aB - pX_aB:.4f}",
            ])
            
    # 6. Write Sorted FASTA (if requested)
    if fasta_out:
        rprint(f"Sorting and writing FASTA to {fasta_out}...")
        # Sort descending based on d_psex_donor
        fasta_records.sort(key=lambda x: x["sort_val"], reverse=True)
        
        with open(fasta_out, "w") as f:
            for rec in fasta_records:
                # Note: mutations strings may contain spaces, which some tools don't like in FASTA headers.
                # But we will adhere strictly to the requested format: walk_idx + "_" + mutations
                f.write(f">{rec['header']}\n{rec['seq']}\n")
            
    total_time = time.time() - t_start_main
    rprint(f"[bold green]Done! (Total time: {total_time:.2f}s)[/bold green]")

if __name__ == "__main__":
    app()