# level2/find_quartet.py
import typer
import csv
from pathlib import Path
from rich.console import Console
from rich import print as rprint

from . import regions, mutate
from .model_backends import openspliceai

app = typer.Typer(add_completion=False)

@app.command()
def main(
    input_x: str = typer.Option(..., "--input-x"),
    model_path: Path = typer.Option(..., "--model-path"),
    out_file: Path = typer.Option("quartet_positions.tsv", "--out"),
    left_flank: str = "",
    right_flank: str = "",
    donor_thresh: float = 0.5,
    acceptor_thresh: float = 0.5,
):
    """
    Identify the target quartet positions and save to TSV.
    """
    X_clean, _ = regions.parse_brackets(input_x)
    
    rprint(f"[bold green]Finding Quartet on X (len={len(X_clean)})...[/bold green]")
    
    # Single Baseline Score
    res0 = openspliceai.score_batch([X_clean], left_flank, right_flank, model_path)[0]
    quartet = mutate.pick_quartet(
        res0.donors, res0.acceptors, 
        donor_thresh, acceptor_thresh, 
        0.05, 0.01, 3
    )
    
    if not quartet:
        rprint("[bold red]No valid quartet found![/bold red]")
        raise typer.Exit(1)

    rprint(f"Quartet Found:\n  dA={quartet.exonA_donor.pos}\n  aP={quartet.psex_acceptor.pos}\n  dP={quartet.psex_donor.pos}\n  aB={quartet.exonB_acceptor.pos}")

    # Save to TSV
    with open(out_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Site", "Position"])
        writer.writerow(["dA", quartet.exonA_donor.pos])
        writer.writerow(["aP", quartet.psex_acceptor.pos])
        writer.writerow(["dP", quartet.psex_donor.pos])
        writer.writerow(["aB", quartet.exonB_acceptor.pos])
    
    rprint(f"[green]Coordinates saved to {out_file}[/green]")

if __name__ == "__main__":
    app()