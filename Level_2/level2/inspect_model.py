# level2/inspect_model.py
import typer
from pathlib import Path
from rich.console import Console
from rich.table import Table
from rich import print as rprint

from . import regions
from .model_backends import openspliceai

app = typer.Typer(add_completion=False)
console = Console()

@app.command()
def main(
    input_x: str = typer.Option(..., "--input-x"),
    model_path: Path = typer.Option(..., "--model-path"),
    left_flank: str = typer.Option("", "--left-flank"),
    right_flank: str = typer.Option("", "--right-flank"),
    min_prob: float = typer.Option(0.01, "--min-prob", help="Minimum probability to display"),
):
    """
    Run a specific model on SeqX and display ALL detected splice sites.
    """
    X_clean, _ = regions.parse_brackets(input_x)
    
    rprint(f"[bold blue]Inspecting Model: {model_path.name}[/bold blue]")
    
    # Run Model
    results = openspliceai.score_batch([X_clean], left_flank, right_flank, model_path)
    res = results[0]
    
    # Display Donors
    table_d = Table(title=f"Donors (p >= {min_prob})")
    table_d.add_column("Position", justify="right")
    table_d.add_column("Probability", justify="right")
    
    for d in res.donors:
        if d.prob >= min_prob:
            table_d.add_row(str(d.pos), f"{d.prob:.4f}")
            
    # Display Acceptors
    table_a = Table(title=f"Acceptors (p >= {min_prob})")
    table_a.add_column("Position", justify="right")
    table_a.add_column("Probability", justify="right")
    
    for a in res.acceptors:
        if a.prob >= min_prob:
            table_a.add_row(str(a.pos), f"{a.prob:.4f}")
            
    console.print(table_d)
    console.print(table_a)

if __name__ == "__main__":
    app()