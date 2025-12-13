# level2/run_council.py
import typer
import subprocess
import glob
import pandas as pd
from pathlib import Path
from rich.console import Console

app = typer.Typer(add_completion=False)
console = Console()

@app.command()
def main(
    input_x: str = typer.Option(..., "--input-x"),
    models_dir: Path = typer.Option(..., "--models-dir", help="Directory containing model_10000nt_rs*.pt files"),
    left_flank: str = typer.Option("", "--left-flank"),
    right_flank: str = typer.Option("", "--right-flank"),
    output_prefix: str = typer.Option("saturation", "--output-prefix"),
):
    model_files = sorted(list(models_dir.glob("*.pt")))
    if not model_files:
        console.print(f"[bold red]No .pt files found in {models_dir}[/bold red]")
        raise typer.Exit(1)
        
    console.print(f"[bold green]Found {len(model_files)} models. Starting Council Run...[/bold green]")
    
    successful_runs = []
    failed_runs = []

    for model_path in model_files:
        model_name = model_path.stem
        console.print(f"\n[bold blue]Running {model_name}...[/bold blue]")
        
        out_tsv = f"{output_prefix}_{model_name}.tsv"
        quartet_out = f"quartet_{model_name}.tsv"
        baseline_out = f"baseline_{model_name}.tsv"
        
        cmd = [
            "python", "-m", "level2.saturation",
            "--input-x", input_x,
            "--model-path", str(model_path),
            "--left-flank", left_flank,
            "--right-flank", right_flank,
            "--out", out_tsv,
            "--quartet-out", quartet_out,
            "--baseline-out", baseline_out,
            "--gain-delta", "0.10",
            "--donor-thresh", "0.10", # Lower threshold for robust detection
            "--acceptor-thresh", "0.10"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            console.print(f"  [green]Success![/green] Output: {out_tsv}")
            successful_runs.append(out_tsv)
        else:
            console.print(f"  [red]Failed![/red]")
            # Extract error reason
            for line in result.stdout.splitlines():
                if "Error" in line or "No valid quartet" in line:
                    console.print(f"    Reason: {line.strip()}")
            failed_runs.append(model_name)

    console.print("\n[bold]Council Run Complete.[/bold]")
    console.print(f"Successful: {len(successful_runs)}")
    console.print(f"Failed: {len(failed_runs)} ({', '.join(failed_runs)})")

if __name__ == "__main__":
    app()