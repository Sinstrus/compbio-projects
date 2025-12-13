# level2/council_saturation.py
import typer
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from rich.console import Console

app = typer.Typer(add_completion=False)
console = Console()

@app.command()
def main(
    input_x: str = typer.Option(..., "--input-x"),
    models_dir: Path = typer.Option(..., "--models-dir"),
    out_file: Path = typer.Option("council_saturation_results.tsv", "--out"),
    
    left_flank: str = typer.Option("", "--left-flank"),
    right_flank: str = typer.Option("", "--right-flank"),
    gain_delta_thresh: float = typer.Option(0.10, "--gain-delta"),
    donor_thresh: float = 0.5,
    acceptor_thresh: float = 0.5,
    
    child_batch_size: int = typer.Option(1000, "--child-batch-size"),
):
    """
    Run saturation mutagenesis on ALL models (recursive) SEQUENTIALLY and aggregate.
    Streams output to terminal so you can see progress.
    """
    model_files = sorted(list(models_dir.rglob("*.pt")))
    if not model_files:
        console.print(f"[bold red]No .pt files found in {models_dir} (recursive)[/bold red]")
        raise typer.Exit(1)
        
    console.print(f"[bold green]Found {len(model_files)} models. Starting Sequential Council...[/bold green]")
    
    successful_files = []
    
    for i, model_path in enumerate(model_files, 1):
        model_label = model_path.stem
        temp_out = f"temp_sat_{model_label}.tsv"
        temp_quartet = f"temp_quartet_{model_label}.tsv"
        
        console.print(f"\n[bold blue]({i}/{len(model_files)}) Running Model: {model_label}[/bold blue]")
        
        cmd = [
            "python", "-m", "level2.saturation",
            "--input-x", input_x,
            "--model-path", str(model_path),
            "--left-flank", left_flank,
            "--right-flank", right_flank,
            "--out", temp_out,
            "--quartet-out", temp_quartet,
            "--gain-delta", str(gain_delta_thresh),
            "--donor-thresh", str(donor_thresh),
            "--acceptor-thresh", str(acceptor_thresh),
            "--batch-size", str(child_batch_size)
        ]
        
        # Stream Output
        try:
            with subprocess.Popen(
                cmd, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.STDOUT, 
                text=True, 
                bufsize=1, 
                universal_newlines=True
            ) as p:
                for line in p.stdout:
                    print(line, end="") 
            
            if p.returncode == 0:
                successful_files.append(temp_out)
                console.print(f"  [green]Success![/green] Output: {temp_out}")
            else:
                console.print(f"  [red]Failed![/red] (See output above)")
        
        except Exception as e:
            console.print(f"[red]Execution error: {e}[/red]")

    if not successful_files:
        console.print("[bold red]All runs failed![/bold red]")
        raise typer.Exit(1)

    # Aggregate results (unchanged logic)
    console.print(f"\n[bold blue]Aggregating {len(successful_files)} runs...[/bold blue]")
    
    dfs = [pd.read_csv(f, sep="\t") for f in successful_files]
    if not dfs: return
    
    base_df = dfs[0].copy()
    prob_cols = ["p_dA", "p_aP", "p_dP", "p_aB", "Max_New_Prob"]
    
    min_len = min(len(d) for d in dfs)
    dfs = [d.iloc[:min_len] for d in dfs]
    base_df = base_df.iloc[:min_len]

    data_stack = np.array([df[prob_cols].values for df in dfs])
    means = np.mean(data_stack, axis=0)
    stds = np.std(data_stack, axis=0)
    
    for i, col in enumerate(prob_cols):
        base_df[col] = means[:, i]
        base_df[f"mean_{col}"] = means[:, i]
        base_df[f"std_{col}"] = stds[:, i]
        
        for j, f in enumerate(successful_files):
            label = Path(f).stem.replace("temp_sat_model_", "") 
            base_df[f"{col}_{label}"] = data_stack[j, :, i]

    base_df.to_csv(out_file, sep="\t", index=False)
    
    for f in successful_files:
        Path(f).unlink()
        q_file = f.replace("temp_sat_", "temp_quartet_")
        if Path(q_file).exists(): Path(q_file).unlink()

    console.print(f"[bold green]Council complete! Results: {out_file}[/bold green]")

if __name__ == "__main__":
    app()