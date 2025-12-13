# level2/council_merge.py
import typer
import pandas as pd
import numpy as np
from pathlib import Path
from rich.console import Console

app = typer.Typer(add_completion=False)
console = Console()

@app.command()
def main(
    tsv_dir: Path = typer.Option(..., "--tsv-dir", help="Directory containing saturation_*.tsv files"),
    out_file: Path = typer.Option("council_merged_results.tsv", "--out"),
    pattern: str = typer.Option("saturation_*.tsv", "--pattern", help="Glob pattern to match TSV files"),
):
    """
    Merge existing saturation TSV files into a single Council result (Mean +/- Std).
    """
    # 1. Find Files
    tsv_files = sorted(list(tsv_dir.glob(pattern)))
    if not tsv_files:
        console.print(f"[bold red]No files matching {pattern} found in {tsv_dir}[/bold red]")
        raise typer.Exit(1)
        
    console.print(f"[bold green]Found {len(tsv_files)} result files. Merging...[/bold green]")
    
    # 2. Load Data
    dfs = []
    valid_files = []
    for f in tsv_files:
        try:
            df = pd.read_csv(f, sep="\t")
            # Basic check to ensure it's a saturation file
            if "p_dA" not in df.columns:
                console.print(f"[yellow]Skipping {f.name}: Missing 'p_dA' column[/yellow]")
                continue
            dfs.append(df)
            valid_files.append(f)
        except Exception as e:
            console.print(f"[yellow]Skipping {f.name}: Read error ({e})[/yellow]")

    if not dfs:
        console.print("[bold red]No valid dataframes loaded![/bold red]")
        raise typer.Exit(1)

    # 3. Check Alignment
    # We assume row counts match. If not, we truncate to the shortest one (risky but functional)
    min_len = min(len(df) for df in dfs)
    dfs = [df.iloc[:min_len] for df in dfs]
    
    base_df = dfs[0].copy()
    prob_cols = ["p_dA", "p_aP", "p_dP", "p_aB", "Max_New_Prob"]
    
    # Handle missing columns gracefully
    final_prob_cols = [c for c in prob_cols if c in base_df.columns]
    
    # 4. Calculate Stats
    data_stack = np.array([df[final_prob_cols].values for df in dfs])
    means = np.mean(data_stack, axis=0)
    stds = np.std(data_stack, axis=0)
    
    for i, col in enumerate(final_prob_cols):
        base_df[col] = means[:, i] # Overwrite main col with mean
        base_df[f"mean_{col}"] = means[:, i]
        base_df[f"std_{col}"] = stds[:, i]
        
        for j, f in enumerate(valid_files):
            # Use regex or simple split to get model name?
            # Try to extract 'rs10' or similar if possible, else use index
            label = f.stem.split("_")[-1] # saturation_model_10000nt_rs10 -> rs10
            base_df[f"{col}_{label}"] = data_stack[j, :, i]

    # 5. Save
    base_df.to_csv(out_file, sep="\t", index=False)
    console.print(f"[bold green]Merged {len(valid_files)} files into {out_file}[/bold green]")

if __name__ == "__main__":
    app()