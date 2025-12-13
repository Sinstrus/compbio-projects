# level2/reconstruct_sequences.py
import typer
import pandas as pd
from pathlib import Path
from rich.console import Console
from . import regions, mutate

app = typer.Typer(add_completion=False)
console = Console()

@app.command()
def main(
    input_x: str = typer.Option(..., "--input-x", help="Original SeqX used for saturation"),
    tsv_file: Path = typer.Argument(..., help="Path to saturation_mutation_results.tsv"),
    out_file: Path = typer.Option(None, "--out", help="Output file (default: saturation_with_seqs.tsv)"),
):
    """
    Reconstruct full DNA sequences from saturation mutagenesis log.
    """
    # Parse X (ignore brackets for seq reconstruction)
    X_clean, _ = regions.parse_brackets(input_x)
    
    console.print(f"Loaded SeqX (len={len(X_clean)})")
    
    df = pd.read_csv(tsv_file, sep="\t")
    
    def get_seq(row):
        mut_str = row["Mutation"]
        if mut_str == "Original":
            return X_clean
        return mutate.apply_mutation_step(X_clean, mut_str)

    console.print(f"Reconstructing sequences for {len(df)} variants...")
    df["Sequence"] = df.apply(get_seq, axis=1)
    
    # Move Sequence column to after Mutation for readability
    cols = list(df.columns)
    cols.insert(2, cols.pop(cols.index("Sequence")))
    df = df[cols]
    
    if not out_file:
        out_file = tsv_file.with_name("saturation_with_seqs.tsv")
        
    df.to_csv(out_file, sep="\t", index=False)
    console.print(f"[green]Done! Wrote {len(df)} sequences to {out_file}[/green]")

if __name__ == "__main__":
    app()