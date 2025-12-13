# level2/plot_evolution.py
import typer
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def main(tsv_file: Path):
    df = pd.read_csv(tsv_file, sep="\t")
    
    rounds = df["Round"].values
    # Plot Deltas
    d_aP = df["delta_aP"].values
    d_dP = df["delta_dP"].values
    
    # Also plot absolute probabilities to see the "real" values
    p_aP = df["pChild_aP"].values
    p_dP = df["pChild_dP"].values

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 1: The Deltas (What you optimized for)
    ax1.plot(rounds, d_aP, 'b-o', label='Delta aP (Acceptor)')
    ax1.plot(rounds, d_dP, 'r-o', label='Delta dP (Donor)')
    ax1.set_title("Phenotype Divergence (Deltas)")
    ax1.set_xlabel("Round")
    ax1.set_ylabel("Change in Probability")
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Plot 2: Absolute Probabilities
    ax2.plot(rounds, p_aP, 'b--', label='Abs aP')
    ax2.plot(rounds, p_dP, 'r--', label='Abs dP')
    ax2.set_title("Absolute Splice Probabilities")
    ax2.set_xlabel("Round")
    ax2.set_ylabel("Probability")
    ax2.set_ylim(0, 1.0)
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    outfile = tsv_file.with_suffix('.png')
    plt.savefig(outfile)
    print(f"Plot saved to {outfile}")

if __name__ == "__main__":
    typer.run(main)