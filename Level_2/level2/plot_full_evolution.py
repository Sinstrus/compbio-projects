# level2/plot_full_evolution.py
import typer
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def main(
    tsv_file: Path = typer.Argument(..., help="Path to evolution_all_variants.tsv"),
    target_site: str = typer.Option("dP", help="Site to plot: dA (exonA-donor), aP (psex-acceptor), dP (psex-donor), aB (exonB-acceptor)"),
):
    """
    Visualize the full evolutionary search space.
    """
    df = pd.read_csv(tsv_file, sep="\t")
    
    # Map shortcodes to CSV column names
    target_map = {
        "dA": "p_dA",
        "aP": "p_aP",
        "dP": "p_dP",
        "aB": "p_aB"
    }
    
    if target_site not in target_map:
        print(f"Error: Unknown target site '{target_site}'. Choose from: dA, aP, dP, aB")
        return

    col_name = target_map[target_site]
    if col_name not in df.columns:
        print(f"Error: Column {col_name} not found in TSV.")
        return

    # Create plot
    plt.figure(figsize=(12, 7))
    
    # Add slight X-axis jitter
    jitter_strength = 0.15
    df['plot_x'] = df['Round'] + np.random.uniform(-jitter_strength, jitter_strength, size=len(df))
    
    # 1. Plot Rejects (Gray)
    rejects = df[df["Status"].isin(["AI_Reject", "Pheno_Reject"])]
    plt.scatter(
        rejects['plot_x'], rejects[col_name], 
        color='gray', alpha=0.2, s=20, label='Rejected', zorder=1
    )
    
    # 2. Plot Passed (Orange)
    passed = df[(df["Status"] == "Passed") & (df["Is_Winner"] == False)]
    plt.scatter(
        passed['plot_x'], passed[col_name], 
        color='orange', alpha=0.5, s=40, label='Passed (Non-Winner)', zorder=2
    )
    
    # 3. Plot Winners (Red) including Round 0
    winners = df[df["Is_Winner"] == True].sort_values("Round")
    
    # Draw line connecting winners
    plt.plot(winners['Round'], winners[col_name], color='red', linewidth=2, alpha=0.8, zorder=3)
    
    # Draw winner dots
    plt.scatter(
        winners['Round'], winners[col_name], 
        color='red', edgecolors='black', alpha=1.0, s=100, marker='*', label='Winner', zorder=4
    )
    
    # Styling
    site_labels = {
        "dA": "Exon A Donor (dA)",
        "aP": "Pseudo-Exon Acceptor (aP)",
        "dP": "Pseudo-Exon Donor (dP)",
        "aB": "Exon B Acceptor (aB)"
    }
    
    plt.title(f"Evolutionary Trajectory: {site_labels.get(target_site, target_site)}", fontsize=14)
    plt.xlabel("Round", fontsize=12)
    plt.ylabel("Probability", fontsize=12)
    
    # Ensure ticks are integers (0, 1, 2...)
    plt.xticks(sorted(df['Round'].unique()))
    plt.ylim(-0.05, 1.05) # Fix y-axis to probability range
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best')
    
    # Save
    outfile = tsv_file.with_name(f"evolution_swarm_{target_site}.png")
    plt.savefig(outfile, dpi=150)
    print(f"Plot saved to {outfile}")

if __name__ == "__main__":
    typer.run(main)