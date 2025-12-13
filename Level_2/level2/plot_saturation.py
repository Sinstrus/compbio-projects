# level2/plot_saturation.py
import typer
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def main(
    tsv_file: Path = typer.Argument(..., help="Path to saturation_mutation_results.tsv"),
    target_site: str = typer.Option("dP", help="Site to plot: dA, aP, dP, aB"),
    quartet_file: Path = typer.Option(None, "--quartet", help="Path to quartet_positions.tsv (optional)"),
):
    df = pd.read_csv(tsv_file, sep="\t")
    
    col_map = {"dA": "p_dA", "aP": "p_aP", "dP": "p_dP", "aB": "p_aB"}
    y_col = col_map.get(target_site)
    
    if not y_col or y_col not in df.columns:
        print(f"Error: Target {target_site} not found.")
        return

    # Read Quartet Positions if available
    q_pos = {}
    if quartet_file and quartet_file.exists():
        q_df = pd.read_csv(quartet_file, sep="\t")
        for _, row in q_df.iterrows():
            q_pos[row["Site"]] = int(row["Position"])
    
    # Extract Baseline Value
    baseline_row = df[df["Position"] == -1]
    baseline_val = baseline_row.iloc[0][y_col] if not baseline_row.empty else 0.0
    
    # Plot data only (ignore Position -1 for scatter)
    plot_df = df[df["Position"] >= 0]

    plt.figure(figsize=(14, 6))
    
    # 1. Base Scatter (Gray)
    plt.scatter(plot_df["Position"], plot_df[y_col], c="lightgray", s=15, alpha=0.5, label="Variant")
    
    # 2. Motif Highlights
    acc_df = plot_df[plot_df["New_Acceptor_Motif"] == True]
    plt.scatter(acc_df["Position"], acc_df[y_col], c="cornflowerblue", s=20, alpha=0.6, label="New YYYNYAG")
    
    gtr_df = plot_df[plot_df["New_GTR_Motif"] == True]
    plt.scatter(gtr_df["Position"], gtr_df[y_col], c="lightcoral", s=20, alpha=0.6, label="New GTR")
    
    # 3. Vertical Lines (Quartet Members)
    # Logic: Exon sites (dA, aB) = Black; Pseudo sites (aP, dP) = Purple
    if q_pos:
        for site, pos in q_pos.items():
            color = "purple" if site in ["aP", "dP"] else "black"
            # Draw line
            plt.axvline(x=pos, color=color, linestyle="-", linewidth=1.5, alpha=0.4)
            # Add Label at top
            plt.text(pos, 1.02, site, ha='center', va='bottom', fontsize=10, color=color, fontweight='bold')

    # 4. Baseline Horizontal
    plt.axhline(y=baseline_val, color='black', linestyle='--', alpha=0.3, label=f"Baseline ({baseline_val:.2f})")
    
    plt.title(f"Saturation Mutagenesis: Effect on {target_site}", fontsize=14)
    plt.xlabel("Nucleotide Position Index", fontsize=12)
    plt.ylabel(f"Probability ({target_site})", fontsize=12)
    plt.ylim(-0.05, 1.1) # Little extra room at top for labels
    plt.legend(loc='upper right')
    
    out_png = tsv_file.with_name(f"saturation_{target_site}.png")
    plt.savefig(out_png, dpi=150)
    print(f"Plot saved to {out_png}")

if __name__ == "__main__":
    typer.run(main)