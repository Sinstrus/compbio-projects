# level2/plot_council_saturation.py
import typer
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

app = typer.Typer(add_completion=False)

@app.command()
def main(
    tsv_file: Path = typer.Argument(..., help="Path to council_saturation_results.tsv"),
    target_site: str = typer.Option("dP", help="Site to plot: dA, aP, dP, aB"),
    quartet_file: Path = typer.Option("quartet_positions.tsv", "--quartet", help="Path to quartet coordinates"),
    out_file: Path = typer.Option(None, "--out", help="Custom output filename (PNG)"),
    delta_thresh: float = typer.Option(0.05, "--delta-thresh"),
    cv_thresh: float = typer.Option(0.20, "--cv-thresh"),
):
    df = pd.read_csv(tsv_file, sep="\t")
    mean_col = f"mean_p_{target_site}"
    std_col = f"std_p_{target_site}"
    if mean_col not in df.columns: return

    baseline_row = df[df["Position"] == -1]
    p_orig = baseline_row.iloc[0][mean_col]
    plot_df = df[df["Position"] >= 0].copy()
    plot_df = plot_df[(plot_df["New_Acceptor_Motif"] == False) & (plot_df["New_GTR_Motif"] == False)]

    plot_df["delta"] = plot_df[mean_col] - p_orig
    plot_df["cv"] = plot_df.apply(lambda row: row[std_col] / row[mean_col] if row[mean_col] > 0.001 else 0.0, axis=1)

    def get_category(row):
        d = row["delta"]; cv = row["cv"]
        if d >= delta_thresh and cv <= cv_thresh: return "robust_up"
        elif d <= -delta_thresh and cv <= cv_thresh: return "robust_down"
        return "uncertain"

    plot_df["category"] = plot_df.apply(get_category, axis=1)

    plt.figure(figsize=(14, 7))
    styles = {
        "uncertain":   {"color": "lightgray", "alpha": 0.3, "label": "Uncertain / Weak"},
        "robust_up":   {"color": "firebrick", "alpha": 0.8, "label": "Robust Increase"},
        "robust_down": {"color": "dodgerblue", "alpha": 0.8, "label": "Robust Decrease"}
    }

    for cat, style in styles.items():
        subset = plot_df[plot_df["category"] == cat]
        if not subset.empty:
            plt.errorbar(subset["Position"], subset[mean_col], yerr=subset[std_col], fmt='none', ecolor=style["color"], alpha=style["alpha"]*0.4, zorder=1)
            plt.scatter(subset["Position"], subset[mean_col], c=style["color"], s=25, alpha=style["alpha"], label=style["label"], zorder=2)

    if quartet_file and quartet_file.exists():
        q_df = pd.read_csv(quartet_file, sep="\t")
        for _, row in q_df.iterrows():
            pos = int(row["Position"])
            color = "purple" if row["Site"] in ["aP", "dP"] else "black"
            plt.scatter(pos, -0.03, marker="^", color=color, s=100, zorder=5, clip_on=False)
            plt.text(pos, -0.08, row["Site"], ha='center', va='top', fontsize=10, color=color, fontweight='bold')

    plt.axhline(y=p_orig, color='black', linestyle='--', alpha=0.5, label=f"Baseline ({p_orig:.2f})")
    plt.title(f"Council Consensus: {target_site}")
    plt.xlabel("Position"); plt.ylabel("Probability")
    plt.ylim(-0.1, 1.1)
    plt.legend()
    
    final_out = out_file if out_file else tsv_file.with_name(f"council_saturation_{target_site}_cv.png")
    plt.savefig(final_out, dpi=150)
    print(f"Plot saved to {final_out}")

if __name__ == "__main__":
    app()