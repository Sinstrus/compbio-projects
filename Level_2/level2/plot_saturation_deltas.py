# level2/plot_saturation_deltas.py
import typer
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

app = typer.Typer(add_completion=False)

@app.command()
def main(
    tsv_file: Path = typer.Argument(..., help="Path to saturation_mutation_results.tsv"),
    target_site: str = typer.Option("dP", help="Site to plot on Y-axis: dA, aP, dP, aB"),
    quartet_file: Path = typer.Option("quartet_positions.tsv", "--quartet", help="Path to quartet coordinates"),
    out_file: Path = typer.Option(None, "--out", help="Custom output filename (PNG)"),
    
    purple_ap_min: float = typer.Option(0.02, "--purple-ap-min", help="Min d_aP to qualify for Purple"),
    purple_dp_max: float = typer.Option(-0.02, "--purple-dp-max", help="Max d_dP to qualify for Purple"),
):
    """
    Plot Saturation Mutagenesis with Phenotype-based Coloring.
    """
    df = pd.read_csv(tsv_file, sep="\t")
    
    col_map = {"dA": "p_dA", "aP": "p_aP", "dP": "p_dP", "aB": "p_aB"}
    y_col = col_map.get(target_site)
    if not y_col or y_col not in df.columns:
        print(f"Error: Target {target_site} not found.")
        return

    baseline_row = df[df["Position"] == -1]
    if baseline_row.empty:
        print("Error: No baseline row (Position -1) found in TSV.")
        return
    
    p_aP_orig = baseline_row.iloc[0]["p_aP"]
    p_dP_orig = baseline_row.iloc[0]["p_dP"]
    y_baseline_val = baseline_row.iloc[0][y_col]

    plot_df = df[df["Position"] >= 0].copy()
    
    # Filter out structural motifs
    initial_count = len(plot_df)
    plot_df = plot_df[
        (plot_df["New_Acceptor_Motif"] == False) & 
        (plot_df["New_GTR_Motif"] == False)
    ]
    filtered_count = len(plot_df)
    print(f"Filtered out {initial_count - filtered_count} variants due to new motifs.")

    def get_color_category(row):
        d_aP = row["p_aP"] - p_aP_orig
        d_dP = row["p_dP"] - p_dP_orig
        
        if d_aP > purple_ap_min and d_dP < purple_dp_max:
            return "purple"
        elif d_aP > 0.02 and d_dP > 0.02:
            return "red"
        elif d_aP < -0.02 and d_dP < -0.02:
            return "blue"
        else:
            return "gray"

    plot_df["category"] = plot_df.apply(get_color_category, axis=1)

    plt.figure(figsize=(14, 7))
    
    style_map = {
        "gray":   {"c": "lightgray",      "alpha": 0.3, "s": 15, "label": "Neutral/Mixed"},
        "blue":   {"c": "cornflowerblue", "alpha": 0.5, "s": 25, "label": "aP↓ dP↓ (Crash)"},
        "red":    {"c": "lightcoral",     "alpha": 0.5, "s": 25, "label": "aP↑ dP↑ (Boost)"},
        "purple": {"c": "purple",         "alpha": 0.9, "s": 40, "label": f"aP>{purple_ap_min}, dP<{purple_dp_max}"}
    }
    
    for cat in ["gray", "blue", "red", "purple"]:
        subset = plot_df[plot_df["category"] == cat]
        if not subset.empty:
            style = style_map[cat]
            plt.scatter(
                subset["Position"], subset[y_col], 
                c=style["c"], 
                alpha=style["alpha"], 
                s=style["s"], 
                label=style["label"]
            )

    if quartet_file and quartet_file.exists():
        q_df = pd.read_csv(quartet_file, sep="\t")
        for _, row in q_df.iterrows():
            site = row["Site"]
            pos = int(row["Position"])
            color = "purple" if site in ["aP", "dP"] else "black"
            width = 2.0 if site in ["aP", "dP"] else 1.5
            plt.axvline(x=pos, color=color, linestyle="-", linewidth=width, alpha=0.5)
            plt.text(pos, 1.02, site, ha='center', va='bottom', fontsize=10, color=color, fontweight='bold')

    plt.axhline(y=y_baseline_val, color='black', linestyle='--', alpha=0.5, label=f"Baseline {target_site}")

    site_labels = {
        "dA": "Exon A Donor (dA)",
        "aP": "Pseudo-Exon Acceptor (aP)",
        "dP": "Pseudo-Exon Donor (dP)",
        "aB": "Exon B Acceptor (aB)"
    }
    plt.title(f"Mutation Effects on {site_labels.get(target_site, target_site)}\n(Purple: aP > {purple_ap_min} & dP < {purple_dp_max})", fontsize=14)
    plt.xlabel("Nucleotide Position Index", fontsize=12)
    plt.ylabel(f"Probability ({target_site})", fontsize=12)
    plt.ylim(-0.05, 1.1)
    plt.legend(loc='upper right')
    
    final_out = out_file if out_file else tsv_file.with_name(f"saturation_deltas_{target_site}.png")
    plt.savefig(final_out, dpi=150)
    print(f"Plot saved to {final_out}")

if __name__ == "__main__":
    app()