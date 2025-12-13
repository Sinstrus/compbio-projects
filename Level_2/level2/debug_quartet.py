# level2/debug_quartet.py
import typer
import sys
from pathlib import Path
from rich import print as rprint
from . import regions
from .model_backends import openspliceai

app = typer.Typer(add_completion=False)

@app.command()
def main(
    input_x: str = typer.Option(..., "--input-x"),
    model_path: Path = typer.Option(..., "--model-path"),
    left_flank: str = typer.Option("", "--left-flank"),
    right_flank: str = typer.Option("", "--right-flank"),
    donor_thresh: float = 0.5, 
    acceptor_thresh: float = 0.5,
):
    """
    Runs the quartet logic in VERBOSE mode with Pre-Flight Checks.
    """
    rprint(f"[bold blue]Diagnostic Run for Model: {model_path.name}[/bold blue]")

    # --- Pre-Flight Checks ---
    if not model_path.exists():
        rprint(f"[bold red]CRITICAL ERROR: Model file not found![/bold red]")
        rprint(f"  Path checked: {model_path}")
        sys.exit(1)
        
    X_clean, _ = regions.parse_brackets(input_x)
    
    rprint(f"  > Input Stats:")
    rprint(f"    - SeqX Length: {len(X_clean)} bp")
    rprint(f"    - Left Flank:  {len(left_flank)} bp")
    rprint(f"    - Right Flank: {len(right_flank)} bp")
    rprint(f"    - Total Input: {len(left_flank) + len(X_clean) + len(right_flank)} bp")

    if len(X_clean) == 0:
        rprint("[bold red]ERROR: SeqX is empty after parsing![/bold red]")
        sys.exit(1)

    # 1. Run Model
    rprint("\n  > Running OpenSpliceAI...")
    try:
        results = openspliceai.score_batch([X_clean], left_flank, right_flank, model_path)
    except Exception as e:
        rprint(f"[bold red]Exception during execution:[/bold red] {e}")
        sys.exit(1)

    res = results[0]
    donors = res.donors
    acceptors = res.acceptors
    
    rprint(f"  > Raw Results: {len(donors)} donors, {len(acceptors)} acceptors.")

    if len(donors) == 0 and len(acceptors) == 0:
        rprint("[bold red]WARNING: No sites found. This usually indicates a bad model path or malformed sequence.[/bold red]")
        rprint("  Check that the model path in your command matches the working saturation command exactly.")
        sys.exit(0)

    # 2. Replicate pick_quartet with logging
    step = 0.05
    floor = 0.01
    min_sep = 3
    
    curr_d = donor_thresh
    curr_a = acceptor_thresh
    
    iteration = 0
    
    while True:
        iteration += 1
        rprint(f"\n[bold yellow]--- Search Iteration {iteration} (Thresh: {curr_d:.3f}) ---[/bold yellow]")
        
        valid_d = [d for d in donors if d.prob >= curr_d]
        valid_a = [a for a in acceptors if a.prob >= curr_a]
        
        if not valid_d: rprint("    [dim]No donors above threshold.[/dim]")
        if not valid_a: rprint("    [dim]No acceptors above threshold.[/dim]")
        
        found_any_inner = False
        
        for dA in valid_d:
            for aB in valid_a:
                if aB.pos <= dA.pos + 3*min_sep: continue
                
                mid_a = [a for a in valid_a if dA.pos + min_sep < a.pos < aB.pos - min_sep]
                mid_d = [d for d in valid_d if dA.pos + min_sep < d.pos < aB.pos - min_sep]
                
                if len(mid_a) > 0 or len(mid_d) > 0:
                    found_any_inner = True
                    rprint(f"    Checking Outer Pair: dA@{dA.pos}({dA.prob:.3f}) ... aB@{aB.pos}({aB.prob:.3f})")
                    rprint(f"      > Inner Candidates: {len(mid_a)} Acc, {len(mid_d)} Don")
                    
                    for aP in mid_a:
                        for dP in mid_d:
                            if aP.pos < dP.pos - min_sep:
                                rprint(f"      [green]SUCCESS![/green] Found Quartet: {dA.pos} -> {aP.pos} -> {dP.pos} -> {aB.pos}")
                                return
                            # else:
                                # rprint(f"      [dim]Reject inner pair {aP.pos}-{dP.pos} (Topology)[/dim]")
        
        if not found_any_inner and valid_d and valid_a:
             rprint("    [dim]Valid outer sites exist, but no inner candidates found between them.[/dim]")

        if curr_d <= floor and curr_a <= floor:
            rprint("\n[bold red]FAILURE: Reached probability floor without finding a quartet.[/bold red]")
            break
            
        curr_d = max(floor, curr_d - step)
        curr_a = max(floor, curr_a - step)

if __name__ == "__main__":
    app()