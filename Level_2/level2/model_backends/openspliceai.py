# level2/model_backends/openspliceai.py
from __future__ import annotations
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict
from rich import print as rprint
import sys

@dataclass
class Site:
    pos: int; prob: float; kind: str

@dataclass
class OsaiResult:
    donors: List[Site]; acceptors: List[Site]

def score_batch(X_sequences, left_flank, right_flank, model_path, flanking_size=10000):
    L_left = len(left_flank)
    N_input = len(X_sequences)
    X_lengths = [len(x) for x in X_sequences]
    results = {i: OsaiResult([], []) for i in range(N_input)}

    with tempfile.TemporaryDirectory(prefix="osai_batch_") as tmpd:
        tmpd_path = Path(tmpd)
        fa = tmpd_path / "input.fa"
        out_dir = tmpd_path / "out"
        
        with open(fa, "w") as f:
            for i, x_seq in enumerate(X_sequences):
                f.write(f">seq_{i}\n{left_flank + x_seq + right_flank}\n")
        
        cmd = ["openspliceai", "predict", "--input-sequence", str(fa), 
               "--flanking-size", str(flanking_size), "--output-dir", str(out_dir), 
               "--model", str(model_path)]
        
        # CHANGED: capture_output=True SILENCES the tool's progress bars
        try: 
            subprocess.run(cmd, check=True, capture_output=True, text=True) 
        except subprocess.CalledProcessError as e:
            # We manually print stderr only if it fails
            rprint(f"[bold red]OpenSpliceAI Error:[/bold red]\n{e.stderr}")
            return [results[i] for i in range(N_input)]
        
        res_dirs = list(out_dir.glob("SpliceAI_*"))
        if res_dirs:
            td = max(res_dirs, key=lambda p: p.stat().st_mtime)
            _parse(td / "donor_predictions.bed", "Donor", results, L_left, X_lengths)
            _parse(td / "acceptor_predictions.bed", "Acceptor", results, L_left, X_lengths)

    return [results[i] for i in range(N_input)]

def _parse(path, kind, res_dict, offset, x_lens):
    if not path.exists(): return
    with path.open() as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 5 or line.startswith("#"): continue
            name = parts[0].split(":")[0]
            if not name.startswith("seq_"): continue
            try:
                idx = int(name.split("_")[1])
                pos = (int(parts[1]) + int(parts[2])) // 2 - offset
                if 0 <= pos < x_lens[idx]:
                    site = Site(pos, float(parts[4]), kind)
                    if kind == "Donor": res_dict[idx].donors.append(site)
                    else: res_dict[idx].acceptors.append(site)
            except: continue
    for r in res_dict.values():
        r.donors.sort(key=lambda s: s.pos)
        r.acceptors.sort(key=lambda s: s.pos)