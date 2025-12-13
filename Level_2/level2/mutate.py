# level2/mutate.py
from __future__ import annotations
import random
from dataclasses import dataclass
from typing import List, Tuple, Optional
from .model_backends.openspliceai import Site

DNA = ["A", "C", "G", "T"]

@dataclass
class Quartet:
    exonA_donor: Site
    psex_acceptor: Site
    psex_donor: Site
    exonB_acceptor: Site

def pick_quartet(
    donors: List[Site],
    acceptors: List[Site],
    d_thr: float, a_thr: float,
    step: float, floor: float, min_sep: int
) -> Optional[Quartet]:
    curr_d = d_thr
    curr_a = a_thr
    while True:
        valid_d = [d for d in donors if d.prob >= curr_d]
        valid_a = [a for a in acceptors if a.prob >= curr_a]
        best_q = None
        best_score = -1.0
        for dA in valid_d:
            for aB in valid_a:
                if aB.pos <= dA.pos + 3*min_sep: continue
                mid_a = [a for a in valid_a if dA.pos + min_sep < a.pos < aB.pos - min_sep]
                mid_d = [d for d in valid_d if dA.pos + min_sep < d.pos < aB.pos - min_sep]
                for aP in mid_a:
                    for dP in mid_d:
                        if aP.pos < dP.pos - min_sep:
                            sc = min(dA.prob, aP.prob, dP.prob, aB.prob)
                            if sc > best_score:
                                best_score = sc
                                best_q = Quartet(dA, aP, dP, aB)
        if best_q: return best_q
        if curr_d <= floor and curr_a <= floor: break
        curr_d = max(floor, curr_d - step)
        curr_a = max(floor, curr_a - step)
    return None

@dataclass
class WalkResult:
    seq: str
    steps: int
    mutations: str

def _mutate_once(seq: str, intervals: List[Tuple[int, int]], rng: random.Random, weights: Tuple[float, float, float]) -> Tuple[str, str]:
    w_sub, w_del, w_ins = weights
    mode = rng.choices(["sub", "del", "ins"], weights=[w_sub, w_del, w_ins])[0]
    allowed_indices = []
    for (start, end) in intervals:
        allowed_indices.extend(range(start, end))
    if not allowed_indices: return seq, ""
    pos = rng.choice(allowed_indices)
    
    if mode == "sub":
        old = seq[pos]
        new_b = rng.choice([b for b in DNA if b != old])
        new_seq = seq[:pos] + new_b + seq[pos+1:]
        return new_seq, f"sub:{pos} {old}>{new_b}"
    elif mode == "del":
        if len(seq) < 2: return seq, ""
        old = seq[pos]
        new_seq = seq[:pos] + seq[pos+1:]
        return new_seq, f"del:{pos} {old}"
    elif mode == "ins":
        new_b = rng.choice(DNA)
        new_seq = seq[:pos] + new_b + seq[pos:]
        return new_seq, f"ins:{pos} {new_b}"
    return seq, ""

def self_avoiding_walks(X: str, editable_intervals: List[Tuple[int,int]], N: int, M: int, weights: Tuple[float,float,float], seed: int) -> List[WalkResult]:
    rng = random.Random(seed)
    results = []
    for _ in range(N):
        curr_seq = X
        curr_muts = []
        for _ in range(M):
            next_seq, desc = _mutate_once(curr_seq, editable_intervals, rng, weights)
            if desc:
                curr_seq = next_seq
                curr_muts.append(desc)
        results.append(WalkResult(curr_seq, len(curr_muts), ";".join(curr_muts)))
    return results

def apply_mutation_step(seq: str, desc: str) -> str:
    if not desc: return seq
    parts = desc.split()
    op_raw = parts[0]
    op, pos_str = op_raw.split(":")
    pos = int(pos_str)
    if pos >= len(seq) and op != "ins": return seq
    if op == "sub":
        _, new_b = parts[1].split(">")
        return seq[:pos] + new_b + seq[pos+1:]
    elif op == "del":
        return seq[:pos] + seq[pos+1:]
    elif op == "ins":
        base = parts[1]
        return seq[:pos] + base + seq[pos:]
    return seq

# --- Motif Checking ---
def is_gtr(seq: str, pos: int) -> bool:
    if pos < 0 or pos + 3 > len(seq): return False
    return seq[pos] == "G" and seq[pos+1] == "T" and seq[pos+2] in "AG"

def is_acceptor_motif(seq: str, pos: int) -> bool:
    if pos < 0 or pos + 7 > len(seq): return False
    for i in [0, 1, 2, 4]:
        if seq[pos+i] not in "CT": return False
    if seq[pos+5] != "A": return False
    if seq[pos+6] != "G": return False
    return True

def check_new_motifs(parent_seq: str, mutations_str: str) -> Tuple[bool, bool]:
    if not mutations_str: return (False, False)
    indices = list(range(len(parent_seq)))
    temp_indices = indices[:]
    steps = mutations_str.split(";")
    for step in steps:
        if not step.strip(): continue
        parts = step.strip().split()
        op_part = parts[0]
        op, pos_s = op_part.split(":")
        pos = int(pos_s)
        if op == "ins": temp_indices.insert(pos, -1)
        elif op == "del":
            if pos < len(temp_indices): temp_indices.pop(pos)

    child_seq = apply_mutation_step(parent_seq, steps[0])
    for s in steps[1:]:
        child_seq = apply_mutation_step(child_seq, s)

    new_gtr = False
    new_acc = False

    for i in range(len(child_seq) - 2):
        if is_gtr(child_seq, i):
            origin = temp_indices[i : i+3]
            if -1 in origin: new_gtr = True 
            elif not (origin[1] == origin[0] + 1 and origin[2] == origin[1] + 1): new_gtr = True
            else:
                parent_idx = origin[0]
                if not is_gtr(parent_seq, parent_idx): new_gtr = True
    
    for i in range(len(child_seq) - 6):
        if is_acceptor_motif(child_seq, i):
            origin = temp_indices[i : i+7]
            if -1 in origin: new_acc = True
            else:
                is_continuous = True
                for k in range(6):
                    if origin[k+1] != origin[k] + 1:
                        is_continuous = False; break
                if not is_continuous: new_acc = True
                else:
                    parent_idx = origin[0]
                    if not is_acceptor_motif(parent_seq, parent_idx): new_acc = True

    return new_gtr, new_acc