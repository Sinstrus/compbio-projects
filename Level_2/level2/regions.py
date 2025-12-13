# level2/regions.py
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple, Iterable

@dataclass(frozen=True)
class Interval:
    """Half-open interval [start, end) in 0-based X coordinates."""
    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 0 or self.end < 0 or self.end < self.start:
            raise ValueError(f"Invalid Interval({self.start}, {self.end})")

    def empty(self) -> bool:
        return self.end <= self.start

    def __len__(self) -> int:
        return self.end - self.start
    
    def intersects(self, other: "Interval") -> bool:
        return not (self.end <= other.start or other.end <= self.start)

def _normalize(seq: str) -> str:
    """Uppercase and validate A/C/G/T plus '[' and ']' characters."""
    s = "".join(seq.split()).upper()
    for ch in s:
        if ch not in "ACGT[]":
            raise ValueError(f"Invalid char '{ch}'. Only A, C, G, T, and [] allowed.")
    return s

def parse_brackets(x_annot: str) -> Tuple[str, List[Interval]]:
    """Returns (clean_sequence, list_of_editable_intervals)."""
    s = _normalize(x_annot)
    out = []
    spans = []
    in_bracket = False
    start_clean = 0

    i = 0
    while i < len(s):
        ch = s[i]
        if ch == "[":
            if in_bracket: raise ValueError("Nested '[' not allowed")
            in_bracket = True
            start_clean = len(out)
        elif ch == "]":
            if not in_bracket: raise ValueError("Unmatched ']'")
            in_bracket = False
            end_clean = len(out)
            if end_clean <= start_clean: raise ValueError("Empty [] block")
            spans.append(Interval(start_clean, end_clean))
        else:
            out.append(ch)
        i += 1
    
    if in_bracket: raise ValueError("Unclosed '['")
    return "".join(out), spans

def _parse_region_pattern(pat: str) -> Tuple[str, List[Interval]]:
    s = _normalize(pat)
    stripped = []
    subspans = []
    i = 0
    running = 0
    while i < len(s):
        ch = s[i]
        if ch == "[":
            j = i + 1
            begin = running
            while j < len(s) and s[j] != "]":
                if s[j] not in "ACGT":
                    raise ValueError(f"Illegal char '{s[j]}' inside [] of pattern")
                stripped.append(s[j])
                running += 1
                j += 1
            if j >= len(s) or s[j] != "]":
                raise ValueError("Unclosed '[' in region pattern")
            end = running
            if end <= begin:
                raise ValueError("Empty [] block in region pattern is not allowed")
            subspans.append(Interval(begin, end))
            i = j + 1
        else:
            if ch not in "ACGT":
                raise ValueError(f"Illegal char '{ch}' in region pattern")
            stripped.append(ch)
            running += 1
            i += 1
    return "".join(stripped), subspans

def _find_all(haystack: str, needle: str) -> List[int]:
    hits = []
    i = 0
    while True:
        j = haystack.find(needle, i)
        if j < 0: break
        hits.append(j)
        i = j + 1 
    return hits

def _merge_intervals(spans: Iterable[Interval]) -> List[Interval]:
    ivs = sorted((iv for iv in spans if not iv.empty()), key=lambda t: (t.start, t.end))
    if not ivs: return []
    out = [ivs[0]]
    for iv in ivs[1:]:
        prev = out[-1]
        if prev.end >= iv.start:
            out[-1] = Interval(prev.start, max(prev.end, iv.end))
        else:
            out.append(iv)
    return out

def apply_region_patterns(X_clean: str, patterns: str) -> List[Interval]:
    XU = _normalize(X_clean)
    merged = []
    for raw in [p for p in (patterns or "").split(",") if p.strip()]:
        stripped, rel_spans = _parse_region_pattern(raw)
        hits = _find_all(XU, stripped)
        if len(hits) == 0:
            raise ValueError(f"Region pattern not found in X: {raw!r}")
        if len(hits) > 1:
            raise ValueError(f"Region pattern is not unique in X (matches={hits}): {raw!r}")
        base = hits[0]
        for iv in rel_spans:
            merged.append(Interval(base + iv.start, base + iv.end))
    return _merge_intervals(merged)