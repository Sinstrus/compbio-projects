# Critical Rules: Insertion Point Calculations

> Part of: [Agent Instructions](../README.md)

---

## Rule: Amino Acid Insertion Points (BUG-004)

**Context:** When designing insertions into a protein sequence at a specific amino acid position (e.g., "insert VHH3 after SKTINGSG, before QNQQTLKF"), you MUST calculate the DNA insertion point programmatically from the amino acid coordinates. NEVER hard-code DNA positions.

**The Problem:**
Hard-coding DNA positions leads to off-by-one-codon errors that produce incorrect protein sequences. This is a CRITICAL error because:
- The synthesized construct will have the wrong sequence
- Flanking sequences will be incorrect
- Protein function may be compromised
- Synthesis costs are wasted (~$200-400)

**The Correct Approach:**

```python
def find_aa_insertion_point(sequence, cds_start, upstream_aa, downstream_aa):
    """
    Calculate DNA insertion point from amino acid boundaries.

    Args:
        sequence: Full DNA sequence
        cds_start: CDS start position (0-indexed)
        upstream_aa: AA sequence before insertion (e.g., "SKTINGSG")
        downstream_aa: AA sequence after insertion (e.g., "QNQQTLKF")

    Returns:
        DNA position for insertion (0-indexed)
    """
    # Translate CDS to amino acids
    cds_seq = sequence[cds_start:]
    cds_aa = str(Seq(cds_seq).translate())

    # Find upstream sequence
    upstream_pos = cds_aa.find(upstream_aa)
    if upstream_pos == -1:
        raise ValueError(f"Upstream '{upstream_aa}' not found")

    # Calculate DNA position after upstream sequence
    aa_end_pos = upstream_pos + len(upstream_aa)
    dna_insertion_point = cds_start + (aa_end_pos * 3)

    # Verify downstream is immediately after
    downstream_pos = cds_aa.find(downstream_aa, aa_end_pos)
    if downstream_pos != aa_end_pos:
        raise ValueError(f"Downstream '{downstream_aa}' not immediately after")

    return dna_insertion_point
```

**Mandatory Verification:**

After calculating the insertion point, ALWAYS verify flanking sequences:

```python
# Verify upstream
before_dna = sequence[insertion_point - len(upstream_aa)*3 : insertion_point]
before_aa = str(Seq(before_dna).translate())
assert before_aa == upstream_aa, f"Upstream is {before_aa}, expected {upstream_aa}"

# Verify downstream
after_dna = sequence[insertion_point : insertion_point + len(downstream_aa)*3]
after_aa = str(Seq(after_dna).translate())
assert after_aa == downstream_aa, f"Downstream is {after_aa}, expected {downstream_aa}"
```

**Example Error (BUG-004):**
```
INTENDED: ...SKTINGSG | INSERT | QNQQTLKF...
ACTUAL:   ...SKTINGSGQ | INSERT | NQQTLKF...  <- Off by 1 codon!

Root cause: Used hard-coded position 3746 instead of calculating from AA coords (3743)
```

**Checklist:**
- [ ] Calculate insertion point from amino acid coordinates (NEVER hard-code)
- [ ] Verify upstream flanking sequence matches expectations
- [ ] Verify downstream flanking sequence matches expectations
- [ ] Include flanking verification in build script
- [ ] Document expected flanking sequences in design report

---

## Rule: Dipeptide Insertions (BUG-006)

**Context:** When inserting into an N-terminal position after a dipeptide (e.g., "MA-ADGYLPD"), you MUST count BOTH amino acids in the dipeptide before calculating the DNA insertion point.

**The Problem:**
Notation like "MA-ADGYLPD" can be misinterpreted:
- WRONG interpretation: "Insert after M, before A-ADGYLPD" (1 amino acid)
- CORRECT interpretation: "Insert after M-A, before ADGYLPD" (2 amino acids)

The dash "-" indicates the insertion point, not a separator between amino acids.

**Example Error (BUG-006):**
```
INTENDED: M-A-[VHH3]-[GGGGS5]-A-D-G-Y-L...
          ^ Insert after M-A dipeptide (2 amino acids)

ACTUAL:   M-[VHH3]-[GGGGS5]-A-A-D-G-Y-L...
          ^ Inserted after M only (1 amino acid) <- Missing first A!

Root cause: Used insertion_point = 2381 (after ATG) instead of 2384 (after ATG GCT)
```

**The Correct Approach:**

```python
def verify_dipeptide_insertion(sequence, cds_start, dipeptide, downstream_aa):
    """
    Calculate DNA insertion point for dipeptide-based N-terminal insertions.
    """
    cds_seq = sequence[cds_start:]
    cds_aa = str(Seq(cds_seq).translate())

    if not cds_aa.startswith(dipeptide):
        raise ValueError(f"CDS does not start with {dipeptide}, got {cds_aa[:2]}")

    # After dipeptide (2 amino acids = 6 bp)
    insertion_point = cds_start + (len(dipeptide) * 3)

    # Verify downstream
    downstream_pos = len(dipeptide)
    if not cds_aa[downstream_pos:].startswith(downstream_aa):
        raise ValueError(f"Downstream mismatch")

    return insertion_point
```

**Notation Clarity:**
```
# AMBIGUOUS:
"MA-ADGYLPD" - Could mean "M" + "A-ADGYLPD" or "MA" + "ADGYLPD"

# CLEAR:
"Insert after M-A dipeptide, before A-D-G-Y-L"
"Notation: M-A | [INSERT] | A-D-G-Y-L"
```

**Checklist:**
- [ ] Count amino acids explicitly (dipeptide = 2 amino acids, not 1)
- [ ] Calculate insertion point from amino acid count (never hard-code)
- [ ] Verify upstream flanking sequence matches dipeptide
- [ ] Verify downstream flanking sequence matches expectations
- [ ] Use clear notation: "M-A | [INSERT] | A-D-G-Y-L"
