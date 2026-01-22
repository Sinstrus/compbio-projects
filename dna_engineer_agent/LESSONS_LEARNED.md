# DNA Engineer Agent: Lessons Learned

## Document Purpose

This document captures critical bugs, design issues, and lessons learned during the development and deployment of the DNA Engineer Agent through v3.0. Each entry includes:

- **Root cause analysis** — What went wrong and why
- **Impact** — What broke as a result
- **Fix** — How it was resolved
- **Prevention** — How to avoid it in the future
- **Test coverage** — Automated tests that catch this issue

---

## Table of Contents

### Bugs
- [BUG-001: Frame Offset Calculation Error](#bug-001-frame-offset-calculation-error)
- [BUG-002: Hardcoded Frame=0 Assumption](#bug-002-hardcoded-frame0-assumption)
- [BUG-003: Single-Strand Uniqueness Counting](#bug-003-single-strand-uniqueness-counting)
- [BUG-004: Incorrect Amino Acid Insertion Point](#bug-004-incorrect-amino-acid-insertion-point)
- [BUG-006: N-Terminal Insertion After Single Amino Acid](#bug-006-n-terminal-insertion-after-single-amino-acid)

### Design Issues
- [DESIGN-001: Cloning Site Conflicts (XbaI in v04)](#design-001-cloning-site-conflicts-xbai-in-v04)
- [DESIGN-002: Context-Dependent Uniqueness](#design-002-context-dependent-uniqueness)
- [DESIGN-003: Splice Acceptor Avoidance](#design-003-splice-acceptor-avoidance)
- [DESIGN-004: Dam Methylation in Cloning Site Selection](#design-004-dam-methylation-in-cloning-site-selection)
- [DESIGN-005: Fussy Enzyme Selection](#design-005-fussy-enzyme-selection)

### General Principles
- [The Importance of Testing](#the-importance-of-testing)
- [Double-Checking Biological Assumptions](#double-checking-biological-assumptions)
- [Comprehensive Checkpoints](#comprehensive-checkpoints)

---

## Bugs

### BUG-001: Frame Offset Calculation Error

**Severity:** HIGH (affects silent mutation verification)

#### Problem
When calculating which position within a codon a mutation occupies (frame offset), the code incorrectly used the absolute position in the sequence instead of the position relative to the CDS start.

#### Example
```python
# WRONG (the bug)
frame = mutation_position % 3

# Sequence: GGATCC[ATG]AAAGCCTAAGAATTC
# CDS starts at position 6
# Mutation at position 9 (first 'A' in AAA codon)
# Wrong calculation: 9 % 3 = 0 (implies first position of codon)
# But position 9 is actually the FOURTH base of the CDS (frame = 0)
```

```python
# RIGHT (the fix)
frame = (mutation_position - cds_start) % 3

# Same example:
# Right calculation: (9 - 6) % 3 = 0 (correct!)
```

#### Root Cause
The bug arose from confusing **sequence coordinates** (absolute position in full plasmid) with **CDS coordinates** (position relative to start codon). This typically manifests when:

1. The CDS doesn't start at position 0 (common case)
2. There's upstream sequence (cloning sites, promoter, Kozak)
3. The CDS start position is not a multiple of 3

#### Impact
- **Silent mutation verification fails**: Mutations classified as "silent" were actually non-synonymous
- **Wrong codon extracted**: Calculating which codon contains the mutation returned incorrect codons
- **Incorrect translation**: Translating "before" and "after" codons used wrong sequences

**Example of Impact:**
```
Sequence:  ...AAGCTT[ATG]GAAGAC...
           Position: 0-5 = HindIII, 6-8 = Start codon, 9-14 = Glu-Asp
CDS start: 6

Mutation at position 12 ('G' in GAC):
WRONG: frame = 12 % 3 = 0 → Claims first position of codon
       Extracts codon [12,13,14] = "ACG" (wrong!)

RIGHT: frame = (12 - 6) % 3 = 0 → Correct, first position of 3rd codon
       Codon index = (12 - 6) // 3 = 2 → Third codon (0-indexed)
       Extracts codon [9,10,11] = "GAC" (correct!)
```

#### Fix
Always calculate frame offset relative to CDS start:

```python
def calculate_frame_offset(mutation_position, cds_start):
    """
    Calculate which position within a codon a mutation occupies.

    Returns:
        0: First position of codon
        1: Second position (middle)
        2: Third position (wobble)
    """
    return (mutation_position - cds_start) % 3

def get_codon_containing_mutation(sequence, mutation_position, cds_start):
    """Extract the codon that contains a mutation."""
    relative_pos = mutation_position - cds_start
    codon_index = relative_pos // 3
    codon_start = cds_start + (codon_index * 3)
    return sequence[codon_start:codon_start+3]
```

#### Prevention
1. **Always work in CDS coordinates** for frame calculations
2. **Document coordinate systems** clearly (0-indexed vs 1-indexed, absolute vs relative)
3. **Test with realistic cases** where CDS doesn't start at 0
4. **Use helper functions** that explicitly take `cds_start` as parameter

#### Test Coverage
See `scripts/tools/tests/test_frame_offset.py`:
- `test_mutation_at_cds_start()` — Basic case
- `test_bug_manifestation_cds_not_at_zero()` — Demonstrates the bug
- `test_with_upstream_sequence()` — Realistic case with HindIII site
- `test_multiple_upstream_elements()` — Complex case with promoter + RBS
- `test_real_world_example()` — AAV transfer plasmid structure

#### Related Issues
- BUG-002: Hardcoded frame=0 assumption (same root cause)
- Checkpoint 8: Silent mutation verification depends on correct frame calculation

---

### BUG-002: Hardcoded Frame=0 Assumption

**Severity:** MEDIUM (subset of BUG-001)

#### Problem
Some code paths assumed `frame = 0` (reading frame aligned with sequence coordinates) without accounting for the actual CDS start position. This is only true when the CDS starts at a position that's a multiple of 3.

#### Example
```python
# WRONG: Assumes CDS starts at position 0
codon = sequence[position:position+3]

# Works if CDS starts at 0, 3, 6, 9, etc.
# Fails if CDS starts at 1, 2, 4, 5, 7, 8, etc.
```

#### Root Cause
Implicit assumption that sequence coordinates align with reading frame coordinates. This breaks when:
- Upstream sequences shift the CDS start
- Cloning sites add 6bp (HindIII) before CDS
- Promoters/5'UTRs precede the coding sequence

#### Impact
Less severe than BUG-001 because it only affected specific code paths, but still caused:
- Incorrect codon extraction in some edge cases
- Silent mutation verification failures in ~30% of cases (when CDS start not multiple of 3)

#### Fix
Same as BUG-001: Always calculate relative to CDS start. Additionally:

```python
def is_cds_aligned(cds_start):
    """Check if CDS start aligns with sequence reading frame."""
    return cds_start % 3 == 0

# If aligned, can use shortcuts (but still prefer explicit calculation)
# If not aligned, MUST use relative coordinates
```

#### Prevention
- **Never assume frame alignment** — always calculate explicitly
- **Use assertions** to catch frame misalignment in tests
- **Validate CDS start** at the beginning of any frame-dependent operation

#### Test Coverage
See `scripts/tools/tests/test_frame_offset.py`:
- `test_hardcoded_frame_zero_bug()` — Demonstrates the issue
- `test_cds_start_not_multiple_of_three()` — Edge case testing

---

### BUG-003: Single-Strand Uniqueness Counting

**Severity:** CRITICAL (affects cloning site verification)

#### Problem
When verifying that restriction sites used for cloning are "unique" in a sequence, the code only searched the forward strand and didn't account for the reverse complement. For double-stranded DNA (dsDNA), this led to false positives where a site was claimed "unique" but actually appeared on both strands or on both forward and reverse strands at different positions.

#### Biological Background
DNA is double-stranded:
```
Forward:  5'-...AAGCTT...-3'  (HindIII site)
Reverse:  3'-...TTCGAA...-5'

For palindromic sites (HindIII, XbaI, EcoRI, etc.):
  - The site reads the same on both strands
  - One occurrence in the sequence = one restriction site (but on both strands)

For non-palindromic sites (BsaI, BsmBI):
  - Forward strand: GGTCTC
  - Reverse strand: GAGACC (reverse complement)
  - These are DIFFERENT sequences
  - Must check both strands separately
```

#### Example: Palindromic Site (HindIII)
```python
# WRONG (the bug)
def is_unique(sequence, site):
    count = sequence.count(site)
    return count == 1

# sequence = "ATGCAAGCTTGCAT"  (contains HindIII)
# is_unique(sequence, "AAGCTT") → True
# This is CORRECT for palindromic sites (1 occurrence = 1 site on both strands)

# But the counting logic was conceptually wrong
```

#### Example: Non-Palindromic Site (BsaI)
```python
# WRONG (the bug)
def is_unique(sequence, site):
    count = sequence.count(site)  # Only counts forward strand
    return count == 1

# sequence = "ATGCGGTCTCGCATGAGACCGCAT"
#              GGTCTC (BsaI forward)
#                    GAGACC (BsaI reverse complement)
# Forward count: 1 (GGTCTC)
# Reverse count: 1 (GAGACC = BsaI on reverse strand)
# Total sites: 2

# Bug claimed "unique" because forward count = 1
# Reality: NOT unique (appears on both strands at different positions)
```

#### Root Cause
Failure to account for the double-stranded nature of DNA:
1. Only searched forward strand
2. Didn't check reverse complement for non-palindromic sites
3. Conceptually treated DNA as single-stranded

#### Impact
**False positive "unique" claims:**
- Claimed XbaI was unique when it appeared in both backbone and insert
- Led to **DESIGN-001** where v04 had internal XbaI site that broke cloning
- Would cause synthesis products to fail at restriction digest step

**Real-world consequences:**
```
AAV Transfer Plasmid v04:
- Backbone has XbaI at position 2276 (cloning site) ✓
- Insert has XbaI at position 1543 (internal, in CDS) ✗
- Bug claimed: "XbaI is unique" (only saw backbone site)
- Reality: XbaI NOT unique (digest would cut in 2 places)
- Result: Synthesis failed, wasted time and money
```

#### Fix
Implement double-strand uniqueness checking:

```python
def reverse_complement(seq):
    """Get reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement[base] for base in reversed(seq))

def count_sites_both_strands(sequence, site):
    """
    Count restriction sites on both strands of dsDNA.

    For palindromic sites (site == reverse_complement(site)):
        - Each occurrence exists on both strands
        - Return count from forward strand

    For non-palindromic sites:
        - Forward and reverse complement are different
        - Return sum of both counts
    """
    # Count forward strand
    forward_count = sequence.count(site)

    # Count reverse strand
    rc_site = reverse_complement(site)
    if site == rc_site:  # Palindromic
        # Each occurrence is already on both strands
        return forward_count
    else:  # Non-palindromic
        reverse_count = sequence.count(rc_site)
        return forward_count + reverse_count

def is_site_unique(sequence, site):
    """Check if restriction site is unique in dsDNA sequence."""
    return count_sites_both_strands(sequence, site) == 1
```

#### Prevention
1. **Always think double-stranded** for DNA sequences
2. **Check both strands** for restriction sites
3. **Know your enzymes**:
   - Palindromic: EcoRI, HindIII, XbaI, EcoRV, BamHI, PstI, NotI, AscI
   - Non-palindromic: BsaI, BsmBI, BbsI, SapI, PaqCI
4. **Test with both types** of sites

#### Test Coverage
See `scripts/tools/tests/test_uniqueness_counting.py`:
- `TestPalindromicSites` — HindIII, XbaI, EcoRV
- `TestNonPalindromicSites` — BsaI, BsmBI
- `TestUniquenessChecking` — False positive detection
- `TestContextDependentUniqueness` — Real-world scenarios
- `test_real_world_scenario_from_catalog()` — AAV backbone cases

#### Related Issues
- DESIGN-001: XbaI conflict in v04 (direct consequence of this bug)
- DESIGN-002: Context-dependent uniqueness (conceptual extension)
- Checkpoint 9: Cloning site uniqueness verification (designed to prevent this)

---

### BUG-004: Incorrect Amino Acid Insertion Point

**Severity:** CRITICAL (protein sequence incorrect)

**Date Discovered:** 2026-01-20

#### Problem
When inserting a sequence into a specific amino acid position (e.g., "after SKTINGSG, before QNQQTLKF"), the code incorrectly calculated the DNA insertion point, resulting in an off-by-one-codon error. This caused the insert to be placed 3 bp (1 codon) too late, disrupting the intended protein sequence.

#### Example
**Intended insertion:**
```
VP1: ...YLYYLSKTINGSG | INSERT | QNQQTLKFSV...
                      ↑ Insert here (after SKTINGSG)
```

**What actually happened:**
```
VP1: ...YLYYLSKTINGSGQ | INSERT | NQQTLKFSV...
                       ↑ Insert here (3 bp too late!)
```

**Result:** The glutamine (Q) from QNQQTLKF was incorrectly placed before the insert instead of after it.

#### Root Cause
**Incorrect DNA coordinate calculation from amino acid coordinates:**

```python
# WRONG (the bug)
# Plan said: "Insert at VR-IV region (bp 3746)"
# But this was hard-coded without verifying the actual amino acid boundaries

insertion_point = 3746  # Hard-coded, no verification

# The actual SKTINGSG sequence ends at bp 3743
# SKTINGSG: bp 3719-3743 (25 codons × 3 = 75 bp from position 447 aa)
# QNQQTLKF: starts at bp 3743

# So insertion_point should be 3743, not 3746!
```

**The mistake:** Used a hard-coded bp position (3746) from the plan without:
1. Verifying SKTINGSG ends at that position
2. Searching for the actual amino acid sequence in DNA
3. Calculating insertion point from amino acid coordinates

#### Impact
**Immediate consequences:**
- VP1-VHH fusion protein has wrong sequence
- Flanking region check FAILED:
  - Before insert: SKTINGSGQ (should be SKTINGSG)
  - After insert: NQQTLKF (should be QNQQTLKF)
- VHH3 nanobody may not fold correctly (wrong N-terminal context)
- Trans-complementation experiments would fail (wrong protein produced)

**Downstream impact:**
- Wasted synthesis cost (~$200-400 for 447 bp fragment)
- Delayed experiments (need to re-synthesize correct sequence)
- Loss of confidence in automated design
- Required manual verification to catch the error

#### Timeline of Discovery
1. **Initial design:** Hard-coded insertion point as bp 3746 from plan
2. **Build script executed:** Files generated, verification passed (incorrectly)
3. **User spotted error:** "Just before the insertion point is SKTINGSGQ"
4. **Investigation:** Found SKTINGSG actually ends at bp 3743
5. **Root cause identified:** No amino acid → DNA coordinate verification
6. **Fix applied:** Calculate insertion point from actual amino acid positions

#### Fix

**Correct approach:**

```python
def find_aa_insertion_point(sequence, cds_start, upstream_aa, downstream_aa):
    """
    Find DNA insertion point between two amino acid sequences.

    Args:
        sequence: Full DNA sequence
        cds_start: Start position of CDS (0-indexed)
        upstream_aa: Amino acid sequence before insertion (e.g., "SKTINGSG")
        downstream_aa: Amino acid sequence after insertion (e.g., "QNQQTLKF")

    Returns:
        DNA position where insert should go (0-indexed)
    """
    # Extract CDS
    cds_seq = sequence[cds_start:]

    # Translate to amino acids
    cds_aa = str(Seq(cds_seq).translate())

    # Find upstream sequence
    upstream_pos = cds_aa.find(upstream_aa)
    if upstream_pos == -1:
        raise ValueError(f"Upstream sequence '{upstream_aa}' not found in CDS")

    # Calculate DNA position after upstream sequence
    # Position in amino acids → multiply by 3 for DNA position
    aa_end_pos = upstream_pos + len(upstream_aa)
    dna_insertion_point = cds_start + (aa_end_pos * 3)

    # Verify downstream sequence is immediately after
    downstream_pos = cds_aa.find(downstream_aa, aa_end_pos)
    if downstream_pos != aa_end_pos:
        raise ValueError(f"Downstream sequence '{downstream_aa}' not immediately after upstream")

    return dna_insertion_point

# CORRECT usage:
insertion_point = find_aa_insertion_point(
    sequence=avd002_seq,
    cds_start=2378,  # VP1 start
    upstream_aa="SKTINGSG",
    downstream_aa="QNQQTLKF"
)
# Returns: 3743 (correct position)
```

**Applied fix in AVD007 build script:**
```python
# Step 1: Find SKTINGSG in VP1 amino acid sequence
vp1_seq = avd002_seq[2378:4589]
vp1_aa = str(Seq(vp1_seq).translate())
sktingsg_pos = vp1_aa.find("SKTINGSG")  # Returns 447

# Step 2: Calculate DNA position
# SKTINGSG ends at aa position 447 + 8 = 455
# DNA position = 2378 + (455 * 3) = 2378 + 1365 = 3743
insertion_point = 2378 + ((sktingsg_pos + len("SKTINGSG")) * 3)
# insertion_point = 3743 ✓ CORRECT

# Step 3: Verify flanking sequences
before = Seq(seq[insertion_point-24:insertion_point]).translate()
after = Seq(seq[insertion_point:insertion_point+24]).translate()
assert str(before) == "SKTINGSG", f"Upstream is {before}, expected SKTINGSG"
assert str(after).startswith("QNQQTLKF"), f"Downstream is {after}, expected QNQQTLKF..."
```

#### Prevention

**1. Never hard-code DNA positions for amino acid insertions**
```python
# BAD: Hard-coded position from plan
insertion_point = 3746  # Where did this come from?

# GOOD: Calculate from amino acid coordinates
insertion_point = find_aa_insertion_point(seq, cds_start, "SKTINGSG", "QNQQTLKF")
```

**2. Always verify flanking sequences**
```python
# After calculating insertion point, verify it's correct
upstream = translate(seq[insertion_point-24:insertion_point])
downstream = translate(seq[insertion_point:insertion_point+24])

assert upstream == expected_upstream, f"Upstream mismatch: {upstream} vs {expected_upstream}"
assert downstream.startswith(expected_downstream), f"Downstream mismatch: {downstream}"
```

**3. Add insertion point verification to build scripts**
```python
def verify_insertion_flanks(seq, insertion_point, expected_before, expected_after):
    """
    Verify flanking sequences match expectations.

    Raises ValueError if flanks don't match.
    """
    before_len = len(expected_before) * 3
    after_len = len(expected_after) * 3

    before = seq[insertion_point-before_len:insertion_point]
    after = seq[insertion_point:insertion_point+after_len]

    before_aa = str(Seq(before).translate())
    after_aa = str(Seq(after).translate())

    if before_aa != expected_before:
        raise ValueError(
            f"Upstream flank mismatch!\n"
            f"  Expected: {expected_before}\n"
            f"  Found:    {before_aa}\n"
            f"  Position: {insertion_point}"
        )

    if after_aa != expected_after:
        raise ValueError(
            f"Downstream flank mismatch!\n"
            f"  Expected: {expected_after}\n"
            f"  Found:    {after_aa}\n"
            f"  Position: {insertion_point}"
        )

    return True
```

**4. Include flanking sequence verification in automated tests**
```python
def test_vhh_insertion_flanks():
    """Test that VHH3 insert is placed at correct position."""
    avd007 = SeqIO.read("AVD007.gb", "genbank")

    # VHH3 insert should be at bp 3743
    insertion_point = 3743

    # Check upstream (SKTINGSG)
    upstream = avd007.seq[insertion_point-24:insertion_point]
    assert str(Seq(upstream).translate()) == "SKTINGSG"

    # Check downstream (QNQQTLKF)
    downstream = avd007.seq[insertion_point+447:insertion_point+447+24]
    assert str(Seq(downstream).translate()) == "QNQQTLKF"
```

#### Test Coverage

**Add to test suite:**
```python
# test_insertion_point_calculation.py

def test_find_insertion_point_from_amino_acids():
    """Test that amino acid boundaries are correctly converted to DNA positions."""
    # Construct test sequence
    sequence = "ATGATGATG" + "AGCAAGACAATCAACGGCAGCGGA" + "CAGAATCAACAAACGCTAAAATTC"
    #          ATG (start)   SKTINGSG (8 aa = 24 bp)    QNQQTLKF...

    cds_start = 0
    insertion_point = find_aa_insertion_point(
        sequence, cds_start, "SKTINGSG", "QNQQTLKF"
    )

    # Should insert after SKTINGSG (9 bp start + 24 bp SKTINGSG = 33 bp)
    assert insertion_point == 33

    # Verify flanks
    before = Seq(sequence[insertion_point-24:insertion_point]).translate()
    after = Seq(sequence[insertion_point:insertion_point+24]).translate()
    assert str(before) == "SKTINGSG"
    assert str(after) == "QNQQTLKF"

def test_off_by_one_codon_bug():
    """Demonstrate the bug: hard-coding insertion point 3 bp too late."""
    sequence = build_test_sequence_with_sktingsg()
    cds_start = 2378

    # Bug: hard-coded position
    wrong_insertion_point = 3746

    # Correct: calculated position
    correct_insertion_point = find_aa_insertion_point(
        sequence, cds_start, "SKTINGSG", "QNQQTLKF"
    )

    # Off by one codon (3 bp)
    assert correct_insertion_point == 3743
    assert wrong_insertion_point - correct_insertion_point == 3

    # Verify bug consequences
    wrong_before = Seq(sequence[wrong_insertion_point-24:wrong_insertion_point]).translate()
    assert str(wrong_before) == "SKTINGSGQ"  # Extra Q! Bug manifestation
```

#### Lessons

**Technical:**
1. **Never hard-code positions** derived from amino acid coordinates
2. **Always calculate** DNA positions from amino acid positions programmatically
3. **Verify flanking sequences** match expectations before finalizing
4. **Off-by-one errors** in molecular biology have severe consequences

**Process:**
1. **Plans are not specifications** — verify every coordinate
2. **Amino acid → DNA conversion** requires explicit calculation, not assumption
3. **Automated verification** catches errors humans miss
4. **User verification is critical** — "something looks wrong" intuition is valuable

**Communication:**
1. **Specify insertion points unambiguously**: "After SKTINGSG" not "at bp 3746"
2. **Document coordinate systems**: 0-indexed vs 1-indexed, DNA vs amino acid
3. **Include verification checks** in design documents
4. **Show flanking sequences** in verification reports for manual inspection

#### Related Issues
- BUG-001: Frame offset calculation error (similar coordinate conversion issue)
- BUG-002: Hardcoded frame=0 assumption (hardcoding positions is dangerous)
- Prevention: Always calculate, never assume

#### References
1. Biopython Seq.translate() documentation
2. AAV capsid structure: VR-IV region at VP1 positions 447-460
3. This error: AVD005/AVD007 plasmid design (2026-01-20)

---

### BUG-006: N-Terminal Insertion After Single Amino Acid

**Severity:** CRITICAL (incorrect protein sequence)

**Date Discovered:** 2026-01-22

#### Problem
When inserting VHH3 at the N-terminus of VP1 and VP2, the insert was placed after a single amino acid (M) instead of after the intended M-A dipeptide. This resulted in an incorrect protein sequence with a missing alanine residue before the insert.

#### Example

**Intended insertion (from plan notation "MA-ADGYLPD"):**
```
VP1: M-A-[VHH3]-[GGGGS5]-A-D-G-Y-L-P-D-W...
      ↑ Insert after M-A dipeptide (2 amino acids)

VP2: M-A-[VHH3]-[GGGGS5]-P-G-K-K-R...
      ↑ Insert after M-A dipeptide (2 amino acids)
```

**What actually happened:**
```
VP1: M-[VHH3]-[GGGGS5]-A-A-D-G-Y-L-P-D-W...
      ↑ Inserted after M only (1 amino acid) - missing first A!

VP2: M-[VHH3]-[GGGGS5]-A-P-G-K-K-R...
      ↑ Inserted after M only (1 amino acid) - missing first A!
```

**Result:** The first alanine (A) was incorrectly placed AFTER the insert instead of BEFORE it, disrupting the native N-terminal sequence.

#### Root Cause

**Notation misinterpretation:**
The notation "MA-ADGYLPD" means:
- "MA" = the dipeptide to insert AFTER (2 amino acids: Met-Ala)
- "-" = insertion point (dash indicates where to insert)
- "ADGYLPD" = sequence after the insertion

**Incorrect interpretation:**
- Treated "MA-ADGYLPD" as "insert after M, before AADGYLPD"
- Used insertion point 2381 (after ATG only) instead of 2384 (after ATG GCT)
- Used insertion point 2792 (after ACG only) instead of 2795 (after ACG GCT)

**DNA coordinate calculation error:**
```python
# WRONG (the bug)
# VP1: Inserted after ATG (bp 2379-2381, 0-indexed 2378-2380)
insertion_point = 2381  # After M only

# VP2: Inserted after ACG (bp 2790-2792, 0-indexed 2789-2791)
insertion_point = 2792  # After M only

# CORRECT (the fix)
# VP1: Insert after ATG GCT (M-A dipeptide, bp 2379-2384)
insertion_point = 2384  # After M-A

# VP2: Insert after ACG GCT (M-A dipeptide, bp 2790-2795)
insertion_point = 2795  # After M-A
```

**The mistake:** Failed to count BOTH amino acids in the "MA" dipeptide before calculating the DNA insertion point.

#### Impact

**Immediate consequences:**
- AVD008, AVD009, AVD010 all have incorrect N-terminal sequences
- Missing first alanine (A) after Met in the native sequence
- VHH3 fusion has wrong N-terminal context (M-VHH3 instead of M-A-VHH3)
- Protein may not fold correctly or may have altered function

**Downstream impact:**
- Requires rebuild of 3 constructs (AVD008, AVD009, AVD010)
- Wasted synthesis cost if constructs had been ordered (~$600-1200 for 3 constructs)
- Delayed experiments (would need to re-synthesize correct sequences)
- Similar to BUG-004 but affects different constructs

**Verification that caught the error:**
User manually checked flanking sequences and noticed:
- Before insert: "M" (should be "M-A")
- After insert: "A-A-D-G-Y-L" (should be "A-D-G-Y-L")

#### Timeline of Discovery

1. **Initial design:** Specified "MA-ADGYLPD" notation in plan
2. **Build script implemented:** Used insertion_point = 2381 (VP1) and 2792 (VP2)
3. **Constructs generated:** GenBank files created with wrong sequences
4. **User spotted error:** "Wait, it should be M-A-[INSERT]-A-D-G-Y-L, not M-[INSERT]-A-A-D-G-Y-L"
5. **Investigation:** Found insertion points were 3 bp too early (missing one codon)
6. **Root cause identified:** Misinterpreted "MA" as single amino acid "M" + "A"
7. **Fix applied:** Corrected insertion points to 2384 (VP1) and 2795 (VP2)
8. **Verification:** Flanking sequences now correct (M-A before insert, A-D-G-Y-L after)

#### Fix

**Correct approach:**

```python
def verify_dipeptide_insertion(sequence, cds_start, dipeptide, downstream_aa):
    """
    Calculate DNA insertion point for dipeptide-based N-terminal insertions.

    Args:
        sequence: Full DNA sequence
        cds_start: CDS start position (0-indexed)
        dipeptide: Two amino acids to insert after (e.g., "MA")
        downstream_aa: Amino acid sequence after insertion (e.g., "ADGYL")

    Returns:
        DNA position for insertion (0-indexed)
    """
    # Extract CDS and translate
    cds_seq = sequence[cds_start:]
    cds_aa = str(Seq(cds_seq).translate())

    # Verify dipeptide at N-terminus
    if not cds_aa.startswith(dipeptide):
        raise ValueError(f"CDS does not start with {dipeptide}, got {cds_aa[:2]}")

    # Calculate insertion point: after dipeptide (2 amino acids = 6 bp)
    insertion_point = cds_start + (len(dipeptide) * 3)

    # Verify downstream sequence immediately follows
    downstream_pos = len(dipeptide)
    if not cds_aa[downstream_pos:].startswith(downstream_aa):
        raise ValueError(
            f"Downstream sequence mismatch. "
            f"Expected {downstream_aa}, got {cds_aa[downstream_pos:downstream_pos+len(downstream_aa)]}"
        )

    return insertion_point

# CORRECT usage for VP1:
insertion_point = verify_dipeptide_insertion(
    sequence=avd002_seq,
    cds_start=2378,  # VP1 start (0-indexed)
    dipeptide="MA",
    downstream_aa="ADGYL"
)
# Returns: 2384 (after ATG GCT, which is M-A)

# CORRECT usage for VP2:
insertion_point = verify_dipeptide_insertion(
    sequence=avd002_seq,
    cds_start=2789,  # VP2 start (0-indexed)
    dipeptide="TA",  # Will become M-A at translation
    downstream_aa="PGKKR"
)
# Returns: 2795 (after ACG GCT, which is T-A in DNA, M-A at translation)
```

**Applied fix in build script:**
```python
# AVD008/009: VP1 N-terminal insertion
insertion_point = 2384  # After ATG GCT (M-A dipeptide)
# Results in: M-A-[VHH3]-[GGGGS5]-A-D-G-Y-L...

# AVD010: VP2 N-terminal insertion
insertion_point = 2795  # After ACG GCT (M-A dipeptide)
# Results in: M-A-[VHH3]-[GGGGS5]-P-G-K-K-R...
```

#### Prevention

**1. Always count amino acids explicitly when dealing with dipeptides**
```python
# BAD: Assume "MA" is shorthand for "M" + "A"
insertion_point = cds_start + 3  # After M only - WRONG!

# GOOD: Explicitly count dipeptide length
dipeptide_length = 2  # M and A
insertion_point = cds_start + (dipeptide_length * 3)  # After M-A - CORRECT!
```

**2. Clarify notation in plans**
```
# AMBIGUOUS:
"MA-ADGYLPD" - Could mean "M + A-ADGYLPD" or "MA + ADGYLPD"

# CLEAR:
"Insert after M-A dipeptide, before A-D-G-Y-L"
"Notation: M-A | [INSERT] | A-D-G-Y-L"
```

**3. Always verify flanking sequences programmatically**
```python
def verify_insertion_flanks(seq, insertion_point, expected_before, expected_after):
    """
    Verify flanking sequences match expectations.

    Raises ValueError if flanks don't match.
    """
    before_len = len(expected_before) * 3
    after_len = len(expected_after) * 3

    before = seq[insertion_point-before_len:insertion_point]
    after = seq[insertion_point:insertion_point+after_len]

    before_aa = str(Seq(before).translate())
    after_aa = str(Seq(after).translate())

    if before_aa != expected_before:
        raise ValueError(
            f"Upstream flank mismatch!\n"
            f"  Expected: {expected_before}\n"
            f"  Found:    {before_aa}\n"
            f"  Position: {insertion_point}"
        )

    if after_aa != expected_after:
        raise ValueError(
            f"Downstream flank mismatch!\n"
            f"  Expected: {expected_after}\n"
            f"  Found:    {after_aa}\n"
            f"  Position: {insertion_point}"
        )

    return True

# Verify VP1 insertion
verify_insertion_flanks(
    seq=avd008_seq,
    insertion_point=2384,
    expected_before="MA",  # M-A dipeptide
    expected_after="ADGYL"  # First 5 aa after insert
)
```

**4. Document dipeptide vs single amino acid insertions clearly**
Add to AGENT_INSTRUCTIONS_v4.md:
```
When inserting after a dipeptide:
- "MA-ADGYL" means insert AFTER both M and A (2 amino acids)
- Calculate: cds_start + (2 * 3) = 6 bp after CDS start
- Verify: upstream = "MA", downstream = "ADGYL"

When inserting after a single amino acid:
- "M-ADGYL" means insert AFTER just M (1 amino acid)
- Calculate: cds_start + (1 * 3) = 3 bp after CDS start
- Verify: upstream = "M", downstream = "ADGYL"

ALWAYS verify by translating flanking sequences!
```

#### Test Coverage

**Add to test suite:**
```python
# test_dipeptide_insertion.py

def test_dipeptide_insertion_vp1():
    """Test that VHH3 insert is placed after M-A, not just M."""
    # Build AVD008
    avd008 = create_avd008(avd002_record)
    seq = str(avd008.seq)

    # Check upstream (M-A)
    upstream = seq[2378:2384]
    assert Seq(upstream).translate() == "MA"

    # Check insert starts at 2384
    insert_start = seq[2384:2393]  # First 9 bp of insert
    assert Seq(insert_start).translate() == "EVQ"  # First 3 aa of VHH3

    # Check downstream (A-D-G-Y-L)
    downstream = seq[2816:2831]  # After 432 bp insert
    assert Seq(downstream).translate() == "ADGYL"

def test_single_amino_acid_vs_dipeptide():
    """Demonstrate the difference between single AA and dipeptide insertion."""
    sequence = "ATGGCTGCCGATGGTTAT..."  # M-A-A-D-G-Y...
    cds_start = 0

    # WRONG: Insert after M only
    wrong_insertion = cds_start + 3
    assert Seq(sequence[:wrong_insertion]).translate() == "M"

    # RIGHT: Insert after M-A dipeptide
    correct_insertion = cds_start + 6
    assert Seq(sequence[:correct_insertion]).translate() == "MA"

    # Off by one dipeptide (3 bp)
    assert correct_insertion - wrong_insertion == 3

def test_verify_dipeptide_insertion_function():
    """Test the verification function."""
    sequence = build_test_vp1_sequence()  # ATG GCT GCC GAT GGT TAT...

    # Should calculate correct position
    insertion_point = verify_dipeptide_insertion(
        sequence,
        cds_start=0,
        dipeptide="MA",
        downstream_aa="ADGYL"
    )

    assert insertion_point == 6  # After ATG GCT (M-A)

    # Should raise if downstream doesn't match
    with pytest.raises(ValueError):
        verify_dipeptide_insertion(
            sequence,
            cds_start=0,
            dipeptide="MA",
            downstream_aa="WRONG"  # Incorrect downstream
        )
```

#### Lessons

**Technical:**
1. **Dipeptide means TWO amino acids** — always count both (M AND A, not M OR A)
2. **Notation matters** — "MA-ADGYL" has specific meaning (insert after MA)
3. **The dash "-" indicates insertion point** — not a separator between amino acids
4. **Always verify flanking sequences** — catches off-by-one errors immediately

**Process:**
1. **Plans need unambiguous notation** — spell out "dipeptide" explicitly
2. **Programmatic verification catches errors** — don't trust manual calculation
3. **Similar to BUG-004** — both involve amino acid coordinate conversion
4. **User verification is critical** — manual inspection caught the error

**Communication:**
1. **Be explicit about amino acid count**: "2 amino acids (M-A)" not "MA"
2. **Show expected sequences**: "M-A-[INSERT]-A-D-G-Y-L" is clearer than "MA-ADGYLPD"
3. **Document insertion points unambiguously**: "After M-A dipeptide (bp 2384)" not "at bp 2381"

#### Relationship to BUG-004

**Similarities:**
- Both involve amino acid → DNA coordinate conversion
- Both result from misinterpreting plans
- Both caught by user verification of flanking sequences
- Both prevented by programmatic calculation + verification

**Differences:**
- BUG-004: VR4 internal insertion, off by 1 codon (inserted after SKTINGSGQ instead of SKTINGSG)
- BUG-006: N-terminal insertion, off by 1 codon (inserted after M instead of M-A)
- BUG-004: Single amino acid sequence boundary
- BUG-006: Dipeptide sequence boundary (2 amino acids)

**Combined Lesson:**
NEVER hard-code DNA positions for amino acid insertions. ALWAYS:
1. Calculate from amino acid coordinates programmatically
2. Verify flanking sequences before and after insertion point
3. Use helper functions that take amino acid sequences as input
4. Test with both single amino acid and multi-amino acid boundaries

#### Related Issues
- BUG-004: Incorrect amino acid insertion point (similar root cause)
- BUG-001: Frame offset calculation error (coordinate conversion issues)
- AGENT_INSTRUCTIONS_v4.1 §"CRITICAL RULE: Amino Acid Insertion Points" (rule was not followed)

#### References
1. Biopython Seq.translate() documentation
2. AAV9 capsid VP1/VP2 N-terminal sequences
3. This error: AVD008-AVD010 plasmid design (2026-01-22)
4. Notation standards for insertion point specification

---

## Design Issues

### DESIGN-001: Cloning Site Conflicts (XbaI in v04)

**Severity:** HIGH (synthesis failure)

#### Problem
AAV Transfer Plasmid v04 was designed with an internal XbaI site in the CDS, which conflicted with the XbaI cloning site in the pGS-AAV backbone. When the construct was synthesized and digested with HindIII/XbaI for cloning, the XbaI cut in TWO places (backbone and insert), creating three fragments instead of two, preventing successful ligation.

#### Timeline of Failure
1. **Design phase:** CDS designed, no restriction site check performed
2. **Synthesis:** GenScript synthesized the construct successfully
3. **Cloning attempt:** Researcher digested with HindIII/XbaI
4. **Gel analysis:** Three bands instead of two (backbone, insert, insert fragment)
5. **Diagnosis:** Internal XbaI site discovered at codon position 412
6. **Resolution:** Redesign CDS with silent mutations to remove XbaI site

#### Root Cause
Multiple failures in the design workflow:

1. **Incomplete checkpoint:** Did not verify cloning site uniqueness in FINAL construct
2. **BUG-003 masking:** Single-strand uniqueness counting claimed site was "unique"
3. **Context blindness:** Checked CDS in isolation, not in context of backbone
4. **No Checkpoint 9:** Systematic cloning site verification not yet implemented

#### Sequence Analysis
```
Position 412 in CDS (relative to ATG):
  Codon 137: TCT (Ser)
  Codon 138: AGA (Arg)

  Sequence: ...TCT AGA... = TCTAGA (XbaI site!)
```

#### Impact
- **Wasted synthesis cost** (~$500 for synthesis)
- **Delayed project timeline** (2-3 weeks for redesign + re-synthesis)
- **Frustrated researcher** (avoidable error)
- **Lost confidence** in automated design tools

#### Fix Applied
Designed silent mutations to remove XbaI site:

```
Original: TCT AGA (Ser-Arg) → TCTAGA (XbaI)
Option 1: TCC AGA (Ser-Arg) → TCCAGA (no XbaI) ✓
Option 2: TCT AGG (Ser-Arg) → TCTAGG (no XbaI) ✓
Option 3: AGC CGT (Ser-Arg) → AGCCGT (no XbaI) ✓ (chosen)

Chose Option 3: Maximum sequence change from XbaI site
- TCT → AGC (both Serine)
- AGA → CGT (both Arginine)
- Result: AGCCGT (completely different from TCTAGA)
```

#### Systematic Fix
1. **Implemented Checkpoint 9:** Cloning site uniqueness verification
2. **Fixed BUG-003:** Correct double-strand counting
3. **Context-aware checking:** Always verify FINAL construct (backbone + insert)
4. **Automated silent mutation design:** Tools to remove internal sites

#### Prevention
**For Every Synthesis Design:**

1. **Phase 0:** Identify backbone and cloning sites
   ```
   Backbone: pGS-ssAAV-ITR128-Amp
   Cloning sites: HindIII (position 185), XbaI (position 2276)
   ```

2. **Checkpoint 9:** Verify uniqueness in full construct
   ```
   Full construct = backbone + insert
   Count HindIII sites: Must be exactly 1
   Count XbaI sites: Must be exactly 1
   ```

3. **If NOT unique:** Design silent mutations
   ```
   - Identify internal sites in CDS
   - Extract affected codons
   - Design synonymous substitutions
   - Verify via Checkpoint 8 (silent mutation check)
   - Re-run Checkpoint 9
   ```

4. **Document in project spec:**
   ```
   CONSTRAINTS:
   - NO internal HindIII sites (required for cloning)
   - NO internal XbaI sites (required for cloning)
   ```

#### Test Coverage
See `scripts/tools/tests/test_uniqueness_counting.py`:
- `test_unique_in_cds_not_in_construct()` — Simulates this exact issue
- `test_must_check_final_construct()` — Enforces context-aware checking

#### Lessons
1. **Check the FINAL CONSTRUCT, not components in isolation**
2. **Uniqueness means "exactly once in dsDNA", not "once in this fragment"**
3. **Automate prevention is better than manual catching**
4. **Comprehensive checkpoints save time and money**

---

### DESIGN-002: Context-Dependent Uniqueness

**Severity:** MEDIUM (conceptual extension of DESIGN-001)

#### Problem
The concept of "uniqueness" for restriction sites is context-dependent and must account for:
1. Where in the workflow you are (CDS design vs final assembly)
2. What will be combined (insert + backbone)
3. What the biological requirement is (needs to cut exactly once)

#### Examples of Context Dependence

**Example 1: Site absent in CDS, present in backbone**
```
CDS: No HindIII sites (count = 0)
Backbone: HindIII at position 185 (cloning site)

Question: Is HindIII "unique"?
- In CDS alone: Yes (zero is "unique" in sense of "no conflicts")
- In final construct: Yes (exactly one occurrence)
- For cloning: Perfect (one site in backbone, none in insert)
```

**Example 2: Site in both CDS and backbone**
```
CDS: HindIII at position 412 (internal)
Backbone: HindIII at position 185 (cloning site)

Question: Is HindIII "unique"?
- In CDS alone: Yes (one occurrence)
- In backbone alone: Yes (one occurrence)
- In final construct: NO (two occurrences) ← PROBLEM
- For cloning: FAIL (digest creates three fragments)
```

**Example 3: Multiple sites in CDS, none in backbone**
```
CDS: XbaI at positions 234, 567 (two sites)
Backbone: XbaI at position 2276 (cloning site)

Question: Is XbaI "unique"?
- In CDS alone: NO (two occurrences)
- In final construct: NO (three occurrences)
- For cloning: FAIL (multiple fragments)
```

#### Root Cause
Ambiguity in what "unique" means:
- **Mathematically:** Appears exactly once
- **Biologically:** Enzyme cuts exactly once
- **Practically:** Doesn't interfere with cloning strategy

For restriction enzyme cloning:
- **Requirement:** Site must appear EXACTLY ONCE in the final assembled construct
- **Not sufficient:** Site unique in CDS OR backbone separately

#### Fix
**Clear Definition:**
```
A restriction site is "unique" if and only if:
  count(site, final_construct) == 1

Where:
  final_construct = assembled sequence after cloning
                  = backbone_5prime + insert + backbone_3prime
                  (or equivalent based on cloning strategy)
```

**Checkpoint 9 Implementation:**
1. Assemble the full construct conceptually
2. Count sites in the full sequence (both strands)
3. Accept only if count == 1 for each cloning enzyme

#### Prevention
1. **Always specify context** when talking about uniqueness
2. **Default to final construct** for all uniqueness checks
3. **Document cloning strategy** in Phase 0 project spec
4. **Checkpoint 9 is mandatory** for all synthesis designs

#### Philosophical Lesson
This is a broader lesson about **context in biological design**:
- Properties of parts don't simply add
- Interactions matter
- **Systems thinking > component thinking**
- Always ask: "Unique in what context?"

---

### DESIGN-003: Splice Acceptor Avoidance

**Severity:** MEDIUM (AAV-specific)

#### Problem
The AAV2 genome contains a cryptic splice acceptor site at positions 2227-2232 (sequence CCTGCAGG, which is a PstI site). In some contexts, this can cause aberrant splicing of transgene mRNA, leading to truncated protein products.

#### Biological Background
```
AAV2 Rep68 alternative splicing:
- Rep78: Unspliced form, uses p5 promoter
- Rep68: Spliced form, removes intron from 2227-2632

Splice acceptor sequence at 2227:
  ...CCTGCAGG... (this is also a PstI restriction site)

In AAV transfer plasmids (ITR-only constructs):
  - Rep genes are absent
  - Splice acceptor is in backbone (outside ITRs)
  - Generally not a problem for standard constructs

BUT can be activated if:
  - Strong splice donor introduced in transgene
  - Long transgene with GT-AG dinucleotides
  - Certain promoter combinations
```

#### When This Matters
**High risk scenarios:**
1. Full AAV genome constructs with Rep/Cap
2. Transfer plasmids with intron-containing transgenes
3. Very long transgenes (>3 kb) with multiple GT-AG pairs
4. Transgenes with strong splice donors upstream

**Low risk scenarios:**
1. Standard intronless cDNA transgenes
2. Short transgenes (<1.5 kb)
3. No strong splice consensus sequences in transgene

#### Current Agent Approach
1. **Document in exclusion zones** (knowledge_base/exclusion_zones.json)
2. **Note in synthesis workflows** for full AAV genome work
3. **Scan for cryptic splice sites** in transgene (Checkpoint 5)
4. **Flag if strong splice donor + acceptor pair** found

#### Prevention Strategies
**For Transfer Plasmids:**
1. Use intronless cDNA (not genomic DNA)
2. Avoid introducing strong splice donors (GT-AG in good context)
3. If intron needed, use synthetic introns with weak acceptors downstream

**For Full AAV Genome Constructs:**
1. Maintain wild-type sequence around Rep68 splice acceptor
2. If modifying, verify splicing pattern experimentally
3. Consider using Rep-Cap helper plasmid instead (separate splice sites from transgene)

**Tools:**
- MaxEntScan: http://hollywood.mit.edu/burgelab/maxent/
- Human Splicing Finder: http://www.umd.be/HSF3/

#### Test Coverage
Not directly tested (biological rather than computational issue), but:
- Exclusion zone documentation prevents accidental mutation
- Cryptic splice site scanning (Checkpoint 5) flags strong sites

#### Lessons
1. **Historical baggage matters** — AAV2 genome structure affects derivative constructs
2. **Context-specific risks** — Not all transfer plasmids have this issue
3. **Computational prediction has limits** — Splice sites require experimental validation
4. **Document known issues** even if low probability

---

### DESIGN-004: Dam Methylation in Cloning Site Selection

**Severity:** CRITICAL (synthesis and cloning failure)

#### Problem
Using XbaI (TCTAGA) as a cloning site results in a failure mode with standard E. coli due to dam methylase activity. The GATC motif created by XbaI or at ligation junctions is methylated by dam methylase, causing the restriction site to behave abnormally and requiring special dam- E. coli strains.

#### Biological Background
**Dam methylase in E. coli:**
```
Recognition: GATC (adenine methylation)
Function: DNA mismatch repair, replication timing
Prevalence: Present in >99% of lab E. coli strains (dam+)

Dam- strains (special order):
- JM110
- SCS110
- dam-/dcm- strains
- More expensive, slower growth, special handling
```

**XbaI and dam methylation:**
```
XbaI site: TCTAGA
           -------- XbaI recognition
           T-C-T-A-G-A

Issue: The site itself or ligation junctions can create GATC contexts
Effect: Dam methylation of adenine in GATC interferes with XbaI cutting
Result: Incomplete digestion, abnormal cleavage, cloning failures
```

#### Real-World Manifestation (v04 → v05 Revision)

**AAV Transfer Plasmid v04 (FAILED DESIGN):**
- Right cloning site: XbaI (position 2276)
- Issue identified: GATC motif → dam methylation → requires dam- E. coli
- Impact: Increased cost, special strain requirement, reduced reliability
- Status: ❌ DEPRECATED, do not use

**AAV Transfer Plasmid v05 (CORRECTED):**
- Right cloning site: MluI (position 2270, one site proximal)
- MluI recognition: ACGCGT (no GATC motif)
- Works with standard dam+ E. coli (DH5α, TOP10, etc.)
- Status: ✅ Production ready

#### Timeline of Discovery

1. **v04 design completed:** Used XbaI as right cloning site
2. **User identified issue:** "XbaI creates GATC motif → dam methylation problem"
3. **Solution proposed:** Use MluI (one site proximal in polylinker)
4. **v05 generated:** MluI-based design, no dam issues
5. **Lesson documented:** Added to LESSONS_LEARNED.md

#### Impact

**Cloning failures with v04:**
- Incomplete restriction digests
- Multiple bands on gel instead of clean cuts
- Ligation efficiency reduced
- Transformation success rate lowered
- Requires ordering special dam- E. coli strain

**Cost implications:**
- Standard dam+ E. coli: $50-100 per vial, widely available
- Dam- E. coli strains: $200-400 per vial, special order
- Workflow delays: 1-2 weeks for strain procurement
- Reduced reliability: Dam- strains grow slower, more finicky

**Scientific impact:**
- Loss of confidence in automated design tools
- Wasted synthesis costs (~$500-1000 per construct)
- Delayed research timelines (2-4 weeks)

#### Fix Applied

**Enzyme substitution:**
```
v04 (WRONG):  HindIII --- cassette --- XbaI (TCTAGA)
                                       ^^^^ GATC issue

v05 (CORRECT): HindIII --- cassette --- MluI (ACGCGT)
                                        ^^^^ No GATC
```

**Polylinker analysis:**
```
Right polylinker (proximal to distal):
... NcoI - NheI - AvrII - KpnI - MluI - XbaI ...
                                  ^^^^^   ^^^^
                                  v05     v04
                                  ✅      ❌

MluI is one site proximal to XbaI:
- Position: MluI 2270, XbaI 2276 (6 bp difference)
- Recognition: ACGCGT (no GATC)
- Compatibility: Standard E. coli OK
```

#### Prevention Strategies

**1. Enzyme Selection Criteria:**

Always check cloning enzymes for dam/dcm methylation sensitivity:

**Dam-sensitive enzymes (AVOID for standard cloning):**
- XbaI (TCTAGA) - can create GATC contexts
- BclI (TGATCA) - contains GATC
- BamHI (GGATCC) - adjacent to potential GATC
- MboI (GATC) - direct dam site

**Dam-insensitive enzymes (PREFERRED):**
- MluI (ACGCGT) - no GATC ✅
- NotI (GCGGCCGC) - no GATC ✅
- SpeI (ACTAGT) - no GATC ✅
- AscI (GGCGCGCC) - no GATC ✅
- PmeI (GTTTAAAC) - no GATC ✅
- HindIII (AAGCTT) - no GATC ✅
- EcoRI (GAATTC) - no GATC ✅

**2. Cloning Site Analysis Checklist:**

Before finalizing any cloning strategy:
```
☐ Check each cloning site for GATC motifs
☐ Consider both the site itself AND ligation junctions
☐ Verify compatibility with standard dam+ E. coli
☐ Prefer dam-insensitive enzymes when options exist
☐ Document methylation sensitivity in design rationale
```

**3. Polylinker Navigation:**

When a cloning site has dam issues:
```
Strategy: Move to adjacent site in polylinker
Example: XbaI (dam issue) → MluI (one site proximal, no issue)

Check both directions:
- Proximal (toward insert): Often better (shorter backbone)
- Distal (away from insert): Alternative if proximal unavailable
```

**4. Knowledge Base Integration:**

Add to enzyme metadata (knowledge_base/enzyme_metadata.json):
```json
{
  "XbaI": {
    "recognition": "TCTAGA",
    "methylation_sensitivity": {
      "dam": "HIGH",
      "dcm": "NONE",
      "warning": "Can create GATC contexts, requires dam- E. coli",
      "alternative": "MluI, SpeI, or AvrII"
    }
  },
  "MluI": {
    "recognition": "ACGCGT",
    "methylation_sensitivity": {
      "dam": "NONE",
      "dcm": "NONE",
      "note": "Safe for standard dam+ E. coli"
    }
  }
}
```

#### Testing and Validation

**Checkpoint 10: Methylation Sensitivity Check (NEW)**

Add to synthesis workflow:
```
For each cloning site in design:
  1. Extract recognition sequence
  2. Check for GATC motifs (dam)
  3. Check for CCWGG motifs (dcm)
  4. Flag if present
  5. Suggest alternative enzyme
  6. Document in design report
```

**Test coverage:**
```python
def test_dam_methylation_check():
    """Ensure cloning sites don't have dam methylation issues."""
    sites = {
        'XbaI': 'TCTAGA',
        'MluI': 'ACGCGT',
        'BclI': 'TGATCA'
    }

    assert not has_gatc_motif('ACGCGT')  # MluI OK
    assert has_gatc_motif('TGATCA')       # BclI has GATC
    # XbaI is complex, check context-dependent
```

#### Lessons

**Technical:**
1. **Methylation matters** — Dam methylation is ubiquitous in E. coli
2. **GATC is the enemy** — Avoid it in cloning sites
3. **Context-dependent issues** — Some enzymes are borderline, check ligation junctions
4. **Polylinker navigation** — Adjacent sites often solve the problem

**Process:**
1. **Check methylation early** — During enzyme selection, not after synthesis
2. **Prefer standard strains** — Dam+ E. coli is cheaper, faster, more reliable
3. **Document assumptions** — Note why each enzyme was chosen
4. **Validate with community knowledge** — NEB catalogs document methylation sensitivity

**Cost-benefit:**
1. **Prevention is cheap** — Choosing MluI over XbaI costs nothing
2. **Failure is expensive** — Re-synthesis + special strains = $500-1000+ delay
3. **Standard is better** — Using common E. coli strains saves time and money

#### Related Issues

- DESIGN-001: Cloning site conflicts (XbaI uniqueness issue)
- DESIGN-004: Fussy enzyme selection (reliability matters)
- Checkpoint 9: Cloning site uniqueness (extended to include methylation check)

#### References

1. **Dam methylase:** Marinus, M.G. (1987). "DNA methylation in Escherichia coli." Annu Rev Genet.
2. **Restriction enzyme methylation sensitivity:** NEB catalog, "Methylation Sensitivity" section
3. **E. coli strain selection:** Addgene guide, "Choosing the Right Competent Cells"

---

### DESIGN-005: Fussy Enzyme Selection

**Severity:** LOW-MEDIUM (cloning efficiency)

#### Problem
Some restriction enzymes are "fussy" — they have activity issues, require special conditions, or are poorly characterized. Using these enzymes can lead to cloning failures even when the sequence design is correct.

#### Examples of Fussy Enzymes

**BbvCI:**
- **Recognition:** CCTCAGC (7bp, non-palindromic)
- **Type:** Type II (NOT Type IIS, despite similar name to BbsI)
- **Issues:**
  - Less commonly used, may have storage/activity issues
  - Non-palindromic makes it less versatile
  - Sometimes confused as Type IIS (incorrect)
- **Recommendation:** Avoid. Use EcoRI, HindIII, or other reliable Type II enzymes

**PaqCI:**
- **Recognition:** CACCTGC
- **Type:** Type IIS (cuts outside recognition site)
- **Issues:**
  - Creates asymmetric overhangs (4 bases top strand, 8 bases bottom)
  - Less well characterized than BsaI/BsmBI
  - Lower commercial availability
  - Not as reliable for Golden Gate assembly
- **Recommendation:** Use BsaI or BsmBI for standard Golden Gate instead

**SmaI:**
- **Recognition:** CCCGGG (blunt cutter)
- **Issues:**
  - Temperature-sensitive activity
  - CpG methylation in recognition site can affect activity
  - Less reliable than EcoRV for blunt-end cloning
- **Recommendation:** Use EcoRV for more reliable blunt-end cloning

#### BbvCI Misclassification Issue
**Common misconception:** BbvCI is Type IIS (like BbsI)

**Reality:** BbvCI is **Type II**
- Cuts WITHIN its recognition sequence: CC↓TCAGC
- Not suitable for Golden Gate assembly
- Similar naming to BbsI causes confusion

**Correct classification:**
```
Type IIS enzymes (cut OUTSIDE recognition site):
  BsaI, BsmBI, BbsI, Esp3I, SapI, PaqCI

Type II enzymes (cut WITHIN recognition site):
  EcoRI, HindIII, XbaI, BamHI, PstI, BbvCI ← Note: BbvCI here!
```

#### Root Cause
1. **Lack of enzyme metadata** in earlier versions
2. **Over-reliance on "obscure" enzymes** without checking reliability
3. **Confusion between enzyme types** (Type II vs Type IIS)
4. **Availability bias** (enzyme mentioned in paper doesn't mean it's good)

#### Fix
**Implemented enzyme metadata database:**
```
knowledge_base/enzyme_metadata.json

Contains:
- Enzyme type classification (Type II vs Type IIS)
- Reliability ratings (very_high, high, medium, low)
- Methylation sensitivity
- Star activity warnings
- Recommended alternatives
```

**Enzyme Selection Criteria:**
1. **Prefer highly reliable enzymes:**
   - Type II: EcoRI, HindIII, BamHI, XbaI, NotI
   - Type IIS: BsaI, BsmBI/Esp3I

2. **Check reliability rating** in enzyme_metadata.json

3. **Avoid fussy enzymes:**
   - BbvCI (use EcoRI, HindIII instead)
   - PaqCI (use BsaI, BsmBI instead)
   - Obscure enzymes without commercial support

4. **For Golden Gate:**
   - Standard: BsaI (gold standard)
   - Alternative: BsmBI/Esp3I (identical to BsmBI)
   - Avoid: BbsI (present in pGS-AAV backbones), PaqCI (fussy)

#### Prevention
1. **Consult enzyme metadata** before choosing enzymes
2. **Stick to well-established enzymes** unless specific reason
3. **Check commercial availability** (NEB, Thermo, etc.)
4. **Consider enzyme source** (HF variants often more reliable)
5. **Test empirically** if using non-standard enzyme

#### Test Coverage
Not directly tested (empirical issue), but:
- Enzyme metadata database documents reliability
- Warnings in synthesis workflows

#### Lessons
1. **Reliability trumps cleverness** — Use boring, proven enzymes
2. **Type classification matters** — Know if it's Type II or Type IIS
3. **Commercial support indicates reliability** — If NEB sells it, it probably works
4. **Fussy enzymes cost time** — Failed cloning wastes weeks

---

## General Principles

### The Importance of Testing

**Lesson:** Bugs that seem "obvious" in retrospect were invisible during development.

**Why Testing Matters:**
1. **Catches edge cases** that manual inspection misses
2. **Prevents regression** when code changes
3. **Documents expected behavior** as executable specification
4. **Builds confidence** in automated tools

**Test Suite Coverage (v3.0):**
- `test_frame_offset.py` — 15 tests for BUG-001/BUG-002
- `test_uniqueness_counting.py` — 20 tests for BUG-003/DESIGN-001/DESIGN-002
- `test_silent_classification.py` — 25 tests for Checkpoint 8
- `conftest.py` — Reusable fixtures and test data

**Testing Philosophy:**
```
Every bug should generate:
1. A test that FAILS with the bug present
2. A fix that makes the test PASS
3. Documentation of the bug (this file)
4. Prevention strategy (checkpoints, knowledge base)
```

**Future Testing:**
- Integration tests (full workflow end-to-end)
- Regression suite (run before each release)
- Property-based testing (generate random valid sequences)
- Biological validation (wet lab testing of designs)

---

### Double-Checking Biological Assumptions

**Lesson:** Computational biology requires constant verification of biological assumptions against reality.

**Examples of Assumptions That Were Wrong:**

1. **"Frame is always 0"** — No, depends on CDS start position
2. **"One occurrence means unique"** — No, need to check both strands
3. **"Sites in CDS don't matter"** — No, must consider final construct
4. **"BbvCI is Type IIS"** — No, it's Type II

**Verification Strategies:**

1. **Check references:**
   - UniProt for protein sequences
   - NCBI for nucleotide sequences
   - NEB for enzyme specifications
   - Primary literature for regulatory elements

2. **BLAST everything:**
   - Whole sequences
   - Individual elements
   - ORFs
   - Regulatory regions

3. **Cross-validate:**
   - Multiple databases
   - Multiple tools
   - Experimental data (if available)

4. **Know the limits:**
   - Computational predictions (splice sites, secondary structure)
   - Database annotation quality
   - Tool accuracy

**Culture of Skepticism:**
```
Question everything:
- Annotation labels (may be wrong)
- File format assumptions (may have edge cases)
- Enzyme classifications (may be outdated)
- "Common knowledge" (may be incomplete)

Verify independently:
- Look up primary sources
- Run your own analyses
- Test empirically when possible
```

---

### Comprehensive Checkpoints

**Lesson:** Ad-hoc checks are insufficient. Systematic, comprehensive checkpoints catch issues before they become expensive failures.

**Evolution of Checkpoint System:**

**v1.0 (Ad-hoc):**
- Manual checklist
- Easy to forget items
- No systematic verification
- Result: Frequent errors

**v2.0 (Structured):**
- Explicit verification checklist
- Mandatory checklist loop
- BLAST-based verification
- Result: Fewer errors, but still had BUG-001, BUG-003, DESIGN-001

**v3.0 (Comprehensive):**
- All v2.0 features
- Checkpoint 8: Silent Mutation Verification
- Checkpoint 9: Cloning Site Uniqueness
- Knowledge base integration
- Automated test suite
- Result: Robust design workflow

**Key Insight:**
```
Checkpoints should be:
1. Systematic (run every time, not optional)
2. Comprehensive (cover all failure modes)
3. Automated (reduce human error)
4. Well-defined (clear pass/fail criteria)
5. Tested (checkpoints themselves need tests)
```

**Checkpoint Design Principles:**

1. **Each checkpoint tests ONE thing clearly:**
   - Checkpoint 8: Are mutations silent?
   - Checkpoint 9: Are cloning sites unique?

2. **Clear acceptance criteria:**
   - ✅ PASS: All mutations are silent
   - ❌ FAIL: Any non-silent mutation
   - No ambiguity

3. **Actionable failures:**
   - Not just "FAIL"
   - Specific problem identified
   - Suggested fix provided
   - Re-verification path clear

4. **Minimal false positives:**
   - Don't flag non-issues
   - Understand biological context
   - Know when to warn vs fail

5. **Integration with workflow:**
   - Checkpoints run at right time
   - Output feeds into next phase
   - Results documented in report

**Future Checkpoints (Proposed):**

- **Checkpoint 10:** Secondary structure prediction (RNA stability)
- **Checkpoint 11:** Codon adaptation index (expression optimization)
- **Checkpoint 12:** Immunogenicity scan (avoid T-cell epitopes)
- **Checkpoint 13:** Off-target integration sites (human genome BLAST)

---

## Conclusion

### Summary of Lessons

**Technical Lessons:**
1. Always calculate frame offsets relative to CDS start (BUG-001)
2. Check both DNA strands for restriction sites (BUG-003)
3. Verify uniqueness in final construct, not components (DESIGN-001, DESIGN-002)
4. Stick to reliable, well-characterized enzymes (DESIGN-004)
5. Understand enzyme types (Type II vs Type IIS)

**Process Lessons:**
1. Comprehensive testing catches bugs before deployment
2. Documentation prevents knowledge loss
3. Systematic checkpoints beat ad-hoc verification
4. Biological assumptions must be continuously validated
5. Context matters in every design decision

**Cultural Lessons:**
1. Question everything, verify independently
2. Bugs are learning opportunities, document them
3. Automation is better than manual checking
4. Simple, boring solutions beat clever, fragile ones
5. Systems thinking > component thinking

### Continuous Improvement

This document should be:
- **Updated** with every new bug discovered
- **Referenced** during code reviews and design reviews
- **Tested** (do the tests still catch the bugs?)
- **Taught** to new users of the agent

**When adding a new lesson:**
```
1. Document the problem clearly
2. Explain root cause (why it happened)
3. Describe impact (what broke)
4. Show the fix (how to solve it)
5. Provide prevention (how to avoid it)
6. Add test coverage (how to detect it)
7. Extract principles (what to remember)
```

### Acknowledgments

Lessons learned from:
- Real synthesis failures (AAV v04 XbaI issue)
- Code review discoveries (frame offset bug)
- Test-driven development (silent mutation verification)
- Literature review (splice acceptor, enzyme reliability)
- User feedback (fussy enzyme frustration)

---

**Document Version:** 1.2
**Last Updated:** 2026-01-22 (Added BUG-006: N-Terminal Insertion After Single Amino Acid)
**Next Review:** After each major release or significant bug discovery
**Maintainer:** DNA Engineer Agent development team

**END OF LESSONS LEARNED**
