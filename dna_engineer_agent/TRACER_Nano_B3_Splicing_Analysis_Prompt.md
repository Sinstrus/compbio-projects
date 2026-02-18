# TRACER-Nano Idea B3: Splicing-Based VHH Display Feasibility Analysis

## Project Context

This analysis supports the TRACER-Nano platform at Voyager Therapeutics. The goal is to display VHH nanobodies on AAV capsids for receptor-mediated tropism redirection. For library screening applications, genotype-phenotype linkage is essential — each capsid must display only the VHH encoded by its packaged genome.

The "B3" concept proposes using inefficient splicing (intron retention) to generate a mixed population of VHH-VP fusion proteins (~10-20%) and normal VP proteins (~80-90%) from a SINGLE genetic construct.

---

## The Proposed Architecture

```
                    ┌─────────────────── INTRON ───────────────────┐
                    │                                              │
5'───[Exon1]───GT......ATG─VHH─(GGGGS)₅─X─branch─polyY......AG───[VP2 coding]───3'
               │                                              │
               splice                                         splice
               donor                                          acceptor
```

### When SPLICED (target: 80-90% of transcripts)

- Intron is removed entirely (ATG, VHH, linker, X region all excised)
- Mature mRNA: `[Exon1]───[VP2 coding]`
- Translation initiates at native VP2 ACG start codon
- Product: **Normal VP2 protein**

### When RETAINED (target: 10-20% of transcripts)

- Intron remains in mature mRNA
- The ATG inside the "intron" provides a strong translation initiation site
- Product: **Met-VHH-(GGGGS)₅-X-VP2 fusion protein**

### Why This Matters for Libraries

| Approach | Genotype-Phenotype Linkage | Stoichiometry Control |
|----------|---------------------------|----------------------|
| Mosaic (two plasmids) | ❌ Broken - capsid can package either plasmid | ✅ Via DNA ratio |
| Splicing (B3, single construct) | ✅ Maintained - one plasmid per cell | ✅ Via splice efficiency |

---

## Analysis Tasks

### Task 1: Splice Site Sequence Constraints

Analyze the consensus sequences for mammalian splice sites and determine the minimum requirements for weak (~10-20% retention) splicing.

#### 1.1 Splice Donor (5' splice site)

Standard consensus: `MAG|GURAGU` (where `|` is the splice junction, M=A/C, R=A/G)

Questions to answer:
- What are the absolutely required nucleotides (the GU dinucleotide)?
- What positions can be mutated to WEAKEN splicing?
- What does the exon side (upstream of junction) need to look like?

Provide a table of splice donor strength vs sequence.

#### 1.2 Splice Acceptor (3' splice site)

Standard consensus: `(Y)ₙ-N-C-A-G|G` where:
- `(Y)ₙ` = polypyrimidine tract (C and U/T only)
- Branch point `YUNAY` is typically 18-40 nt upstream
- `AG` dinucleotide is essential
- First nucleotide of downstream exon is preferentially G

Questions to answer:
- What is the MINIMUM polypyrimidine tract length for:
  - Strong splicing (>95% spliced)?
  - Medium splicing (~50% spliced)?
  - Weak splicing (~80-90% spliced, 10-20% retained)?
- How does branch point sequence affect retention rate?
- Can we design for ~10-20% retention specifically?

#### 1.3 Tools and References

Use or reference:
- MaxEntScan (Burge lab) for splice site scoring
- Human Splicing Finder
- SpliceAI predictions
- Literature: Mount SM 1982, Burge & Karlin 1997

---

### Task 2: Translation of Splice Site Sequences When Retained

When the intron is RETAINED, the splice site sequences become part of the mRNA and will be translated. This creates amino acid constraints.

#### 2.1 Splice Donor Translation

If we place the splice donor immediately after the (GGGGS)₅ linker:

```
...(GGGGS)₅───[splice_donor_sequence]───VHH_continues_or_ends...
```

Wait - I need to reconsider the architecture. The splice donor should be BEFORE the VHH (at the 5' end of the intron), not after. Let me reconsider:

```
[Exon1]──GT─ATG─VHH─(GGGGS)₅─X─AG──[Exon2/VP2]
         │                      │
         donor                  acceptor
```

So the splice donor consensus `GT` comes immediately after the upstream exon. When retained, this `GT` is translated as part of the reading frame.

Question: What codon context does `GT` create? Depends on the upstream exon's last nucleotide.

#### 2.2 Polypyrimidine Tract Translation (Critical Constraint)

The "X" region must contain a polypyrimidine tract (only C and T/U nucleotides). When retained and translated, what amino acids can these encode?

Create a complete table:

| Codon | Nucleotides | Amino Acid | Acceptable in Linker? |
|-------|-------------|------------|----------------------|
| TTT | pyrimidine only | Phe | ⚠️ Bulky, hydrophobic |
| TTC | pyrimidine only | Phe | ⚠️ Bulky, hydrophobic |
| TCT | pyrimidine only | Ser | ✅ Small, flexible |
| TCC | pyrimidine only | Ser | ✅ Small, flexible |
| TTA | NOT pure pyrimidine | - | N/A |
| TTG | NOT pure pyrimidine | - | N/A |
| CTT | pyrimidine only | Leu | ⚠️ Hydrophobic |
| CTC | pyrimidine only | Leu | ⚠️ Hydrophobic |
| CTA | NOT pure pyrimidine | - | N/A |
| CTG | NOT pure pyrimidine | - | N/A |
| CCT | pyrimidine only | Pro | ❌ Causes kinks |
| CCC | pyrimidine only | Pro | ❌ Causes kinks |
| CCA | NOT pure pyrimidine | - | N/A |
| CCG | NOT pure pyrimidine | - | N/A |

Complete this analysis for ALL codons and identify:
1. Which amino acids are ONLY encodable by pure pyrimidine codons?
2. Which amino acids are IMPOSSIBLE to encode with pure pyrimidines?
3. What is the "best case" linker sequence from a polyY tract?

#### 2.3 The CAG|G Junction

The 3' splice site ends with `CAG|G` (where `|` is the junction, G is first nt of exon).

When retained, `CAG` encodes Glutamine (Q). 
- Is Gln acceptable in a linker? (Yes, polar, flexible)
- Does the downstream G create issues? (Depends on next two nucleotides)

#### 2.4 Branch Point Translation

The branch point consensus is `YUNAY` (typically `CURAY` with the A being the branch adenosine).
- This is NOT pure pyrimidine
- It breaks the polyY tract
- What amino acids does `CURAY` encode in different frames?

---

### Task 3: Reading Frame Preservation

This is the critical feasibility question. Can we maintain the reading frame across the entire construct when the intron is retained?

#### 3.1 Component Lengths

Calculate and verify:

| Component | Amino Acids | Nucleotides | Divisible by 3? |
|-----------|-------------|-------------|-----------------|
| Start codon (ATG) | 1 (Met) | 3 | ✅ |
| VHH3 (anti-ALPL) | 118 | 354 | ✅ |
| (GGGGS)₅ linker | 25 | 75 | ✅ |
| X region (splice acceptor) | ? | ? | Must be ✅ |
| **Total before VP2** | 144 + X/3 | 432 + X | X mod 3 = 0 |

#### 3.2 Frame Constraint Equation

For the fusion protein to be in-frame with VP2:

```
(upstream_exon_contribution) + (ATG) + (VHH) + (GGGGS)₅ + (X) ≡ 0 (mod 3)
```

If VHH = 354 bp and (GGGGS)₅ = 75 bp, then:
- 354 + 75 = 429 bp = 143 codons ✅

So the ATG-VHH-linker portion is already frame-neutral (divisible by 3).

**Therefore: X must also be divisible by 3 to maintain frame.**

#### 3.3 Minimum X Length Analysis

What is the minimum biologically functional X region?

Components needed in X:
1. Branch point sequence: ~7 nt (`YUNAY` core, may need flanking)
2. Spacer between branch point and polyY: variable (can be 0)
3. Polypyrimidine tract: minimum ~8-15 nt for weak splicing
4. Spacer before AG: typically 0-3 nt
5. AG dinucleotide: 2 nt

Minimum estimate: 7 + 0 + 12 + 0 + 2 = **21 nt minimum**

But we need X mod 3 = 0, so minimum is **21 nt (7 codons)** or **24 nt (8 codons)**.

#### 3.4 Design X to Be Frame-Compatible

Design the shortest X sequence that:
- Contains functional (but weak) splice acceptor
- Is exactly divisible by 3
- Translates to acceptable amino acids

---

### Task 4: Detailed X Region Design

#### 4.1 Design Constraints Summary

| Constraint | Requirement | Flexibility |
|------------|-------------|-------------|
| Branch point | YUNAY motif, A is branch point | Position 18-40 nt upstream of AG |
| PolyY tract | Minimum ~10 nt of C/T only | Can be longer for stronger splicing |
| 3' splice site | AG dinucleotide essential | Must be exactly AG |
| Frame | Total length mod 3 = 0 | Can pad with acceptable codons |
| Amino acids | No stops, minimal Pro | Prefer Ser, Phe, Leu |

#### 4.2 Proposed X Sequence (First Attempt)

Design a 24 nt (8 codon) X region:

```
Position:  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24
Sequence:  C   T   A   A   C   T   T   C   T   T   T   C   T   T   T   C   C   T   T   C   C   A   G   G
           └─────branch─────┘   └──────────polypyrimidine tract──────────┘           └─AG─┘└─exon
           
Codons:    |  CTA  |  ACT  |  TCT  |  TTC  |  TTT  |  CCT  |  TCC  |  AGG  |
AAs:       |  Leu  |  Thr  |  Ser  |  Phe  |  Phe  |  Pro  |  Ser  |  Arg  |
                                                      ❌Pro
```

Problem: CCT encodes Pro (causes kinks). Need to redesign.

#### 4.3 Iterate on X Design

Try alternative sequences that avoid Pro codons (CCN) while maintaining splice function.

The challenge: polypyrimidine tracts that are pure C/T will inevitably encode some combination of:
- Phe (TTT, TTC)
- Leu (CTT, CTC) 
- Ser (TCT, TCC)
- Pro (CCT, CCC) ← want to avoid

To avoid Pro, we need to avoid `CC` dinucleotides followed by `T` or `C`.

Design strategy: Use `TCT`, `TTC`, `TTT`, `CTT` but not `CCT` or `CCC`.

```
Revised polyY: TCTTTCTTTCTTCC
               Ser-Phe-Phe-Leu-?
```

Continue this analysis to find an optimal X sequence.

#### 4.4 Final X Sequence Recommendation

Provide:
- Complete nucleotide sequence
- Codon breakdown
- Amino acid translation
- Predicted splice efficiency (use MaxEntScan or similar)
- Any caveats or risks

---

### Task 5: VP1 vs VP2 Considerations

The B3 architecture could be applied to either VP1 or VP2 N-terminus. Analyze differences:

#### 5.1 VP2 N-terminal Insertion

- VP2 starts at ACG (weak start codon, encodes Thr)
- Native VP2 N-terminus: M-A-A-D-G-Y-L-P-D-W...
- After B3: M-VHH-(GGGGS)₅-X-A-A-D-G-Y-L-P-D-W... (when retained)

Question: Does the ACG→ATG(Met-VHH...) override the native VP2 start?

#### 5.2 VP1 N-terminal Insertion

- VP1 starts at ATG (strong start codon)
- Native VP1 N-terminus: M-A-A-D-G-Y-L-P-D-W...
- Would need to insert the "intron" after the native ATG

Question: How does this affect VP2 and VP3 expression from the same transcript?

#### 5.3 Recommendation

Which VP is better suited for the B3 approach and why?

---

### Task 6: Failure Mode Analysis

#### 6.1 Nonsense-Mediated Decay (NMD)

NMD degrades transcripts with premature stop codons. Could intron retention trigger NMD?

Check:
- Does the retained "intron" introduce any in-frame stop codons?
- Is there an exon-exon junction >50 nt downstream of any stop?
- What is the predicted NMD susceptibility?

#### 6.2 Cryptic Splice Sites in VHH

The VHH sequence (~350 bp) may contain cryptic GT...AG pairs that could cause mis-splicing.

Scan VHH3 (anti-ALPL) sequence for:
- Cryptic splice donors (GT followed by splice-like context)
- Cryptic splice acceptors (polypyrimidine + AG)
- Use MaxEntScan to score any potential sites

If cryptic sites exist, can they be removed by silent mutations?

#### 6.3 Start Codon Competition

Multiple ATGs in a transcript can compete for ribosome initiation.

Analyze:
- Is the ATG in the retained intron in good Kozak context? (want: GCC**ATG**G)
- Are there upstream ATGs that could interfere?
- Could leaky scanning bypass the intended start?

#### 6.4 Alternative Splicing Outcomes

What if splicing occurs at unexpected sites?

Map all possible splice outcomes:
1. Intended: Full intron removal → Normal VP2
2. Intended: Full intron retention → VHH-VP2 fusion
3. Unintended: Partial splicing → Truncated products?
4. Unintended: Cryptic site usage → Wrong junctions?

---

### Task 7: Experimental Validation Plan

If the design appears feasible, outline experiments to test it:

#### 7.1 Splice Reporter Assay

Design a minigene reporter to test splicing in isolation:
- Exon1 - [test intron] - Exon2-GFP
- Measure GFP expression (requires correct splicing)
- RT-PCR to quantify spliced vs retained transcripts

#### 7.2 Western Blot Validation

Expected bands for VP2 constructs:
- Normal VP2: ~65 kDa
- VHH-VP2 fusion: ~65 + 15 = ~80 kDa

If B3 works, should see BOTH bands from a single construct, with ratio ~80:20.

#### 7.3 Functional Transduction Assay

Test if VHH-VP2 capsids retain:
- Capsid assembly capability
- VHH binding function (ELISA with ALPL)
- Transduction of ALPL+ cells

---

## Deliverables

Please provide:

1. **Splice site constraint analysis** - Table of requirements for weak vs strong splicing

2. **Polypyrimidine codon table** - Complete analysis of what amino acids can/cannot be encoded

3. **Optimized X region sequence** - Nucleotide sequence with:
   - Annotation of functional elements
   - Codon/amino acid breakdown
   - Predicted splice efficiency score

4. **Frame verification calculation** - Explicit math showing frame is preserved

5. **Risk matrix** - Table of potential failure modes with likelihood and mitigation

6. **GO/NO-GO recommendation** - Is B3 worth pursuing experimentally?

---

## Success Criteria

The B3 design is considered FEASIBLE if:

- [ ] X region can be ≤30 amino acids while maintaining splice acceptor function
- [ ] All translated amino acids in X are acceptable for a linker (no Pro, no stops)
- [ ] Reading frame is preserved (verified by in silico translation)
- [ ] No high-confidence cryptic splice sites in VHH sequence
- [ ] No obvious NMD triggers
- [ ] Predicted retention rate is tunable in the 5-30% range
- [ ] Design is compatible with both VP1 and VP2 insertion

---

## Reference Information

### AAV9 Cap Gene Coordinates (for context)

- VP1 start (ATG): position 2379 (1-indexed)
- VP2 start (ACG): position 2790 (1-indexed)  
- VP3 start (ATG): position 2986 (1-indexed)

### Key Literature

- Mount SM (1982) NAR - Splice site consensus sequences
- Burge & Karlin (1997) - MaxEntScan basis
- Braunschweig et al (2014) Genome Research - Intron retention
- PMC12385487 - Referenced in TRACER-Nano slides for IR patterns

### VHH3 Sequence

If available in project files, use the actual VHH3 (anti-ALPL) sequence. Otherwise, use a representative VHH sequence (~118 aa).

---

## Notes for the Agent

- This is exploratory/feasibility analysis, not final design
- Be explicit about uncertainties and assumptions
- If the concept is fundamentally flawed, say so clearly
- Suggest alternative approaches if B3 appears infeasible
- Consider that ~10-20% retention is the TARGET, but 5-30% would be acceptable
