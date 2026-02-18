# Synthetic Fragment Boundary Selection (DESIGN-005)

> Part of: [Agent Instructions](../README.md)

**Purpose:** Select optimal restriction sites for synthetic fragment ordering.

**Critical Rule:** Fragment boundaries must use sites that FLANK the modified region, NOT sites within it.

## Procedure

1. **Map the modified region:**
   - Compare parent vs child sequence
   - Find first and last nucleotide that differs
   - This defines the "modified region"

2. **Identify all restriction sites:**
   - Map unique sites in the final construct
   - Classify each as INSIDE or OUTSIDE the modified region

3. **Select flanking sites:**
   - UPSTREAM: Closest unique site BEFORE modified region
   - DOWNSTREAM: Closest unique site AFTER modified region

4. **Verify:**
   - Fragment includes ALL changes between parent and child
   - Both restriction sites are OUTSIDE the modified region
   - Sites are unique in the construct (appear exactly once)

## Example

```
Modified region: bp 1520-4148 (VP1-VHH in AVD005)

Sites INSIDE (cannot use as boundaries):
- BsrGI at 3353 — 1833 bp INTO VP1-VHH
- BmtI at 3510 — 1990 bp INTO VP1-VHH

Sites OUTSIDE (can use as boundaries):
- AvrII at 1676 — 156 bp BEFORE VP1 start
- BsrGI at 3358 — AFTER VP1-VHH end

Selected fragment: AvrII (1676) to BsrGI (3358) = 1,683 bp
```
