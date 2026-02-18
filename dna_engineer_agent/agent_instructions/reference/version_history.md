# Version History

> Part of: [Agent Instructions](../README.md)

## v4.1 (2026-01-20)
- Added CRITICAL RULE: Amino acid insertion point calculation (prevents BUG-004)
- Mandatory programmatic DNA position calculation from amino acid coordinates
- Flanking sequence verification requirements for protein insertions
- Documented BUG-004: Incorrect amino acid insertion point (off-by-one-codon error)

## v4.0 (2026-01-14)
- Added Checkpoint 10: Parent-Child Sequence Verification
- Added critical rule for multi-construct builds (use each parent's own sequence)
- Added DESIGN-005: Synthetic fragment boundary selection guidelines
- Added DESIGN-006: VP1/VP2/VP3 knockout design patterns
- Added DESIGN-007: VHH insertion documentation
- Documented BUG-005: Parent sequence mismatch in multi-construct builds

## v3.0 (2025-12-19)
- Added Phase 0.2: Synthesis workflow selection
- Added Checkpoint 8: Silent Mutation Verification
- Added Checkpoint 9: Cloning Site Uniqueness Verification
- Integrated knowledge base files (workflows, exclusion zones, enzyme metadata)
- Added comprehensive test suite
- Documented lessons learned and bug prevention strategies

## v2.0
- Goal-driven requirements derivation
- Mandatory checklist execution loop
- BLAST-based verification workflow
- Structural rule validation

## v1.0 (archived)
- Basic annotation verification
- Manual checklist approach
