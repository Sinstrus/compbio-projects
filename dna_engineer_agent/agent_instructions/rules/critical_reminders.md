# Critical Reminders

> Part of: [Agent Instructions](../README.md)

1. **Never Skip Verification** — Even if the user only asks about one element, verify ALL elements required by the system. Use the checklist.

2. **References Are Truth** — Alignment results override file annotations. Report discrepancies clearly.

3. **Be Explicit About Uncertainty** — If a reference is unavailable, mark as UNVERIFIED and document why.

4. **Follow Acceptance Criteria Rigorously** — Don't rationalize borderline results. >=85% = verified, <85% = flag it.

5. **Generate Actionable Reports** — Don't just list problems — suggest solutions with specific coordinates and proposed fixes.

6. **Frame Offset Calculation (BUG-001)** — ALWAYS calculate relative to CDS start. WRONG: `position % 3`. RIGHT: `(position - cds_start) % 3`. See: `scripts/tools/tests/test_frame_offset.py`

7. **Double-Strand Uniqueness (BUG-003)** — For palindromic restriction sites, count on both strands. See: `scripts/tools/tests/test_uniqueness_counting.py`

8. **Exclusion Zone Awareness** — Load and respect exclusion zones from `knowledge_base/exclusion_zones.json`. Never propose mutations in critical regions (ITR boundaries, Kozak, polyA signal, promoter core).

9. **Type IIS vs Type II Classification** — BbvCI is Type II (NOT Type IIS). PaqCI IS Type IIS. See: `knowledge_base/enzyme_metadata.json`

10. **Checkpoint 8 & 9 Integration** — These often run in a loop: CP9 finds non-unique sites → design silent mutations → CP8 verifies silence → apply → re-run CP9 → repeat until both pass.

11. **Checkpoint 10: Parent-Child (BUG-005)** — ALWAYS run after building any construct from a parent plasmid. Verifies only intentional changes exist.

12. **Multi-Construct Builds (BUG-005)** — Use each parent's OWN sequence. WRONG: Extract VP1 from plasmid A for both constructs. RIGHT: Extract VP1 from each construct's own parent.

13. **Synthetic Fragment Boundaries (DESIGN-005)** — Fragment boundaries must FLANK the modified region, not be inside it. Find closest unique sites OUTSIDE the region of change.
