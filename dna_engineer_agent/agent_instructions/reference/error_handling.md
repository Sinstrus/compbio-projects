# Error Handling and Edge Cases

> Part of: [Agent Instructions](../README.md)

## Unknown Systems
1. Ask the user directly
2. Provide a list of supported systems
3. Offer to analyze as "generic plasmid" (reduced verification)

## Missing References
1. Mark element as UNVERIFIED
2. Document the limitation clearly
3. Suggest where the user can find the reference
4. Continue with other verifications

## Conflicting Annotations
1. Trust the BLAST alignment
2. Report the discrepancy explicitly
3. Provide both sets of coordinates
4. Flag for user review

## Borderline Identity Scores (65-75%)
1. Report as PARTIAL
2. Provide alignment details
3. Note whether this could be a serotype variant or true mismatch
4. Suggest targeted Sanger sequencing if experimental plasmid

## Checkpoint Failures During Design
1. DO NOT proceed with synthesis
2. Document the specific failure (which site, which codon, why)
3. Propose specific mutations to fix
4. Re-verify after proposed fixes
5. Iterate until checkpoints pass
