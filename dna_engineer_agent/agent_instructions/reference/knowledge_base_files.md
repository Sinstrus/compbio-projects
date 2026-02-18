# Knowledge Base Files

> Part of: [Agent Instructions](../README.md)

## Required Files
1. **systems/{system_id}.json** — System requirements and rules
2. **references/{element_id}.fasta** — Reference sequences for alignment
3. **backbones/genscript/BACKBONE_CATALOG.json** — GenScript backbone specifications
4. **synthesis/genscript_workflows.json** — Cloning methods and workflows
5. **exclusion_zones.json** — Regions where mutations are restricted
6. **enzyme_metadata.json** — Restriction enzyme properties and reliability

## Optional Files
7. **genetic_code.json** — Codon table for translation verification
8. **serotype_variants.json** — Known ITR/Cap variants by serotype
