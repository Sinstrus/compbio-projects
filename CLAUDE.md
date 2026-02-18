# Workspace Structure

## DNA Construct Engineering

### Construct Bank
All DNA constructs (.gb, .dna files) are stored in a single location:
`dna_engineer_agent/constructs/`

The full inventory with narrative descriptions is in:
`dna_engineer_agent/CONSTRUCT_REGISTRY.md`

### Toolset
`dna_engineer_agent/` is a shared toolset — NOT a project directory.
Key tools: construct_verifier.py, design_reasoner.py, splice_donor_hbonds.py
Reference data: knowledge_base/, parts_library/, backbones/, cis_elements/
See `dna_engineer_agent/CLAUDE.md` for full toolset documentation.

### Rules for New DNA Constructs
When building new DNA constructs:
1. Output ALL construct files (.gb, .dna) to `dna_engineer_agent/constructs/`
2. Use the next available AVD number (check CONSTRUCT_REGISTRY.md for the latest)
3. Name files as: AVD{NNN}-{backbone}-{modification}-{site}.gb
4. Update CONSTRUCT_REGISTRY.md with ID, filename, and 2-3 sentence description
5. Use construct_verifier.py to verify before finalizing (see dna_engineer_agent/CLAUDE.md)
6. NEVER store construct files inside project directories

### Projects
Individual projects live directly under ~/projects/:
- AVD_VHH_Display_ALPL/ — Anti-ALPL/TfR1 VHH capsid display
- splice_donor_pipeline/ — Splice donor strength characterization
- aav_transfer_plasmid/ — Transfer plasmid engineering
- aav_vhh_v1/ — AAV VHH v1 (early prototype)
- nocap_restriction_sites/ — NoCap restriction site analysis
- executive_assistant/ — Meeting/task assistant
- ppt_maker/ — Voyager slide builder
