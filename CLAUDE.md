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
3. **AVD numbers MUST be unique to a single DNA sequence. Never assign an existing AVD number to a different sequence. Skipping numbers is fine; duplication is not. Multiple file formats or annotation variants of the same sequence share one AVD number.**
4. Name files as: AVD{NNN}-{backbone}-{modification}-{site}.gb
5. Update CONSTRUCT_REGISTRY.md with ID, filename, and 2-3 sentence description
6. Use construct_verifier.py to verify before finalizing (see dna_engineer_agent/CLAUDE.md)
7. NEVER store construct files inside project directories

### Projects
Individual projects live directly under ~/projects/:
- AVD_VHH_Display_ALPL/ — Anti-ALPL/TfR1 VHH capsid display
- splice_donor_pipeline/ — Splice donor strength characterization
- aav_transfer_plasmid/ — Transfer plasmid engineering
- aav_vhh_v1/ — AAV VHH v1 (early prototype)
- nocap_restriction_sites/ — NoCap restriction site analysis
- executive_assistant/ — Meeting/task assistant
- ppt_maker/ — Voyager slide builder (DEPRECATED — decks migrated to voyager_slide_templates)
- voyager_slide_templates/ — Voyager PPTX-compatible slide authoring

### Slide Decks (PPPxxx Barcoding)
All slide decks are tracked by PPP number and stored in `voyager_slide_templates/decks/PPPxxx-name/`.
Registry: `voyager_slide_templates/PRESENTATION_REGISTRY.md`
Authoritative guide: `voyager_slide_templates/CLAUDE.md`

Rules for new decks:
1. Use the next available PPP number (check PRESENTATION_REGISTRY.md)
2. **PPP numbers are unique and permanent.** Never reassign a PPP number to a different deck.
3. Name folders as: `PPPxxx-descriptive-name/`
4. Update PRESENTATION_REGISTRY.md with PPP number, name, slide count, date, and status
5. Build: `cd voyager_slide_templates/decks/PPPxxx-name/ && npm install && node build.js`

Key conventions:
- Slides are bare `<body>` HTML with inline styles, exported via html2pptx pipeline
- All text must be in `<p>`, `<h1>`-`<h6>`, `<ul>`, or `<ol>` tags — bare `<span>` elements are silently dropped
- Every CSS gradient needs `data-rasterize="gradient"`; branding via `<img>` tags only
- Shared html2pptx library at `voyager_slide_templates/html2pptx/`
- Shared Voyager branding PNGs at `voyager_slide_templates/assets/`
- Legacy decks in `ppt_maker/` are deprecated — all decks now live in voyager_slide_templates
