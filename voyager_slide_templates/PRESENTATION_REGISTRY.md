# Presentation Registry

All slide decks are stored in `voyager_slide_templates/decks/` and indexed by PPP number.
Each PPPxxx folder is self-contained with HTML slides, a `build.js`, and a `package.json`.

## Registry

| PPP | Folder | Name | Slides | Date | Parent | Status | Notes |
|-----|--------|------|--------|------|--------|--------|-------|
| PPP001 | PPP001-aav-vhh-poc | AAV VHH Display POC | 13 | 2025-01 | — | archived | First POC deck, no Voyager branding |
| PPP002 | PPP002-aav-vhh-voyager | AAV VHH Display Voyager | 13 | 2025-01 | PPP001 | archived | Voyager branding applied |
| PPP003 | PPP003-aav-vhh-voyager-v2 | AAV VHH Display Voyager v2 | 13 | 2025-01 | PPP002 | archived | Updated version |
| PPP004 | PPP004-tracer-b3-feasibility | TRACER-Nano B3 Feasibility | 13 | 2025-01 | — | active | B3 insertion site feasibility analysis |
| PPP005 | PPP005-tracer-vr4-splicing | TRACER-Nano VR4 Splicing | 12 | 2025-01 | — | archived | Superseded by PPP011 |
| PPP006 | PPP006-intron-retention-library | VHH3 Intron Retention Library | 14 | 2025-02 | — | archived | Original intron retention library deck |
| PPP007 | PPP007-intron-retention-voyager | VHH3 Intron Retention Voyager | 14 | 2025-02 | PPP006 | archived | Voyager-branded version of PPP006 |
| PPP008 | PPP008-first-principles-vr4 | First-Principles VR4 Donor Library | 10 | 2025-02 | — | active | AVD100-131 synthetic splice donor constructs |
| PPP009 | PPP009-alpl-vhh-display | ALPL VHH Capsid Display | 13 | 2025-01 | — | active | Anti-ALPL VHH display on AAV9 |
| PPP010 | PPP010-nterm-intron-retention | N-Term Intron Retention VP1/VP2 | 10 | 2026-02 | — | active | N-terminal VP1/VP2 intron retention extensions |
| PPP011 | PPP011-intron-retention-unified | Intron Retention Unified | 15 | 2026-02 | PPP007,PPP010 | active | Unified platform deck: VR4 library + N-terminal positions |

## Conventions

- **PPP numbers are unique and permanent.** Never reassign a PPP number to a different deck. Skipping numbers is fine; duplication is not.
- **Status values:** `active` (current, maintained), `archived` (superseded or historical, still buildable).
- **Parent:** PPP number(s) of the deck this one was derived from.
- **Build:** `cd decks/PPPxxx-name/ && npm install && node build.js`

## Shared Infrastructure

- **html2pptx/**: Shared build library at `voyager_slide_templates/html2pptx/`. All `build.js` files import from `../../html2pptx/index.cjs`.
- **assets/**: Shared Voyager branding PNGs at `voyager_slide_templates/assets/`. Decks access these via an `assets` symlink pointing to `../../assets`.
- **deck-assets/**: Per-deck folder for deck-specific images (diagrams, custom backgrounds). Not all decks have one.
