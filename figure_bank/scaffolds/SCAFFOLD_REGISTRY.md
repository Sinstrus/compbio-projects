# Scaffold Registry (SCF###)

Locked figure-layout templates. Every new FIG must be rendered through one of these scaffolds.
See `figure_bank/CLAUDE.md` for the scaffold mandate.

All scaffolds use the **960 × 465** viewBox (10 in × 4.84 in @ 96 dpi — Voyager slide content area).

When adding a new scaffold:
1. Use the next available SCF number (`get_next_number(...)` works on this file).
2. SCF numbers are unique and permanent. **Skipped numbers are fine; reuse is not.** SCF004 was deprecated and removed during wave 2; the number is retired.
3. Add a row below.
4. Implement the module and `preview()` function.
5. Run `python -m figure_bank.scaffolds._previews regenerate SCF###` and visually inspect the output.

## Generic layout scaffolds (Wave 1)

| SCF    | Module                  | Slots                                                                          | Preview                       | Description                                                       |
|--------|-------------------------|--------------------------------------------------------------------------------|-------------------------------|-------------------------------------------------------------------|
| SCF001 | SCF001_single_panel.py  | title, subtitle?, body OR body_image, legend?                                  | previews/SCF001-preview.svg   | Single panel — universal fallback for diagrams, charts, one-offs  |
| SCF002 | SCF002_two_panel.py     | title, subtitle?, panel_{left,right}{header, body, caption}, legend?           | previews/SCF002-preview.svg   | Two side-by-side panels (locked symmetric geometry)               |
| SCF003 | SCF003_three_panel.py   | title, subtitle?, panel_{left,center,right}{header, body, caption}, legend?    | previews/SCF003-preview.svg   | Three side-by-side panels (hard cap at 3)                         |

> SCF004 (parameterized 3-5-stage card workflow) was retired. Use the per-step-count scaffolds SCF008-SCF012 instead.

## Domain scaffolds (Wave 2)

| SCF    | Module                                  | Slots                                                                                        | Preview                       | Description                                                          |
|--------|-----------------------------------------|----------------------------------------------------------------------------------------------|-------------------------------|----------------------------------------------------------------------|
| SCF005 | SCF005_repcap_three_zone.py             | 33 slots: title, subtitle?, key_pill, 2× plasmid box × {header + 3× VP{label/bar/icon/tag}}, result{header, body}, 4× legend{swatch, label} | previews/SCF005-preview.svg   | RepCap 3-zone: WT plasmid + VHH plasmid → resulting capsid (mosaic / trans-comp modes) |
| SCF006 | SCF006_architecture_matrix.py           | 6 figure-level + 7 per-row × ≤8 rows: title, subtitle?, 4 column headers, rows[i].{border, name, xref, col1, col2, col3a, col3b} | previews/SCF006-preview.svg   | Architecture comparison matrix (4–8 rows, 4 columns, color-coded categories) |
| SCF007 | SCF007_pulldown_pipeline_5stage.py      | 24 slots: title.bar.{title, badge}, 5× stages[i]{icon, label, sublabel}, readout.{label, body} | previews/SCF007-preview.svg   | Pulldown / protocol pipeline — dark title bar + 5 circle stages + readout strip |

## Workflow card scaffolds (Wave 2 — per step count)

| SCF    | Module                                  | Layout      | Card width | Card slots                                | Preview                       |
|--------|-----------------------------------------|-------------|------------|-------------------------------------------|-------------------------------|
| SCF008 | SCF008_card_workflow_3step.py           | 3 single row | 290 px    | label (max 22) + badge (max 12) + icon + body (max 130 / 5 lines) | previews/SCF008-preview.svg |
| SCF009 | SCF009_card_workflow_4step.py           | 4 single row | 211 px    | label (max 18) + badge (max 10) + icon + body (max 90 / 5 lines)  | previews/SCF009-preview.svg |
| SCF010 | SCF010_card_workflow_5step.py           | 5 single row | 164 px    | label (max 14) + badge (max 8)  + icon + body (max 70 / 5 lines)  | previews/SCF010-preview.svg |
| SCF011 | SCF011_card_workflow_6step_wrap.py      | 6 (3+3 wrap) | 271 px    | label (max 18) + badge (max 8)  + icon + body (max 60 / 2 lines)  | previews/SCF011-preview.svg |
| SCF012 | SCF012_card_workflow_7step_wrap.py      | 7 (4+3 wrap, row 2 centered) | 197 px | label (max 14) + badge (max 6) + icon + body (max 50 / 2 lines) | previews/SCF012-preview.svg |

## Plot-grid scaffolds (Wave 3 — per panel count)

The "money slide" family — for matplotlib dose-response and signal-vs-dose
data plots. Each scaffold owns ALL text via locked slots: figure title +
subtitle + per-panel { id, title, x-axis label, y-axis label (rotated -90) }.

**Caller responsibility:** matplotlib output must be **text-clean** —
suppress `ax.set_title('')`, `fig.suptitle('')`, `ax.set_xlabel('')`,
`ax.set_ylabel('')`. The image embedded into `panels[i].body` should be
just the axes + data + legend. The scaffold handles all surrounding text.

Panel-title and axis-label fonts scale by panel count: 1-3 panels use 13-16pt
(large-axis tier, axis labels readable against larger plots), 4-9 panels use
10-13pt, 10-12 panels use 9-11pt.

| SCF    | Module                                  | Panels | Layout            | Per-panel size | Preview                       |
|--------|-----------------------------------------|--------|-------------------|----------------|-------------------------------|
| SCF013 | SCF013_plot_grid_1panel.py              | 1      | 1×1 full canvas   | 920 × 372      | previews/SCF013-preview.svg   |
| SCF014 | SCF014_plot_grid_2panel.py              | 2      | 1×2 side-by-side  | 454 × 372      | previews/SCF014-preview.svg   |
| SCF015 | SCF015_plot_grid_3panel.py              | 3      | 1×3 row           | 298 × 372      | previews/SCF015-preview.svg   |
| SCF016 | SCF016_plot_grid_4panel.py              | 4      | 2×2 grid          | 454 × 180      | previews/SCF016-preview.svg   |
| SCF017 | SCF017_plot_grid_5panel_ribbon.py       | 5      | 1×5 ribbon        | 176 × 372      | previews/SCF017-preview.svg   |
| SCF018 | SCF018_plot_grid_6panel.py              | 6      | 2×3 grid          | 300 × 181      | previews/SCF018-preview.svg   |
| SCF019 | SCF019_plot_grid_8panel.py              | 8      | 2×4 grid          | 222 × 181      | previews/SCF019-preview.svg   |
| SCF020 | SCF020_plot_grid_9panel.py              | 9      | 2×5 minus 1 (last cell empty) | 176 × 181 | previews/SCF020-preview.svg |
| SCF021 | SCF021_plot_grid_10panel.py             | 10     | 2×5 grid          | 176 × 181      | previews/SCF021-preview.svg   |
| SCF022 | SCF022_plot_grid_11panel.py             | 11     | 2×6 minus 1 (last cell empty) | 145 × 181 | previews/SCF022-preview.svg |
| SCF023 | SCF023_plot_grid_12panel.py             | 12     | 2×6 grid          | 145 × 181      | previews/SCF023-preview.svg   |

## Tall-canvas variants (Wave 3+)

For shared-axis figures where Y-axis resolution matters more than fitting the
standard Voyager 960×465 content area. **ViewBox is non-standard** — these
SCFs are intended for full-slide figure use, not standard slide content
embedding. Use only when the visual story requires more vertical room (log
Y range spanning 3+ decades, dose-response across many orders of magnitude).

| SCF    | Module                                  | Panels | Layout | viewBox | Per-panel size | Body Y resolution | Preview |
|--------|-----------------------------------------|--------|--------|---------|----------------|-------------------|---------|
| SCF024 | SCF024_plot_grid_4panel_tall.py         | 4      | 2×2    | 960×800 | 454 × 347      | ~298 px → ~100 px/decade @ 3 dec | previews/SCF024-preview.svg |

> Note: 7-panel skipped intentionally — no clean grid arrangement; use SCF019 (8) or SCF018 (6) instead.

> **No plot-grid scaffold uses more than 2 rows.** For dose-response and signal-vs-dose data, taller plots (longer Y axis) win over wider plots (longer X axis): line separation in magnitude space is harder to read by eye than dose spacing on a log X axis. Above 8 panels, add columns rather than rows — 9/10/11/12 panel layouts are (2×5)−1 / 2×5 / (2×6)−1 / 2×6 respectively, never 3+ rows.

### Picking between scaffolds — Y-axis vs X-axis tradeoff

When two scaffolds could both fit your panel count (e.g. SCF018 6-panel 2×3 vs a hypothetical 1×6 ribbon), ask:

- **Does the data benefit more from a longer Y axis or a longer X axis?**
- **Long Y wins** when: line separation in magnitude space matters (dose-response curves at different titers, signal-vs-dose, comparing transduction efficiency between constructs). Lines that look identical on a squashed Y become readable when the Y axis stretches.
- **Long X wins** when: many distinct X values need clear spacing (e.g. timecourse with 20+ time points, sequence position scans). Rare in this codebase — most "money slides" are dose-response.
- **Default to long Y.** Log-spaced X axes (dose, dilution) compress gracefully; eye reads log spacing fine even when squished. Y-axis magnitude differences do not.

Each workflow card has the same shape: dark header band (label + badge) → drawing-only icon area → body caption.
Wrap-around variants (SCF011/012) use shorter cards (170 tall vs 282) to fit two rows on the canvas, with an
overshoot-hook wrap arrow that approaches the destination card's left edge perpendicularly.

## Slot conventions

- **TextSlot** locks font size, weight, alignment, color, and `max_chars`. Overflow → `SlotOverflowError`.
- **WrappedTextSlot** is multi-line; locks `max_chars`, `max_lines`, and `line_height` in addition.
- **BodySlot** is a free-form SVG region. Bounds-check elements; `allow_text=False` rejects any `<text>`.
- **ImageSlot** embeds `<image href>` at locked dimensions; no content validation.
- All viewBoxes are 960 × 465 (16:9 widescreen). Coordinates are in SVG user units (≈ pixels at 96 dpi).

## Hard rules

1. **No inline SVG generation for FIGs.** Use a scaffold or stop and propose a new SCF.
2. **Never widen a slot box "just for this figure."** Either shorten the content or pick a different scaffold.
3. **Every drawing in a BodySlot, every text in a Text/WrappedTextSlot — no mixing.**
4. **All structured-input dataclass fields are required.** No optional defaults that create per-row visual irregularity. If a row genuinely has no second line, pass a deliberate placeholder (e.g. `"—"`).
5. **`show_outlines=True` is for previews only.** Real figures have no slot outlines.

## Hard-fail summary

The runtime raises:
- `SlotOverflowError` — text exceeds `max_chars`, rendered width exceeds slot width, or wrapped content needs more lines than `max_lines`.
- `SlotTypeError` — wrong dataclass type, missing required field, bad enum value, wrong list length, `<text>` in a no-text BodySlot, etc.
- `BodyOverflowError` — SVG element with explicit coords falls outside a BodySlot's box (skipped inside `<g transform>` blocks; clipPath enforces visual containment).

There is no soft-fail path. Callers must shorten content, pick a different scaffold, or propose a new one.

## Structured-table scaffolds (Wave 4)

Variable-height figures for design panel comparisons. **ViewBox width is 960; height is computed at render time** from the number of rows and their text-line counts.

| SCF    | Module                                  | ViewBox              | Description                                                                         |
|--------|-----------------------------------------|----------------------|-------------------------------------------------------------------------------------|
| SCF025 | SCF025_design_panel_overview.py         | 960 × (variable)     | 8-column design panel table: pill strip, outcome range bars, target-zone overlay    |

SCF025 row-height rule: `total_lines = 1 + len(description) + (1 if annotation else 0)` → 4+ lines: 60 px | 3 lines: 50 px | 2 lines: 40 px. Bar axis ticks locked at 0/20/40/60/80/90/100% (3.27 px/%).

## Roadmap (Wave 5+ — not yet built)

Bar for adding a scaffold: at least 3 concrete figures sharing the same structural pattern.

- **Tier-Library Table** — ranked rows with embedded mini bar chart per row (FIG125)
- **Display-Library Triplet** — RepCap + ITR cassette + result (FIG145–147; need cleaner cassette pattern first)
- **Dual-State Mechanism** — left/right state diagram (FIG060; only 1 example so far)
- **Pulldown Conditions Matrix** — close to SCF006 but distinct shape (FIG115)
- **Western blot lane diagrams** — wrap existing `gel_render.py` as an SCF
- **Pipelines at non-5 stage counts** — wait for examples
- **Card workflows beyond 7 steps** — wait for examples
- **7-panel plot grid** — skipped in wave 3, no clean grid arrangement; revisit if needed
