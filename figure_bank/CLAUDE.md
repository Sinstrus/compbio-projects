# Figure Bank — Design Rules & Verification

Authoritative design rules for all figures in the bank. Loaded as system instructions when working in `figure_bank/`.

## Mandatory Checklist

Before finalizing ANY figure:

1. [ ] FIG number obtained from `get_next_number("figure_bank/FIGURE_REGISTRY.md", "FIG")`
2. [ ] File deposited in `figure_bank/`, NEVER in project directories
3. [ ] FIGURE_REGISTRY.md updated with FIG number, filename, format, project, related constructs/decks, and description
4. [ ] `validate_bank()` passes (no orphans, no missing files, no duplicate numbers)
5. [ ] SVG entities use numeric references: `&#215;` not `&times;`, `&#8594;` not `&rarr;`, `&#8212;` not `&mdash;`, `&#160;` not `&nbsp;`
6. [ ] **DESIGN-011 alignment**: content center-justified in grid cells; headers centered over columns; legends/footnotes left-aligned
7. [ ] Lineage noted ("Modified from **FIG###**") in registry description if derived from existing figure
8. [ ] Minimum font size: 14pt for annotations, 16pt for primary text
9. [ ] `validate_svg()` passes (entity + viewBox + XML well-formedness checks)

## DESIGN-011 Quick Reference — Hierarchical Alignment

| Element | Alignment | Why |
|---------|-----------|-----|
| Content inside grid cells (IDs, badges, icons) | **Center** within cell | Cell is a visual container |
| Column/row headers | **Center** over the cell span they label | Headers "own" their column |
| Legends, footnotes, keys | **Left-align** | Read like prose, left-to-right |
| Titles | **Left-align** | Anchored to top-left reading origin |

**SVG centering formula** for paired text (e.g. "AVD006 + AVD002"):
```xml
<!-- Cell center = (cell_x + cell_width/2) -->
<text x="{center-7}" text-anchor="end">AVD006</text>
<text x="{center}" text-anchor="middle">+</text>
<text x="{center+7}">AVD002</text>
```

**Centering a badge group** (N items):
```
group_width = N * item_w + (N-1) * gap
start_x = center - group_width / 2
```

This applies to ALL grid/matrix figures: experimental conditions, plate layouts, comparison tables, dose grids.

## SVG Gotchas

1. **Entity references**: XML does not define HTML entities. Always use numeric character references:
   - `&#215;` not `&times;` (multiplication/cross)
   - `&#8594;` not `&rarr;` (right arrow)
   - `&#8212;` not `&mdash;` (em dash)
   - `&#160;` not `&nbsp;` (non-breaking space)

2. **Arrow double-rendering**: Use EITHER `<line>` with `marker-end` OR a text arrow character — never both (FIG001 had doubled `>` + arrowhead).

3. **ViewBox required**: Every SVG must have a `viewBox` attribute with 4 numeric values (e.g. `viewBox="0 0 900 420"`).

4. **Slide-optimized canvas**: Use `viewBox="0 0 900 420"` (16:9 widescreen) for figures intended for slide decks.

## Font Hierarchy (all Arial/Helvetica)

| Level | Size | Weight | Use |
|-------|------|--------|-----|
| Title | 24pt | bold | Figure title |
| Pill label | 20pt | bold | Key label (orange/green) |
| Box header | 18pt | bold | Construct names |
| VP labels | 16pt | bold | VP1/VP2/VP3 |
| Legend text | 16pt | normal | Legend items |
| Annotations | 14pt | normal | Secondary annotations (minimum size) |

## Programmatic Verification

```python
from asset_registry import get_next_number, validate_bank, validate_svg

# Get next FIG number
next_fig = get_next_number("../figure_bank/FIGURE_REGISTRY.md", "FIG")

# Audit registry/bank consistency
issues = validate_bank("../figure_bank/FIGURE_REGISTRY.md", "../figure_bank/", "FIG")

# Check SVG for common anti-patterns
issues = validate_svg("../figure_bank/FIG001-example.svg")
```

## Detailed Template Coordinates

See `memory/figure_design.md` for full slide-optimized layout coordinates (VP bars, box positions, icosahedron geometry, VHH markers).
