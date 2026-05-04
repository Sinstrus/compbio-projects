# Figure Bank — Scaffold-Based Figure Generation

Authoritative rules for all figures in this bank. Loaded as system instructions when working in `figure_bank/`.

## The mandate

**New FIGs MUST be rendered through a scaffold (SCF###).** Inline SVG generation
— hand-writing `<text>` tags with hand-calculated `x`/`y` coordinates — is
forbidden. Layout decisions belong to scaffolds, not to per-figure code.

The scaffold system enforces:
- **Locked viewBox** — every scaffold uses 960 × 465 (Voyager slide content area: 10 in × 4.84 in @ 96 dpi)
- **Locked slot geometry** — bounds, fonts, weights, alignments, `max_chars`
- **No mixing** — text in TextSlot/WrappedTextSlot, drawing in BodySlot. BodySlot can be `allow_text=False` to enforce object-only regions
- **Hard-fail on overflow** — `SlotOverflowError`, `SlotTypeError`, `BodyOverflowError` raised at render time

There is no soft-fail path. There is no post-hoc visual audit (the old
`audit_figs.sh` / Vision-API approach has been removed).

## Scaffold catalog

See `scaffolds/SCAFFOLD_REGISTRY.md` for the full registry. Quick map:

### Generic layout scaffolds (Wave 1)

| SCF | Use case |
|---|---|
| **SCF001** | Single panel — universal fallback (schematics, single charts via image embed, one-offs) |
| **SCF002** | Two side-by-side panels — comparisons, before/after |
| **SCF003** | Three side-by-side panels — three-condition contrasts (hard cap at 3) |

> SCF004 (parameterized 3-5 stage card workflow) was deprecated and removed during wave 2. Use SCF008-SCF012 instead.

### Domain scaffolds (Wave 2)

| SCF | Use case |
|---|---|
| **SCF005** | RepCap 3-zone — WT plasmid + VHH plasmid → resulting capsid (mosaic / trans-comp modes). Structured VPBar input. |
| **SCF006** | Architecture comparison matrix — 4–8 rows, 4 columns, color-coded category borders. |
| **SCF007** | Pulldown / protocol pipeline — dark title bar + 5 circle stages with icons + readout strip. |

### Workflow card scaffolds (Wave 2 — per step count)

| SCF | Layout |
|---|---|
| **SCF008** | 3-step single row card workflow |
| **SCF009** | 4-step single row card workflow |
| **SCF010** | 5-step single row card workflow (tightest single-row variant) |
| **SCF011** | 6-step wrap-around (3+3, both rows aligned) |
| **SCF012** | 7-step wrap-around (4+3, row 2 centered) |

Each card has the shape: dark header band (label + badge) → drawing-only icon area → body caption.

### Plot-grid scaffolds (Wave 3 — per panel count, "money slides")

For matplotlib dose-response and signal-vs-dose plots. Each scaffold owns ALL
text (figure title/subtitle, per-panel id/title/x-axis/y-axis labels) in
locked slots; matplotlib output must be **text-clean** (suppress
`set_title`, `set_xlabel`, `set_ylabel`).

| SCF | Panels | Layout |
|---|---|---|
| **SCF013** | 1 | full canvas |
| **SCF014** | 2 | 1×2 side-by-side (large axis fonts) |
| **SCF015** | 3 | 1×3 row (large axis fonts) |
| **SCF016** | 4 | 2×2 grid |
| **SCF017** | 5 | 1×5 ribbon (narrow tall) |
| **SCF018** | 6 | 2×3 grid |
| **SCF019** | 8 | 2×4 grid |
| **SCF020** | 9 | 2×5 with bottom-right empty |
| **SCF021** | 10 | 2×5 grid |
| **SCF022** | 11 | 2×6 with bottom-right empty |
| **SCF023** | 12 | 2×6 grid |
| **SCF024** | 4 | 2×2 (TALL — 960×800 viewBox, ~100 px/decade for log-Y) |

> 7-panel skipped (no clean grid). **No plot-grid scaffold uses >2 rows.** Above 8 panels, add columns rather than rows — taller plots (longer Y) beat wider plots (longer X) for dose-response data, since line separation in magnitude space is harder to read by eye than log-dose spacing. 1-3 panel scaffolds use 13-16pt axis fonts; 4-9 panel use 10-13pt; 10-12 panel use 9-11pt.

> **SCF024 is non-standard** — its 960×800 viewBox does not fit a Voyager content slide; use it for full-slide figures or shared-axis log-Y figures where SCF016's compressed 142-px body is too short. FIG178 (EXP26000387 absolute RLU 2×2) is the canonical SCF024 use case.

**Scaffold selection — Y-axis vs X-axis tradeoff:** when two scaffolds could fit your panel count, default to the layout with a longer Y axis. Log-spaced X axes (dose, dilution) compress gracefully; magnitude differences on Y do not.

**Bar chart orientation — use horizontal bars when there are ≥6 named categories.** The scaffold controls the outer chrome but the matplotlib body is the caller's responsibility. Vertical bars with rotated tick labels crowd badly at 9+ conditions inside a ~430px wide body. Horizontal bars (`ax.barh`) put the condition names on the Y axis where each row has ~35 px of space at `fontsize=8` — no rotation needed. Reverse the data array so the first condition appears at the top (matplotlib y=0 is bottom). FIG193 (EXP26000510 9-condition LiCOR western) is the canonical example.

## Matplotlib output via PanelRenderer

**Mandate: matplotlib panels embedded in plot-grid scaffolds (SCF013–SCF024) must
use `figure_bank.matplotlib_panel.PanelRenderer`.** Hand-tuned `fontsize=N,
rotation=R` was the source of every text-clipping/overflow bug we hit during
wave 3 retrofits (FIG174, FIG177, FIG178, FIG193). The renderer eliminates
that whole failure mode by sizing labels deterministically and detecting
collisions before save.

What it does, in one paragraph:

`PanelRenderer` owns the per-panel figsize/dpi/margins config; from those it
derives a per-label width budget and a perpendicular clearance budget (so
rot=90 labels can't clip into the bottom margin — the FIG177 bug). Calling
`renderer.size_x_labels(ax, labels=...)` measures the labels via
`text.get_window_extent()` (rotation-aware), runs the strategy ladder
(default → shrink → rotate 30° → rotate 90° → flip to barh, gated by
`axis_priority`), and applies the chosen `fontsize`/`rotation`/`ha` to `ax`.
`renderer.check_collisions(ax)` does a pairwise AABB check on tick labels +
legend frame + `ax.texts` (value annotations) and returns a list of detected
overlaps. `renderer.save_as_data_uri(fig, path)` applies margins + saves +
returns a base64 data URI for embedding in the scaffold's `ImageSlot`.

Pattern (FIG178 canonical, see `dna_engineer_agent/scripts/build_retrofit_fig178.py`):

```python
from figure_bank.matplotlib_panel import PanelRenderer

# ONE renderer per panel-group → shared figsize/margins is the shared-axis lock.
renderer = PanelRenderer(
    body_w_px=429, body_h_px=298,           # scaffold body box
    figsize_inches=(2.86, 1.99),            # match body aspect
    dpi=200,
    margins={"left": 0.16, "right": 0.97, "top": 0.95, "bottom": 0.13},
    axis_priority="y",                      # log-Y is the story; never barh-flip
)

def _render_panel(plot_fn, x_labels, out_path):
    fig, ax = renderer.fig_ax()
    plot_fn(ax)                             # plots data + sets xticks (positions only)
    ax.set_title(''); ax.set_xlabel(''); ax.set_ylabel('')   # scaffold owns text
    x_sizing = renderer.size_x_labels(ax, labels=x_labels)
    y_sizing = renderer.size_y_labels(ax, labels=Y_TICK_LABELS)
    for w in x_sizing.warnings + y_sizing.warnings:
        print(f"  ! sizing: {w}")
    for c in renderer.check_collisions(ax):
        print(f"  ! collision: {c}")
    return renderer.save_as_data_uri(fig, out_path)
```

**Shared-axis figures** (SCF024 et al): use **one** `PanelRenderer` instance
for all N panels. The shared-figsize/shared-margins lock that makes panels
pixel-align is now expressed as one config object, not N copies of the same
hand-tune.

**`axis_priority="y"`** (the default) blocks the horizontal-bar flip strategy
— preferred for dose-response and signal-vs-dose figures where Y-axis
magnitude resolution wins. Pass `"x"` only when the data benefits from a
longer X axis (rare; manual label-density situations).

**Do not pass `bbox_inches='tight'` to `savefig`.** It conflicts with
`subplots_adjust` and breaks the shared-axis pixel alignment across panels.
The renderer's locked margins are the right answer; if labels seem cropped,
inspect the sizer warnings (the clearance gate may be telling you to
shorten labels or use a smaller body scaffold).

### Auto-placement: legends and corner text vs data overlap

`renderer.check_collisions(ax)` (default `include_data=True`) flags legends
and `ax.texts` that cover plotted lines/bars/scatter — not just inter-text
overlaps. Two helpers act on what it detects:

- **`renderer.auto_place_legend(ax, **legend_kwargs) -> PlacementResult`** —
  use this instead of `ax.legend(loc=...)`. Tries 8 compass positions
  (`upper left`, `upper center`, `upper right`, `center right`, `lower right`,
  `lower center`, `lower left`, `center left`); accepts the first with zero
  data overlap. Falls back through half-decade then full-decade headroom on
  log axes (or +25%/+50% on linear) if every interior position overlaps.
  When even the headroom ladder can't clear, returns the **least-bad**
  position with `overlapped=True`.

- **`renderer.place_corner_text(ax, text, prefer='upper right', **text_kwargs)
  -> PlacementResult`** — replaces hand-coded
  `ax.text(0.96, 0.96, ..., transform=ax.transAxes)` corner badges. Same
  ladder, but the candidate space is 8 anchor points in `transAxes` (4
  corners + 4 edge midpoints).

**Shared-axis grids: use the two-pass workflow instead of `allow_headroom=False`.**
In shared-axis-locked figures (FIG178, FIG181-187), each panel is a separate
matplotlib figure, so growing one panel's ylim desyncs the grid. The correct
pattern is:

1. **Pass 1 — probe:** For each legend-bearing panel, render a throw-away figure,
   then call `probe_legend_headroom(ax, **legend_kwargs) -> str`. This is
   non-destructive (restores `ylim` and removes the test legend). It returns the
   headroom label needed: `''`, `'half_decade'`, `'full_decade'`, or `'exhausted'`.

2. **Reconcile:** `max_headroom_label([label_a, label_b, ...])` returns the most
   demanding label across all panels.

3. **Compute shared ylim:** `headroom_ylim_top(base_top, label, yscale='log')` converts
   the label to a new numeric top. Apply `(base_bot, new_top)` uniformly to all panels.

4. **Pass 2 — render:** Render all panels with `ylim_override=shared_ylim` and
   `allow_headroom=False` (no further escalation — ylim is already extended).

See `build_retrofit_fig178.py` (`_probe_shared_ylim()` + `build_fig178()`) for the
canonical example. FIG178 needed `full_decade` headroom for panel A's 5-line legend,
extending the shared top from `2e6 → 2e7`.

`allow_headroom=False` without a probe pass is now the **fallback** for grids where
headroom would cause more harm than the residual overlap (e.g. a 9-panel figure where
extending ylim compresses all bars). In that case the helper still picks the least-bad
compass position and returns `overlapped=True` so the build script logs it honestly.

Full design: `figure_bank/SPEC_matplotlib_panel_sizer.md`. Tests:
`figure_bank/test_matplotlib_panel.py` (runnable as plain Python).

## Building a new FIG

1. **Pick a scaffold.** Read `scaffolds/SCAFFOLD_REGISTRY.md` and the relevant
   `SCF###_*.py` module. If nothing fits, **stop and ask the user** — propose
   a new SCF### shape. Do not fall back to inline SVG.

2. **Get the next FIG number:**
   ```python
   from asset_registry import get_next_number
   fig = get_next_number("../figure_bank/FIGURE_REGISTRY.md", "FIG")
   ```

3. **Write a small build script** (in `dna_engineer_agent/scripts/`) that
   imports the scaffold and calls `render(...)`:
   ```python
   from figure_bank.scaffolds import SCF005_repcap_three_zone as scf
   from figure_bank.scaffolds.SCF005_repcap_three_zone import (
       VPBar, VHH_X_VR4, sample_capsid_body, MOSAIC_VHH_COLOR,
   )

   svg = scf.render(
       title="VR4 Loop Mosaic",
       subtitle="DNA Ratio 1:1",
       mode="mosaic",
       wt={"header": "AVD002 — WT", "bars": [
           VPBar("VP1", "ATG", 80.1),
           VPBar("VP2", "ACG", 65.5),
           VPBar("VP3", "ATG", 59.0),
       ]},
       vhh={"header": "AVD006 — VR4-VHH", "bars": [
           VPBar("VP1", "ATG", 95.6, vhh_x=VHH_X_VR4),
           VPBar("VP2", "ACG", 80.7, vhh_x=VHH_X_VR4),
           VPBar("VP3", "ATG", 74.3, vhh_x=VHH_X_VR4),
       ]},
       result={
           "header": "Resulting Capsid",
           "body_svg": sample_capsid_body(n_vhh_markers=5, vhh_color=MOSAIC_VHH_COLOR),
       },
   )
   Path(f"figure_bank/{fig}-vr4-loop-mosaic.svg").write_text(svg)
   ```

4. **Handle exceptions.** If `render()` raises:
   - `SlotOverflowError` — shorten the offending content, OR pick a different scaffold with a larger slot.
   - `SlotTypeError` — usually a missing required field, wrong list length, or text in a no-text BodySlot. Fix the call site.
   - `BodyOverflowError` — body SVG has elements outside the body box. Move them in (coords are body-relative).
   - **Never** widen a slot box for one figure — that defeats the lockdown.

5. **Update registries:**
   - `FIGURE_REGISTRY.md` — add the FIG row
   - Run `validate_bank("FIGURE_REGISTRY.md", "figure_bank/", "FIG")` to confirm consistency

6. **Pre-flight SVG sanity check** (entity references, viewBox, XML well-formedness):
   ```python
   from asset_registry import validate_svg
   validate_svg("figure_bank/FIG###-*.svg")
   ```
   This catches stale entity issues (`&times;` instead of `&#215;`) but is
   not the layout enforcer — the scaffold already handled that.

## When a scaffold doesn't fit

Two valid responses, in this order:

1. **Re-frame the figure** to fit an existing scaffold. Most figures that
   "need" a custom layout can be expressed as one of SCF001–SCF012 with
   different slot content. SCF001 (single panel free-form body) is the
   universal fallback.

2. **Propose a new SCF.** Write a one-paragraph description of the new
   scaffold's shape (slots, locks, viewBox is always 960×465) and ask the
   user before implementing. **Bar for adding a scaffold: at least 3
   concrete figures sharing the same structural pattern.** New scaffolds
   get the next available SCF number from `scaffolds/SCAFFOLD_REGISTRY.md`.

**Not valid:** writing inline SVG, copy-pasting an existing FIG's layout
manually, or "just this once" overrides.

## Scaffold authoring (for new SCF###)

Each scaffold lives in `figure_bank/scaffolds/SCF###_descriptive_name.py` and
exposes:

- `BARCODE` — string `"SCF###"`
- `DESCRIPTION` — one-line description for the registry
- `VIEWBOX` — `(0, 0, 960, 465)` (locked across all scaffolds)
- Module-level slot constants (`TextSlot` / `WrappedTextSlot` / `BodySlot` /
  `ImageSlot`) — these encode the locked geometry. For per-row/per-card
  slots, write a factory function that derives the slot from the row index.
- One or more structured-input `@dataclass(frozen=True)` types (e.g. `Card`,
  `Stage`, `VPBar`, `Row`). **All fields required, no optional defaults that
  create per-instance visual irregularity.** If a field can be empty, use
  a deliberate placeholder (e.g. `"—"`).
- `render(*, ...)` — emitter function (kw-only args). Validates structured
  input shape, then composes the SVG by calling slot.render() / slot.embed().
- `preview(show_outlines=True)` — fills slots with filler content. **Every
  slot's preview content should declare its lock** (font + max_chars, e.g.
  `"name — 15pt bold (max 28)"`).
- Optional convenience helpers like `sample_capsid_body()` /
  `sample_card_icon()` for callers who want a placeholder while they figure
  out the real domain drawing.

After editing or adding a scaffold:
```bash
python -m figure_bank.scaffolds._previews regenerate         # all
python -m figure_bank.scaffolds._previews regenerate SCF005  # just one
```
Then visually inspect `scaffolds/previews/SCF###-preview.svg` (outlined,
shows every slot bound) and `SCF###-preview-clean.svg` (what real figures
look like).

### Content validation — pixel-width, not character count

**Character-count checks (`len(text) > max_chars`) are forbidden in scaffold
slot validators.** A 45-char limit in a 183px column at 8.5pt allows ~225px
of text — wider than the column — and raises no error. Character width varies
by glyph; counting chars is a meaningless proxy for what actually renders.

Use `_check_px()` with real Liberation Sans font metrics instead:

```python
from ._runtime import measure_text_width, SlotOverflowError

def _check_px(slot: str, text: str, font_size: float, font_weight: str, max_px: float) -> None:
    w = measure_text_width(text, font_size, font_weight)  # type: ignore[arg-type]
    if w > max_px:
        raise SlotOverflowError(slot, w, max_px, unit=" px")
```

- `font_size` must match the SVG `font-size` attribute exactly (float OK — PIL accepts it)
- `font_weight` must be `"bold"` or `"normal"` to select the right font file
- `max_px` = actual column pixel width (from locked geometry constants) minus a small safety margin (4–12 px)

**SCF025 is the canonical example:** col 1 enforced at 179 px (10pt bold / 8.5pt normal / 8pt normal), col 4 at 113 px (8pt bold / 7.5pt normal). The error message reports the exact measured width and limit in px.

**When `SlotOverflowError` fires: shorten the build script content. Never raise
the pixel limit to accommodate content.** The limit is the authority. If you need
more room, use a different scaffold or redesign the figure layout.

## Bookkeeping

- `get_next_number(registry, prefix)` — next available barcode (works for FIG and SCF)
- `validate_bank(registry, dir, prefix)` — orphans/missing/duplicates
- `validate_svg(path)` — XML well-formedness, viewBox presence, entity sanity
- File naming: `FIG###-descriptive-slug.svg`

## Hard rules summary

1. No inline SVG. Use a scaffold.
2. No widening slots for "just this figure."
3. Symmetry is mandatory — multi-panel scaffolds enforce identical slot geometry.
4. `SlotOverflowError` means shorten content in the build script. **Never raise
   the pixel limit or character limit to accommodate content.** The limit is the
   authority. If the content won't compress, use a different scaffold.
5. FIG numbers are unique and permanent. SCF numbers are unique and permanent.
   SCF004 was retired during wave 2 and is not reusable.
6. Existing pre-scaffold FIGs (FIG001–FIG190) are legacy and untouched. Only
   new FIGs go through scaffolds.
7. **Scaffold slot validators must use pixel-width measurement, not character
   counting.** Use `_check_px()` with `measure_text_width()` from `_runtime.py`.
   Character limits like `len(text) > 45` are meaningless — glyph widths vary.
   See "Content validation" section above for the canonical pattern.

## SVG entity reminders

These still matter inside `body_svg` fragments — the scaffold doesn't rewrite them:

- `&#215;` not `&times;`
- `&#8594;` not `&rarr;`
- `&#8212;` not `&mdash;`
- `&#160;` not `&nbsp;`
- Alpha as `fill-opacity="0.75"`, not 8-digit hex `#ffffffcc`
