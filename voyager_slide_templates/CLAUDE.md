# Voyager Slide Templates — PPTX-Compatible Authoring Guide

This directory contains Voyager-branded HTML slide decks that export to PowerPoint via the html2pptx pipeline.

## PPPxxx Barcoding System

All decks are tracked by a unique PPP number, analogous to the AVDxxx system for DNA constructs.

**Registry:** `PRESENTATION_REGISTRY.md` — the authoritative manifest with PPP IDs, names, and status.

### Rules for New Decks
1. Use the next available PPP number (check `PRESENTATION_REGISTRY.md`)
2. **PPP numbers are unique and permanent.** Never reassign a PPP number to a different deck. Skipping numbers is fine; duplication is not.
3. Name folders as: `PPPxxx-descriptive-name/`
4. Update `PRESENTATION_REGISTRY.md` with PPP number, folder name, slide count, date, status, and description
5. Build with: `cd decks/PPPxxx-name/ && npm install && node build.js`

## Critical Rule

**DO NOT use the CSS gallery files** (`css/`, `templates/`) for authoring new decks. Those are visual reference only — CSS `background-image` properties and CSS variables break the html2pptx conversion pipeline.

All new decks MUST use the **inline html2pptx format** documented below.

## Slide Dimensions

All slides are 960 x 540 px (16:9 ratio). The `<body>` tag defines the slide boundary.

## Title Slide Template

```html
<body style="width:960px;height:540px;margin:0;position:relative;font-family:Calibri,Arial,sans-serif;background:#fff;">
  <!-- Neuron bg (right side, above footer) -->
  <div style="position:absolute;top:0;right:0;width:100%;height:480px;">
    <img src="assets/Title_Slide_Background_Neuron_Art1.png" style="position:absolute;right:0;top:0;height:100%;"/>
  </div>
  <!-- White fade overlay -->
  <div style="position:absolute;left:0;top:0;width:60%;height:480px;background:linear-gradient(90deg,#fff 30%,rgba(255,255,255,0.85) 55%,transparent);" data-rasterize="gradient"></div>
  <!-- Logo top-left -->
  <img src="assets/Title_Slide_Voyager_Logo1.png" style="position:absolute;left:36px;top:28px;width:200px;"/>
  <!-- Title block -->
  <div style="position:absolute;left:36px;top:140px;max-width:55%;">
    <h1 style="font-family:'Century Gothic',CenturyGothic,AppleGothic,sans-serif;font-size:24pt;color:#001446;font-weight:bold;margin:0 0 6px 0;">Title Text</h1>
    <p style="font-size:14pt;color:#003482;margin:0 0 10px 0;font-weight:bold;">Subtitle</p>
    <p style="font-size:9pt;color:#64707E;margin:0;">Author | Date</p>
  </div>
  <!-- Footer banner (60px) -->
  <div style="position:absolute;bottom:0;left:0;width:100%;height:60px;">
    <img src="assets/Title_Slide_Footer_Blue_Banner1.png" style="position:absolute;top:0;left:0;width:100%;height:100%;object-fit:fill;"/>
    <img src="assets/Title_Slide_Footer_Neuron_Art1.png" style="position:absolute;right:0;bottom:0;height:100%;"/>
    <p style="position:absolute;right:30px;top:20px;margin:0;font-size:8pt;color:#fff;font-family:Calibri,Arial,sans-serif;letter-spacing:3px;">C O N F I D E N T I A L</p>
  </div>
</body>
```

## Content Slide Template

```html
<body class="col" style="width:960px;height:540px;margin:0;background:#fff;font-family:Calibri,Arial,sans-serif;">
  <!-- Header: 75px with img-based banner + absolute logo -->
  <div style="position:relative;height:75px;flex-shrink:0;">
    <img src="assets/Body_Slide_Blue_Title_Background_Medium_Level_Banner1.png" style="position:absolute;top:0;left:0;width:100%;height:100%;"/>
    <div class="row items-center" style="position:absolute;top:0;left:0;width:100%;height:100%;padding:0 20px 0 30px;">
      <h1 style="font-family:'Century Gothic',CenturyGothic,AppleGothic,sans-serif;font-size:18pt;color:#fff;font-weight:bold;margin:0;flex:1;">Slide Title</h1>
    </div>
    <img src="assets/Body_Slide_Blue_Title_Voyager_Logo1.png" style="position:absolute;right:12px;top:0;height:75px;width:128px;"/>
  </div>
  <!-- Content area -->
  <div class="col" style="padding:10px 30px 0 30px;flex:1;gap:8px;">
    <!-- Slide content here -->
  </div>
  <!-- Footer -->
  <div style="height:44px;padding:4px 30px 0 30px;">
    <p style="font-size:7pt;color:#64707E;margin:0;">N&ensp;|&ensp;<span style="letter-spacing:3px;">CONFIDENTIAL</span></p>
  </div>
</body>
```

## How Assets Get Into PPTX

| Asset | HTML Pattern | PPTX Result |
|---|---|---|
| Blue header banner | `<img src="assets/Body_Slide_Blue_Title_Background_Medium_Level_Banner1.png" style="position:absolute;..."/>` | Positioned image at top |
| Header logo | `<img src="assets/Body_Slide_Blue_Title_Voyager_Logo1.png" style="position:absolute;right:12px;top:0;..."/>` | Right-aligned in header |
| Title neuron art | `<img src="assets/Title_Slide_Background_Neuron_Art1.png" style="position:absolute;right:0;top:0;..."/>` | Background image |
| White fade overlay | `<div style="background:linear-gradient(...)" data-rasterize="gradient">` | Playwright screenshots gradient to PNG |
| Footer banner | `<img>` tags with absolute positioning | Positioned images at bottom |
| Navy gradient boxes | `<div style="background:linear-gradient(90deg,#001446,#003482);..." data-rasterize="gradient">` | Rasterized to PNG |

**Key**: Every visual element is either an `<img>` tag (directly extracted) or a `data-rasterize="gradient"` element (screenshot-captured by Playwright). Nothing relies on external CSS.

## Mandatory Rules

1. **No `<!DOCTYPE>`, `<html>`, `<head>`, or `<link>` tags** — slides are bare `<body>` elements
2. **No CSS variables** — use inline hex colors (see table below)
3. **No CSS `background-image`** — use `<img>` tags with `position:absolute`
4. **`data-rasterize="gradient"` on EVERY CSS gradient** — required for Playwright to capture them
5. **Logo must be absolute-positioned** outside the flex row, not inside it
6. **Assets via `<img>` tags only** — referenced from `assets/` directory
7. **All visible text MUST be in `<p>`, `<h1>`-`<h6>`, `<ul>`, or `<ol>` tags** — html2pptx only extracts text from these supported tags. `<span>` is ONLY for inline formatting within those tags. Bare `<span>` elements (e.g., as flex children) are silently dropped in PPTX output.

## Color Reference (Inline Hex, Not CSS Variables)

| Name | Hex | Usage |
|---|---|---|
| Voyager Navy | `#001446` | Gradient start, primary dark |
| Voyager Blue | `#003482` | Gradient end, headings, accents |
| Bright Orange | `#FF7400` | Secondary accent, highlights |
| Mid Blue | `#4E8EC9` | Tertiary accent |
| Slate | `#64707E` | Footer text, labels |
| Green | `#4CAF50` | Status badges, positive indicators |
| Red | `#c62828` | Warnings, KO indicators |

## CSS Utility Classes

The html2pptx engine supports these utility classes:

| Class | Effect |
|---|---|
| `class="col"` | `display:flex; flex-direction:column` |
| `class="row"` | `display:flex; flex-direction:row` |
| `class="row items-center"` | Row with `align-items:center` |

## Content Patterns

### Navy Gradient Box (Banner/Callout)
```html
<div style="background:linear-gradient(90deg,#001446,#003482);padding:8px 14px;border-radius:6px;" data-rasterize="gradient">
  <p style="font-size:8pt;color:rgba(255,255,255,0.7);margin:0 0 2px 0;">LABEL</p>
  <p style="font-size:9pt;color:#fff;margin:0;font-weight:bold;">Content text</p>
</div>
```

### Div-Based Table Header Row
```html
<div class="row" style="gap:3px;margin-bottom:3px;">
  <div style="flex:1;background:linear-gradient(90deg,#001446,#003482);padding:4px 6px;border-radius:3px;" data-rasterize="gradient">
    <p style="font-size:7pt;color:#fff;margin:0;font-weight:bold;">Column Header</p>
  </div>
  <!-- repeat for each column -->
</div>
```

### Metric Row
```html
<div class="row" style="gap:8px;">
  <div style="flex:1;padding:6px;text-align:center;">
    <p style="font-family:'Century Gothic',...;font-size:14pt;color:#003482;margin:0;font-weight:bold;">42</p>
    <p style="font-size:7pt;color:#666;margin:1px 0 0 0;">Label</p>
  </div>
  <!-- repeat -->
</div>
```

### Status Badge
```html
<div style="background:#4CAF50;padding:3px 10px;border-radius:4px;">
  <p style="font-size:8pt;color:#fff;margin:0;font-weight:bold;">PASS</p>
</div>
```

## Build Tooling

### Shared Infrastructure
- **html2pptx/**: Shared build library at `voyager_slide_templates/html2pptx/`. All `build.js` files import from `../../html2pptx/index.cjs`.
- **assets/**: Shared Voyager branding PNGs at `voyager_slide_templates/assets/`. Decks access these via an `assets` symlink pointing to `../../assets`.
- **deck-assets/**: Per-deck folder for deck-specific images (diagrams, custom backgrounds). Not all decks have one.

### Deck Structure
```
voyager_slide_templates/
├── CLAUDE.md                      ← this file (authoritative rules)
├── PRESENTATION_REGISTRY.md       ← PPPxxx manifest
├── html2pptx/                     ← shared build library
├── assets/                        ← shared Voyager branding PNGs
├── css/                           ← CSS gallery (visual reference ONLY, do NOT use)
├── templates/                     ← HTML template gallery (visual reference ONLY)
└── decks/
    └── PPPxxx-descriptive-name/
        ├── build.js               ← requires('../../html2pptx/index.cjs')
        ├── package.json
        ├── assets/                ← symlink → ../../assets/
        ├── deck-assets/           ← deck-specific images (optional)
        ├── slide01-title.html
        ├── slide02-*.html
        ├── ...
        ├── index.html             ← browser preview gallery (optional)
        └── *.pptx                 ← built output
```

### Build Commands
```bash
cd decks/PPPxxx-name/
npm install          # first time only
node build.js        # produces .pptx file
```

### Pre-Build Checklist
- [ ] Every slide is a bare `<body>` element (no DOCTYPE/html/head)
- [ ] All colors are inline hex (no CSS variables)
- [ ] All branding assets are `<img>` tags (no CSS background-image)
- [ ] Every CSS gradient has `data-rasterize="gradient"`
- [ ] Logo is absolute-positioned outside the header flex row
- [ ] All visible text is in `<p>`, `<h1>`-`<h6>`, `<ul>`, or `<ol>` — no bare `<span>` elements
- [ ] Footer has slide number and CONFIDENTIAL text
- [ ] Assets symlink points to `../../assets/`
- [ ] `build.js` imports from `../../html2pptx/index.cjs`
