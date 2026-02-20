"""Unified registry utility for FIG###, PPP###, and AVD### barcoded asset banks.

Parses markdown-table registries, assigns next available numbers, and audits
bank directories for consistency. Zero external dependencies (stdlib only).
"""

import re
import xml.etree.ElementTree as ET
from pathlib import Path

# Matches a table row whose first data cell starts with PREFIX + digits.
# Example: "| FIG001 | FIG001-vr4-loop-mosaic.svg | SVG | ..."
_ROW_RE = re.compile(
    r"^\|\s*(?P<barcode>[A-Z]+\d+)\s*\|"  # first cell: barcode
    r"\s*(?P<rest>.+)"                      # remaining cells
)


def parse_registry(registry_path, prefix):
    """Parse a markdown registry, returning {number_int: {barcode, filename, raw_row}}.

    Scans every line for table rows whose first cell matches *prefix* followed
    by digits (e.g. ``FIG001``, ``AVD133``, ``PPP011``).  The second cell is
    taken as the filename/folder.

    Parameters
    ----------
    registry_path : str or Path
        Path to the registry markdown file.
    prefix : str
        Barcode prefix to match (e.g. "FIG", "AVD", "PPP").

    Returns
    -------
    dict[int, dict]
        Keyed by the numeric part of the barcode.  Each value contains:
        ``barcode`` (str), ``filename`` (str), ``raw_row`` (str).
    """
    path = Path(registry_path)
    if not path.exists():
        return {}

    entries = {}
    prefix_upper = prefix.upper()

    for line in path.read_text().splitlines():
        m = _ROW_RE.match(line.strip())
        if not m:
            continue
        barcode = m.group("barcode").strip()
        if not barcode.upper().startswith(prefix_upper):
            continue
        # Extract numeric part
        num_str = barcode[len(prefix_upper):]
        if not num_str.isdigit():
            continue
        num = int(num_str)

        # Second cell = filename/folder
        rest = m.group("rest")
        cells = [c.strip() for c in rest.split("|")]
        filename = cells[0] if cells else ""

        # Allow multiple files per barcode (e.g. AVD002 has 3 format variants)
        # Keep the first occurrence as the canonical entry
        if num not in entries:
            entries[num] = {
                "barcode": barcode,
                "filename": filename,
                "raw_row": line.strip(),
            }

    return entries


def get_next_number(registry_path, prefix, zero_pad=3):
    """Return the next available barcode string (e.g. ``FIG007``).

    Parameters
    ----------
    registry_path : str or Path
        Path to the registry markdown file.
    prefix : str
        Barcode prefix (e.g. "FIG").
    zero_pad : int
        Minimum digit width (default 3 → ``FIG007``).

    Returns
    -------
    str
        Next barcode, e.g. ``"FIG007"``.
    """
    entries = parse_registry(registry_path, prefix)
    if not entries:
        return f"{prefix.upper()}{1:0{zero_pad}d}"
    highest = max(entries.keys())
    return f"{prefix.upper()}{highest + 1:0{zero_pad}d}"


def validate_bank(registry_path, bank_dir, prefix, file_glob="*"):
    """Audit an asset bank for consistency.

    Returns a list of issue dicts, each with ``type``, ``detail``, and
    ``severity`` keys.  An empty list means no issues found.

    Issue types
    -----------
    - ``orphan_file``: file in *bank_dir* matching *prefix* but not in registry
    - ``missing_file``: entry in registry whose filename is absent from *bank_dir*
    - ``duplicate_number``: same barcode number appears on multiple distinct rows
      (only detectable when the registry has been hand-edited; ``parse_registry``
      deduplicates so this check re-scans the raw file)

    Parameters
    ----------
    registry_path : str or Path
        Path to the registry markdown file.
    bank_dir : str or Path
        Directory containing the asset files.
    prefix : str
        Barcode prefix (e.g. "FIG", "AVD", "PPP").
    file_glob : str
        Glob pattern for files to consider (default ``"*"``).
    """
    issues = []
    reg_path = Path(registry_path)
    bank = Path(bank_dir)
    prefix_upper = prefix.upper()

    entries = parse_registry(registry_path, prefix)

    # --- Collect all filenames mentioned in registry ---
    registry_filenames = set()
    for entry in entries.values():
        if entry["filename"]:
            registry_filenames.add(entry["filename"])
    # Also collect filenames from duplicate rows (parse_registry keeps first only)
    if reg_path.exists():
        for line in reg_path.read_text().splitlines():
            m = _ROW_RE.match(line.strip())
            if not m:
                continue
            barcode = m.group("barcode").strip()
            if not barcode.upper().startswith(prefix_upper):
                continue
            rest = m.group("rest")
            cells = [c.strip() for c in rest.split("|")]
            if cells and cells[0]:
                registry_filenames.add(cells[0])

    # --- Orphan files: in bank but not in registry ---
    if bank.is_dir():
        for f in sorted(bank.glob(file_glob)):
            if f.is_file() and f.name.upper().startswith(prefix_upper) and f.name not in registry_filenames:
                # Skip the registry file itself
                if f.suffix.lower() == ".md":
                    continue
                issues.append({
                    "type": "orphan_file",
                    "detail": f"File '{f.name}' is in bank but not in registry",
                    "severity": "WARNING",
                })

    # --- Missing files: in registry but not in bank ---
    if bank.is_dir():
        bank_files = {f.name for f in bank.iterdir() if f.is_file()}
        for entry in entries.values():
            fn = entry["filename"]
            if fn and fn not in bank_files:
                issues.append({
                    "type": "missing_file",
                    "detail": f"Registry entry '{entry['barcode']}' references '{fn}' but file not found in bank",
                    "severity": "ERROR",
                })

    # --- Duplicate numbers: same prefix+number on multiple rows ---
    if reg_path.exists():
        seen_numbers = {}
        for line in reg_path.read_text().splitlines():
            m = _ROW_RE.match(line.strip())
            if not m:
                continue
            barcode = m.group("barcode").strip()
            if not barcode.upper().startswith(prefix_upper):
                continue
            num_str = barcode[len(prefix_upper):]
            if not num_str.isdigit():
                continue
            num = int(num_str)
            # Extract filename from second cell
            rest = m.group("rest")
            cells = [c.strip() for c in rest.split("|")]
            filename = cells[0] if cells else ""
            if num in seen_numbers:
                # Duplicate number — but same AVD can have multiple format variants
                # Only flag if the filenames differ AND it's not an expected multi-format case
                prev_fn = seen_numbers[num]
                if filename != prev_fn:
                    # Check if they share the same barcode prefix (e.g. AVD002 with .dna and .gb)
                    # This is allowed per CLAUDE.md rules
                    prev_base = prev_fn.split("-", 1)[0] if prev_fn else ""
                    curr_base = filename.split("-", 1)[0] if filename else ""
                    if prev_base == curr_base:
                        # Same construct, different format — allowed
                        continue
                issues.append({
                    "type": "duplicate_number",
                    "detail": f"Number {barcode} appears multiple times with different names",
                    "severity": "ERROR",
                })
            else:
                seen_numbers[num] = filename

    return issues


# ---------------------------------------------------------------------------
# SVG validation
# ---------------------------------------------------------------------------

# HTML entities that are invalid in SVG/XML — must use numeric references instead.
_BAD_ENTITIES = {
    "&times;": "&#215;",
    "&mdash;": "&#8212;",
    "&rarr;": "&#8594;",
    "&nbsp;": "&#160;",
    "&larr;": "&#8592;",
    "&ndash;": "&#8211;",
    "&bull;": "&#8226;",
}

_BAD_ENTITY_RE = re.compile(
    r"&(?:" + "|".join(k[1:-1] for k in _BAD_ENTITIES) + r");"
)

_FONT_SIZE_RE = re.compile(
    r'font-size\s*[:=]\s*["\']?\s*(\d+(?:\.\d+)?)\s*(?:px|pt)?',
    re.IGNORECASE,
)

_VIEWBOX_RE = re.compile(
    r'viewBox\s*=\s*["\']'
    r'\s*(-?\d+(?:\.\d+)?)\s+'
    r'(-?\d+(?:\.\d+)?)\s+'
    r'(-?\d+(?:\.\d+)?)\s+'
    r'(-?\d+(?:\.\d+)?)\s*'
    r'["\']',
)


def validate_svg(svg_path):
    """Check an SVG file for common anti-patterns.

    Returns a list of issue dicts, each with ``type``, ``detail``, and
    ``severity`` keys.  An empty list means no issues found.

    Issue types
    -----------
    - ``bad_entity``: HTML entity that is invalid in SVG XML
    - ``missing_viewbox``: No ``viewBox`` attribute on root ``<svg>``
    - ``xml_parse_error``: File is not well-formed XML
    - ``small_font``: Font size below 12 (likely too small for slides)

    Parameters
    ----------
    svg_path : str or Path
        Path to the SVG file to validate.
    """
    issues = []
    path = Path(svg_path)

    if not path.exists():
        issues.append({
            "type": "file_not_found",
            "detail": f"SVG file not found: {path}",
            "severity": "ERROR",
        })
        return issues

    raw = path.read_text(encoding="utf-8")

    # --- Bad entities (check raw text BEFORE XML parse) ---
    for m in _BAD_ENTITY_RE.finditer(raw):
        entity = m.group(0)
        replacement = _BAD_ENTITIES.get(entity, "numeric reference")
        issues.append({
            "type": "bad_entity",
            "detail": f"Invalid SVG entity '{entity}' — use '{replacement}' instead",
            "severity": "ERROR",
        })

    # --- XML well-formedness ---
    try:
        root = ET.fromstring(raw)
    except ET.ParseError as e:
        issues.append({
            "type": "xml_parse_error",
            "detail": f"SVG is not well-formed XML: {e}",
            "severity": "ERROR",
        })
        # Can't do further checks if XML is broken
        return issues

    # --- viewBox check ---
    viewbox = root.get("viewBox")
    if viewbox is None:
        issues.append({
            "type": "missing_viewbox",
            "detail": "Root <svg> has no viewBox attribute",
            "severity": "WARNING",
        })
    else:
        if not _VIEWBOX_RE.search(f'viewBox="{viewbox}"'):
            issues.append({
                "type": "missing_viewbox",
                "detail": f"viewBox value '{viewbox}' does not have 4 numeric values",
                "severity": "WARNING",
            })

    # --- Font size check (scan raw text for both attribute and CSS forms) ---
    for m in _FONT_SIZE_RE.finditer(raw):
        size = float(m.group(1))
        if size < 12:
            issues.append({
                "type": "small_font",
                "detail": f"Font size {m.group(1)} is below minimum (12); consider >=14 for annotations, >=16 for primary text",
                "severity": "WARNING",
            })

    return issues
