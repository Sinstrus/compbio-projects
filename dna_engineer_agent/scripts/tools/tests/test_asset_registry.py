"""Tests for asset_registry.py — unified registry utility."""

import os
import tempfile
from pathlib import Path

import pytest

# Allow imports from parent directory
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from asset_registry import parse_registry, get_next_number, validate_bank, validate_svg


# ---------------------------------------------------------------------------
# Paths to live registries
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent  # dna_engineer_agent/
FIGURE_REGISTRY = REPO_ROOT.parent / "figure_bank" / "FIGURE_REGISTRY.md"
CONSTRUCT_REGISTRY = REPO_ROOT / "CONSTRUCT_REGISTRY.md"
PRESENTATION_REGISTRY = REPO_ROOT.parent / "voyager_slide_decks" / "PRESENTATION_REGISTRY.md"


# ---------------------------------------------------------------------------
# parse_registry
# ---------------------------------------------------------------------------
class TestParseRegistry:
    def test_parse_figure_registry(self):
        """Parses real FIGURE_REGISTRY.md and returns FIG entries."""
        if not FIGURE_REGISTRY.exists():
            pytest.skip("FIGURE_REGISTRY.md not found")
        entries = parse_registry(FIGURE_REGISTRY, "FIG")
        assert len(entries) >= 6, f"Expected >= 6 FIG entries, got {len(entries)}"
        assert 1 in entries
        assert entries[1]["barcode"] == "FIG001"
        assert "FIG001" in entries[1]["filename"]

    def test_parse_construct_registry(self):
        """Parses real CONSTRUCT_REGISTRY.md with multiple table sections."""
        if not CONSTRUCT_REGISTRY.exists():
            pytest.skip("CONSTRUCT_REGISTRY.md not found")
        entries = parse_registry(CONSTRUCT_REGISTRY, "AVD")
        # Should find AVD001 through at least AVD133
        assert len(entries) >= 100, f"Expected >= 100 AVD entries, got {len(entries)}"
        assert 1 in entries
        assert entries[1]["barcode"] == "AVD001"
        assert 133 in entries

    def test_parse_presentation_registry(self):
        """Parses real PRESENTATION_REGISTRY.md."""
        if not PRESENTATION_REGISTRY.exists():
            pytest.skip("PRESENTATION_REGISTRY.md not found")
        entries = parse_registry(PRESENTATION_REGISTRY, "PPP")
        assert len(entries) >= 11, f"Expected >= 11 PPP entries, got {len(entries)}"
        assert 1 in entries
        assert entries[1]["barcode"] == "PPP001"

    def test_parse_nonexistent_file(self):
        """Returns empty dict for missing file."""
        entries = parse_registry("/nonexistent/path/REGISTRY.md", "FIG")
        assert entries == {}

    def test_parse_ignores_header_rows(self):
        """Table header/separator rows are not parsed as entries."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".md", delete=False) as f:
            f.write("| FIG | Filename |\n")
            f.write("|-----|----------|\n")
            f.write("| FIG001 | test.svg |\n")
            f.name
        try:
            entries = parse_registry(f.name, "FIG")
            assert len(entries) == 1
            assert 1 in entries
        finally:
            os.unlink(f.name)


# ---------------------------------------------------------------------------
# get_next_number
# ---------------------------------------------------------------------------
class TestGetNextNumber:
    def test_get_next_fig_number(self):
        """Next FIG number is one past the highest entry in the registry."""
        if not FIGURE_REGISTRY.exists():
            pytest.skip("FIGURE_REGISTRY.md not found")
        entries = parse_registry(FIGURE_REGISTRY, "FIG")
        highest = max(entries.keys())
        next_num = get_next_number(FIGURE_REGISTRY, "FIG")
        assert next_num == f"FIG{highest + 1:03d}"

    def test_get_next_number_empty_registry(self):
        """Returns PREFIX001 for empty registry."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".md", delete=False) as f:
            f.write("# Empty Registry\n\n| FIG | Filename |\n|-----|----------|\n")
            f.name
        try:
            assert get_next_number(f.name, "FIG") == "FIG001"
        finally:
            os.unlink(f.name)

    def test_get_next_number_nonexistent_file(self):
        """Returns PREFIX001 when registry file doesn't exist."""
        assert get_next_number("/nonexistent/REGISTRY.md", "FIG") == "FIG001"

    def test_get_next_number_custom_padding(self):
        """Respects zero_pad parameter."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".md", delete=False) as f:
            f.write("| PPP01 | folder1 |\n")
            f.name
        try:
            assert get_next_number(f.name, "PPP", zero_pad=2) == "PPP02"
        finally:
            os.unlink(f.name)

    def test_get_next_avd_number(self):
        """Returns correct next AVD number from real construct registry."""
        if not CONSTRUCT_REGISTRY.exists():
            pytest.skip("CONSTRUCT_REGISTRY.md not found")
        next_num = get_next_number(CONSTRUCT_REGISTRY, "AVD")
        # Should be AVD134 or higher (133 entries exist)
        num = int(next_num[3:])
        assert num >= 134


# ---------------------------------------------------------------------------
# validate_bank
# ---------------------------------------------------------------------------
class TestValidateBank:
    def test_validate_bank_clean(self):
        """No issues on current figure_bank (golden path)."""
        if not FIGURE_REGISTRY.exists():
            pytest.skip("FIGURE_REGISTRY.md not found")
        bank_dir = FIGURE_REGISTRY.parent
        issues = validate_bank(FIGURE_REGISTRY, bank_dir, "FIG")
        assert issues == [], f"Expected no issues, got: {issues}"

    def test_validate_bank_finds_orphans(self):
        """Detects files in bank that aren't in the registry."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a registry with one entry
            reg = Path(tmpdir) / "REGISTRY.md"
            reg.write_text(
                "| FIG | Filename |\n"
                "|-----|----------|\n"
                "| FIG001 | FIG001-test.svg |\n"
            )
            # Create the registered file plus an orphan
            (Path(tmpdir) / "FIG001-test.svg").write_text("<svg/>")
            (Path(tmpdir) / "FIG002-orphan.svg").write_text("<svg/>")

            issues = validate_bank(reg, tmpdir, "FIG")
            orphans = [i for i in issues if i["type"] == "orphan_file"]
            assert len(orphans) == 1
            assert "FIG002-orphan.svg" in orphans[0]["detail"]

    def test_validate_bank_finds_missing(self):
        """Detects registry entries whose files are absent from bank."""
        with tempfile.TemporaryDirectory() as tmpdir:
            reg = Path(tmpdir) / "REGISTRY.md"
            reg.write_text(
                "| FIG | Filename |\n"
                "|-----|----------|\n"
                "| FIG001 | FIG001-test.svg |\n"
                "| FIG002 | FIG002-missing.svg |\n"
            )
            # Only create FIG001
            (Path(tmpdir) / "FIG001-test.svg").write_text("<svg/>")

            issues = validate_bank(reg, tmpdir, "FIG")
            missing = [i for i in issues if i["type"] == "missing_file"]
            assert len(missing) == 1
            assert "FIG002-missing.svg" in missing[0]["detail"]

    def test_duplicate_number_detected(self):
        """Catches duplicate FIG numbers with different base names in registry."""
        with tempfile.TemporaryDirectory() as tmpdir:
            reg = Path(tmpdir) / "REGISTRY.md"
            reg.write_text(
                "| FIG | Filename |\n"
                "|-----|----------|\n"
                "| FIG001 | FIG001-original.svg |\n"
                "| FIG001 | totally-different-name.svg |\n"
            )
            (Path(tmpdir) / "FIG001-original.svg").write_text("<svg/>")
            (Path(tmpdir) / "totally-different-name.svg").write_text("<svg/>")

            issues = validate_bank(reg, tmpdir, "FIG")
            dupes = [i for i in issues if i["type"] == "duplicate_number"]
            assert len(dupes) == 1
            assert "FIG001" in dupes[0]["detail"]

    def test_avd_multi_format_not_flagged(self):
        """AVD entries with multiple format variants (same base name) are allowed."""
        with tempfile.TemporaryDirectory() as tmpdir:
            reg = Path(tmpdir) / "REGISTRY.md"
            reg.write_text(
                "| AVD | Filename |\n"
                "|-----|----------|\n"
                "| AVD002 | AVD002-Rep2Mut2Cap9-6R-wt.dna |\n"
                "| AVD002 | AVD002-Rep2Mut2Cap9-6R-wt-5fold.gb |\n"
            )
            (Path(tmpdir) / "AVD002-Rep2Mut2Cap9-6R-wt.dna").write_text("")
            (Path(tmpdir) / "AVD002-Rep2Mut2Cap9-6R-wt-5fold.gb").write_text("")

            issues = validate_bank(reg, tmpdir, "AVD")
            dupes = [i for i in issues if i["type"] == "duplicate_number"]
            assert len(dupes) == 0, f"Multi-format variants should not be flagged: {dupes}"


# ---------------------------------------------------------------------------
# validate_svg
# ---------------------------------------------------------------------------
class TestValidateSvg:
    """Tests for SVG anti-pattern detection."""

    VALID_SVG = (
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 900 420">'
        '<rect width="100" height="100" fill="red"/>'
        '<text font-size="16">Hello</text>'
        '</svg>'
    )

    def test_valid_svg_no_issues(self):
        """A well-formed SVG with valid entities and viewBox produces no issues."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".svg", delete=False) as f:
            f.write(self.VALID_SVG)
        try:
            issues = validate_svg(f.name)
            assert issues == [], f"Expected no issues, got: {issues}"
        finally:
            os.unlink(f.name)

    def test_bad_entity_times(self):
        """Detects &times; (should be &#215;)."""
        svg = (
            '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100">'
            '<text>&times;</text></svg>'
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".svg", delete=False) as f:
            f.write(svg)
        try:
            issues = validate_svg(f.name)
            bad = [i for i in issues if i["type"] == "bad_entity"]
            assert len(bad) >= 1
            assert "&#215;" in bad[0]["detail"]
        finally:
            os.unlink(f.name)

    def test_bad_entity_mdash(self):
        """Detects &mdash; (should be &#8212;)."""
        svg = (
            '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100">'
            '<text>&mdash;</text></svg>'
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".svg", delete=False) as f:
            f.write(svg)
        try:
            issues = validate_svg(f.name)
            bad = [i for i in issues if i["type"] == "bad_entity"]
            assert len(bad) >= 1
            assert "&#8212;" in bad[0]["detail"]
        finally:
            os.unlink(f.name)

    def test_bad_entity_rarr(self):
        """Detects &rarr; (should be &#8594;)."""
        svg = (
            '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100">'
            '<text>&rarr;</text></svg>'
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".svg", delete=False) as f:
            f.write(svg)
        try:
            issues = validate_svg(f.name)
            bad = [i for i in issues if i["type"] == "bad_entity"]
            assert len(bad) >= 1
            assert "&#8594;" in bad[0]["detail"]
        finally:
            os.unlink(f.name)

    def test_missing_viewbox(self):
        """Warns when root <svg> has no viewBox attribute."""
        svg = (
            '<svg xmlns="http://www.w3.org/2000/svg">'
            '<rect width="100" height="100"/></svg>'
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".svg", delete=False) as f:
            f.write(svg)
        try:
            issues = validate_svg(f.name)
            vb = [i for i in issues if i["type"] == "missing_viewbox"]
            assert len(vb) == 1
            assert vb[0]["severity"] == "WARNING"
        finally:
            os.unlink(f.name)

    def test_xml_parse_error(self):
        """Detects malformed XML (unclosed tag)."""
        svg = '<svg xmlns="http://www.w3.org/2000/svg"><rect></svg>'
        with tempfile.NamedTemporaryFile(mode="w", suffix=".svg", delete=False) as f:
            f.write(svg)
        try:
            issues = validate_svg(f.name)
            parse_err = [i for i in issues if i["type"] == "xml_parse_error"]
            assert len(parse_err) == 1
            assert parse_err[0]["severity"] == "ERROR"
        finally:
            os.unlink(f.name)

    def test_small_font_detected(self):
        """Warns on font-size below 12."""
        svg = (
            '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100">'
            '<text font-size="10">Tiny</text>'
            '<text font-size="16">Normal</text>'
            '</svg>'
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".svg", delete=False) as f:
            f.write(svg)
        try:
            issues = validate_svg(f.name)
            small = [i for i in issues if i["type"] == "small_font"]
            assert len(small) == 1
            assert "10" in small[0]["detail"]
        finally:
            os.unlink(f.name)

    def test_file_not_found(self):
        """Returns error for nonexistent file."""
        issues = validate_svg("/nonexistent/path/test.svg")
        assert len(issues) == 1
        assert issues[0]["type"] == "file_not_found"
        assert issues[0]["severity"] == "ERROR"

    def test_multiple_bad_entities(self):
        """Detects multiple different bad entities in one file."""
        svg = (
            '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100">'
            '<text>&times; &rarr; &mdash;</text></svg>'
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".svg", delete=False) as f:
            f.write(svg)
        try:
            issues = validate_svg(f.name)
            bad = [i for i in issues if i["type"] == "bad_entity"]
            assert len(bad) == 3
        finally:
            os.unlink(f.name)

    def test_numeric_entities_ok(self):
        """Numeric character references (&#215; etc.) should NOT be flagged."""
        svg = (
            '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100">'
            '<text>&#215; &#8594; &#8212;</text></svg>'
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".svg", delete=False) as f:
            f.write(svg)
        try:
            issues = validate_svg(f.name)
            bad = [i for i in issues if i["type"] == "bad_entity"]
            assert len(bad) == 0
        finally:
            os.unlink(f.name)
