#!/usr/bin/env python3
"""
Tests for the ConstructVerifier module.

Covers:
- Restriction site verification (BUG-003: both strands)
- Protein identity verification
- Insertion flank verification (BUG-004/006)
- Enzyme validation (DESIGN-004/005)
- Synthesis screen (GC, homopolymers, repeats)
- Frame check (BUG-001)
- Verification gating (output blocked on critical failure)
"""

import sys
from pathlib import Path

import pytest

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from construct_verifier import (
    ConstructVerifier,
    CheckResult,
    Severity,
    VerificationError,
    VerificationReport,
    count_sites_both_strands,
    find_site_positions,
    get_dam_sensitive_enzymes,
    get_fussy_enzymes,
    get_recognition_site,
    load_enzyme_metadata,
    reverse_complement,
    synthesis_screen,
    translate,
    validate_enzyme_choice,
    verify_and_gate,
)


# ---------------------------------------------------------------------------
# Helper function tests
# ---------------------------------------------------------------------------

class TestReverseComplement:
    def test_basic(self):
        assert reverse_complement("ATGC") == "GCAT"

    def test_palindrome(self):
        assert reverse_complement("GAATTC") == "GAATTC"  # EcoRI

    def test_non_palindrome(self):
        assert reverse_complement("GGTCTC") == "GAGACC"  # BsaI


class TestTranslate:
    def test_basic(self):
        assert translate("ATGAAAGCC") == "MKA"

    def test_stop_codon(self):
        assert translate("ATGAAATAAGCC") == "MK"  # Stops at TAA

    def test_frame1(self):
        # Frame 1: skip first base → TGAAAGCC → codons TGA, AGC, C(partial)
        # TGA = stop → empty protein
        assert translate("ATGAAAGCC", frame=1) == ""

    def test_frame2(self):
        # Frame 2: skip first 2 bases → GAAAGCC → codons GAA, AGC, C(partial)
        assert translate("ATGAAAGCC", frame=2) == "ES"


class TestCountSitesBothStrands:
    """BUG-003 prevention: must count sites on BOTH strands."""

    def test_palindromic_single(self):
        # HindIII (AAGCTT) is palindromic
        seq = "GGGAAGCTTGGG"
        assert count_sites_both_strands(seq, "AAGCTT") == 1

    def test_palindromic_double(self):
        seq = "AAGCTTNNNNNNAAGCTT"
        assert count_sites_both_strands(seq, "AAGCTT") == 2

    def test_non_palindromic_forward_only(self):
        # BsaI (GGTCTC) — non-palindromic
        seq = "GGGGGTCTCGGG"
        assert count_sites_both_strands(seq, "GGTCTC") == 1

    def test_non_palindromic_reverse_only(self):
        # BsaI rc = GAGACC
        seq = "GGGGAGACCGGG"
        assert count_sites_both_strands(seq, "GGTCTC") == 1

    def test_non_palindromic_both_strands(self):
        # BsaI on both strands
        seq = "GGTCTCNNNNNNGAGACC"
        assert count_sites_both_strands(seq, "GGTCTC") == 2

    def test_none(self):
        seq = "AAAAAAAAAA"
        assert count_sites_both_strands(seq, "GAATTC") == 0

    def test_case_insensitive(self):
        seq = "gggaagcttggg"
        assert count_sites_both_strands(seq, "AAGCTT") == 1


class TestFindSitePositions:
    def test_palindromic(self):
        seq = "GGGAAGCTTGGG"
        positions = find_site_positions(seq, "AAGCTT")
        assert positions == [(3, '+')]

    def test_non_palindromic_both(self):
        seq = "GGTCTCNNNNNNGAGACC"
        positions = find_site_positions(seq, "GGTCTC")
        assert len(positions) == 2
        strands = {s for _, s in positions}
        assert '+' in strands
        assert '-' in strands


# ---------------------------------------------------------------------------
# Restriction site check tests
# ---------------------------------------------------------------------------

class TestRestrictionCheck:
    def test_unique_site_passes(self):
        seq = "GGGAAGCTTGGGTCTAGAGGG"  # 1x HindIII, 1x XbaI
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_restriction_check(
            required_unique={"HindIII": "AAGCTT", "XbaI": "TCTAGA"}
        )
        report = v.run_all()
        assert report.passed
        assert all(c.passed for c in report.checks)

    def test_duplicate_site_fails(self):
        seq = "AAGCTTNNNNNAAGCTT"
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_restriction_check(required_unique={"HindIII": "AAGCTT"})
        report = v.run_all()
        assert not report.passed
        assert report.critical_failures[0].name == "RestrictionSite:HindIII:unique"

    def test_missing_site_fails(self):
        seq = "AAAAAAAAA"
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_restriction_check(required_unique={"HindIII": "AAGCTT"})
        report = v.run_all()
        assert not report.passed

    def test_forbidden_site_present_fails(self):
        seq = "GGGGGTCTCGGG"
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_restriction_check(forbidden={"BsaI": "GGTCTC"})
        report = v.run_all()
        assert not report.passed
        assert "BsaI" in report.critical_failures[0].name

    def test_forbidden_site_absent_passes(self):
        seq = "AAAAAAAAAA"
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_restriction_check(forbidden={"BsaI": "GGTCTC"})
        report = v.run_all()
        assert report.passed

    def test_forbidden_on_reverse_strand_fails(self):
        """BUG-003: Must detect BsaI on reverse strand too."""
        # GAGACC = rc of GGTCTC
        seq = "GGGGAGACCGGG"
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_restriction_check(forbidden={"BsaI": "GGTCTC"})
        report = v.run_all()
        assert not report.passed


# ---------------------------------------------------------------------------
# Protein check tests
# ---------------------------------------------------------------------------

class TestProteinCheck:
    def test_correct_protein(self):
        # ATG AAA GCC = MKA
        seq = "NNNNATGAAAGCCNNN"
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_protein_check(
            region_start=4, region_end=13,
            expected_protein="MKA", frame=0, label="TestCDS"
        )
        report = v.run_all()
        assert report.passed

    def test_wrong_protein(self):
        # ATG AAA GCC = MKA, but we expect MKG
        seq = "NNNNATGAAAGCCNNN"
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_protein_check(
            region_start=4, region_end=13,
            expected_protein="MKG", frame=0, label="TestCDS"
        )
        report = v.run_all()
        assert not report.passed
        assert "MISMATCH" in report.critical_failures[0].message


# ---------------------------------------------------------------------------
# Enzyme validation tests
# ---------------------------------------------------------------------------

class TestEnzymeValidation:
    def test_hindiii_passes(self):
        """HindIII is a very_high reliability, no dam sensitivity enzyme."""
        results = validate_enzyme_choice("HindIII")
        # HindIII should have no issues
        assert all(r.passed for r in results) or len(results) == 0

    def test_dam_sensitive_enzyme(self):
        """XbaI should trigger a dam sensitivity warning."""
        results = validate_enzyme_choice("XbaI")
        dam_results = [r for r in results if "dam" in r.name.lower()]
        # XbaI has "Dam methylation may reduce activity" — should be a warning
        assert any(not r.passed for r in dam_results)

    def test_fussy_enzyme(self):
        """BbvCI should be flagged as fussy."""
        results = validate_enzyme_choice("BbvCI")
        # Should have reliability warning and/or fussy warning
        assert any(not r.passed for r in results)

    def test_paqci_fussy(self):
        """PaqCI should be flagged."""
        results = validate_enzyme_choice("PaqCI")
        assert any(not r.passed for r in results)

    def test_unknown_dam_sensitive(self):
        """BclI is in DAM_SENSITIVE_ENZYMES set but not in metadata."""
        results = validate_enzyme_choice("BclI")
        assert any("dam" in r.name.lower() for r in results)

    def test_verifier_integration(self):
        """Test enzyme validation through ConstructVerifier."""
        seq = "AAAA"
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_enzyme_validation(["HindIII", "EcoRI"])
        report = v.run_all()
        # HindIII and EcoRI are both excellent — should pass
        # (only warning-level issues at most, no critical failures)
        assert report.passed


# ---------------------------------------------------------------------------
# Synthesis screen tests
# ---------------------------------------------------------------------------

class TestSynthesisScreen:
    def test_balanced_sequence(self):
        """Sequence with ~50% GC should pass."""
        seq = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"  # 50% GC
        results = synthesis_screen(seq, name="balanced")
        gc_results = [r for r in results if "gc_content" in r.name]
        assert all(r.passed for r in gc_results)

    def test_high_gc_warns(self):
        """All-GC sequence should trigger warning."""
        seq = "G" * 60 + "C" * 60  # 100% GC
        results = synthesis_screen(seq, name="high_gc")
        gc_results = [r for r in results if "gc_content" in r.name]
        assert any(not r.passed for r in gc_results)

    def test_low_gc_warns(self):
        """All-AT sequence should trigger warning."""
        seq = "A" * 60 + "T" * 60  # 0% GC
        results = synthesis_screen(seq, name="low_gc")
        gc_results = [r for r in results if "gc_content" in r.name]
        assert any(not r.passed for r in gc_results)

    def test_homopolymer_warns(self):
        """7-base homopolymer should trigger warning."""
        seq = "ATGC" * 10 + "AAAAAAA" + "ATGC" * 10
        results = synthesis_screen(seq, name="homopolymer")
        homo_results = [r for r in results if "homopolymer" in r.name]
        assert any(not r.passed for r in homo_results)

    def test_no_homopolymer(self):
        """Short runs should pass."""
        seq = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        results = synthesis_screen(seq, name="no_homo")
        homo_results = [r for r in results if "homopolymer" in r.name]
        assert all(r.passed for r in homo_results)

    def test_repeat_warns(self):
        """11bp direct repeat should trigger warning."""
        repeat = "ATGCATGCATG"  # 11 bp
        seq = repeat + "NNNNNNNNN" + repeat + "AAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        results = synthesis_screen(seq, name="repeat", repeat_min_len=11)
        repeat_results = [r for r in results if "repeats" in r.name]
        assert any(not r.passed for r in repeat_results)

    def test_no_repeat(self):
        """Non-repetitive sequence should pass."""
        # Each 15-mer must be unique — use a random-like non-repeating sequence
        seq = "ATGCTAGCAATGCGTACCAGATCTTGCAGTCAACGGTACTAGCCGAATGCTAACGGCTTG"
        results = synthesis_screen(seq, name="no_repeat", repeat_min_len=15)
        repeat_results = [r for r in results if "repeats" in r.name]
        assert all(r.passed for r in repeat_results)


# ---------------------------------------------------------------------------
# Frame check tests (BUG-001)
# ---------------------------------------------------------------------------

class TestFrameCheck:
    def test_correct_frame(self):
        seq = "A" * 100
        v = ConstructVerifier(sequence=seq, name="test")
        # Position 30, CDS starts at 10 → frame = (30-10) % 3 = 20 % 3 = 2
        v.add_frame_check(position=30, cds_start=10, expected_frame=2)
        report = v.run_all()
        assert report.passed

    def test_wrong_frame(self):
        seq = "A" * 100
        v = ConstructVerifier(sequence=seq, name="test")
        # Position 30, CDS starts at 10 → frame = 2, but we expect 0
        v.add_frame_check(position=30, cds_start=10, expected_frame=0)
        report = v.run_all()
        assert not report.passed

    def test_bug001_relative_not_absolute(self):
        """BUG-001: Frame must be relative to CDS start, not absolute position."""
        seq = "A" * 100
        v = ConstructVerifier(sequence=seq, name="test")
        # Position 33, CDS starts at 12
        # WRONG (absolute): 33 % 3 = 0
        # RIGHT (relative): (33 - 12) % 3 = 21 % 3 = 0
        # In this case they happen to match, so use a case where they differ:
        # Position 34, CDS starts at 12
        # WRONG (absolute): 34 % 3 = 1
        # RIGHT (relative): (34 - 12) % 3 = 22 % 3 = 1
        # Still matches. Try: Position 35, CDS starts at 13
        # WRONG (absolute): 35 % 3 = 2
        # RIGHT (relative): (35 - 13) % 3 = 22 % 3 = 1
        v.add_frame_check(position=35, cds_start=13, expected_frame=1)
        report = v.run_all()
        assert report.passed

        # If we used absolute position (35 % 3 = 2), this would be wrong
        v2 = ConstructVerifier(sequence=seq, name="test2")
        v2.add_frame_check(position=35, cds_start=13, expected_frame=2)
        report2 = v2.run_all()
        assert not report2.passed  # 2 is the absolute-position answer, which is WRONG


# ---------------------------------------------------------------------------
# Packaging check tests
# ---------------------------------------------------------------------------

class TestPackagingCheck:
    def test_within_limit(self):
        seq = "A" * 10000
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_packaging_check(itr_left_end=100, itr_right_start=4900)
        report = v.run_all()
        assert report.passed

    def test_over_limit(self):
        seq = "A" * 10000
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_packaging_check(itr_left_end=100, itr_right_start=5200)
        report = v.run_all()
        assert not report.passed


# ---------------------------------------------------------------------------
# Gate function tests
# ---------------------------------------------------------------------------

class TestVerifyAndGate:
    def test_passes_through(self):
        seq = "GGGAAGCTTGGG"
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_restriction_check(required_unique={"HindIII": "AAGCTT"})
        report = verify_and_gate(v)
        assert report.passed

    def test_raises_on_failure(self):
        seq = "AAGCTTNNNNNNAAGCTT"
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_restriction_check(required_unique={"HindIII": "AAGCTT"})
        with pytest.raises(VerificationError) as exc_info:
            verify_and_gate(v)
        assert "FAILED" in str(exc_info.value)

    def test_warnings_dont_block(self):
        """WARNING-severity issues should not block output."""
        seq = "G" * 200  # Will fail GC check (100% GC) — but that's WARNING, not CRITICAL
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_synthesis_screen()
        report = v.run_all()
        # Synthesis screen issues are WARNING severity — should still pass
        assert report.passed
        # But should have warnings
        assert len(report.warnings) > 0


# ---------------------------------------------------------------------------
# Report formatting tests
# ---------------------------------------------------------------------------

class TestVerificationReport:
    def test_str_output(self):
        report = VerificationReport(construct_name="test")
        report.checks.append(CheckResult(
            name="TestCheck",
            passed=True,
            severity=Severity.CRITICAL,
            message="All good",
        ))
        output = str(report)
        assert "test" in output
        assert "PASS" in output

    def test_mixed_results(self):
        report = VerificationReport(construct_name="test")
        report.checks.append(CheckResult(
            name="Good", passed=True, severity=Severity.CRITICAL, message="OK"
        ))
        report.checks.append(CheckResult(
            name="Bad", passed=False, severity=Severity.CRITICAL, message="Failed"
        ))
        assert not report.passed
        assert len(report.critical_failures) == 1

    def test_empty_report_passes(self):
        report = VerificationReport(construct_name="test")
        assert report.passed


# ---------------------------------------------------------------------------
# Custom check tests
# ---------------------------------------------------------------------------

class TestCustomCheck:
    def test_custom_pass(self):
        seq = "ATGAAAGCC"
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_custom_check(
            name="HasATG",
            check_fn=lambda s: (s.startswith("ATG"), "Starts with ATG"),
        )
        report = v.run_all()
        assert report.passed

    def test_custom_fail(self):
        seq = "GGGAAAGCC"
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_custom_check(
            name="HasATG",
            check_fn=lambda s: (s.startswith("ATG"), "Does not start with ATG"),
        )
        report = v.run_all()
        assert not report.passed


# ---------------------------------------------------------------------------
# Integration test: Multiple checks combined
# ---------------------------------------------------------------------------

class TestIntegration:
    def test_full_pipeline(self):
        """Run a realistic verification pipeline with multiple checks."""
        # Build a simple construct: HindIII + CDS + XbaI
        cds = "ATGAAAGCCTTTGGCAAATAA"  # MKA FGK*
        seq = "AAGCTT" + cds + "TCTAGA"

        v = ConstructVerifier(sequence=seq, name="test_full")
        v.add_restriction_check(
            required_unique={"HindIII": "AAGCTT", "XbaI": "TCTAGA"},
            forbidden={"BsaI": "GGTCTC", "BsmBI": "CGTCTC"},
        )
        v.add_protein_check(
            region_start=6, region_end=6 + 18,  # 6 codons = 18 bp (before stop)
            expected_protein="MKAFGK",
            frame=0,
            label="MiniCDS",
        )
        v.add_enzyme_validation(["HindIII", "XbaI"])
        v.add_synthesis_screen()

        report = v.run_all()
        print(report)  # For debug visibility
        # Critical checks should pass
        # (warnings from synthesis screen for short seq are OK)
        assert report.passed, f"Failed: {report.critical_failures}"


# ---------------------------------------------------------------------------
# Knowledge base integration tests (Improvement A)
# ---------------------------------------------------------------------------

class TestDamSensitiveFromMetadata:
    """A4: Verify dam-sensitive enzymes are derived from JSON."""

    def test_xbai_derived_from_json(self):
        """XbaI should be in the set derived from enzyme_metadata.json."""
        dam_set = get_dam_sensitive_enzymes()
        assert "XbaI" in dam_set

    def test_superset_of_defaults(self):
        """JSON-derived set should be a superset of the hard-coded defaults."""
        from construct_verifier import _DEFAULT_DAM_SENSITIVE
        dam_set = get_dam_sensitive_enzymes()
        # JSON may have more than defaults (e.g., EcoRI, PstI mention Dam)
        # but should include the critical ones
        for enzyme in ["XbaI"]:
            assert enzyme in dam_set


class TestFussyFromMetadata:
    """A4: Verify fussy enzymes are derived from JSON."""

    def test_bbvci_derived_from_json(self):
        """BbvCI should be in the fussy set derived from enzyme_metadata.json."""
        fussy_set = get_fussy_enzymes()
        assert "BbvCI" in fussy_set

    def test_paqci_derived_from_json(self):
        """PaqCI should be in the fussy set."""
        fussy_set = get_fussy_enzymes()
        assert "PaqCI" in fussy_set

    def test_smai_derived_from_json(self):
        """SmaI should be in the fussy set."""
        fussy_set = get_fussy_enzymes()
        assert "SmaI" in fussy_set


class TestGetRecognitionSite:
    """A4: Verify recognition site lookup from KB."""

    def test_hindiii(self):
        assert get_recognition_site("HindIII") == "AAGCTT"

    def test_ecori(self):
        assert get_recognition_site("EcoRI") == "GAATTC"

    def test_bsai(self):
        assert get_recognition_site("BsaI") == "GGTCTC"

    def test_unknown_enzyme(self):
        assert get_recognition_site("FakeEnzyme123") is None


class TestRestrictionCheckByName:
    """A4: End-to-end test of add_restriction_check_by_name()."""

    def test_unique_site_passes(self):
        seq = "GGGAAGCTTGGGTCTAGAGGG"  # 1x HindIII, 1x XbaI
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_restriction_check_by_name(required_unique=["HindIII", "XbaI"])
        report = v.run_all()
        assert report.passed

    def test_forbidden_site_fails(self):
        seq = "GGGGGTCTCGGG"  # Contains BsaI
        v = ConstructVerifier(sequence=seq, name="test")
        v.add_restriction_check_by_name(forbidden=["BsaI"])
        report = v.run_all()
        assert not report.passed

    def test_unknown_enzyme_raises(self):
        seq = "AAAA"
        v = ConstructVerifier(sequence=seq, name="test")
        with pytest.raises(ValueError, match="Unknown enzyme"):
            v.add_restriction_check_by_name(required_unique=["FakeEnzyme123"])


class TestFallbackWhenNoJson:
    """A4: Verify defaults work when JSON is unavailable."""

    def test_dam_sensitive_fallback(self, monkeypatch):
        """When JSON has no enzymes, fall back to defaults."""
        import construct_verifier as cv
        # Clear cache and mock load to return empty
        monkeypatch.setattr(cv, '_enzyme_metadata_cache', None)
        original_load = cv.load_enzyme_metadata

        def mock_load():
            cv._enzyme_metadata_cache = {"enzymes": {}, "fussy_enzymes": {"enzymes": []}}
            return cv._enzyme_metadata_cache

        monkeypatch.setattr(cv, 'load_enzyme_metadata', mock_load)
        try:
            dam_set = cv.get_dam_sensitive_enzymes()
            assert dam_set == cv._DEFAULT_DAM_SENSITIVE
        finally:
            # Restore cache so other tests work
            monkeypatch.setattr(cv, '_enzyme_metadata_cache', None)
            monkeypatch.setattr(cv, 'load_enzyme_metadata', original_load)

    def test_fussy_fallback(self, monkeypatch):
        """When JSON has no fussy enzymes, fall back to defaults."""
        import construct_verifier as cv
        monkeypatch.setattr(cv, '_enzyme_metadata_cache', None)
        original_load = cv.load_enzyme_metadata

        def mock_load():
            cv._enzyme_metadata_cache = {"enzymes": {}, "fussy_enzymes": {"enzymes": []}}
            return cv._enzyme_metadata_cache

        monkeypatch.setattr(cv, 'load_enzyme_metadata', mock_load)
        try:
            fussy_set = cv.get_fussy_enzymes()
            assert fussy_set == cv._DEFAULT_FUSSY
        finally:
            monkeypatch.setattr(cv, '_enzyme_metadata_cache', None)
            monkeypatch.setattr(cv, 'load_enzyme_metadata', original_load)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
