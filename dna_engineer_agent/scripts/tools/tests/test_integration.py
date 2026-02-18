#!/usr/bin/env python3
"""
Integration tests for full construct verification (Plan item 3b).

These tests load actual GenBank plasmid files and verify end-to-end correctness:
1. GenBank roundtrip: read → verify features → verify sequence properties
2. Restriction site uniqueness in assembled constructs
3. Known construct properties (AVD002 base, AVD006 with VHH3 insertion)
4. ConstructVerifier integration with real construct data

These tests require the actual plasmid files to be present.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from construct_verifier import (
    ConstructVerifier,
    count_sites_both_strands,
    reverse_complement,
    translate,
    verify_and_gate,
)


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

PLASMID_DIR = Path(__file__).parent.parent.parent.parent / "constructs"

# Skip all tests in this file if plasmid files are not present
pytestmark = pytest.mark.skipif(
    not PLASMID_DIR.exists(),
    reason="Plasmid directory not found"
)


def _try_import_biopython():
    """Try to import Biopython; skip tests if not available."""
    try:
        from Bio import SeqIO
        return SeqIO
    except ImportError:
        pytest.skip("Biopython not installed")


# ---------------------------------------------------------------------------
# Known reference values for AVD002 (base plasmid)
# ---------------------------------------------------------------------------

AVD002_EXPECTED = {
    "length": 7104,
    "name": "AVD002",
    "cloning_sites": {
        "HindIII": {"site": "AAGCTT", "count": 1},
        "XbaI": {"site": "TCTAGA", "count": 1},
    },
}

# ---------------------------------------------------------------------------
# Known reference values for AVD006 (VHH3 at VR-IV)
# ---------------------------------------------------------------------------

AVD006_EXPECTED = {
    "length": 7551,
    "name": "AVD006",
    "insertion_size": 447,  # bp added vs AVD002
    "cloning_sites": {
        "HindIII": {"site": "AAGCTT", "count": 1},
        "XbaI": {"site": "TCTAGA", "count": 1},
    },
}


# ---------------------------------------------------------------------------
# GenBank roundtrip tests
# ---------------------------------------------------------------------------

class TestGenBankRoundtrip:
    """Verify that GenBank files are parseable and contain expected data."""

    def test_avd002_parseable(self):
        """AVD002 base plasmid can be parsed."""
        SeqIO = _try_import_biopython()
        gb_path = PLASMID_DIR / "AVD002-Rep2Mut2Cap9-6R-wt-5fold.gb"
        if not gb_path.exists():
            pytest.skip(f"{gb_path} not found")

        record = SeqIO.read(str(gb_path), "genbank")
        assert len(record.seq) == AVD002_EXPECTED["length"]
        assert len(record.features) > 0

    def test_avd006_parseable(self):
        """AVD006 VHH3-VR4 plasmid can be parsed."""
        SeqIO = _try_import_biopython()
        gb_path = PLASMID_DIR / "AVD006-Rep2Mut2Cap9-VP1-VHH3-VR4.gb"
        if not gb_path.exists():
            pytest.skip(f"{gb_path} not found")

        record = SeqIO.read(str(gb_path), "genbank")
        assert len(record.seq) == AVD006_EXPECTED["length"]

    def test_avd006_size_difference(self):
        """AVD006 should be exactly 447 bp larger than AVD002."""
        SeqIO = _try_import_biopython()
        avd002_path = PLASMID_DIR / "AVD002-Rep2Mut2Cap9-6R-wt-5fold.gb"
        avd006_path = PLASMID_DIR / "AVD006-Rep2Mut2Cap9-VP1-VHH3-VR4.gb"

        if not avd002_path.exists() or not avd006_path.exists():
            pytest.skip("Required plasmid files not found")

        avd002 = SeqIO.read(str(avd002_path), "genbank")
        avd006 = SeqIO.read(str(avd006_path), "genbank")

        size_diff = len(avd006.seq) - len(avd002.seq)
        assert size_diff == AVD006_EXPECTED["insertion_size"], \
            f"Expected {AVD006_EXPECTED['insertion_size']} bp insertion, got {size_diff}"


# ---------------------------------------------------------------------------
# Restriction site verification on real constructs
# ---------------------------------------------------------------------------

class TestRestrictionSitesReal:
    """Verify cloning site uniqueness in actual constructs (Checkpoint 9)."""

    def _load_sequence(self, filename):
        SeqIO = _try_import_biopython()
        gb_path = PLASMID_DIR / filename
        if not gb_path.exists():
            pytest.skip(f"{gb_path} not found")
        record = SeqIO.read(str(gb_path), "genbank")
        return str(record.seq).upper()

    def test_avd002_hindiii_unique(self):
        """AVD002 should have exactly 1 HindIII site."""
        seq = self._load_sequence("AVD002-Rep2Mut2Cap9-6R-wt-5fold.gb")
        assert count_sites_both_strands(seq, "AAGCTT") == 1

    def test_avd002_xbai_unique(self):
        """AVD002 should have exactly 1 XbaI site."""
        seq = self._load_sequence("AVD002-Rep2Mut2Cap9-6R-wt-5fold.gb")
        assert count_sites_both_strands(seq, "TCTAGA") == 1

    def test_avd006_hindiii_unique(self):
        """AVD006 should have exactly 1 HindIII site (insertion didn't create new one)."""
        seq = self._load_sequence("AVD006-Rep2Mut2Cap9-VP1-VHH3-VR4.gb")
        assert count_sites_both_strands(seq, "AAGCTT") == 1

    def test_avd006_xbai_unique(self):
        """AVD006 should have exactly 1 XbaI site."""
        seq = self._load_sequence("AVD006-Rep2Mut2Cap9-VP1-VHH3-VR4.gb")
        assert count_sites_both_strands(seq, "TCTAGA") == 1


# ---------------------------------------------------------------------------
# ConstructVerifier on real constructs
# ---------------------------------------------------------------------------

class TestConstructVerifierReal:
    """Run the full ConstructVerifier pipeline on actual constructs."""

    def _load_sequence(self, filename):
        SeqIO = _try_import_biopython()
        gb_path = PLASMID_DIR / filename
        if not gb_path.exists():
            pytest.skip(f"{gb_path} not found")
        record = SeqIO.read(str(gb_path), "genbank")
        return str(record.seq).upper()

    def test_avd002_full_verification(self):
        """Run full verification on AVD002 base plasmid."""
        seq = self._load_sequence("AVD002-Rep2Mut2Cap9-6R-wt-5fold.gb")

        verifier = ConstructVerifier(sequence=seq, name="AVD002")
        verifier.add_restriction_check(
            required_unique={"HindIII": "AAGCTT", "XbaI": "TCTAGA"},
        )
        verifier.add_enzyme_validation(["HindIII", "XbaI"])
        verifier.add_synthesis_screen()

        report = verifier.run_all()
        # Cloning sites should be unique — this is CRITICAL
        restriction_checks = [c for c in report.checks if "RestrictionSite" in c.name]
        assert all(c.passed for c in restriction_checks), \
            f"Restriction site check failed: {[c for c in restriction_checks if not c.passed]}"

    def test_avd006_full_verification(self):
        """Run full verification on AVD006 (VHH3 insertion at VR-IV)."""
        seq = self._load_sequence("AVD006-Rep2Mut2Cap9-VP1-VHH3-VR4.gb")

        verifier = ConstructVerifier(sequence=seq, name="AVD006")
        verifier.add_restriction_check(
            required_unique={"HindIII": "AAGCTT", "XbaI": "TCTAGA"},
        )
        verifier.add_enzyme_validation(["HindIII", "XbaI"])
        verifier.add_synthesis_screen()

        report = verifier.run_all()
        # Cloning sites must remain unique after VHH3 insertion
        restriction_checks = [c for c in report.checks if "RestrictionSite" in c.name]
        assert all(c.passed for c in restriction_checks), \
            f"Restriction site check failed: {[c for c in restriction_checks if not c.passed]}"

    def test_avd006_parent_child(self):
        """Verify AVD006 vs AVD002 parent-child relationship."""
        avd002_seq = self._load_sequence("AVD002-Rep2Mut2Cap9-6R-wt-5fold.gb")
        avd006_seq = self._load_sequence("AVD006-Rep2Mut2Cap9-VP1-VHH3-VR4.gb")

        verifier = ConstructVerifier(sequence=avd006_seq, name="AVD006")
        verifier.add_parent_child_check(
            parent_sequence=avd002_seq,
            label="AVD006-vs-AVD002",
        )

        report = verifier.run_all()
        # Should have size difference info
        size_checks = [c for c in report.checks if "size" in c.name.lower()]
        assert len(size_checks) > 0


# ---------------------------------------------------------------------------
# Feature presence tests
# ---------------------------------------------------------------------------

class TestFeaturePresence:
    """Verify that key features are annotated in GenBank files."""

    def test_avd006_has_vhh3_feature(self):
        """AVD006 should have a VHH3-related feature annotation."""
        SeqIO = _try_import_biopython()
        gb_path = PLASMID_DIR / "AVD006-Rep2Mut2Cap9-VP1-VHH3-VR4.gb"
        if not gb_path.exists():
            pytest.skip(f"{gb_path} not found")

        record = SeqIO.read(str(gb_path), "genbank")
        feature_labels = []
        for feat in record.features:
            label = feat.qualifiers.get("label", [""])[0]
            feature_labels.append(label)

        # Should have VHH3 or VHH-related feature
        vhh_features = [l for l in feature_labels if "VHH" in l.upper()]
        assert len(vhh_features) > 0, \
            f"No VHH feature found. Features: {feature_labels[:20]}"

    def test_avd006_has_rep_and_cap(self):
        """AVD006 should have Rep and Cap CDS features."""
        SeqIO = _try_import_biopython()
        gb_path = PLASMID_DIR / "AVD006-Rep2Mut2Cap9-VP1-VHH3-VR4.gb"
        if not gb_path.exists():
            pytest.skip(f"{gb_path} not found")

        record = SeqIO.read(str(gb_path), "genbank")
        cds_labels = []
        for feat in record.features:
            if feat.type == "CDS":
                label = feat.qualifiers.get("label", [""])[0]
                cds_labels.append(label)

        # Should have Rep and Cap/VP features
        has_rep = any("Rep" in l for l in cds_labels)
        has_cap = any("VP" in l or "Cap" in l for l in cds_labels)
        assert has_rep, f"No Rep CDS found. CDS labels: {cds_labels}"
        assert has_cap, f"No Cap/VP CDS found. CDS labels: {cds_labels}"


# ---------------------------------------------------------------------------
# All-constructs sweep
# ---------------------------------------------------------------------------

class TestAllConstructsSweep:
    """Quick sanity checks across all available plasmid GenBank files."""

    def test_all_plasmids_parseable(self):
        """Every .gb file in the plasmid directory should be parseable."""
        SeqIO = _try_import_biopython()
        gb_files = list(PLASMID_DIR.glob("AVD*.gb"))
        if not gb_files:
            pytest.skip("No AVD*.gb files found")

        failures = []
        for gb_path in gb_files:
            try:
                record = SeqIO.read(str(gb_path), "genbank")
                assert len(record.seq) > 0
            except Exception as e:
                failures.append(f"{gb_path.name}: {e}")

        assert not failures, f"Failed to parse: {failures}"

    def test_all_plasmids_have_hindiii(self):
        """All Rep2Mut2 AVD plasmids should have exactly 1 HindIII site (cloning site)."""
        SeqIO = _try_import_biopython()
        gb_files = list(PLASMID_DIR.glob("AVD*.gb"))
        if not gb_files:
            pytest.skip("No AVD*.gb files found")

        failures = []
        for gb_path in gb_files:
            # Skip temp files and splicing reporters (different backbone, no HindIII)
            if "temp" in gb_path.name.lower():
                continue
            if "SplicingReporter" in gb_path.name:
                continue
            try:
                record = SeqIO.read(str(gb_path), "genbank")
                seq = str(record.seq).upper()
                count = count_sites_both_strands(seq, "AAGCTT")
                if count != 1:
                    failures.append(f"{gb_path.name}: {count} HindIII sites")
            except Exception:
                pass  # Skip files that can't be parsed

        assert not failures, f"HindIII not unique in: {failures}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
