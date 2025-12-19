"""
Pytest configuration and fixtures for DNA Engineer Agent tests.
"""

import pytest
from pathlib import Path


@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return Path(__file__).parent / "test_data"


@pytest.fixture
def sample_cds():
    """Sample CDS sequence for testing (ATG...TAA)."""
    return "ATGAAAGCCTAA"  # ATG-AAA-GCC-TAA (Met-Lys-Ala-Stop)


@pytest.fixture
def sample_cds_with_context():
    """Sample CDS with upstream/downstream context."""
    return {
        'upstream': 'GGATCC',  # BamHI site
        'cds': 'ATGAAAGCCTAA',
        'downstream': 'GAATTC',  # EcoRI site
        'full_seq': 'GGATCCATGAAAGCCTAAGAATTC'
    }


@pytest.fixture
def codon_table():
    """Standard genetic code codon table."""
    return {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }


@pytest.fixture
def enzyme_sites():
    """Common restriction enzyme recognition sequences."""
    return {
        'HindIII': 'AAGCTT',
        'XbaI': 'TCTAGA',
        'EcoRI': 'GAATTC',
        'BamHI': 'GGATCC',
        'BsaI': 'GGTCTC',
        'BsmBI': 'CGTCTC',
        'BbsI': 'GAAGAC',
    }


@pytest.fixture
def backbone_catalog_path():
    """Path to the backbone catalog JSON file."""
    return Path(__file__).parent.parent.parent.parent / 'backbones' / 'genscript' / 'BACKBONE_CATALOG.json'
