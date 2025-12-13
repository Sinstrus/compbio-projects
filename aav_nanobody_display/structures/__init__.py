"""
Structure handling module for AAV nanobody display pipeline.

Submodules:
- fetch: Download structures from RCSB and AlphaFold DB
- parse: Parse and analyze capsid structures
- assemble: Build 60-mer capsids from monomers
"""

from .fetch import StructureFetcher, fetch_capsid
from .parse import CapsidParser, CapsidMonomer, VRLoop, parse_capsid
from .assemble import assemble_capsid_from_monomer, extract_asymmetric_unit

__all__ = [
    "StructureFetcher",
    "fetch_capsid",
    "CapsidParser", 
    "CapsidMonomer",
    "VRLoop",
    "parse_capsid",
    "assemble_capsid_from_monomer",
    "extract_asymmetric_unit",
]
