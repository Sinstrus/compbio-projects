#!/usr/bin/env python3
"""
DNA Engineer Agent - Environment Validation Script

Run this script to verify that the agent environment is correctly configured.
It checks:
1. Python dependencies (Biopython)
2. Project file structure
3. JSON schema validity
4. Basic sequence manipulation capability
"""

import sys
import json
from pathlib import Path

def check_dependencies():
    """Verify required Python packages are installed."""
    print("=" * 60)
    print("DEPENDENCY CHECK")
    print("=" * 60)
    
    required = {
        'Bio': 'biopython',
        'Bio.Seq': 'biopython',
        'Bio.SeqIO': 'biopython',
        'Bio.SeqFeature': 'biopython',
    }
    
    missing = []
    for module, package in required.items():
        try:
            __import__(module)
            print(f"  ✓ {module}")
        except ImportError:
            print(f"  ✗ {module} (install: pip install {package})")
            missing.append(package)
    
    if missing:
        print(f"\n  MISSING PACKAGES: pip install {' '.join(set(missing))}")
        return False
    
    print("\n  All dependencies satisfied.\n")
    return True


def check_file_structure(base_path: Path):
    """Verify the expected project structure exists."""
    print("=" * 60)
    print("FILE STRUCTURE CHECK")
    print("=" * 60)
    
    expected_files = [
        "AGENT_INSTRUCTIONS_v2.md",
        "VERSION.json",
        "knowledge_base/SCHEMA.json",
        "knowledge_base/systems/recombinant_aav_production.json",
        "knowledge_base/references/aav_references.json",
        "cis_elements/manifest.json",
        "parts_library/schema.json",
        "parts_library/linkers.json",
    ]
    
    expected_dirs = [
        "knowledge_base",
        "knowledge_base/systems",
        "knowledge_base/references",
        "cis_elements",
        "parts_library", 
        "projects",
        "reports",
        "scripts",
        "test_data",
    ]
    
    all_ok = True
    
    for dir_name in expected_dirs:
        dir_path = base_path / dir_name
        if dir_path.is_dir():
            print(f"  ✓ {dir_name}/")
        else:
            print(f"  ✗ {dir_name}/ (missing directory)")
            all_ok = False
    
    for file_name in expected_files:
        file_path = base_path / file_name
        if file_path.is_file():
            print(f"  ✓ {file_name}")
        else:
            print(f"  ✗ {file_name} (missing file)")
            all_ok = False
    
    if all_ok:
        print("\n  File structure OK.\n")
    else:
        print("\n  Some files/directories missing.\n")
    
    return all_ok


def validate_json_files(base_path: Path):
    """Check that JSON files are valid and parseable."""
    print("=" * 60)
    print("JSON VALIDATION")
    print("=" * 60)
    
    json_files = [
        "knowledge_base/SCHEMA.json",
        "knowledge_base/systems/recombinant_aav_production.json",
        "knowledge_base/references/aav_references.json",
        "cis_elements/manifest.json",
        "parts_library/schema.json",
        "parts_library/linkers.json",
    ]
    
    all_ok = True
    for file_name in json_files:
        file_path = base_path / file_name
        if not file_path.exists():
            print(f"  ✗ {file_name} (file not found)")
            all_ok = False
            continue
            
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            print(f"  ✓ {file_name} (valid JSON)")
        except json.JSONDecodeError as e:
            print(f"  ✗ {file_name} (invalid JSON: {e})")
            all_ok = False
    
    if all_ok:
        print("\n  All JSON files valid.\n")
    
    return all_ok


def test_biopython_operations():
    """Run a simple Biopython test to verify functionality."""
    print("=" * 60)
    print("BIOPYTHON FUNCTIONALITY TEST")
    print("=" * 60)
    
    try:
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio import SeqIO
        import io
        
        # Create a test sequence
        test_seq = Seq("ATGAAACCCGGGTTTTAA")
        record = SeqRecord(
            test_seq,
            id="TEST001",
            name="TestPlasmid",
            description="Test sequence for environment validation",
            annotations={"molecule_type": "DNA"}
        )
        
        # Add an annotation
        feature = SeqFeature(
            FeatureLocation(0, 18, strand=1),
            type="CDS",
            qualifiers={
                "label": ["TestGene"],
                "note": ["[Action]: Created | [Reason]: Validation | [Risk]: None"]
            }
        )
        record.features.append(feature)
        
        # Test translation
        protein = test_seq.translate()
        expected_protein = "MKPGF*"
        assert str(protein) == expected_protein, f"Translation failed: got {protein}"
        print(f"  ✓ Translation: {test_seq} → {protein}")
        
        # Test GenBank writing
        output = io.StringIO()
        SeqIO.write(record, output, "genbank")
        gb_content = output.getvalue()
        assert "TestGene" in gb_content, "Feature annotation not in output"
        print(f"  ✓ GenBank output: {len(gb_content)} bytes")
        
        # Test reverse complement
        rc = test_seq.reverse_complement()
        print(f"  ✓ Reverse complement: {test_seq} → {rc}")
        
        print("\n  Biopython operations verified.\n")
        return True
        
    except Exception as e:
        print(f"  ✗ Biopython test failed: {e}\n")
        return False


def list_available_parts(base_path: Path):
    """Display available parts libraries."""
    print("=" * 60)
    print("AVAILABLE PARTS LIBRARIES")
    print("=" * 60)
    
    parts_dir = base_path / "parts_library"
    if not parts_dir.exists():
        print("  Parts library directory not found.")
        return
    
    for json_file in parts_dir.glob("*.json"):
        if json_file.name == "schema.json":
            continue
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            lib_name = data.get("library_name", json_file.stem)
            parts = data.get("parts", [])
            print(f"\n  [{json_file.name}] {lib_name}")
            print(f"  {'─' * 40}")
            for part in parts:
                part_id = part.get("id", "?")
                part_name = part.get("name", "?")
                seq_type = part.get("sequence_type", "?")
                print(f"    • {part_id}: {part_name} ({seq_type})")
        except Exception as e:
            print(f"  Error reading {json_file.name}: {e}")
    
    print()


def list_cis_elements(base_path: Path):
    """Display loaded cis-elements."""
    print("=" * 60)
    print("LOADED CIS-ELEMENTS (Fallback)")
    print("=" * 60)
    
    manifest_path = base_path / "cis_elements" / "manifest.json"
    if not manifest_path.exists():
        print("  Cis-element manifest not found.")
        return
    
    try:
        with open(manifest_path, 'r') as f:
            data = json.load(f)
        
        elements = data.get("elements", [])
        print(f"\n  {len(elements)} elements loaded (for fallback use)\n")
            
    except Exception as e:
        print(f"  Error reading manifest: {e}")
    
    print()


def list_biological_systems(base_path: Path):
    """Display available biological systems from knowledge base."""
    print("=" * 60)
    print("AVAILABLE BIOLOGICAL SYSTEMS")
    print("=" * 60)
    
    systems_dir = base_path / "knowledge_base" / "systems"
    if not systems_dir.exists():
        print("  Systems directory not found.")
        return
    
    for json_file in systems_dir.glob("*.json"):
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            system_id = data.get("system_id", json_file.stem)
            name = data.get("name", "Unknown")
            goal = data.get("goal", "No goal defined")
            
            # Count requirements
            proteins = data.get("proteins", {})
            cis_elements = data.get("cis_elements", {})
            construct_types = data.get("construct_types", {})
            
            print(f"\n  [{system_id}]")
            print(f"  Name: {name}")
            print(f"  Goal: {goal}")
            print(f"  Proteins defined: {len(proteins)}")
            print(f"  Cis-elements defined: {len(cis_elements)}")
            print(f"  Construct types: {len(construct_types)}")
            
        except Exception as e:
            print(f"  Error reading {json_file.name}: {e}")
    
    print()


def main():
    """Run all validation checks."""
    print("\n" + "=" * 60)
    print("  DNA ENGINEER AGENT - ENVIRONMENT VALIDATION")
    print("=" * 60 + "\n")
    
    # Determine base path (script is in scripts/, so go up one level)
    script_path = Path(__file__).parent.parent  # Go from scripts/ to root
    if (script_path / "AGENT_INSTRUCTIONS.md").exists():
        base_path = script_path
    else:
        base_path = Path.cwd()
    
    print(f"  Base path: {base_path}\n")
    
    # Show agent version
    version_file = base_path / "VERSION.json"
    if version_file.exists():
        try:
            with open(version_file, 'r') as f:
                version_data = json.load(f)
            print(f"  Agent Version: {version_data.get('agent_version', 'unknown')}")
            print(f"  Release Date:  {version_data.get('release_date', 'unknown')}\n")
        except:
            pass
    
    # Run checks
    results = {
        "dependencies": check_dependencies(),
        "file_structure": check_file_structure(base_path),
        "json_valid": validate_json_files(base_path),
        "biopython_ops": test_biopython_operations(),
    }
    
    # Show available resources
    list_biological_systems(base_path)
    list_available_parts(base_path)
    list_cis_elements(base_path)
    
    # Summary
    print("=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)
    
    all_passed = all(results.values())
    for check, passed in results.items():
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {status}: {check}")
    
    print()
    if all_passed:
        print("  🟢 AGENT READY: All checks passed.")
        print("     You can now provide a .gb file and design request.")
    else:
        print("  🔴 AGENT NOT READY: Some checks failed.")
        print("     Please resolve issues above before proceeding.")
    
    print("\n" + "=" * 60 + "\n")
    
    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
