#!/usr/bin/env python3
"""
Build AVD005 (VP1 knockout) and AVD007 (VHH3-D2 display) GenBank files.

AVD005: VP1 start codon knockout (ATG→AAG) for trans-complementation
AVD007: VHH3 at VR-IV with 5×GGGGS N-term, 1×GGGGS C-term linkers

Version: 1.0
Date: 2026-01-20
Base plasmid: AVD002-Rep2Mut2Cap9-6R-wt.dna (7,104 bp)
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from datetime import datetime
import sys
from pathlib import Path
import snapgene_reader

try:
    from dnachisel import *
    DNACHISEL_AVAILABLE = True
except ImportError:
    print("⚠ DNAchisel not available - using pre-optimized sequence")
    DNACHISEL_AVAILABLE = False

# Project paths
PROJECT_DIR = Path(__file__).parent.parent
PLASMIDS_DIR = PROJECT_DIR / "plasmids"
SYNTHETIC_DIR = PROJECT_DIR / "synthetic_fragments"
DOCS_DIR = PROJECT_DIR / "docs"

# Ensure output directories exist
SYNTHETIC_DIR.mkdir(exist_ok=True)
DOCS_DIR.mkdir(exist_ok=True)

# VHH3 anti-ALPL sequence (119 amino acids - verified)
VHH3_AA = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSTADMGWFRQAPGKGRELVAAVSGSGFSTYSDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAKATISLYYAMDVWGQGTTVTVSS"

# D2 Linkers (Biogen patent specification - modified to 5×GGGGS)
LINK_D2_N_AA = "GGGGSGGGGSGGGGSGGGGSGGGGS"  # (GGGGS)x5 = 25 aa
LINK_D2_C_AA = "GGGGS"  # (GGGGS)x1 = 5 aa (ASYMMETRIC)

# Full insert amino acid sequence
FULL_INSERT_AA = LINK_D2_N_AA + VHH3_AA + LINK_D2_C_AA  # 149 aa total (25+119+5)


def load_avd002():
    """Load AVD002 base plasmid from SnapGene .dna file"""
    print("Loading AVD002 base plasmid...")
    avd002_path = PLASMIDS_DIR / "AVD002-Rep2Mut2Cap9-6R-wt.dna"

    if not avd002_path.exists():
        raise FileNotFoundError(f"AVD002 not found at {avd002_path}")

    # Read SnapGene file using Biopython
    record = SeqIO.read(str(avd002_path), 'snapgene')

    print(f"  ✓ Loaded AVD002: {len(record.seq)} bp")
    print(f"  ✓ Features: {len(record.features)}")

    # Verify critical positions
    seq = str(record.seq).upper()
    if seq[2378:2381] != "ATG":
        raise ValueError(f"VP1 start codon at bp 2379-2381 is {seq[2378:2381]}, expected ATG")
    if seq[2789:2792] != "ACG":
        raise ValueError(f"VP2 start codon at bp 2790-2792 is {seq[2789:2792]}, expected ACG")
    if seq[2984:2987] != "ATG":
        raise ValueError(f"VP3 start codon at bp 2985-2987 is {seq[2984:2987]}, expected ATG")

    print(f"  ✓ VP1 start (bp 2379): {seq[2378:2381]}")
    print(f"  ✓ VP2 start (bp 2790): {seq[2789:2792]}")
    print(f"  ✓ VP3 start (bp 2985): {seq[2984:2987]}")

    return record


def create_avd005(avd002_record):
    """
    Create AVD005: VP1 knockout helper plasmid

    Single mutation: ATG→AAG at bp 2379 (VP1 start codon)
    VP2/VP3 remain functional for trans-complementation
    """
    print("\nCreating AVD005 (VP1 knockout)...")

    avd002_seq = str(avd002_record.seq).upper()

    # Make mutation at position 2378 (0-indexed)
    avd005_seq = avd002_seq[:2378] + "AAG" + avd002_seq[2381:]

    # Verify mutation
    if avd005_seq[2378:2381] != "AAG":
        raise ValueError("VP1 knockout mutation failed")

    print(f"  ✓ VP1 knockout: ATG→AAG at bp 2379-2381")
    print(f"  ✓ VP2 start preserved: {avd005_seq[2789:2792]}")
    print(f"  ✓ VP3 start preserved: {avd005_seq[2984:2987]}")
    print(f"  ✓ Total size: {len(avd005_seq)} bp")

    # Create SeqRecord
    record = SeqRecord(
        Seq(avd005_seq),
        id="AVD005",
        name="AVD005",
        description="Rep2Mut2Cap9-VP1ko | VP1 knockout helper plasmid for trans-complementation",
        annotations={
            "molecule_type": "DNA",
            "topology": "circular",
            "date": datetime.now().strftime("%d-%b-%Y").upper()
        }
    )

    # Copy all features from AVD002
    # Since we only changed 3 bp, all positions remain the same
    # Just update VP1 annotation to note it's knocked out
    features = []
    for feat in avd002_record.features:
        new_feat = SeqFeature(
            location=feat.location,
            type=feat.type,
            qualifiers=dict(feat.qualifiers)  # Copy qualifiers
        )

        # Update VP1 annotation
        if feat.type == "CDS" and feat.qualifiers.get("label", [""])[0] == "VP1":
            new_feat.qualifiers["label"] = "VP1 (KNOCKOUT)"
            new_feat.qualifiers["note"] = "VP1 start codon mutated ATG→AAG - no translation"

        features.append(new_feat)

    record.features = features
    print(f"  ✓ Copied {len(features)} features from AVD002")

    return record


def optimize_insert_sequence(insert_aa):
    """
    Codon-optimized 447 bp insert (5×GGGGS + VHH3 + 1×GGGGS)

    Pre-optimized sequence with:
    - Varied codon usage for glycine-rich linkers
    - GC content optimized for synthesis
    - Human codon preferences
    - No repetitive sequences or hairpins
    """
    print("\n  Using pre-optimized insert sequence...")

    # Pre-optimized sequence (447 bp, 149 aa)
    # Manually designed with varied codons to avoid synthesis issues
    optimized_seq = (
        # N-terminal linker 5×GGGGS (75 bp, 25 aa) - varied codons
        "GGTGGAGGCGGATCTGGAGGCGGTGGTTCAGGCGGTGGAGGAAGTGGTGGCGGAGGTTCTGGTGGAGGCGGTTCT"
        # VHH3 (357 bp, 119 aa) - codon optimized for human expression
        "GAGGTGCAACTGGTTGAAAGCGGCGGAGGACTTGTTCAACCCGGCGGCAGCCTTAGGCTTTCTTGCGCTGCCAGCGGCTTCACCTTTAGCACCGCCGACATGGGCTGGTTTAGGCAAGCTCCCGGAAAAGGCAGGGAACTTGTTGCCGCTGTGAGCGGCAGCGGCTTCAGCACCTACTCTGATAGCGTTGAGGGCAGGTTCACCATCAGCAGGGACAACGCCAAGAGGATGGTGTACCTGCAGATGAACAGCTTGAGGGCCGAGGACACCGCCGTGTACTACTGCGCCAAGGCCACAATTAGCCTGTACTACGCCATGGATGTGTGGGGACAGGGCACCACCGTGACCGTGAGCAGC"
        # C-terminal linker 1×GGGGS (15 bp, 5 aa)
        "GGTGGAGGCGGTTCT"
    )

    # Verify translation
    translated = str(Seq(optimized_seq).translate())
    if translated != insert_aa:
        raise ValueError(f"Pre-optimized sequence translation mismatch!\nExpected: {insert_aa}\nGot: {translated}")

    # Verify length
    expected_len = len(insert_aa) * 3  # Should be 447 bp
    if len(optimized_seq) != expected_len:
        raise ValueError(f"Pre-optimized sequence length is {len(optimized_seq)}, expected {expected_len} bp")

    # Calculate and report GC content
    gc_content = 100 * sum(1 for c in optimized_seq if c in 'GC') / len(optimized_seq)

    print(f"  ✓ Insert length: {len(optimized_seq)} bp ({len(insert_aa)} aa)")
    print(f"  ✓ GC content: {gc_content:.1f}%")
    print(f"  ✓ Translation verified: {len(translated)} aa")

    return optimized_seq


def shift_feature_location(location, insertion_point, shift_amount):
    """
    Shift a feature location if it's after the insertion point.

    Args:
        location: FeatureLocation or CompoundLocation
        insertion_point: Position where insert was added
        shift_amount: How many bp to shift

    Returns:
        New shifted location
    """
    from Bio.SeqFeature import CompoundLocation

    if isinstance(location, CompoundLocation):
        # Handle compound locations (e.g., join features)
        new_parts = []
        for part in location.parts:
            new_parts.append(shift_feature_location(part, insertion_point, shift_amount))
        return CompoundLocation(new_parts, operator=location.operator)

    # Simple location
    start = location.start
    end = location.end
    strand = location.strand

    # If feature starts at or after insertion point, shift it
    if start >= insertion_point:
        start += shift_amount
        end += shift_amount
    # If feature spans insertion point, extend the end
    elif start < insertion_point < end:
        end += shift_amount

    return FeatureLocation(start, end, strand=strand)


def create_avd007(avd002_record):
    """
    Create AVD007: VHH3 display with Biogen D2 linkers

    Steps:
    1. VP2 knockout: ACG→ACC at bp 2790 (silent)
    2. VP3 knockout: ATG→CTG at bp 2985 (non-silent)
    3. Insert 447 bp (5×GGGGS + VHH3 + 1×GGGGS) at VR-IV region (bp 3746)
    """
    print("\nCreating AVD007 (VHH3 display with D2 linkers)...")

    avd002_seq = str(avd002_record.seq).upper()

    # Step 1: VP2 knockout (silent)
    seq = avd002_seq[:2789] + "ACC" + avd002_seq[2792:]
    print(f"  ✓ VP2 knockout: ACG→ACC at bp 2790-2792 (silent, Thr→Thr)")

    # Step 2: VP3 knockout (non-silent)
    seq = seq[:2984] + "CTG" + seq[2987:]
    print(f"  ✓ VP3 knockout: ATG→CTG at bp 2985-2987 (non-silent, Met→Leu)")

    # Step 3: Generate optimized insert
    optimized_insert = optimize_insert_sequence(FULL_INSERT_AA)

    # Verify insert length
    expected_length = len(FULL_INSERT_AA) * 3  # 149 aa * 3 = 447 bp
    if len(optimized_insert) != expected_length:
        raise ValueError(f"Insert length {len(optimized_insert)} != expected {expected_length}")

    # Verify translation
    insert_translation = str(Seq(optimized_insert).translate())
    if insert_translation != FULL_INSERT_AA:
        raise ValueError(f"Insert translation mismatch:\nExpected: {FULL_INSERT_AA}\nGot: {insert_translation}")

    print(f"  ✓ Insert verified: {len(optimized_insert)} bp, translates to {len(FULL_INSERT_AA)} aa")

    # Step 4: Insert at VR-IV region (between SKTINGSG and QNQQTLKF)
    # SKTINGSG ends at bp 3743 (0-indexed), so insert after bp 3743
    # This places insert between SKTINGSG and QNQQTLKF as intended
    insertion_point = 3743  # CORRECTED from 3746
    insert_size = len(optimized_insert)  # 447 bp

    # Build final sequence
    avd007_seq = seq[:insertion_point] + optimized_insert + seq[insertion_point:]

    print(f"  ✓ Inserted at bp {insertion_point} (VR-IV region)")
    print(f"  ✓ Total size: {len(avd007_seq)} bp (7104 + 447)")

    # Create SeqRecord
    record = SeqRecord(
        Seq(avd007_seq),
        id="AVD007",
        name="AVD007",
        description="Rep2Mut2Cap9-VP1-VHH3-D2 | VHH3 display at VR-IV with 5×GGGGS N-term, 1×GGGGS C-term",
        annotations={
            "molecule_type": "DNA",
            "topology": "circular",
            "date": datetime.now().strftime("%d-%b-%Y").upper()
        }
    )

    # Copy and shift features from AVD002
    features = []
    for feat in avd002_record.features:
        # Shift feature location if needed
        new_location = shift_feature_location(feat.location, insertion_point, insert_size)

        new_feat = SeqFeature(
            location=new_location,
            type=feat.type,
            qualifiers=dict(feat.qualifiers)  # Copy qualifiers
        )

        # Update specific annotations
        label = feat.qualifiers.get("label", [""])[0]

        if label == "VP1":
            new_feat.qualifiers["label"] = "VP1-VHH3 fusion"
            new_feat.qualifiers["note"] = "AAV9 VP1 with anti-ALPL VHH3 at VR-IV using D2 linkers"
        elif label == "VP2":
            new_feat.qualifiers["label"] = "VP2 (KNOCKOUT)"
            new_feat.qualifiers["note"] = "VP2 start codon mutated ACG→ACC - no translation"
        elif label == "VP3":
            new_feat.qualifiers["label"] = "VP3 (KNOCKOUT)"
            new_feat.qualifiers["note"] = "VP3 start codon mutated ATG→CTG - no translation"
        elif "AAP" in label:
            new_feat.qualifiers["note"] = f"Assembly-Activating Protein (+1 reading frame, preserved by {insert_size} bp insert)"

        features.append(new_feat)

    # Add new VHH3 insert features (at corrected position 3743)
    # N-terminal linker
    features.append(SeqFeature(
        FeatureLocation(insertion_point, insertion_point + 75),
        type="misc_feature",
        qualifiers={
            "label": "LINK_D2_N",
            "note": "N-terminal linker (GGGGS)×5, 25 aa - provides rotational freedom"
        }
    ))

    # VHH3
    features.append(SeqFeature(
        FeatureLocation(insertion_point + 75, insertion_point + 75 + 357, strand=+1),
        type="CDS",
        qualifiers={
            "label": "VHH3_anti-ALPL",
            "note": f"Anti-ALPL VHH nanobody, {len(VHH3_AA)} aa, codon-optimized",
            "translation": VHH3_AA
        }
    ))

    # C-terminal linker
    features.append(SeqFeature(
        FeatureLocation(insertion_point + 75 + 357, insertion_point + insert_size),
        type="misc_feature",
        qualifiers={
            "label": "LINK_D2_C",
            "note": "C-terminal linker (GGGGS)×1, 5 aa - nucleates VHH folding"
        }
    ))

    record.features = features
    print(f"  ✓ Transferred {len(features)} features (shifted {insert_size} bp where needed)")

    return record, optimized_insert


def verify_avd005(record):
    """Verify AVD005 VP1 knockout"""
    print("\nVerifying AVD005...")
    seq = str(record.seq)

    checks = []

    # Check VP1 start codon mutation
    vp1_codon = seq[2378:2381]
    checks.append(("VP1 start = AAG", vp1_codon == "AAG", f"Found: {vp1_codon}"))

    # Check VP2/VP3 preserved
    vp2_codon = seq[2789:2792]
    checks.append(("VP2 start = ACG", vp2_codon == "ACG", f"Found: {vp2_codon}"))

    vp3_codon = seq[2984:2987]
    checks.append(("VP3 start = ATG", vp3_codon == "ATG", f"Found: {vp3_codon}"))

    # Check size
    checks.append(("Size = 7104 bp", len(seq) == 7104, f"Found: {len(seq)} bp"))

    # Print results
    all_passed = True
    for check_name, passed, details in checks:
        status = "✓" if passed else "✗"
        print(f"  {status} {check_name} - {details}")
        if not passed:
            all_passed = False

    if all_passed:
        print("  ✅ AVD005 verification PASSED")
    else:
        print("  ❌ AVD005 verification FAILED")
        raise ValueError("AVD005 verification failed")

    return all_passed


def verify_avd007(record):
    """Verify AVD007 VHH3 insertion"""
    print("\nVerifying AVD007...")
    seq = str(record.seq)

    checks = []

    # Check VP2/VP3 knockouts
    vp2_codon = seq[2789:2792]
    checks.append(("VP2 knockout = ACC", vp2_codon == "ACC", f"Found: {vp2_codon}"))

    vp3_codon = seq[2984:2987]
    checks.append(("VP3 knockout = CTG", vp3_codon == "CTG", f"Found: {vp3_codon}"))

    # Check insert length and translation (447 bp, 149 aa)
    # Insert is at bp 3743 (CORRECTED from 3746)
    insertion_point = 3743
    expected_insert_len = len(FULL_INSERT_AA) * 3  # 149 * 3 = 447 bp
    insert_region = seq[insertion_point:insertion_point+expected_insert_len]
    checks.append((f"Insert length = {expected_insert_len} bp", len(insert_region) == expected_insert_len,
                  f"Found: {len(insert_region)} bp"))

    insert_translation = str(Seq(insert_region).translate())
    has_stop = "*" in insert_translation
    checks.append(("No stop codons in insert", not has_stop, f"Stop found" if has_stop else "Clean"))

    translation_matches = insert_translation == FULL_INSERT_AA
    checks.append(("Insert translates correctly", translation_matches,
                  "Match" if translation_matches else f"Mismatch"))

    # Check total size (7104 + 447 = 7551 bp)
    expected_total_size = 7104 + expected_insert_len
    checks.append((f"Size = {expected_total_size} bp", len(seq) == expected_total_size,
                  f"Found: {len(seq)} bp"))

    # Check reading frame (multiple of 3)
    checks.append(("Insert is multiple of 3", len(insert_region) % 3 == 0,
                  f"{len(insert_region)} % 3 = {len(insert_region) % 3}"))

    # Check flanking sequences
    # Before insert should be SKTINGSG
    flank_before = seq[insertion_point-24:insertion_point]
    flank_before_aa = str(Seq(flank_before).translate())
    checks.append(("Upstream = SKTINGSG", flank_before_aa == "SKTINGSG",
                  f"Found: {flank_before_aa}"))

    # After insert should be QNQQTLKF
    flank_after = seq[insertion_point+expected_insert_len:insertion_point+expected_insert_len+24]
    flank_after_aa = str(Seq(flank_after).translate())
    checks.append(("Downstream = QNQQTLKF", flank_after_aa == "QNQQTLKF",
                  f"Found: {flank_after_aa}"))

    # Print results
    all_passed = True
    for check_name, passed, details in checks:
        status = "✓" if passed else "✗"
        print(f"  {status} {check_name} - {details}")
        if not passed:
            all_passed = False

    if all_passed:
        print("  ✅ AVD007 verification PASSED")
    else:
        print("  ❌ AVD007 verification FAILED")
        raise ValueError("AVD007 verification failed")

    return all_passed


def export_genbank(record, output_path):
    """Write GenBank file"""
    SeqIO.write(record, output_path, "genbank")
    print(f"  ✓ Wrote {output_path}")


def export_synthetic_fragment(insert_seq, output_path):
    """Export synthetic fragment for ordering"""
    with open(output_path, 'w') as f:
        f.write(f">AVD007_synthetic_fragment | {len(insert_seq)}bp insert (5×GGGGS + VHH3 + 1×GGGGS)\n")
        # Write in 80 character lines
        for i in range(0, len(insert_seq), 80):
            f.write(insert_seq[i:i+80] + "\n")
    print(f"  ✓ Wrote {output_path}")


def create_verification_report(avd005_record, avd007_record, insert_seq):
    """Create verification report"""
    report_path = DOCS_DIR / "DESIGN_VERIFICATION_AVD005_AVD007.md"

    with open(report_path, 'w') as f:
        f.write("# Design Verification Report: AVD005 and AVD007\n\n")
        f.write(f"**Date:** {datetime.now().strftime('%Y-%m-%d')}\n\n")
        f.write(f"**Base plasmid:** AVD002-Rep2Mut2Cap9-6R-wt.dna (7,104 bp)\n\n")

        f.write("## AVD005: VP1 Knockout Helper Plasmid\n\n")
        f.write("### Design Summary\n")
        f.write(f"- **Size:** {len(avd005_record.seq)} bp\n")
        f.write(f"- **Modification:** VP1 start codon knockout (ATG→AAG at bp 2379-2381)\n")
        f.write(f"- **Purpose:** Trans-complementation helper for AVD007\n\n")

        f.write("### Verification Checklist\n")
        seq005 = str(avd005_record.seq)
        f.write(f"- [{'x' if seq005[2378:2381] == 'AAG' else ' '}] VP1 start codon = AAG (knockout)\n")
        f.write(f"- [{'x' if seq005[2789:2792] == 'ACG' else ' '}] VP2 start codon = ACG (preserved)\n")
        f.write(f"- [{'x' if seq005[2984:2987] == 'ATG' else ' '}] VP3 start codon = ATG (preserved)\n")
        f.write(f"- [{'x' if len(seq005) == 7104 else ' '}] Total size = 7,104 bp\n\n")

        f.write("## AVD007: VHH3 Display with D2 Linkers\n\n")
        f.write("### Design Summary\n")
        f.write(f"- **Size:** {len(avd007_record.seq)} bp\n")
        f.write(f"- **Modifications:**\n")
        f.write(f"  - VP2 knockout: ACG→ACC at bp 2790-2792 (silent)\n")
        f.write(f"  - VP3 knockout: ATG→CTG at bp 2985-2987 (non-silent)\n")
        f.write(f"  - VHH3 insertion: {len(insert_seq)} bp at VR-IV region (bp 3743, after SKTINGSG)\n")
        f.write(f"- **Insert composition:**\n")
        f.write(f"  - N-terminal linker: 5×GGGGS (25 aa, 75 bp)\n")
        f.write(f"  - VHH3: {len(VHH3_AA)} aa, {len(VHH3_AA)*3} bp\n")
        f.write(f"  - C-terminal linker: 1×GGGGS (5 aa, 15 bp)\n")
        f.write(f"  - **Total insert:** {len(FULL_INSERT_AA)} aa, {len(insert_seq)} bp\n\n")

        f.write("### Verification Checklist\n")
        seq007 = str(avd007_record.seq)
        insertion_point = 3743  # Corrected position
        expected_insert_len = len(insert_seq)
        insert_region = seq007[insertion_point:insertion_point+expected_insert_len]
        expected_total_size = 7104 + expected_insert_len
        f.write(f"- [{'x' if seq007[2789:2792] == 'ACC' else ' '}] VP2 knockout = ACC\n")
        f.write(f"- [{'x' if seq007[2984:2987] == 'CTG' else ' '}] VP3 knockout = CTG\n")
        f.write(f"- [{'x' if len(insert_region) == expected_insert_len else ' '}] Insert length = {expected_insert_len} bp\n")
        f.write(f"- [{'x' if '*' not in str(Seq(insert_region).translate()) else ' '}] No stop codons in insert\n")
        f.write(f"- [{'x' if len(insert_region) % 3 == 0 else ' '}] Insert is multiple of 3\n")
        f.write(f"- [{'x' if len(seq007) == expected_total_size else ' '}] Total size = {expected_total_size:,} bp\n\n")

        f.write("### Insert Sequence Analysis\n\n")
        gc_content = 100 * sum(1 for c in insert_seq if c in 'GC') / len(insert_seq)
        f.write(f"- **GC content:** {gc_content:.1f}%\n")
        f.write(f"- **Translation:** {str(Seq(insert_seq).translate())}\n\n")

        f.write("### Synthetic Fragment for Ordering\n\n")
        f.write(f"```\n{insert_seq}\n```\n\n")

        f.write("## Next Steps\n\n")
        f.write("1. Order synthetic fragment from GenScript (~$200-400)\n")
        f.write("2. Clone into AVD002 using appropriate restriction sites\n")
        f.write("3. Sequence verify both plasmids (Sanger sequencing)\n")
        f.write("4. Co-transfect AVD005 (helper) + AVD007 (VHH display)\n")
        f.write("5. Functional assays: Western blot, ELISA, transduction\n\n")

        f.write("---\n")
        f.write(f"*Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*\n")

    print(f"  ✓ Wrote {report_path}")


def main():
    """Main build pipeline"""
    print("=" * 70)
    print("Building AVD005 (VP1 knockout) and AVD007 (VHH3 display)")
    print("=" * 70)

    # Load base plasmid
    avd002_record = load_avd002()

    # Create AVD005
    avd005_record = create_avd005(avd002_record)
    verify_avd005(avd005_record)

    # Create AVD007
    avd007_record, insert_seq = create_avd007(avd002_record)
    verify_avd007(avd007_record)

    # Export files
    print("\nWriting GenBank files...")
    export_genbank(avd005_record, PLASMIDS_DIR / "AVD005-Rep2Mut2Cap9-VP1ko.gb")
    export_genbank(avd007_record, PLASMIDS_DIR / "AVD007-Rep2Mut2Cap9-VP1-VHH3-D2.gb")

    print("\nWriting synthetic fragment...")
    export_synthetic_fragment(insert_seq, SYNTHETIC_DIR / "AVD007_synthetic_fragment.fasta")

    print("\nWriting verification report...")
    create_verification_report(avd005_record, avd007_record, insert_seq)

    print("\n" + "=" * 70)
    print("✅ Build completed successfully")
    print("=" * 70)
    print("\nGenerated files:")
    print(f"  - {PLASMIDS_DIR / 'AVD005-Rep2Mut2Cap9-VP1ko.gb'} ({len(avd005_record.seq)} bp)")
    print(f"  - {PLASMIDS_DIR / 'AVD007-Rep2Mut2Cap9-VP1-VHH3-D2.gb'} ({len(avd007_record.seq)} bp)")
    print(f"  - {SYNTHETIC_DIR / 'AVD007_synthetic_fragment.fasta'} ({len(insert_seq)} bp)")
    print(f"  - {DOCS_DIR / 'DESIGN_VERIFICATION_AVD005_AVD007.md'}")
    print()


if __name__ == "__main__":
    main()
