#!/usr/bin/env python3
"""
Build AVD005-AVD011 AAV plasmid constructs with VHH3 display variants.

Constructs:
- AVD005: VP1 knockout only
- AVD006: VHH3 in VR4, VP2/3 native
- AVD007: VHH3 in VR4, VP2/3 knockouts
- AVD008: VHH3 at VP1 N-term, VP2/3 native
- AVD009: VHH3 at VP1 N-term, VP2/3 knockouts
- AVD010: VHH3 at VP2 N-term, VP1/3 native
- AVD011: VP2 knockout only

Version: 1.0
Date: 2026-01-22
Base plasmid: AVD002-Rep2Mut2Cap9-6R-wt.dna (7,104 bp)
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from datetime import datetime
import sys
from pathlib import Path

try:
    from dnachisel import *
    DNACHISEL_AVAILABLE = True
except ImportError:
    print("⚠ DNAchisel not available - using pre-optimized sequences")
    DNACHISEL_AVAILABLE = False

# Project paths
PROJECT_DIR = Path(__file__).parent.parent
PLASMIDS_DIR = PROJECT_DIR / "plasmids"
SYNTHETIC_DIR = PROJECT_DIR / "synthetic_fragments"
DOCS_DIR = PROJECT_DIR / "docs"

# Ensure output directories exist
SYNTHETIC_DIR.mkdir(exist_ok=True)
DOCS_DIR.mkdir(exist_ok=True)

# VHH3 anti-ALPL sequence (119 amino acids)
VHH3_AA = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSTADMGWFRQAPGKGRELVAAVSGSGFSTYSDSVEGRFTISRDNAKRMVYLQMNSLRAEDTAVYYCAKATISLYYAMDVWGQGTTVTVSS"

# Linker sequences
GGGGS5_AA = "GGGGSGGGGSGGGGSGGGGSGGGGS"  # 25 aa
GGGGS1_AA = "GGGGS"  # 5 aa

# Pre-optimized DNA sequences (from previous work)
VHH3_DNA = (
    "GAGGTGCAACTGGTTGAAAGCGGCGGAGGACTTGTTCAACCCGGCGGCAGCCTTAGGCTTTCTTGCGCTGCCAGCGGC"
    "TTCACCTTTAGCACCGCCGACATGGGCTGGTTTAGGCAAGCTCCCGGAAAAGGCAGGGAACTTGTTGCCGCTGTGAGC"
    "GGCAGCGGCTTCAGCACCTACTCTGATAGCGTTGAGGGCAGGTTCACCATCAGCAGGGACAACGCCAAGAGGATGGTG"
    "TACCTGCAGATGAACAGCTTGAGGGCCGAGGACACCGCCGTGTACTACTGCGCCAAGGCCACAATTAGCCTGTACTAC"
    "GCCATGGATGTGTGGGGACAGGGCACCACCGTGACCGTGAGCAGC"
)

GGGGS5_DNA = (
    "GGTGGAGGCGGATCTGGAGGCGGTGGTTCAGGCGGTGGAGGAAGTGGTGGCGGAGGTTCTGGTGGAGGCGGTTCT"
)

GGGGS1_DNA = "GGTGGAGGCGGTTCT"


def load_avd002():
    """Load AVD002 base plasmid from SnapGene .dna file"""
    print("Loading AVD002 base plasmid...")
    avd002_path = PLASMIDS_DIR / "AVD002-Rep2Mut2Cap9-6R-wt.dna"

    if not avd002_path.exists():
        raise FileNotFoundError(f"AVD002 not found at {avd002_path}")

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


def shift_feature_location(location, insertion_point, shift_amount):
    """Shift a feature location if it's after the insertion point."""
    from Bio.SeqFeature import CompoundLocation

    if isinstance(location, CompoundLocation):
        new_parts = []
        for part in location.parts:
            new_parts.append(shift_feature_location(part, insertion_point, shift_amount))
        return CompoundLocation(new_parts, operator=location.operator)

    start = location.start
    end = location.end
    strand = location.strand

    if start >= insertion_point:
        start += shift_amount
        end += shift_amount
    elif start < insertion_point < end:
        end += shift_amount

    return FeatureLocation(start, end, strand=strand)


def create_avd005(avd002_record):
    """
    AVD005: VP1 knockout only

    Single mutation: ATG→AAG at bp 2379 (VP1 start codon)
    VP2/VP3 remain functional
    """
    print("\nCreating AVD005 (VP1 knockout)...")

    avd002_seq = str(avd002_record.seq).upper()
    avd005_seq = avd002_seq[:2378] + "AAG" + avd002_seq[2381:]

    # Verify mutation
    if avd005_seq[2378:2381] != "AAG":
        raise ValueError("VP1 knockout mutation failed")

    print(f"  ✓ VP1 knockout: ATG→AAG at bp 2379-2381")
    print(f"  ✓ Total size: {len(avd005_seq)} bp")

    record = SeqRecord(
        Seq(avd005_seq),
        id="AVD005",
        name="AVD005",
        description="Rep2Mut2Cap9-VP1ko | VP1 knockout helper plasmid",
        annotations={
            "molecule_type": "DNA",
            "topology": "circular",
            "date": datetime.now().strftime("%d-%b-%Y").upper()
        }
    )

    # Copy features and update VP1 annotation
    features = []
    for feat in avd002_record.features:
        new_feat = SeqFeature(
            location=feat.location,
            type=feat.type,
            qualifiers=dict(feat.qualifiers)
        )

        if feat.type == "CDS" and feat.qualifiers.get("label", [""])[0] == "VP1":
            new_feat.qualifiers["label"] = "VP1 (KNOCKOUT)"
            new_feat.qualifiers["note"] = "VP1 start codon mutated ATG→AAG - no translation"

        features.append(new_feat)

    record.features = features
    return record


def create_avd006(avd002_record):
    """
    AVD006: VHH3 in VR4, VP2/3 native

    Insert [GGGGS5]-[VHH3]-[GGGGS1] at bp 3743 (VR4)
    VP2 and VP3 remain functional
    """
    print("\nCreating AVD006 (VHH3 in VR4, VP2/3 native)...")

    avd002_seq = str(avd002_record.seq).upper()

    # Create insert: GGGGS5 + VHH3 + GGGGS1
    insert_dna = GGGGS5_DNA + VHH3_DNA + GGGGS1_DNA
    insertion_point = 3743

    # Verify insert translation
    insert_aa = str(Seq(insert_dna).translate())
    expected_aa = GGGGS5_AA + VHH3_AA + GGGGS1_AA
    if insert_aa != expected_aa:
        raise ValueError(f"Insert translation mismatch:\nExpected: {expected_aa}\nGot: {insert_aa}")

    # Build final sequence
    avd006_seq = avd002_seq[:insertion_point] + insert_dna + avd002_seq[insertion_point:]

    print(f"  ✓ Inserted {len(insert_dna)} bp at bp {insertion_point} (VR4)")
    print(f"  ✓ Total size: {len(avd006_seq)} bp (7104 + 447)")

    record = SeqRecord(
        Seq(avd006_seq),
        id="AVD006",
        name="AVD006",
        description="Rep2Mut2Cap9-VP1-VHH3-VR4 | VHH3 at VR4, VP2/3 native",
        annotations={
            "molecule_type": "DNA",
            "topology": "circular",
            "date": datetime.now().strftime("%d-%b-%Y").upper()
        }
    )

    # Copy and shift features
    features = []
    for feat in avd002_record.features:
        new_location = shift_feature_location(feat.location, insertion_point, len(insert_dna))
        new_feat = SeqFeature(
            location=new_location,
            type=feat.type,
            qualifiers=dict(feat.qualifiers)
        )

        label = feat.qualifiers.get("label", [""])[0]
        if label == "VP1":
            new_feat.qualifiers["label"] = "VP1-VHH3 fusion"
            new_feat.qualifiers["note"] = "AAV9 VP1 with anti-ALPL VHH3 at VR-IV"

        features.append(new_feat)

    # Add VHH3 insert features
    features.append(SeqFeature(
        FeatureLocation(insertion_point, insertion_point + 75),
        type="misc_feature",
        qualifiers={"label": "LINK_D2_N", "note": "N-terminal linker (GGGGS)×5, 25 aa"}
    ))
    features.append(SeqFeature(
        FeatureLocation(insertion_point + 75, insertion_point + 75 + 357, strand=+1),
        type="CDS",
        qualifiers={"label": "VHH3_anti-ALPL", "note": f"Anti-ALPL VHH nanobody, {len(VHH3_AA)} aa", "translation": VHH3_AA}
    ))
    features.append(SeqFeature(
        FeatureLocation(insertion_point + 75 + 357, insertion_point + len(insert_dna)),
        type="misc_feature",
        qualifiers={"label": "LINK_D2_C", "note": "C-terminal linker (GGGGS)×1, 5 aa"}
    ))

    record.features = features
    return record


def create_avd007(avd002_record):
    """
    AVD007: VHH3 in VR4, VP2/3 knockouts

    Modifications:
    1. Insert [GGGGS5]-[VHH3]-[GGGGS1] at bp 3743 (VR4)
    2. VP2 knockout: ACG→ACC at bp 2790 (silent)
    3. VP3 knockout: ATG→CTG at bp 2985 (non-silent)
    """
    print("\nCreating AVD007 (VHH3 in VR4, VP2/3 KO)...")

    avd002_seq = str(avd002_record.seq).upper()

    # Apply knockouts first
    seq = avd002_seq[:2789] + "ACC" + avd002_seq[2792:]
    seq = seq[:2984] + "CTG" + seq[2987:]

    print(f"  ✓ VP2 knockout: ACG→ACC at bp 2790-2792 (silent)")
    print(f"  ✓ VP3 knockout: ATG→CTG at bp 2985-2987 (non-silent)")

    # Create insert
    insert_dna = GGGGS5_DNA + VHH3_DNA + GGGGS1_DNA
    insertion_point = 3743

    # Build final sequence
    avd007_seq = seq[:insertion_point] + insert_dna + seq[insertion_point:]

    print(f"  ✓ Inserted {len(insert_dna)} bp at bp {insertion_point} (VR4)")
    print(f"  ✓ Total size: {len(avd007_seq)} bp (7104 + 447)")

    record = SeqRecord(
        Seq(avd007_seq),
        id="AVD007",
        name="AVD007",
        description="Rep2Mut2Cap9-VP1-VHH3-D2 | VHH3 at VR4, VP2/3 KO",
        annotations={
            "molecule_type": "DNA",
            "topology": "circular",
            "date": datetime.now().strftime("%d-%b-%Y").upper()
        }
    )

    # Copy and shift features
    features = []
    for feat in avd002_record.features:
        new_location = shift_feature_location(feat.location, insertion_point, len(insert_dna))
        new_feat = SeqFeature(
            location=new_location,
            type=feat.type,
            qualifiers=dict(feat.qualifiers)
        )

        label = feat.qualifiers.get("label", [""])[0]
        if label == "VP1":
            new_feat.qualifiers["label"] = "VP1-VHH3 fusion"
            new_feat.qualifiers["note"] = "AAV9 VP1 with anti-ALPL VHH3 at VR-IV"
        elif label == "VP2":
            new_feat.qualifiers["label"] = "VP2 (KNOCKOUT)"
            new_feat.qualifiers["note"] = "VP2 start codon mutated ACG→ACC - no translation"
        elif label == "VP3":
            new_feat.qualifiers["label"] = "VP3 (KNOCKOUT)"
            new_feat.qualifiers["note"] = "VP3 start codon mutated ATG→CTG - no translation"

        features.append(new_feat)

    # Add VHH3 insert features
    features.append(SeqFeature(
        FeatureLocation(insertion_point, insertion_point + 75),
        type="misc_feature",
        qualifiers={"label": "LINK_D2_N", "note": "N-terminal linker (GGGGS)×5, 25 aa"}
    ))
    features.append(SeqFeature(
        FeatureLocation(insertion_point + 75, insertion_point + 75 + 357, strand=+1),
        type="CDS",
        qualifiers={"label": "VHH3_anti-ALPL", "note": f"Anti-ALPL VHH nanobody, {len(VHH3_AA)} aa", "translation": VHH3_AA}
    ))
    features.append(SeqFeature(
        FeatureLocation(insertion_point + 75 + 357, insertion_point + len(insert_dna)),
        type="misc_feature",
        qualifiers={"label": "LINK_D2_C", "note": "C-terminal linker (GGGGS)×1, 5 aa"}
    ))

    record.features = features
    return record


def create_avd008(avd002_record):
    """
    AVD008: VHH3 at VP1 N-terminus, VP2/3 native

    Insert [VHH3]-[GGGGS5] after M-A dipeptide at bp 2384 (after ATG GCT)
    Results in: M-A-[VHH3]-[GGGGS5]-A-D-G-Y-L...
    VP2 and VP3 remain functional
    """
    print("\nCreating AVD008 (VHH3 at VP1 N-term, VP2/3 native)...")

    avd002_seq = str(avd002_record.seq).upper()

    # Create insert: VHH3 + GGGGS5
    insert_dna = VHH3_DNA + GGGGS5_DNA
    insertion_point = 2384  # After ATG GCT (M-A dipeptide)

    # Verify insert translation
    insert_aa = str(Seq(insert_dna).translate())
    expected_aa = VHH3_AA + GGGGS5_AA
    if insert_aa != expected_aa:
        raise ValueError(f"Insert translation mismatch:\nExpected: {expected_aa}\nGot: {insert_aa}")

    # Build final sequence
    avd008_seq = avd002_seq[:insertion_point] + insert_dna + avd002_seq[insertion_point:]

    print(f"  ✓ Inserted {len(insert_dna)} bp at bp {insertion_point} (VP1 N-term)")
    print(f"  ✓ Total size: {len(avd008_seq)} bp (7104 + 432)")

    record = SeqRecord(
        Seq(avd008_seq),
        id="AVD008",
        name="AVD008",
        description="Rep2Mut2Cap9-VP1n-VHH3 | VHH3 at VP1 N-term, VP2/3 native",
        annotations={
            "molecule_type": "DNA",
            "topology": "circular",
            "date": datetime.now().strftime("%d-%b-%Y").upper()
        }
    )

    # Copy and shift features
    features = []
    for feat in avd002_record.features:
        new_location = shift_feature_location(feat.location, insertion_point, len(insert_dna))
        new_feat = SeqFeature(
            location=new_location,
            type=feat.type,
            qualifiers=dict(feat.qualifiers)
        )

        label = feat.qualifiers.get("label", [""])[0]
        if label == "VP1":
            new_feat.qualifiers["label"] = "VP1-VHH3 N-term fusion"
            new_feat.qualifiers["note"] = "AAV9 VP1 with anti-ALPL VHH3 at N-terminus"

        features.append(new_feat)

    # Add VHH3 insert features
    features.append(SeqFeature(
        FeatureLocation(insertion_point, insertion_point + 357, strand=+1),
        type="CDS",
        qualifiers={"label": "VHH3_anti-ALPL", "note": f"Anti-ALPL VHH nanobody, {len(VHH3_AA)} aa", "translation": VHH3_AA}
    ))
    features.append(SeqFeature(
        FeatureLocation(insertion_point + 357, insertion_point + len(insert_dna)),
        type="misc_feature",
        qualifiers={"label": "LINK_GGGGS5", "note": "Linker (GGGGS)×5, 25 aa"}
    ))

    record.features = features
    return record


def create_avd009(avd002_record):
    """
    AVD009: VHH3 at VP1 N-terminus, VP2/3 knockouts

    Modifications:
    1. Insert [VHH3]-[GGGGS5] after M-A dipeptide at bp 2384 (after ATG GCT)
       Results in: M-A-[VHH3]-[GGGGS5]-A-D-G-Y-L...
    2. VP2 knockout: ACG→ACC at bp 2790 (silent)
    3. VP3 knockout: ATG→CTG at bp 2985 (non-silent)
    """
    print("\nCreating AVD009 (VHH3 at VP1 N-term, VP2/3 KO)...")

    avd002_seq = str(avd002_record.seq).upper()

    # Apply knockouts first (BEFORE insertion, so positions don't shift)
    seq = avd002_seq[:2789] + "ACC" + avd002_seq[2792:]
    seq = seq[:2984] + "CTG" + seq[2987:]

    print(f"  ✓ VP2 knockout: ACG→ACC at bp 2790-2792 (silent)")
    print(f"  ✓ VP3 knockout: ATG→CTG at bp 2985-2987 (non-silent)")

    # Create insert
    insert_dna = VHH3_DNA + GGGGS5_DNA
    insertion_point = 2384  # After ATG GCT (M-A dipeptide)

    # Build final sequence
    avd009_seq = seq[:insertion_point] + insert_dna + seq[insertion_point:]

    print(f"  ✓ Inserted {len(insert_dna)} bp at bp {insertion_point} (VP1 N-term)")
    print(f"  ✓ Total size: {len(avd009_seq)} bp (7104 + 432)")

    record = SeqRecord(
        Seq(avd009_seq),
        id="AVD009",
        name="AVD009",
        description="Rep2Mut2Cap9-VP1n-VHH3-D2 | VHH3 at VP1 N-term, VP2/3 KO",
        annotations={
            "molecule_type": "DNA",
            "topology": "circular",
            "date": datetime.now().strftime("%d-%b-%Y").upper()
        }
    )

    # Copy and shift features
    features = []
    for feat in avd002_record.features:
        new_location = shift_feature_location(feat.location, insertion_point, len(insert_dna))
        new_feat = SeqFeature(
            location=new_location,
            type=feat.type,
            qualifiers=dict(feat.qualifiers)
        )

        label = feat.qualifiers.get("label", [""])[0]
        if label == "VP1":
            new_feat.qualifiers["label"] = "VP1-VHH3 N-term fusion"
            new_feat.qualifiers["note"] = "AAV9 VP1 with anti-ALPL VHH3 at N-terminus"
        elif label == "VP2":
            new_feat.qualifiers["label"] = "VP2 (KNOCKOUT)"
            new_feat.qualifiers["note"] = "VP2 start codon mutated ACG→ACC - no translation"
        elif label == "VP3":
            new_feat.qualifiers["label"] = "VP3 (KNOCKOUT)"
            new_feat.qualifiers["note"] = "VP3 start codon mutated ATG→CTG - no translation"

        features.append(new_feat)

    # Add VHH3 insert features
    features.append(SeqFeature(
        FeatureLocation(insertion_point, insertion_point + 357, strand=+1),
        type="CDS",
        qualifiers={"label": "VHH3_anti-ALPL", "note": f"Anti-ALPL VHH nanobody, {len(VHH3_AA)} aa", "translation": VHH3_AA}
    ))
    features.append(SeqFeature(
        FeatureLocation(insertion_point + 357, insertion_point + len(insert_dna)),
        type="misc_feature",
        qualifiers={"label": "LINK_GGGGS5", "note": "Linker (GGGGS)×5, 25 aa"}
    ))

    record.features = features
    return record


def create_avd010(avd002_record):
    """
    AVD010: VHH3 at VP2 N-terminus, VP1/3 native

    Insert [VHH3]-[GGGGS5] after M-A dipeptide at bp 2795 (after ACG GCT)
    Results in: M-A-[VHH3]-[GGGGS5]-P-G-K-K-R...
    Note: ACG translates to M at translation initiation
    VP1 and VP3 remain functional
    """
    print("\nCreating AVD010 (VHH3 at VP2 N-term, VP1/3 native)...")

    avd002_seq = str(avd002_record.seq).upper()

    # Create insert: VHH3 + GGGGS5
    insert_dna = VHH3_DNA + GGGGS5_DNA
    insertion_point = 2795  # After ACG GCT (M-A dipeptide)

    # Verify insert translation
    insert_aa = str(Seq(insert_dna).translate())
    expected_aa = VHH3_AA + GGGGS5_AA
    if insert_aa != expected_aa:
        raise ValueError(f"Insert translation mismatch:\nExpected: {expected_aa}\nGot: {insert_aa}")

    # Build final sequence
    avd010_seq = avd002_seq[:insertion_point] + insert_dna + avd002_seq[insertion_point:]

    print(f"  ✓ Inserted {len(insert_dna)} bp at bp {insertion_point} (VP2 N-term)")
    print(f"  ✓ Total size: {len(avd010_seq)} bp (7104 + 432)")

    record = SeqRecord(
        Seq(avd010_seq),
        id="AVD010",
        name="AVD010",
        description="Rep2Mut2Cap9-VP2n-VHH3 | VHH3 at VP2 N-term, VP1/3 native",
        annotations={
            "molecule_type": "DNA",
            "topology": "circular",
            "date": datetime.now().strftime("%d-%b-%Y").upper()
        }
    )

    # Copy and shift features
    features = []
    for feat in avd002_record.features:
        new_location = shift_feature_location(feat.location, insertion_point, len(insert_dna))
        new_feat = SeqFeature(
            location=new_location,
            type=feat.type,
            qualifiers=dict(feat.qualifiers)
        )

        label = feat.qualifiers.get("label", [""])[0]
        if label == "VP2":
            new_feat.qualifiers["label"] = "VP2-VHH3 N-term fusion"
            new_feat.qualifiers["note"] = "AAV9 VP2 with anti-ALPL VHH3 at N-terminus"

        features.append(new_feat)

    # Add VHH3 insert features
    features.append(SeqFeature(
        FeatureLocation(insertion_point, insertion_point + 357, strand=+1),
        type="CDS",
        qualifiers={"label": "VHH3_anti-ALPL", "note": f"Anti-ALPL VHH nanobody, {len(VHH3_AA)} aa", "translation": VHH3_AA}
    ))
    features.append(SeqFeature(
        FeatureLocation(insertion_point + 357, insertion_point + len(insert_dna)),
        type="misc_feature",
        qualifiers={"label": "LINK_GGGGS5", "note": "Linker (GGGGS)×5, 25 aa"}
    ))

    record.features = features
    return record


def create_avd011(avd002_record):
    """
    AVD011: VP2 knockout only

    Single mutation: ACG→ACC at bp 2790 (VP2 start codon, silent T→T)
    VP1 and VP3 remain functional
    """
    print("\nCreating AVD011 (VP2 knockout)...")

    avd002_seq = str(avd002_record.seq).upper()
    avd011_seq = avd002_seq[:2789] + "ACC" + avd002_seq[2792:]

    # Verify mutation
    if avd011_seq[2789:2792] != "ACC":
        raise ValueError("VP2 knockout mutation failed")

    print(f"  ✓ VP2 knockout: ACG→ACC at bp 2790-2792 (silent, Thr→Thr)")
    print(f"  ✓ Total size: {len(avd011_seq)} bp")

    record = SeqRecord(
        Seq(avd011_seq),
        id="AVD011",
        name="AVD011",
        description="Rep2Mut2Cap9-VP2ko | VP2 knockout, VP1/3 functional",
        annotations={
            "molecule_type": "DNA",
            "topology": "circular",
            "date": datetime.now().strftime("%d-%b-%Y").upper()
        }
    )

    # Copy features and update VP2 annotation
    features = []
    for feat in avd002_record.features:
        new_feat = SeqFeature(
            location=feat.location,
            type=feat.type,
            qualifiers=dict(feat.qualifiers)
        )

        if feat.type == "CDS" and feat.qualifiers.get("label", [""])[0] == "VP2":
            new_feat.qualifiers["label"] = "VP2 (KNOCKOUT)"
            new_feat.qualifiers["note"] = "VP2 start codon mutated ACG→ACC - no translation"

        features.append(new_feat)

    record.features = features
    return record


def verify_construct(record, construct_id, expected_mutations, expected_inserts):
    """Verify construct has correct mutations and insertions"""
    print(f"\nVerifying {construct_id}...")
    seq = str(record.seq)

    checks = []

    # Check mutations
    for mutation_name, position, expected_codon in expected_mutations:
        actual_codon = seq[position:position+3]
        checks.append((mutation_name, actual_codon == expected_codon, f"Expected {expected_codon}, got {actual_codon}"))

    # Check inserts
    for insert_name, position, expected_length in expected_inserts:
        insert_region = seq[position:position+expected_length]
        insert_translation = str(Seq(insert_region).translate())
        has_stop = "*" in insert_translation
        checks.append((f"{insert_name} length", len(insert_region) == expected_length, f"{len(insert_region)} bp"))
        checks.append((f"{insert_name} no stops", not has_stop, "Clean" if not has_stop else "Stop found"))

    # Check total size
    expected_size = 7104 + sum(length for _, _, length in expected_inserts)
    checks.append(("Total size", len(seq) == expected_size, f"{len(seq)} bp (expected {expected_size})"))

    # Print results
    all_passed = True
    for check_name, passed, details in checks:
        status = "✓" if passed else "✗"
        print(f"  {status} {check_name}: {details}")
        if not passed:
            all_passed = False

    if all_passed:
        print(f"  ✅ {construct_id} verification PASSED")
    else:
        print(f"  ❌ {construct_id} verification FAILED")
        raise ValueError(f"{construct_id} verification failed")

    return all_passed


def export_genbank(record, output_path):
    """Write GenBank file"""
    SeqIO.write(record, output_path, "genbank")
    print(f"  ✓ Wrote {output_path}")


def create_verification_report(records_dict):
    """Create comprehensive verification report for all constructs"""
    report_path = DOCS_DIR / "DESIGN_VERIFICATION_AVD005_AVD011.md"

    with open(report_path, 'w') as f:
        f.write("# Design Verification Report: AVD005-AVD011\n\n")
        f.write(f"**Date:** {datetime.now().strftime('%Y-%m-%d')}\n\n")
        f.write(f"**Base plasmid:** AVD002-Rep2Mut2Cap9-6R-wt.dna (7,104 bp)\n\n")

        # Summary table
        f.write("## Summary of Constructs\n\n")
        f.write("| ID | Size (bp) | VP1 | VP2 | VP3 | Description |\n")
        f.write("|----|-----------|-----|-----|-----|-------------|\n")

        construct_info = {
            "AVD005": (7104, "KO (AAG)", "Native (ACG)", "Native (ATG)", "VP1 knockout helper"),
            "AVD006": (7551, "VHH3-VR4", "Native (ACG)", "Native (ATG)", "VHH3 at VR4, VP2/3 functional"),
            "AVD007": (7551, "VHH3-VR4", "KO (ACC)", "KO (CTG)", "VHH3 at VR4, VP2/3 KO"),
            "AVD008": (7536, "VHH3-Nterm", "Native (ACG)", "Native (ATG)", "VHH3 at VP1 N-term, VP2/3 functional"),
            "AVD009": (7536, "VHH3-Nterm", "KO (ACC)", "KO (CTG)", "VHH3 at VP1 N-term, VP2/3 KO"),
            "AVD010": (7536, "Native (ATG)", "VHH3-Nterm", "Native (ATG)", "VHH3 at VP2 N-term, VP1/3 functional"),
            "AVD011": (7104, "Native (ATG)", "KO (ACC)", "Native (ATG)", "VP2 knockout helper")
        }

        for construct_id in sorted(records_dict.keys()):
            size, vp1, vp2, vp3, desc = construct_info[construct_id]
            f.write(f"| {construct_id} | {size:,} | {vp1} | {vp2} | {vp3} | {desc} |\n")

        f.write("\n---\n\n")

        # Detailed sections for each construct
        for construct_id in sorted(records_dict.keys()):
            record = records_dict[construct_id]
            seq = str(record.seq)

            f.write(f"## {construct_id}: {construct_info[construct_id][4]}\n\n")
            f.write(f"**Size:** {len(seq):,} bp\n\n")

            f.write("**Modifications:**\n")
            if "005" in construct_id:
                f.write("- VP1 knockout: ATG→AAG at bp 2379-2381\n")
            elif "006" in construct_id:
                f.write("- VHH3 insert at VR4: 447 bp (GGGGS5-VHH3-GGGGS1) at bp 3743\n")
            elif "007" in construct_id:
                f.write("- VP2 knockout: ACG→ACC at bp 2790-2792 (silent)\n")
                f.write("- VP3 knockout: ATG→CTG at bp 2985-2987 (non-silent)\n")
                f.write("- VHH3 insert at VR4: 447 bp (GGGGS5-VHH3-GGGGS1) at bp 3743\n")
            elif "008" in construct_id:
                f.write("- VHH3 insert at VP1 N-term: 432 bp (VHH3-GGGGS5) at bp 2384\n")
                f.write("- Sequence: M-A-[VHH3]-[GGGGS5]-A-D-G-Y-L...\n")
            elif "009" in construct_id:
                f.write("- VP2 knockout: ACG→ACC at bp 2790-2792 (silent)\n")
                f.write("- VP3 knockout: ATG→CTG at bp 2985-2987 (non-silent)\n")
                f.write("- VHH3 insert at VP1 N-term: 432 bp (VHH3-GGGGS5) at bp 2384\n")
                f.write("- Sequence: M-A-[VHH3]-[GGGGS5]-A-D-G-Y-L...\n")
            elif "010" in construct_id:
                f.write("- VHH3 insert at VP2 N-term: 432 bp (VHH3-GGGGS5) at bp 2795\n")
                f.write("- Sequence: M-A-[VHH3]-[GGGGS5]-P-G-K-K-R... (ACG→M at translation)\n")
            elif "011" in construct_id:
                f.write("- VP2 knockout: ACG→ACC at bp 2790-2792 (silent)\n")

            f.write("\n")

            # Verification checklist
            f.write("**Verification:**\n")

            # Check start codons
            if "005" in construct_id:
                f.write(f"- [{'x' if seq[2378:2381] == 'AAG' else ' '}] VP1 = AAG (knockout)\n")
                f.write(f"- [{'x' if seq[2789:2792] == 'ACG' else ' '}] VP2 = ACG (native)\n")
                f.write(f"- [{'x' if seq[2984:2987] == 'ATG' else ' '}] VP3 = ATG (native)\n")
            elif "011" in construct_id:
                f.write(f"- [{'x' if seq[2378:2381] == 'ATG' else ' '}] VP1 = ATG (native)\n")
                f.write(f"- [{'x' if seq[2789:2792] == 'ACC' else ' '}] VP2 = ACC (knockout)\n")
                f.write(f"- [{'x' if seq[2984:2987] == 'ATG' else ' '}] VP3 = ATG (native)\n")
            elif "006" in construct_id:
                f.write(f"- [{'x' if seq[2378:2381] == 'ATG' else ' '}] VP1 = ATG (native, has VHH3 insert)\n")
                f.write(f"- [{'x' if seq[2789:2792] == 'ACG' else ' '}] VP2 = ACG (native)\n")
                f.write(f"- [{'x' if seq[2984:2987] == 'ATG' else ' '}] VP3 = ATG (native)\n")
                f.write(f"- [x] VHH3 insert at bp 3743: 447 bp\n")
            elif "007" in construct_id:
                f.write(f"- [{'x' if seq[2378:2381] == 'ATG' else ' '}] VP1 = ATG (native, has VHH3 insert)\n")
                f.write(f"- [{'x' if seq[2789:2792] == 'ACC' else ' '}] VP2 = ACC (knockout)\n")
                f.write(f"- [{'x' if seq[2984:2987] == 'CTG' else ' '}] VP3 = CTG (knockout)\n")
                f.write(f"- [x] VHH3 insert at bp 3743: 447 bp\n")
            elif "008" in construct_id:
                f.write(f"- [x] VP1 N-term VHH3 insert at bp 2384 (after M-A): 432 bp\n")
                f.write(f"- [{'x' if seq[2789+432:2792+432] == 'ACG' else ' '}] VP2 = ACG (native)\n")
                f.write(f"- [{'x' if seq[2984+432:2987+432] == 'ATG' else ' '}] VP3 = ATG (native)\n")
            elif "009" in construct_id:
                f.write(f"- [x] VP1 N-term VHH3 insert at bp 2384 (after M-A): 432 bp\n")
                f.write(f"- [{'x' if seq[2789+432:2792+432] == 'ACC' else ' '}] VP2 = ACC (knockout)\n")
                f.write(f"- [{'x' if seq[2984+432:2987+432] == 'CTG' else ' '}] VP3 = CTG (knockout)\n")
            elif "010" in construct_id:
                f.write(f"- [{'x' if seq[2378:2381] == 'ATG' else ' '}] VP1 = ATG (native)\n")
                f.write(f"- [x] VP2 N-term VHH3 insert at bp 2795 (after M-A): 432 bp\n")
                f.write(f"- [{'x' if seq[2984+432:2987+432] == 'ATG' else ' '}] VP3 = ATG (native)\n")

            f.write(f"- [x] Total size: {len(seq):,} bp\n")
            f.write("\n---\n\n")

        # VHH3 sequences
        f.write("## VHH3 Insert Sequences\n\n")
        f.write(f"**VHH3 protein sequence** (119 aa):\n```\n{VHH3_AA}\n```\n\n")
        f.write(f"**VHH3 DNA sequence** (357 bp):\n```\n{VHH3_DNA}\n```\n\n")
        f.write(f"**GGGGS5 linker** (25 aa, 75 bp):\n```\n{GGGGS5_DNA}\n```\n\n")
        f.write(f"**GGGGS1 linker** (5 aa, 15 bp):\n```\n{GGGGS1_DNA}\n```\n\n")

        # Next steps
        f.write("## Next Steps\n\n")
        f.write("1. Sequence verify all plasmids (Sanger sequencing)\n")
        f.write("2. Production screening:\n")
        f.write("   - AVD005 + AVD007: VP1-VHH3 only display\n")
        f.write("   - AVD005 + AVD009: VP1-VHH3 N-term only display\n")
        f.write("   - Compare display efficiency across constructs\n")
        f.write("3. Functional assays:\n")
        f.write("   - Western blot (VP expression, VHH3 detection)\n")
        f.write("   - ELISA (anti-ALPL binding)\n")
        f.write("   - Transduction efficiency\n")
        f.write("   - ALPL-targeted cell binding\n\n")

        f.write("---\n")
        f.write(f"*Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*\n")

    print(f"  ✓ Wrote {report_path}")


def main():
    """Main build pipeline"""
    print("=" * 70)
    print("Building AVD005-AVD011 AAV Plasmid Constructs")
    print("=" * 70)

    # Load base plasmid
    avd002_record = load_avd002()

    records_dict = {}

    # Create all constructs
    print("\n" + "=" * 70)
    print("Creating constructs...")
    print("=" * 70)

    # AVD005
    avd005 = create_avd005(avd002_record)
    verify_construct(avd005, "AVD005",
                    [("VP1 knockout", 2378, "AAG"),
                     ("VP2 preserved", 2789, "ACG"),
                     ("VP3 preserved", 2984, "ATG")],
                    [])
    records_dict["AVD005"] = avd005

    # AVD006
    avd006 = create_avd006(avd002_record)
    verify_construct(avd006, "AVD006",
                    [("VP1 preserved", 2378, "ATG"),
                     ("VP2 preserved", 2789, "ACG"),
                     ("VP3 preserved", 2984, "ATG")],
                    [("VR4 insert", 3743, 447)])
    records_dict["AVD006"] = avd006

    # AVD007
    avd007 = create_avd007(avd002_record)
    verify_construct(avd007, "AVD007",
                    [("VP1 preserved", 2378, "ATG"),
                     ("VP2 knockout", 2789, "ACC"),
                     ("VP3 knockout", 2984, "CTG")],
                    [("VR4 insert", 3743, 447)])
    records_dict["AVD007"] = avd007

    # AVD008
    avd008 = create_avd008(avd002_record)
    verify_construct(avd008, "AVD008",
                    [("VP1 preserved", 2378, "ATG"),
                     ("VP2 preserved", 2789 + 432, "ACG"),
                     ("VP3 preserved", 2984 + 432, "ATG")],
                    [("VP1 N-term insert", 2384, 432)])
    records_dict["AVD008"] = avd008

    # AVD009
    avd009 = create_avd009(avd002_record)
    verify_construct(avd009, "AVD009",
                    [("VP1 preserved", 2378, "ATG"),
                     ("VP2 knockout", 2789 + 432, "ACC"),
                     ("VP3 knockout", 2984 + 432, "CTG")],
                    [("VP1 N-term insert", 2384, 432)])
    records_dict["AVD009"] = avd009

    # AVD010
    avd010 = create_avd010(avd002_record)
    verify_construct(avd010, "AVD010",
                    [("VP1 preserved", 2378, "ATG"),
                     ("VP2 preserved", 2789, "ACG"),
                     ("VP3 preserved", 2984 + 432, "ATG")],
                    [("VP2 N-term insert", 2795, 432)])
    records_dict["AVD010"] = avd010

    # AVD011
    avd011 = create_avd011(avd002_record)
    verify_construct(avd011, "AVD011",
                    [("VP1 preserved", 2378, "ATG"),
                     ("VP2 knockout", 2789, "ACC"),
                     ("VP3 preserved", 2984, "ATG")],
                    [])
    records_dict["AVD011"] = avd011

    # Export all GenBank files
    print("\n" + "=" * 70)
    print("Writing GenBank files...")
    print("=" * 70)

    export_genbank(avd005, PLASMIDS_DIR / "AVD005-Rep2Mut2Cap9-VP1ko.gb")
    export_genbank(avd006, PLASMIDS_DIR / "AVD006-Rep2Mut2Cap9-VP1-VHH3-VR4.gb")
    export_genbank(avd007, PLASMIDS_DIR / "AVD007-Rep2Mut2Cap9-VP1-VHH3-D2.gb")
    export_genbank(avd008, PLASMIDS_DIR / "AVD008-Rep2Mut2Cap9-VP1n-VHH3.gb")
    export_genbank(avd009, PLASMIDS_DIR / "AVD009-Rep2Mut2Cap9-VP1n-VHH3-D2.gb")
    export_genbank(avd010, PLASMIDS_DIR / "AVD010-Rep2Mut2Cap9-VP2n-VHH3.gb")
    export_genbank(avd011, PLASMIDS_DIR / "AVD011-Rep2Mut2Cap9-VP2ko.gb")

    # Create verification report
    print("\n" + "=" * 70)
    print("Creating verification report...")
    print("=" * 70)
    create_verification_report(records_dict)

    print("\n" + "=" * 70)
    print("✅ Build completed successfully")
    print("=" * 70)
    print("\nGenerated files:")
    for construct_id, record in sorted(records_dict.items()):
        print(f"  - {construct_id}: {len(record.seq):,} bp")
    print(f"\nDocumentation:")
    print(f"  - {DOCS_DIR / 'DESIGN_VERIFICATION_AVD005_AVD011.md'}")
    print()


if __name__ == "__main__":
    main()
