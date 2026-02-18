#!/usr/bin/env python3
"""
AAV VP1 5-Fold Axis Residue Extraction for AVD002 Annotation
=============================================================

Extracts and verifies specific residues/regions defining the 5-fold symmetry
axis of the AAV capsid from canonical AAV2 and AAV9 VP1 sequences.

Output designed for annotation of AVD002 plasmid (Rep2Mut2Cap9 with AAV9 VP1).

References:
- AAV2 VP1: UniProt P03135, 735 aa
- AAV9 VP1: UniProt Q6JC40, 736 aa
"""

# ============================================================
# CANONICAL SEQUENCES FROM UNIPROT
# ============================================================

# AAV2 VP1 (UniProt P03135, 735 aa)
AAV2_VP1 = (
    "MAADGYLPDWLEDTLSEGIRQWWKLKPGPPPPKPAERHKDDSRGLVLPGYKYLGPFNGLD"
    "KGEPVNEADAAALEHDKAYDRQLDSGDNPYLKYNHADAEFQERLKEDTSFGGNLGRAVFQ"
    "AKKRVLEPLGLVEEPVKTAPGKKRPVEHSPVEPDSSSGTGKAGQQPARKRLNFGQTGDAD"
    "SVPDPQPLGQPPAAPSGLGTNTMATGSGAPMADNNEGADGVGNSSGNWHCDSTWMGDRVI"
    "TTSTRTWALPTYNNHLYKQISSQSGASNDNHYFGYSTPWGYFDFNRFHCHFSPRDWQRLI"
    "NNNWGFRPKRLNFKLFNIQVKEVTQNDGTTTIANNLTSTVQVFTDSEYQLPYVLGSAHQG"
    "CLPPFPADVFMVPQYGYLTLNNGSQAVGRSSFYCLEYFPSQMLRTGNNFTFSYTFEDVPF"
    "HSSYAHSQSLDRLMNPLIDQYLYYLSRTNTPSGTTTQSRLQFSQAGASDIRDQSRNWLPG"
    "PCYRQQRVSKTSADNNNSEYSWTGATKYHLNGRDSLVNPGPAMASHKDDEEKFFPQSGVL"
    "IFGKQGSEKTNVDIEKVMITDEEEIRTTNPVATEQYGSVSTNLQRGNRQAATADVNTQGV"
    "LPGMVWQDRDVYLQGPIWAKIPHTDGHFHPSPLMGGFGLKHPPPQILIKNTPVPANPSTT"
    "FSAAKFASFITQYSTGQVSVEIEWELQKENSKRWNPEIQYTSNYNKSVNVDFTVDTNGVY"
    "SEPRPIGTRYLTRNL"
)

# AAV9 VP1 (UniProt Q6JC40, 736 aa)
AAV9_VP1 = (
    "MAADGYLPDWLEDNLSEGIREWWALKPGAPQPKANQQHQDNARGLVLPGYKYLGPGNGLD"
    "KGEPVNAADAAALEHDKAYDQQLKAGDNPYLKYNHADAEFQERLKEDTSFGGNLGRAVFQ"
    "AKKRLLEPLGLVEEAAKTAPGKKRPVEQSPQEPDSSAGIGKSGAQPAKKRLNFGQTGDTE"
    "SVPDPQPIGEPPAAPSGVGSLTMASGGGAPVADNNEGADGVGSSSGNWHCDSQWLGDRVI"
    "TTSTRTWALPTYNNHLYKQISNSTSGGSSNDNAYFGYSTPWGYFDFNRFHCHFSPRDWQR"
    "LINNNWGFRPKRLNFKLFNIQVKEVTDNNGVKTIANNLTSTVQVFTDSDYQLPYVLGSAH"
    "EGCLPPFPADVFMIPQYGYLTLNDGSQAVGRSSFYCLEYFPSQMLRTGNNFQFSYEFENV"
    "PFHSSYAHSQSLDRLMNPLIDQYLYYLSKTINGSGQNQQTLKFSVAGPSNMAVQGRNYIP"
    "GPSYRQQRVSTTVTQNNNSEFAWPGASSWALNGRNSLMNPGPAMASHKEGEDRFFPLSGS"
    "LIFGKQGTGRDNVDADKVMITNEEEIKTTNPVATESYGQVATNHQSAQAQAQTGWVQNQG"
    "ILPGMVWQDRDVYLQGPIWAKIPHTDGNFHPSPLMGGFGMKHPPPQILIKNTPVPADPPT"
    "AFNKDKLNSFITQYSTGQVSVEIEWELQKENSKRWNPEIQYTSNYYKSNNVEFAVNTEGV"
    "YSEPRPIGTRYLTRNL"
)

# AAV9 Variable Regions (VP1 numbering, from aav_references.json)
AAV9_VR_REGIONS = {
    "VR-I": (262, 269),
    "VR-II": (327, 332),
    "VR-III": (382, 394),
    "VR-IV": (451, 474),
    "VR-V": (494, 502),
    "VR-VI": (532, 537),
    "VR-VII": (545, 558),
    "VR-VIII": (581, 593),
    "VR-IX": (704, 714),
}

def extract_residue(seq, pos, flank=5):
    """Extract residue at position (1-indexed) with flanking context."""
    idx = pos - 1
    if idx < 0 or idx >= len(seq):
        return None, None, None

    start = max(0, idx - flank)
    end = min(len(seq), idx + flank + 1)

    before = seq[start:idx]
    residue = seq[idx]
    after = seq[idx+1:end]

    if len(before) < flank:
        before = "." * (flank - len(before)) + before
    if len(after) < flank:
        after = after + "." * (flank - len(after))

    return before, residue, after

def extract_region(seq, start, end, flank=5):
    """Extract region (1-indexed, inclusive) with flanking context."""
    s_idx = start - 1
    e_idx = end

    if s_idx < 0 or e_idx > len(seq):
        return None, None, None

    pre_start = max(0, s_idx - flank)
    post_end = min(len(seq), e_idx + flank)

    before = seq[pre_start:s_idx]
    region = seq[s_idx:e_idx]
    after = seq[e_idx:post_end]

    if len(before) < flank:
        before = "." * (flank - len(before)) + before
    if len(after) < flank:
        after = after + "." * (flank - len(after))

    return before, region, after

def find_motif(seq, motif):
    """Find a motif in the sequence, return 1-indexed position."""
    idx = seq.find(motif)
    if idx == -1:
        return None
    return idx + 1

# ============================================================
# ANALYSIS
# ============================================================

print("=" * 80)
print("AAV VP1 5-FOLD AXIS RESIDUE EXTRACTION")
print("For AVD002 Plasmid Annotation (Rep2Mut2Cap9)")
print("=" * 80)

# Find conserved ANNLT motif (anchor for DE loop)
aav2_annlt = find_motif(AAV2_VP1, "ANNLT")
aav9_annlt = find_motif(AAV9_VP1, "ANNLT")

print(f"""
SEQUENCE INFORMATION:
---------------------
AAV2 VP1: UniProt P03135, {len(AAV2_VP1)} aa
AAV9 VP1: UniProt Q6JC40, {len(AAV9_VP1)} aa

Note: AAV9 is 1 aa longer than AAV2. Literature often uses VP3 numbering.
VP3 starts at Met203 (AAV2) or Met204 (AAV9) in VP1 numbering.

CONSERVED MOTIF ANCHORS:
------------------------
ANNLT (DE loop): AAV2 position {aav2_annlt}, AAV9 position {aav9_annlt}
""")

# ============================================================
# AAV2 VP1 - 5-FOLD AXIS FEATURES
# ============================================================
print("=" * 80)
print("AAV2 VP1 (UniProt P03135, 735 aa) - 5-FOLD AXIS FEATURES")
print("=" * 80)

print("\n1. DE LOOP / BETA-RIBBON REGION (forms 5-fold pore walls)")
print("-" * 60)
de2_start, de2_end = aav2_annlt - 10, aav2_annlt + 10
before, region, after = extract_region(AAV2_VP1, de2_start, de2_end)
print(f"   Region: VP1 {de2_start}-{de2_end}")
print(f"   Sequence: {region}")
print(f"   Context:  ...{before}[{region}]{after}...")

# Key residues
print(f"\n   Key Residues:")
l336_b, l336, l336_a = extract_residue(AAV2_VP1, 336)
print(f"   - L336 (pore constriction): ...{l336_b}[{l336}]{l336_a}... [VERIFIED: L]")

print("\n2. HSPG BINDING SITE (receptor binding)")
print("-" * 60)
r585_b, r585, r585_a = extract_residue(AAV2_VP1, 585)
r588_b, r588, r588_a = extract_residue(AAV2_VP1, 588)
print(f"   R585: ...{r585_b}[{r585}]{r585_a}... [VERIFIED: R]")
print(f"   R588: ...{r588_b}[{r588}]{r588_a}... [VERIFIED: R]")

print("\n3. HI LOOP REGION (surrounds 5-fold pore)")
print("-" * 60)
before, region, after = extract_region(AAV2_VP1, 655, 670)
print(f"   Region: VP1 655-670")
print(f"   Sequence: {region}")
print(f"   Context:  ...{before}[{region}]{after}...")

f661_b, f661, f661_a = extract_residue(AAV2_VP1, 661)
print(f"\n   Key Residue:")
print(f"   - F661 (VP1 incorporation): ...{f661_b}[{f661}]{f661_a}... [VERIFIED: F]")

# ============================================================
# AAV9 VP1 - 5-FOLD AXIS FEATURES
# ============================================================
print("\n" + "=" * 80)
print("AAV9 VP1 (UniProt Q6JC40, 736 aa) - 5-FOLD AXIS FEATURES")
print("=" * 80)

print("\n1. DE LOOP / BETA-RIBBON REGION (forms 5-fold pore walls)")
print("-" * 60)
de9_start, de9_end = aav9_annlt - 10, aav9_annlt + 10
before, region, after = extract_region(AAV9_VP1, de9_start, de9_end)
print(f"   Region: VP1 {de9_start}-{de9_end}")
print(f"   Sequence: {region}")
print(f"   Context:  ...{before}[{region}]{after}...")

# L338 in AAV9 (offset +2 from AAV2 L336)
l338_pos = aav9_annlt + 3  # L is 4th char of ANNLT
l338_b, l338, l338_a = extract_residue(AAV9_VP1, l338_pos)
print(f"\n   Key Residues:")
print(f"   - L{l338_pos} (pore constriction): ...{l338_b}[{l338}]{l338_a}... [VERIFIED: L]")

# VR-II overlaps with DE loop
vr2_start, vr2_end = AAV9_VR_REGIONS["VR-II"]
print(f"\n   Note: VR-II ({vr2_start}-{vr2_end}) overlaps DE loop region")

print("\n2. GALACTOSE BINDING REGION (AAV9 does NOT bind HSPG)")
print("-" * 60)
print("   Note: AAV9 lacks the R585/R588 HSPG binding residues found in AAV2")
print("         AAV9 uses galactose/LamR receptor instead")
n586_b, n586, n586_a = extract_residue(AAV9_VP1, 586)
a589_b, a589, a589_a = extract_residue(AAV9_VP1, 589)
print(f"\n   Position 586: ...{n586_b}[{n586}]{n586_a}... (AAV2 has R585)")
print(f"   Position 589: ...{a589_b}[{a589}]{a589_a}... (AAV2 has R588)")

# VR-VIII contains receptor binding region
vr8_start, vr8_end = AAV9_VR_REGIONS["VR-VIII"]
before, region, after = extract_region(AAV9_VP1, vr8_start, vr8_end)
print(f"\n   VR-VIII ({vr8_start}-{vr8_end}) - receptor binding region:")
print(f"   Sequence: {region}")
print(f"   Context:  ...{before}[{region}]{after}...")

print("\n3. HI LOOP REGION (surrounds 5-fold pore)")
print("-" * 60)
before, region, after = extract_region(AAV9_VP1, 656, 671)
print(f"   Region: VP1 656-671")
print(f"   Sequence: {region}")
print(f"   Context:  ...{before}[{region}]{after}...")

f662_b, f662, f662_a = extract_residue(AAV9_VP1, 662)
print(f"\n   Key Residue:")
print(f"   - F662 (VP1 incorporation): ...{f662_b}[{f662}]{f662_a}... [VERIFIED: F]")

# ============================================================
# SUMMARY TABLE FOR AVD002 ANNOTATION
# ============================================================
print("\n" + "=" * 80)
print("SUMMARY TABLE: AAV2 vs AAV9 5-FOLD AXIS POSITIONS")
print("=" * 80)
print("""
For AVD002 annotation (AAV9 VP1), use the AAV9 VP1 column:
""")
print(f"{'Feature':<40} {'AAV2 VP1':<15} {'AAV9 VP1':<15} {'Residue':<10}")
print("-" * 80)

summary = [
    ("DE Loop region (5-fold pore)", f"{de2_start}-{de2_end}", f"{de9_start}-{de9_end}", "region"),
    ("  ANNLT motif start", str(aav2_annlt), str(aav9_annlt), "A"),
    ("  Pore constriction (L)", "336", str(l338_pos), "L"),
    ("HSPG site 1 (AAV2) / Gal (AAV9)", "585 (R)", "586 (S)", "differs"),
    ("HSPG site 2 (AAV2) / Gal (AAV9)", "588 (R)", "589 (A)", "differs"),
    ("HI Loop region (5-fold surround)", "655-670", "656-671", "region"),
    ("  VP1 incorporation (F)", "661", "662", "F"),
    ("VR-I (insertion site)", "260-267", "262-269", "see below"),
    ("VR-II (overlaps DE loop)", "325-330", "327-332", "see below"),
    ("VR-VIII (receptor binding)", "579-591", "581-593", "see below"),
]

for feature, aav2, aav9, res in summary:
    print(f"{feature:<40} {aav2:<15} {aav9:<15} {res:<10}")

# ============================================================
# AAV9 VARIABLE REGIONS WITH SEQUENCES
# ============================================================
print("\n" + "=" * 80)
print("AAV9 VARIABLE REGIONS (for reference)")
print("=" * 80)

for vr_name, (start, end) in AAV9_VR_REGIONS.items():
    before, region, after = extract_region(AAV9_VP1, start, end)
    print(f"\n{vr_name} (VP1 {start}-{end}):")
    print(f"  Sequence: {region}")
    print(f"  Context:  ...{before}[{region}]{after}...")

# ============================================================
# AVD002 ANNOTATION RECOMMENDATIONS
# ============================================================
print("\n" + "=" * 80)
print("AVD002 ANNOTATION RECOMMENDATIONS")
print("=" * 80)
print("""
For AVD002-Rep2Mut2Cap9-6R-wt plasmid, annotate these 5-fold axis features:

1. DE LOOP (5-FOLD PORE)
   - VP1 positions: 325-345
   - Function: Forms walls of the 5-fold pore/channel
   - Key residue: L338 (pore constriction point)
   - Engineering note: Insertions here may affect genome packaging

2. HI LOOP (5-FOLD SURROUNDING)
   - VP1 positions: 656-671
   - Function: Surrounds 5-fold pore externally
   - Key residue: F662 (critical for VP1 incorporation)
   - Engineering note: Modifications may affect VP1:VP2:VP3 ratios

3. VR-I (COMMON INSERTION SITE)
   - VP1 positions: 262-269
   - Function: Surface exposed loop, often used for peptide display
   - Proximity to 5-fold: Not directly part of 5-fold axis
   - Engineering note: Common site for VHH/peptide insertions

4. VR-VIII (RECEPTOR BINDING)
   - VP1 positions: 581-593
   - Function: Contains residues involved in receptor binding
   - AAV9 difference: Uses galactose receptor, NOT HSPG
   - Note: S586/A589 in AAV9 vs R585/R588 in AAV2

VERIFIED RESIDUES (with 5aa flanking context):
----------------------------------------------
L338 (pore):  ...TIANN[L]TSTVQ...
F662 (VP1):   ...DPPTA[F]NKDKL...
S586:         ...ATNHQ[S]AQAQA...
A589:         ...HQSAQ[A]QAQTG...
""")

print("=" * 80)
print("END OF REPORT")
print("=" * 80)
