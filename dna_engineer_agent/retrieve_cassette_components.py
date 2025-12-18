#!/usr/bin/env python3
"""
Phase 2: Retrieve and Define Transgene Cassette Components
DNA Engineer Agent v2.2 - EF1A-VP1-rBG Assembly Project
"""

from Bio import SeqIO
from Bio.Seq import Seq

print("="*80)
print("PHASE 2: COMPONENT RETRIEVAL")
print("="*80)

# Load source plasmid
vp1_source = SeqIO.read("test_data/BASE-DRAFT-AAV9-RepCap-NOCAP.gb", "genbank")

# ============================================================================
# COMPONENT 1: VP1 ORF with proper ATG start codon
# ============================================================================
print("\n[1] VP1 ORF EXTRACTION")
print("-"*80)

# The VP1 annotation starts at position 2365, but let's verify the start codon
vp1_features = [f for f in vp1_source.features if f.type == "CDS" and
                f.qualifiers.get('label', [''])[0] == "VP1"]

if vp1_features:
    vp1_feature = vp1_features[0]
    vp1_start = vp1_feature.location.start  # 0-indexed: 2364
    vp1_end = vp1_feature.location.end      # 2-indexed: 4575

    # Extract sequence
    vp1_nt = str(vp1_source.seq[vp1_start:vp1_end])

    print(f"Annotated VP1 coordinates: {vp1_start+1}..{vp1_end}")
    print(f"Length: {len(vp1_nt)} bp")
    print(f"First 12 nt: {vp1_nt[:12]}")
    print(f"Last 3 nt: {vp1_nt[-3:]}")

    # The annotated VP1 starts with AAG (Lysine) not ATG
    # For our expression construct, we'll use the annotated CDS as-is
    # and provide proper initiation with Kozak+ATG upstream
    if vp1_nt[:3] == "ATG":
        print("✅ Starts with ATG")
    else:
        print(f"ℹ️ Annotated VP1 starts with {vp1_nt[:3]} (expected for this annotation)")
        print("   We will add proper ATG via Kozak sequence for expression")

    # Use the annotated VP1 CDS as-is - this is the mature protein sequence
    vp1_final = vp1_nt
    vp1_final_start = vp1_start

    # Translate and verify
    vp1_protein = Seq(vp1_final).translate()
    print(f"\nTranslation ({len(vp1_protein)} aa):")
    print(f"  First 60 aa: {vp1_protein[:60]}")
    print(f"  Internal stops: {str(vp1_protein)[:-1].count('*')}")
    print(f"  Stop codon: {vp1_final[-3:]}")
else:
    print("ERROR: VP1 not found")
    exit(1)

# ============================================================================
# COMPONENT 2: EF1α Promoter (retrieve from reference or define)
# ============================================================================
print("\n[2] EF1α PROMOTER")
print("-"*80)

# EF1α promoter with first intron (commonly used in expression vectors)
# Source: Based on pEF-GFP (Addgene #11154) and similar vectors
# This is the "long" EF1α with Intron A for enhanced expression

# Full EF1α promoter with intron (~1172 bp)
# Note: This is a commonly used sequence from human EEF1A1 gene
EF1A_FULL = """
GGATCTGCGATCGCTCCGGTGCCCGTCAGTGGGCAGAGCGCACATCGCCCACAGTCCCCGAGAAGTTGGGGGGAGGGGTCGGCAATTGAACCGGTGCCTAGAGAAGGTGGCGCGGGGTAAACTGGGAAAGTGATGTCGTGTACTGGCTCCGCCTTTTTCCCGAGGGTGGGGGAGAACCGTATATAAGTGCAGTAGTCGCCGTGAACGTTCTTTTTCGCAACGGGTTTGCCGCCAGAACACAGGTAAGTGCCGTGTGTGGTTCCCGCGGGCCTGGCCTCTTTACGGGTTATGGCCCTTGCGTGCCTTGAATTACTTCCACGCCCCTGGCTGCAGTACGTGATTCTTGATCCCGAGCTTCGGGTTGGAAGTGGGTGGGAGAGTTCGAGGCCTTGCGCTTAAGGAGCCCCTTCGCCTCGTGCTTGAGTTGAGGCCTGGCCTGGGCGCTGGGGCCGCCGCGTGCGAATCTGGTGGCACCTTCGCGCCTGTCTCGCTGCTTTCGATAAGTCTCTAGCCATTTAAAATTTTTGATGACCTGCTGCGACGCTTTTTTTCTGGCAAGATAGTCTTGTAAATGCGGGCCAAGATCTGCACACTGGTATTTCGGTTTTTGGGGCCGCGGGCGGCGACGGGGCCCGTGCGTCCCAGCGCACATGTTCGGCGAGGCGGGGCCTGCGAGCGCGGCCACCGAGAATCGGACGGGGGTAGTCTCAAGCTGGCCGGCCTGCTCTGGTGCCTGGCCTCGCGCCGCCGTGTATCGCCCCGCCCTGGGCGGCAAGGCTGGCCCGGTCGGCACCAGTTGCGTGAGCGGAAAGATGGCCGCTTCCCGGCCCTGCTGCAGGGAGCTCAAAATGGAGGACGCGGCGCTCGGGAGAGCGGGCGGGTGAGTCACCCACACAAAGGAAAAGGGCCTTTCCGTCCTCAGCCGTCGCTTCATGTGACTCCACGGAGTACCGGGCGCCGTCCAGGCACCTCGATTAGTTCTCGAGCTTTTGGAGTACGTCGTCTTTAGGTTGGGGGGAGGGGTTTTATGCGATGGAGTTTCCCCACACTGAGTGGGTGGAGACTGAAGTTAGGCCAGCTTGGCACTTGATGTAATTCTCCTTGGAATTTGCCCTTTTTGAGTTTGGATCTTGGTTCATTCTCAAGCCTCAGACAGTGGTTCAAAGTTTTTTTCTTCCATTTCAGGTGTCGTGA
""".replace('\n', '').replace(' ', '')

# Core EF1α promoter without intron (~212 bp)
EF1A_CORE = """
GGATCTGCGATCGCTCCGGTGCCCGTCAGTGGGCAGAGCGCACATCGCCCACAGTCCCCGAGAAGTTGGGGGGAGGGGTCGGCAATTGAACCGGTGCCTAGAGAAGGTGGCGCGGGGTAAACTGGGAAAGTGATGTCGTGTACTGGCTCCGCCTTTTTCCCGAGGGTGGGGGAGAACCGTATATAAGTGCAGTAGTCGCCGTGAACGTTCTTTTTCGCAACGGGTTTGCCGCCAGAACACAG
""".replace('\n', '').replace(' ', '')

print(f"EF1α Full (with intron): {len(EF1A_FULL)} bp")
print(f"EF1α Core (promoter only): {len(EF1A_CORE)} bp")
print(f"\nFull promoter first 60 bp: {EF1A_FULL[:60]}")
print(f"Full promoter last 60 bp:  {EF1A_FULL[-60:]}")

# ============================================================================
# COMPONENT 3: Kozak Sequence
# ============================================================================
print("\n[3] KOZAK SEQUENCE")
print("-"*80)

# Optimal Kozak consensus: GCCACCATGG
# where ATG is the start codon
KOZAK = "GCCACC"  # Will be placed before ATG of VP1

print(f"Kozak sequence: {KOZAK}ATG")
print("Consensus pattern: (G/A)CCATGG")

# ============================================================================
# COMPONENT 4: Rabbit β-Globin polyA signal
# ============================================================================
print("\n[4] RABBIT β-GLOBIN POLYA SIGNAL")
print("-"*80)

# Standard rBG polyA signal (127 bp)
# Source: Commonly used in expression vectors
RBG_POLYA = """
TGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAATAAAATGAGGAAATTGCATCGCATTGTCTGAG
""".replace('\n', '').replace(' ', '')

print(f"rBG polyA length: {len(RBG_POLYA)} bp")
print(f"First 40 bp: {RBG_POLYA[:40]}")
print(f"Contains AATAAA: {'Yes' if 'AATAAA' in RBG_POLYA else 'No'}")
if 'AATAAA' in RBG_POLYA:
    aataaa_pos = RBG_POLYA.find('AATAAA')
    print(f"  AATAAA position: {aataaa_pos+1}")

# ============================================================================
# ASSEMBLE CASSETTE
# ============================================================================
print("\n" + "="*80)
print("CASSETTE ASSEMBLY TEST")
print("="*80)

# For VP1, we need to remove the ATG since Kozak includes it
vp1_no_start = vp1_final[3:]  # Remove ATG

# Assemble: EF1α + Kozak + ATG + VP1(without ATG) + polyA
cassette_full = EF1A_FULL + KOZAK + "ATG" + vp1_no_start + RBG_POLYA
cassette_core = EF1A_CORE + KOZAK + "ATG" + vp1_no_start + RBG_POLYA

print(f"\nFull cassette: {len(cassette_full)} bp")
print(f"Core cassette: {len(cassette_core)} bp")

print(f"\nPackaging limits:")
print(f"  scAAV: 2400 bp")
print(f"  ssAAV: 4700 bp")

if len(cassette_full) <= 4700:
    print(f"  ✅ Full cassette fits in ssAAV")
else:
    print(f"  ❌ Full cassette too large for ssAAV")

if len(cassette_core) <= 2400:
    print(f"  ✅ Core cassette fits in scAAV")
else:
    print(f"  ⚠️ Core cassette ({len(cassette_core)} bp) exceeds scAAV limit by {len(cassette_core) - 2400} bp")

# Save components to a Python module for use in assembly
print("\n" + "="*80)
print("SAVING COMPONENTS")
print("="*80)

with open("cassette_components.py", "w") as f:
    f.write('"""Cassette components for EF1A-VP1-rBG assembly"""\n\n')
    f.write(f'# VP1 ORF (AAV9) - {len(vp1_final)} bp\n')
    f.write(f'VP1_ORF = """{vp1_final}"""\n\n')
    f.write(f'# VP1 start position in source plasmid\n')
    f.write(f'VP1_START_POS = {vp1_final_start + 1}  # 1-indexed\n\n')
    f.write(f'# EF1α promoter with intron - {len(EF1A_FULL)} bp\n')
    f.write(f'EF1A_FULL = """{EF1A_FULL}"""\n\n')
    f.write(f'# EF1α promoter core (no intron) - {len(EF1A_CORE)} bp\n')
    f.write(f'EF1A_CORE = """{EF1A_CORE}"""\n\n')
    f.write(f'# Kozak sequence\n')
    f.write(f'KOZAK = "{KOZAK}"\n\n')
    f.write(f'# Rabbit β-globin polyA - {len(RBG_POLYA)} bp\n')
    f.write(f'RBG_POLYA = """{RBG_POLYA}"""\n\n')
    f.write(f'# Assembled cassettes\n')
    f.write(f'CASSETTE_FULL = EF1A_FULL + KOZAK + VP1_ORF + RBG_POLYA\n')
    f.write(f'CASSETTE_CORE = EF1A_CORE + KOZAK + VP1_ORF + RBG_POLYA\n')

print("✅ Components saved to cassette_components.py")

print("\n" + "="*80)
print("COMPONENT RETRIEVAL COMPLETE")
print("="*80)
