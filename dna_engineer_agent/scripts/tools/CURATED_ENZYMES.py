#!/usr/bin/env python3
"""
Curated list of reliable, non-fussy restriction enzymes

Criteria for inclusion:
1. Commonly used workhorse enzymes
2. Stable at standard storage conditions (-20°C)
3. Work well at 37°C with minimal star activity
4. 6-8 bp recognition sites (avoid 4-5 bp - too frequent)
5. Minimal ambiguity codes (avoid excessive N's)
6. Single-site cutters (no double-cutting like BaeI)
7. Available from all major suppliers
8. Well-characterized cutting behavior

Excluded:
- Nicking enzymes (Nb.*, Nt.*)
- Mega-enzymes (I-SceI, I-CeuI, etc.)
- Highly ambiguous sites (>4 N's)
- Temperature-sensitive enzymes
- Rare or uncommon enzymes
"""

RELIABLE_ENZYMES = [
    # Standard 6-bp cutters (most common)
    ("AfeI", "AGCGCT"),
    ("AgeI", "ACCGGT"),
    ("ApaI", "GGGCCC"),
    ("ApaLI", "GTGCAC"),
    ("AvrII", "CCTAGG"),
    ("BamHI", "GGATCC"),
    ("BglII", "AGATCT"),
    ("BmtI", "GCTAGC"),      # Same as NheI
    ("BspEI", "TCCGGA"),
    ("BsrGI", "TGTACA"),
    ("BstBI", "TTCGAA"),
    ("ClaI", "ATCGAT"),
    ("EagI", "CGGCCG"),
    ("EcoRI", "GAATTC"),
    ("EcoRV", "GATATC"),
    ("HindIII", "AAGCTT"),
    ("HpaI", "GTTAAC"),
    ("KpnI", "GGTACC"),
    ("MfeI", "CAATTG"),
    ("MluI", "ACGCGT"),
    ("MscI", "TGGCCA"),
    ("NcoI", "CCATGG"),
    ("NdeI", "CATATG"),
    ("NheI", "GCTAGC"),      # Same as BmtI
    ("NruI", "TCGCGA"),
    ("NsiI", "ATGCAT"),
    ("PciI", "ACATGT"),
    ("PmlI", "CACGTG"),
    ("PstI", "CTGCAG"),
    ("PvuI", "CGATCG"),
    ("PvuII", "CAGCTG"),
    ("SacI", "GAGCTC"),
    ("SacII", "CCGCGG"),
    ("SalI", "GTCGAC"),
    ("ScaI", "AGTACT"),
    ("SmaI", "CCCGGG"),
    ("SnaBI", "TACGTA"),
    ("SpeI", "ACTAGT"),
    ("SphI", "GCATGC"),
    ("StuI", "AGGCCT"),
    ("XbaI", "TCTAGA"),
    ("XhoI", "CTCGAG"),
    ("XmaI", "CCCGGG"),      # Same as SmaI

    # 8-bp cutters (very specific, low star activity)
    ("AscI", "GGCGCGCC"),
    ("AsiSI", "GCGATCGC"),
    ("FseI", "GGCCGGCC"),    # Can be fussy but very useful
    ("NotI", "GCGGCCGC"),
    ("PacI", "TTAATTAA"),
    ("PmeI", "GTTTAAAC"),
    ("SbfI", "CCTGCAGG"),
    ("SrfI", "GCCCGGGC"),
    ("SwaI", "ATTTAAAT"),

    # Useful 7-bp cutters
    ("AflII", "CTTAAG"),
    ("BstZ17I", "GTATAC"),

    # A few common Type IIS enzymes (cut outside recognition site)
    # These are very commonly used despite being Type IIS
    ("BsaI", "GGTCTC"),
    ("BsmBI", "CGTCTC"),
    ("Esp3I", "CGTCTC"),    # Same as BsmBI

    # Some with minimal ambiguity (1-2 positions)
    ("AflIII", "ACRYGT"),   # R = A/G, Y = C/T
    ("BsaWI", "WCCGGW"),    # W = A/T
    ("HincII", "GTYRAC"),   # Y = C/T, R = A/G
]

# Alternative: If you want a super-conservative list, use only these most common ones:
SUPER_CONSERVATIVE = [
    ("AgeI", "ACCGGT"),
    ("ApaI", "GGGCCC"),
    ("AvrII", "CCTAGG"),
    ("BamHI", "GGATCC"),
    ("BglII", "AGATCT"),
    ("EcoRI", "GAATTC"),
    ("EcoRV", "GATATC"),
    ("HindIII", "AAGCTT"),
    ("KpnI", "GGTACC"),
    ("MluI", "ACGCGT"),
    ("NcoI", "CCATGG"),
    ("NdeI", "CATATG"),
    ("NheI", "GCTAGC"),
    ("NotI", "GCGGCCGC"),
    ("PstI", "CTGCAG"),
    ("SacI", "GAGCTC"),
    ("SalI", "GTCGAC"),
    ("SmaI", "CCCGGG"),
    ("SpeI", "ACTAGT"),
    ("XbaI", "TCTAGA"),
    ("XhoI", "CTCGAG"),
]

if __name__ == "__main__":
    print(f"Reliable enzymes: {len(RELIABLE_ENZYMES)}")
    print(f"Super conservative: {len(SUPER_CONSERVATIVE)}")

    # Show which were removed from original list
    from silent_sites import RESTRICTION_ENZYMES
    original_names = {name for name, _ in RESTRICTION_ENZYMES}
    curated_names = {name for name, _ in RELIABLE_ENZYMES}
    removed = original_names - curated_names

    print(f"\nRemoved {len(removed)} enzymes from original list of {len(RESTRICTION_ENZYMES)}")
    print("\nExamples of removed enzymes:")
    for i, enzyme in enumerate(sorted(removed)[:20]):
        print(f"  {enzyme}")
    if len(removed) > 20:
        print(f"  ... and {len(removed) - 20} more")
