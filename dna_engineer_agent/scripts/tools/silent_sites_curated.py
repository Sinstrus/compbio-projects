#!/usr/bin/env python3
"""
Silent restriction site finder with CURATED enzyme list
Wrapper around silent_sites.py that only uses reliable, non-fussy enzymes
"""

# Import everything from silent_sites
from silent_sites import *

# Override RESTRICTION_ENZYMES with curated list
RESTRICTION_ENZYMES = [
    # Standard 6-bp cutters (most common, reliable enzymes)
    ("AfeI", "AGCGCT"),
    ("AflII", "CTTAAG"),
    ("AgeI", "ACCGGT"),
    ("ApaI", "GGGCCC"),
    ("ApaLI", "GTGCAC"),
    ("AvrII", "CCTAGG"),
    ("BamHI", "GGATCC"),
    ("BglII", "AGATCT"),
    ("BmtI", "GCTAGC"),      # Same as NheI
    ("BbvCI", "CCTCAGC"),    # Type IIS, for Golden Gate
    ("BsaI", "GGTCTC"),      # Type IIS, commonly used
    ("BsmBI", "CGTCTC"),     # Type IIS
    ("BspEI", "TCCGGA"),
    ("BsrGI", "TGTACA"),
    ("BstBI", "TTCGAA"),
    ("BstZ17I", "GTATAC"),
    ("ClaI", "ATCGAT"),
    ("EagI", "CGGCCG"),
    ("EcoRI", "GAATTC"),
    ("EcoRV", "GATATC"),
    ("Esp3I", "CGTCTC"),     # Same as BsmBI, Type IIS
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
    ("FseI", "GGCCGGCC"),
    ("NotI", "GCGGCCGC"),
    ("PacI", "TTAATTAA"),
    ("PmeI", "GTTTAAAC"),
    ("SbfI", "CCTGCAGG"),
    ("SrfI", "GCCCGGGC"),
    ("SwaI", "ATTTAAAT"),

    # Minimal ambiguity (1-2 variable positions, commonly used)
    ("AflIII", "ACRYGT"),
    ("BsaWI", "WCCGGW"),
    ("HincII", "GTYRAC"),
]

# Override the find_candidates function to use our curated list
def find_candidates_curated(
    dna_seq: str,
    protein_seq: str,
    max_mutations: int,
    min_length: int,
    roi_seq: Optional[str] = None
) -> List[Candidate]:
    """
    Find all candidate restriction sites using CURATED enzyme list

    This version only searches reliable, non-fussy enzymes
    """
    global RESTRICTION_ENZYMES

    # Temporarily replace the enzyme list
    import silent_sites
    original_enzymes = silent_sites.RESTRICTION_ENZYMES
    silent_sites.RESTRICTION_ENZYMES = RESTRICTION_ENZYMES

    # Call original function
    result = silent_sites.find_candidates(
        dna_seq=dna_seq,
        protein_seq=protein_seq,
        max_mutations=max_mutations,
        min_length=min_length,
        roi_seq=roi_seq
    )

    # Restore original list
    silent_sites.RESTRICTION_ENZYMES = original_enzymes

    return result

if __name__ == "__main__":
    print(f"Curated enzyme list: {len(RESTRICTION_ENZYMES)} reliable enzymes")
    print("\nExcluded:")
    print("  - Nicking enzymes (Nb.*, Nt.*)")
    print("  - Mega-enzymes (I-SceI, etc.)")
    print("  - Highly ambiguous sites (>2 N's)")
    print("  - Double-cutters (BaeI, etc.)")
    print("  - Temperature-sensitive enzymes")
    print("\nIncluded:")
    print("  - Common 6-bp and 8-bp cutters")
    print("  - A few reliable Type IIS enzymes")
    print("  - Enzymes with minimal ambiguity")
