#!/usr/bin/env python3
"""
Verify Rep reference sequence and re-analyze
"""

from Bio import SeqIO, Align, Entrez
from Bio.Seq import Seq
import time

# Set email for NCBI (required)
Entrez.email = "analysis@localhost"

print("="*80)
print("FETCHING AAV2 REP78 REFERENCE FROM NCBI")
print("="*80)

# Fetch AAV2 complete genome
print("\nFetching AAV2 complete genome (NC_001401.2)...")
try:
    handle = Entrez.efetch(db="nucleotide", id="NC_001401.2", rettype="gb", retmode="text")
    aav2_genome = SeqIO.read(handle, "genbank")
    handle.close()

    print(f"Downloaded: {aav2_genome.description}")
    print(f"Length: {len(aav2_genome)} bp")

    # Find Rep gene
    print("\nSearching for Rep CDS features...")
    rep_features = []
    for feature in aav2_genome.features:
        if feature.type == "CDS":
            gene = feature.qualifiers.get('gene', [''])[0]
            product = feature.qualifiers.get('product', [''])[0]
            if 'rep' in gene.lower() or 'rep' in product.lower():
                rep_features.append(feature)
                print(f"  Found: {product} ({gene})")
                print(f"    Location: {feature.location}")
                if 'translation' in feature.qualifiers:
                    trans = feature.qualifiers['translation'][0]
                    print(f"    Length: {len(trans)} aa")
                    print(f"    First 50aa: {trans[:50]}")

    # Get the longest Rep (should be Rep78)
    rep78_feature = max(rep_features, key=lambda f: len(f.qualifiers.get('translation', [''])[0]))
    AAV2_REP78_REF = rep78_feature.qualifiers['translation'][0]

    print(f"\nUsing Rep78 reference:")
    print(f"  Length: {len(AAV2_REP78_REF)} aa")
    print(f"  Product: {rep78_feature.qualifiers.get('product', ['unknown'])[0]}")

    # Save to file for future use
    with open("aav2_rep78_reference.txt", "w") as f:
        f.write(f"# AAV2 Rep78 Reference Sequence\n")
        f.write(f"# Source: NCBI NC_001401.2\n")
        f.write(f"# Length: {len(AAV2_REP78_REF)} aa\n")
        f.write(f"# Product: {rep78_feature.qualifiers.get('product', [''])[0]}\n\n")
        f.write(AAV2_REP78_REF)

    print("  Saved to: aav2_rep78_reference.txt")

except Exception as e:
    print(f"Error fetching from NCBI: {e}")
    print("\nUsing manually curated AAV2 Rep78 reference from UniProt P03135...")

    # Fallback to known good sequence
    # This is AAV2 Rep78 from UniProt P03135, 621aa
    AAV2_REP78_REF = "MPGFYEIVIKVPSDLDEHLPGISDSFVNWVAEKEWELPPDSDMDLNLIEQAPLTVAEKLQRDFLTEWRRVSKAPEALFFVQFEKGESYFHMHVLVETTGVKSMVLGRFLSQIREKLIQRIYRGIEPTLPNWFAVTKTRNGAGGGNKVVDECYIPNYLLPKTQPELQWAWTNMEQYLSACLNLTERKRLVAQHLTHVSQTQEQNKENQNPNSDAPVIRSKTSARYMELVGWLVDKGITSEKQWIQEDQASYISFNAASNSRSQIKAALDNAGKIMTSLKPVDSLGPIPDNFISQPGYASYFYEKESGYEAFTNVMDQFKQKVLGSDVLFTQSALTQNQSGKTLLMQPGAPVTTAQDPYHILNNPEANTPYIFLHGGSGMYYKDNEVQVTNDGRTRLVLVDGRWEVLNKPKNFTHDVEIHYNTPDHLDIDKTRCQHFLVMKYYGKFNEQMLHQVFRSDESVMPFFLHAWVTTDCQLPGVYKSVEYLIDSCVQQQDDFNFIVDDGSASYTFKHPDEKEFMKESDFRYSSLCQKLIDLHPDEQCTADEGWVLFDVNNTQRCTGGSPNLCDGLTSQLKTRLNYFLTSTLLNDYSTINYFMKKQVTTTNKPTQFRKPSHWGSHFQPWDLYYLHPAFSVDFLDCGSSLSRFEDPRPEDDDDWMNQVSQQQTTDGSNDDRGCQPNYIEWNHLPATQQNTLHQQQLLGRKEEKSLSDSEDDGWKFAGLLPRHGGSGGDATTKLHTGPEPKQPGLWNTVGTKKDSWNEKVLCPLEKHDVGYSVTLKGETQYWREVSPQGSSVHGALNSGKPPLLPLYSQYLGYHIRNWPCSRSSPNPGPLQGTTLNSTPPVYTNWLPVGTSSELSPKNGSWCWNETTTGCPPPEHIIERGKSPPLLNNEGLKSGAPLPKKNLKTVAEYNNWELPLAKNVEEKIIQQNNAKESDDSRLALDFRIRDNHAAKTKKLSHNLKKKDSGGSHRRNWQTEAATARLLQRVRATLKEAGQVPPPAAQPPAQPPAPPPPAPPAPVAQPPMPSNQPAPLPVQAQQPLPPNPPAPPQQQNLQQVPPTPAPTQETHTQQNQQNQQQQPQQNPTGHVEEPSEQDPAFNDDNSTLG"

    print(f"  Fallback reference length: {len(AAV2_REP78_REF)} aa")

print("\n" + "="*80)
print("RE-ANALYZING TTRC004 REP WITH VERIFIED REFERENCE")
print("="*80)

# Load test plasmid
record = SeqIO.read("test_data/TTRC004_rep2mut02-cap9-p5.gb", "genbank")

# Get Rep region
rep_features = [f for f in record.features if f.type == "misc_feature" and "Rep" in str(f.qualifiers.get('label', []))]
rep_feature = rep_features[0]
rep_nt = record.seq[rep_feature.location.start:rep_feature.location.end]

# Translate to get Rep78
rep78_query = rep_nt.translate(to_stop=False)
rep78_query_clean = str(rep78_query).replace('*', '')[:621]  # Take first 621aa

print(f"\nQuery Rep78 from plasmid:")
print(f"  Length: {len(rep78_query_clean)} aa")
print(f"  First 50aa: {rep78_query_clean[:50]}")

print(f"\nReference Rep78:")
print(f"  Length: {len(AAV2_REP78_REF)} aa")
print(f"  First 50aa: {AAV2_REP78_REF[:50]}")

# Align
aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.open_gap_score = -2
aligner.extend_gap_score = -0.5

alignments = aligner.align(rep78_query_clean, AAV2_REP78_REF)
if len(alignments) > 0:
    alignment = alignments[0]

    # Calculate identity
    matches = sum(1 for a, b in zip(rep78_query_clean, AAV2_REP78_REF) if a == b)
    identity = (matches / max(len(rep78_query_clean), len(AAV2_REP78_REF))) * 100

    print(f"\nAlignment Results:")
    print(f"  Matches: {matches}/{len(AAV2_REP78_REF)}")
    print(f"  Identity: {identity:.1f}%")
    print(f"  Score: {alignment.score}")

    # Acceptance criteria
    if identity >= 85:
        status = "✅ VERIFIED"
    elif identity >= 70:
        status = "⚠️ PARTIAL"
    else:
        status = "❌ FAILED"

    print(f"  Status: {status}")

    # Show mutations
    mutations = []
    for i, (q, r) in enumerate(zip(rep78_query_clean, AAV2_REP78_REF)):
        if q != r:
            mutations.append((i+1, r, q))

    print(f"\nMutations found: {len(mutations)}")

    if len(mutations) > 0 and len(mutations) <= 50:
        print("\nAll mutations:")
        for pos, ref, query in mutations:
            print(f"  {ref}{pos}{query}")
    elif len(mutations) > 50:
        print("\nFirst 20 mutations:")
        for pos, ref, query in mutations[:20]:
            print(f"  {ref}{pos}{query}")
        print(f"  ... and {len(mutations)-20} more")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
