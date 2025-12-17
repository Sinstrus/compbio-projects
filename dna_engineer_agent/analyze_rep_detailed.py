#!/usr/bin/env python3
"""
Detailed Rep protein analysis for TTRC004
"""

from Bio import SeqIO, Align
from Bio.Seq import Seq

# Load the plasmid
record = SeqIO.read("test_data/TTRC004_rep2mut02-cap9-p5.gb", "genbank")

# AAV2 Rep78 reference (YP_680422.1)
AAV2_REP78_REF = """MPGFYEIVIKVPSDLDEHLPGISDSFVNWVAEKEWELPPDSDMDLNLIEQAPLTVAEKLQRDFLTEWRRVSKAPEALFFVQFEKGESYFHMHVLVETTGVKSMVLGRFLSQIREKLIQRIYRGIEPTLPNWFAVTKTRNGAGGGNKVVDECYIPNYLLPKTQPELQWAWTNMEQYLSACLNLTERKRLVAQHLTHVSQTQEQNKENQNPNSDAPVIRSKTSARYMELVGWLVDKGITSEKQWIQEDQASYISFNAASNSRSQIKAALDNAGKIMTSLKPVDSLGPIPDNFISQPGYASYFYEKESGYEAFTNVMDQFKQKVLGSDVLFTQSALTQNQSGKTLLMQPGAPVTTAQDPYHILNNPEANTPYIFLHGGSGMYYKDAEVQVTNDGRTRLVLVDGRWEVLNKPKNFTHDVEIHYNTPDHLDIDKTRCQHFLVMKYYGKFNEQMLHQVFRSDESVMPFFLHAWVTTDCQLPGVYKSVEYLIDSCVQQQDDFNFIVDDGSASYTFKHPDEKEFMKESDFRYSSLCQKLIDLHPDEQCTADEGWVLFDVNNTQRCTGGSPNLCDGLTSQLKTRLNYFLTSTLLNDYSTINYFMKKQVTTTNKPVPFRPNQPDTITQFFPKNNKEWWQGGLTVTAKQTDNCWTNRFLIEGEKNKFSVCKPLRVFKKHLQHVFENNEECQNGGVEQGIQQGIKDQKAKIELMEELQSASAAAPVKKASVKVEEVTEDMFDFQVGTQVDKSVCKDFTDWQFLKTTTVKGVTVEDLKNTFIINSSVPSVQRHSSWLSQADAAFNKVEAFLSESQLPGMVWQDRDVYLQPFIPHLNPAMSSLYSALNSLSKVADILAEVRSRSRAELLKQFYEETDISKTNFSLQSSRTPPAGSPQSSLDLSTFQHQNLAQEDKSTFQWDIKDDIRLPVGAQGLGERIVQVMNGFYKGLWGDSRLLQGPDSWLKPNNVGGGLGLGTFLGNGLAGVLGNKLGLQAESGGGLGAVDPVVLSNSVLVNGGGLQTQLVFNKGVAVVNASSGGLALAEALVEVNTPPGSGARAVSNLVPWYQVDVIAGNLALNLTSGGGLALGVEQVAVSGGSIGVAGGMLGKGSGTLALDAVEVGGGFGAVDVGSGALALGKSRVNKDDGSTQYWDINTGNSSTQYWDLNKQGGSTVSWDIKSGGSTVQWDIKNGGGSTVSWDLKNAGGKIPWDVKDQGGDIVQWDLKNAGGQIPWDLKNGGGNISWDLNNQGGQINWTLDKQGGTIKWDLKDQGGTINWDLSIKGGTKAKWDITDQGGTIKWDLKNQGGNIKWDVNDQGGKVKWDLKDQGGQIKWDLKNQGGNIKWDITDQGGNIKWDIKDQGGSIKWDLKNQGGKVKWDLKDQGGKIKWNIKDQGGKVKWNIKDQGG""".replace(" ", "").replace("\n", "")

print("="*80)
print("REP PROTEIN DETAILED ANALYSIS")
print("="*80)

# Get Rep region
rep_features = [f for f in record.features if f.type == "misc_feature" and "Rep" in str(f.qualifiers.get('label', []))]
if not rep_features:
    print("ERROR: No Rep feature found")
    exit(1)

rep_feature = rep_features[0]
rep_start = rep_feature.location.start
rep_end = rep_feature.location.end
rep_coords = f"{rep_start+1}..{rep_end}"

print(f"\nRep region: {rep_coords}")
print(f"Rep region length: {rep_end - rep_start} bp ({(rep_end - rep_start)/3:.1f} aa)")

# Extract nucleotide sequence
rep_nt_seq = record.seq[rep_start:rep_end]

# Translate in all 3 frames and find ORFs
print("\nSearching for ORFs in all 3 reading frames...")

all_orfs = []
for frame in range(3):
    frame_seq = rep_nt_seq[frame:]
    frame_seq = frame_seq[:len(frame_seq)//3*3]
    protein = frame_seq.translate(to_stop=False)

    # Find all ORFs (ATG to stop or end)
    import re
    for match in re.finditer(r'M[^*]*', str(protein)):
        orf_start_aa = match.start()
        orf_end_aa = match.end()
        orf_seq = match.group()
        orf_len = len(orf_seq.rstrip('*'))

        if orf_len >= 300:  # Only consider ORFs >= 300aa for Rep proteins
            nt_start = rep_start + frame + orf_start_aa * 3
            nt_end = rep_start + frame + orf_end_aa * 3
            all_orfs.append({
                'frame': frame,
                'orf_seq': orf_seq,
                'length': orf_len,
                'nt_start': nt_start + 1,  # 1-indexed
                'nt_end': nt_end,
                'aa_start_in_feature': orf_start_aa,
            })

print(f"Found {len(all_orfs)} ORFs >= 300aa")

# Sort by length (descending)
all_orfs.sort(key=lambda x: x['length'], reverse=True)

# Align each ORF to Rep78 reference
def align_sequences(seq1, seq2):
    """Simple alignment"""
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq1, seq2)
    if len(alignments) > 0:
        alignment = alignments[0]
        # Count matches
        matches = sum(1 for i, (a, b) in enumerate(zip(seq1, seq2)) if a == b)
        identity = (matches / max(len(seq1), len(seq2))) * 100
        return identity, alignment.score
    return 0, 0

print("\n" + "="*80)
print("REP78/68/52/40 IDENTIFICATION")
print("="*80)

# The longest ORF should be Rep78 (or Rep52 if p19 is used)
# In AAV2, Rep78 is 621aa, Rep52 is 397aa (starts at aa 225 of Rep78)

for i, orf in enumerate(all_orfs[:5]):  # Check top 5 ORFs
    orf_seq_clean = orf['orf_seq'].rstrip('*')
    identity, score = align_sequences(orf_seq_clean, AAV2_REP78_REF[:621])

    print(f"\nORF {i+1}:")
    print(f"  Frame: {orf['frame']}")
    print(f"  Coordinates: {orf['nt_start']}..{orf['nt_end']}")
    print(f"  Length: {orf['length']} aa")
    print(f"  Identity vs AAV2 Rep78 (621aa): {identity:.1f}%")

    # Determine which Rep protein this is
    if orf['length'] >= 600:
        # Likely Rep78 or Rep68
        if identity >= 70:
            print(f"  Candidate: Rep78 (unspliced)")
            # Check for Rep68 splice variant
            # Rep68 is 1-536 of Rep78
            rep68_identity, _ = align_sequences(orf_seq_clean[:536], AAV2_REP78_REF[:536])
            print(f"  Rep68 region (1-536) identity: {rep68_identity:.1f}%")
        else:
            print(f"  Candidate: Unknown (low identity)")

    elif 380 <= orf['length'] <= 420:
        # Likely Rep52 or Rep40
        # Rep52 corresponds to aa 225-621 of Rep78
        rep52_ref = AAV2_REP78_REF[224:621]  # 0-indexed
        identity_rep52, _ = align_sequences(orf_seq_clean, rep52_ref)
        print(f"  Identity vs Rep52 region (225-621): {identity_rep52:.1f}%")

        if identity_rep52 >= 70:
            print(f"  Candidate: Rep52 (unspliced, p19-driven)")
            # Check for Rep40 splice variant (225-536)
            rep40_ref = AAV2_REP78_REF[224:536]
            rep40_identity, _ = align_sequences(orf_seq_clean[:312], rep40_ref)
            print(f"  Rep40 region (225-536) identity: {rep40_identity:.1f}%")

    # Status
    if identity >= 85:
        print(f"  Status: ✅ VERIFIED")
    elif identity >= 70:
        print(f"  Status: ⚠️ PARTIAL")
    else:
        print(f"  Status: ❌ FAILED or not a Rep protein")

# Check for mutations
print("\n" + "="*80)
print("MUTATION ANALYSIS (Rep2mut02)")
print("="*80)

# The top ORF should be Rep78/52
if all_orfs:
    top_orf = all_orfs[0]
    top_orf_seq = top_orf['orf_seq'].rstrip('*')

    # Align to Rep78 reference and find differences
    print(f"\nComparing top ORF ({top_orf['length']}aa) to AAV2 Rep78 reference...")

    # Simple comparison for mutations
    if len(top_orf_seq) == 621:
        # Full Rep78
        ref_seq = AAV2_REP78_REF[:621]
        mutations = []
        for i, (query, ref) in enumerate(zip(top_orf_seq, ref_seq)):
            if query != ref:
                mutations.append((i+1, ref, query))

        print(f"Found {len(mutations)} amino acid differences")

        if len(mutations) <= 20:
            print("\nMutations:")
            for pos, ref, query in mutations:
                print(f"  Position {pos}: {ref} → {query}")
        else:
            print(f"\nShowing first 20 mutations:")
            for pos, ref, query in mutations[:20]:
                print(f"  Position {pos}: {ref} → {query}")

        # Calculate identity
        identity_pct = ((621 - len(mutations)) / 621) * 100
        print(f"\nOverall identity: {identity_pct:.1f}%")

    elif 380 <= len(top_orf_seq) <= 420:
        # Likely Rep52
        ref_seq = AAV2_REP78_REF[224:621]
        mutations = []
        min_len = min(len(top_orf_seq), len(ref_seq))
        for i in range(min_len):
            if top_orf_seq[i] != ref_seq[i]:
                mutations.append((i+1, ref_seq[i], top_orf_seq[i]))

        print(f"Found {len(mutations)} amino acid differences from Rep52 reference")

        if len(mutations) <= 20:
            print("\nMutations:")
            for pos, ref, query in mutations:
                print(f"  Position {pos} (of Rep52): {ref} → {query}")
                print(f"           (Position {pos+224} in Rep78 numbering)")
        else:
            print(f"\nShowing first 20 mutations:")
            for pos, ref, query in mutations[:20]:
                print(f"  Position {pos} (of Rep52): {ref} → {query}")

print("\n" + "="*80)
print("PROMOTER-DRIVEN EXPRESSION PREDICTION")
print("="*80)

# Check which promoters are present
promoters = [f for f in record.features if f.type == "regulatory" and "promoter" in str(f.qualifiers.get('regulatory_class', []))]

print("\nPromoters found:")
for prom in promoters:
    label = prom.qualifiers.get('label', ['unknown'])[0]
    coords = f"{prom.location.start+1}..{prom.location.end}"
    print(f"  {label}: {coords}")

print("\nExpected Rep expression:")
print("  p5 promoter → drives Rep78/68 (if present)")
print("  p19 promoter → drives Rep52/40 (if present)")
print("  p40 promoter → drives VP1/2/3 and AAP")

# Check spatial relationships
p5_present = any('p5' in str(f.qualifiers.get('label', [''])).lower() for f in promoters)
p19_present = any('p19' in str(f.qualifiers.get('label', [''])).lower() for f in promoters)

print("\nIn this construct:")
if p5_present:
    print("  ✅ p5 promoter present → should express Rep78/68")
else:
    print("  ⚠️ p5 promoter not found")

if p19_present:
    print("  ✅ p19 promoter present → should express Rep52/40")
else:
    print("  ⚠️ p19 promoter not found")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
