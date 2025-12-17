#!/usr/bin/env python3
"""
Debug VP1/VP2/VP3 nesting issue
"""

from Bio import SeqIO

# Load the plasmid
record = SeqIO.read("test_data/TTRC004_rep2mut02-cap9-p5.gb", "genbank")

# Get the features
vp1_features = [f for f in record.features if f.type == "CDS" and "VP1" in str(f.qualifiers.get('label', []))]
vp2_features = [f for f in record.features if f.type == "CDS" and "VP2" in str(f.qualifiers.get('label', []))]
vp3_features = [f for f in record.features if f.type == "CDS" and "VP3" in str(f.qualifiers.get('label', []))]

vp1_feature = vp1_features[0]
vp2_feature = vp2_features[0]
vp3_feature = vp3_features[0]

vp1_translation = vp1_feature.qualifiers.get('translation', [''])[0]
vp2_translation = vp2_feature.qualifiers.get('translation', [''])[0]
vp3_translation = vp3_feature.qualifiers.get('translation', [''])[0]

print("VP1 Details:")
print(f"  Coordinates: {vp1_feature.location.start+1}..{vp1_feature.location.end}")
print(f"  Length: {len(vp1_translation)} aa")
print(f"  First 50aa: {vp1_translation[:50]}")
print(f"  Last 50aa: {vp1_translation[-50:]}")

print("\nVP2 Details:")
print(f"  Coordinates: {vp2_feature.location.start+1}..{vp2_feature.location.end}")
print(f"  Length: {len(vp2_translation)} aa")
print(f"  First 50aa: {vp2_translation[:50]}")
print(f"  Last 50aa: {vp2_translation[-50:]}")

print("\nVP3 Details:")
print(f"  Coordinates: {vp3_feature.location.start+1}..{vp3_feature.location.end}")
print(f"  Length: {len(vp3_translation)} aa")
print(f"  First 50aa: {vp3_translation[:50]}")
print(f"  Last 50aa: {vp3_translation[-50:]}")

print("\n" + "="*80)
print("NESTING CHECK")
print("="*80)

# Check if they share the same C-terminus
print("\nC-terminus comparison:")
print(f"VP1 ends with: ...{vp1_translation[-30:]}")
print(f"VP2 ends with: ...{vp2_translation[-30:]}")
print(f"VP3 ends with: ...{vp3_translation[-30:]}")

if vp1_translation[-30:] == vp2_translation[-30:]:
    print("✅ VP1 and VP2 share C-terminus")
else:
    print("❌ VP1 and VP2 DO NOT share C-terminus")

if vp1_translation[-30:] == vp3_translation[-30:]:
    print("✅ VP1 and VP3 share C-terminus")
else:
    print("❌ VP1 and VP3 DO NOT share C-terminus")

# Check if VP2 and VP3 are substrings of VP1
print("\nSubstring check:")
if vp2_translation in vp1_translation:
    pos = vp1_translation.find(vp2_translation)
    print(f"✅ VP2 is a substring of VP1, starting at position {pos+1}")
else:
    print("❌ VP2 is NOT a substring of VP1")

if vp3_translation in vp1_translation:
    pos = vp1_translation.find(vp3_translation)
    print(f"✅ VP3 is a substring of VP1, starting at position {pos+1}")
else:
    print("❌ VP3 is NOT a substring of VP1")

# Check reading frame consistency
print("\n" + "="*80)
print("READING FRAME CHECK")
print("="*80)

vp1_start = vp1_feature.location.start
vp2_start = vp2_feature.location.start
vp3_start = vp3_feature.location.start

print(f"VP1 starts at nt {vp1_start+1}")
print(f"VP2 starts at nt {vp2_start+1}")
print(f"VP3 starts at nt {vp3_start+1}")

vp2_offset = (vp2_start - vp1_start) % 3
vp3_offset = (vp3_start - vp1_start) % 3

print(f"\nVP2 frame offset from VP1: {vp2_offset}")
print(f"VP3 frame offset from VP1: {vp3_offset}")

if vp2_offset == 0:
    print("✅ VP2 is in the same reading frame as VP1")
else:
    print(f"❌ VP2 is in frame +{vp2_offset} relative to VP1")

if vp3_offset == 0:
    print("✅ VP3 is in the same reading frame as VP1")
else:
    print(f"❌ VP3 is in frame +{vp3_offset} relative to VP1")

# Now extract and translate directly from nucleotide sequence to verify
print("\n" + "="*80)
print("DIRECT TRANSLATION VERIFICATION")
print("="*80)

vp1_nt = record.seq[vp1_start:vp1_feature.location.end]
vp2_nt = record.seq[vp2_start:vp2_feature.location.end]
vp3_nt = record.seq[vp3_start:vp3_feature.location.end]

vp1_direct = str(vp1_nt.translate())
vp2_direct = str(vp2_nt.translate())
vp3_direct = str(vp3_nt.translate())

print(f"VP1 direct translation matches annotation: {vp1_direct.rstrip('*') == vp1_translation}")
print(f"VP2 direct translation matches annotation: {vp2_direct.rstrip('*') == vp2_translation}")
print(f"VP3 direct translation matches annotation: {vp3_direct.rstrip('*') == vp3_translation}")

# Check if VP2 and VP3 translations are C-terminal fragments of VP1 direct translation
print(f"\nVP2 is C-terminal fragment of VP1: {vp1_direct.endswith(vp2_direct)}")
print(f"VP3 is C-terminal fragment of VP1: {vp1_direct.endswith(vp3_direct)}")

# Show the expected relationships for AAV9
print("\n" + "="*80)
print("EXPECTED STRUCTURE (AAV9)")
print("="*80)

AAV9_VP1_REF = """MAADGYLPDWLEDNLSEGIREWWALKPGAPQPKANQQHQDNARGLVLPGYKYLGPGNGLDKGEPVNAADAAALEHDKAYDQQLKAGDNPYLKYNHADAEFQERLKEDTSFGGNLGRAVFQAKKRLLEPLGLVEEAAKTAPGKKRPVEQSPQEPDSSAGIGKSGAQPAKKRLNFGQTGDTESVPDPQPIGEPPAAPSGVGSLTMASGGGAPVADNNEGADGVGSSSGNWHCDSQWLGDRVITTSTRTWALPTYNNHLYKQISNSTSGGSSNDNAYFGYSTPWGYFDFNRFHCHFSPRDWQRLINNNWGFRPKRLNFKLFNIQVKEVTDNNGVKTIANNLTSTVQVFTDSDYQLPYVLGSAHEGCLPPFPADVFMIPQYGYLTLNDGSQAVGRSSFYCLEYFPSQMLRTGNNFQFSYEFENVPFHSSYAHSQSLDRLMNPLIDQYLYYLSKTINGSGQNQQTLKFSVAGPSNMAVQGRNYIPGPSYRQQRVSTTVTQNNNSEFAWPGASSWALNGRNSLMNPGPAMASHKEGEDRFFPLSGSLIFGKQGTGRDNVDADKVMITNEEEIKTTNPVATESYGQVATNHQSAQAQAQTGWVQNQGILPGMVWQDRDVYLQGPIWAKIPHTDGNFHPSPLMGGFGMKHPPPQILIKNTPVPADPPTAFNKDKLNSFITQYSTGQVSVEIEWELQKENSKRWNPEIQYTSNYYKSNNVEFAVNTEGVYSEPRPIGTRYLTRNL"""

print(f"Reference VP1 length: {len(AAV9_VP1_REF)} aa")
print(f"VP2 should start ~138aa into VP1")
print(f"VP3 should start ~204aa into VP1")

expected_vp2 = AAV9_VP1_REF[137:]
expected_vp3 = AAV9_VP1_REF[203:]

print(f"\nExpected VP2 length: {len(expected_vp2)} aa")
print(f"Actual VP2 length: {len(vp2_translation)} aa")

print(f"\nExpected VP3 length: {len(expected_vp3)} aa")
print(f"Actual VP3 length: {len(vp3_translation)} aa")

# Check if the plasmid VP2/VP3 match the expected positions
vp2_calc_start = len(vp1_translation) - len(vp2_translation)
vp3_calc_start = len(vp1_translation) - len(vp3_translation)

print(f"\nCalculated VP2 start position in VP1: aa {vp2_calc_start + 1}")
print(f"Expected VP2 start position in VP1: aa 138")

print(f"\nCalculated VP3 start position in VP1: aa {vp3_calc_start + 1}")
print(f"Expected VP3 start position in VP1: aa 204")
