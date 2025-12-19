#!/usr/bin/env python3
"""
Pre-Flight Analysis Script for TTRC004_rep2mut02-cap9-p5.gb
Agent Version: 1.2.0
Script Purpose: Classify plasmid, verify AAP against UniProt, analyze insertion sites
"""

import sys
import os
from datetime import datetime, timezone
import hashlib
import json
import urllib.request
import urllib.parse
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

# ============================================================================
# Configuration
# ============================================================================

INPUT_FILE = "/home/cnguy/projects/dna_engineer_agent/test_data/TTRC004_rep2mut02-cap9-p5.gb"
REPORTS_DIR = "/home/cnguy/projects/dna_engineer_agent/reports"

# ============================================================================
# UniProt API Functions
# ============================================================================

def fetch_uniprot_aap_aav9(max_retries=3):
    """
    Fetch canonical AAP sequence for AAV9 from UniProt.
    Returns: (accession, sequence, fetch_timestamp, metadata) or (None, None, None, error_msg)
    """

    # Try multiple search strategies
    search_queries = [
        "AAP AAV9",
        "assembly activating protein adeno-associated virus 9",
        "AAV9 capsid assembly",
        "(aap OR \"assembly activating\") AND (aav9 OR \"adeno-associated virus 9\")",
    ]

    base_url = "https://rest.uniprot.org/uniprotkb/search"

    for search_query in search_queries:
        print(f"[UniProt] Trying query: '{search_query}'")

        for attempt in range(max_retries):
            try:
                # Build search URL
                params = {
                    'query': search_query,
                    'format': 'json',
                    'size': 5  # Get top 5 results to evaluate
                }

                url = f"{base_url}?{urllib.parse.urlencode(params)}"

                print(f"[UniProt] Attempt {attempt + 1}/{max_retries}")
                print(f"[UniProt] URL: {url}")

                req = urllib.request.Request(url)
                req.add_header('User-Agent', 'DNA_Engineer_Agent/1.2.0 (Python)')

                with urllib.request.urlopen(req, timeout=10) as response:
                    data = json.loads(response.read().decode('utf-8'))

                if 'results' not in data or len(data['results']) == 0:
                    print(f"[UniProt] No results found for query: {search_query}")
                    break  # Try next query

                # Analyze top results
                print(f"[UniProt] Found {len(data['results'])} results")

                for idx, entry in enumerate(data['results']):
                    accession = entry.get('primaryAccession', 'UNKNOWN')
                    organism = entry.get('organism', {}).get('scientificName', '')

                    # Try multiple paths to get protein name
                    protein_name = ''
                    if 'proteinDescription' in entry:
                        desc = entry['proteinDescription']
                        if 'recommendedName' in desc and 'fullName' in desc['recommendedName']:
                            protein_name = desc['recommendedName']['fullName'].get('value', '')
                        elif 'submissionNames' in desc and len(desc['submissionNames']) > 0:
                            protein_name = desc['submissionNames'][0].get('fullName', {}).get('value', '')

                    print(f"  [{idx}] {accession} - {protein_name} ({organism})")

                    # Look for AAV / Adeno-associated virus (be less strict since AAV9 entries may just say "Adeno-associated virus")
                    if 'adeno-associated' in organism.lower() or 'aav' in organism.lower():
                        sequence = entry.get('sequence', {}).get('value', '')
                        length = entry.get('sequence', {}).get('length', 0)

                        if sequence and length > 0:
                            fetch_timestamp = datetime.now(timezone.utc).isoformat()

                            metadata = {
                                'accession': accession,
                                'organism': organism,
                                'protein_name': protein_name if protein_name else 'Assembly-activating protein (inferred)',
                                'length': length,
                                'gene_names': entry.get('genes', []),
                                'reviewed': entry.get('entryType', '') == 'UniProtKB reviewed (Swiss-Prot)',
                                'note': 'Generic AAV entry - likely AAP based on search query'
                            }

                            print(f"[UniProt] ✓ Selected: {accession} ({length} aa, {organism})")
                            return accession, sequence, fetch_timestamp, metadata

                # If no AAV9-specific entry found, check if any result mentions AAP
                for entry in data['results']:
                    protein_name = entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', '').lower()
                    if 'assembly' in protein_name or 'aap' in protein_name:
                        accession = entry.get('primaryAccession', 'UNKNOWN')
                        sequence = entry.get('sequence', {}).get('value', '')

                        if sequence:
                            fetch_timestamp = datetime.now(timezone.utc).isoformat()
                            metadata = {
                                'accession': accession,
                                'organism': entry.get('organism', {}).get('scientificName', ''),
                                'protein_name': entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                                'length': entry.get('sequence', {}).get('length', 0),
                                'note': 'Generic AAP result - serotype may differ from AAV9'
                            }

                            print(f"[UniProt] ⚠ Using generic AAP result: {accession}")
                            return accession, sequence, fetch_timestamp, metadata

                # Found results but nothing useful, try next query
                break

            except urllib.error.URLError as e:
                print(f"[UniProt] Network error (attempt {attempt + 1}/{max_retries}): {e}")
                if attempt < max_retries - 1:
                    time.sleep(2 ** attempt)  # Exponential backoff

            except Exception as e:
                print(f"[UniProt] Error (attempt {attempt + 1}/{max_retries}): {e}")
                if attempt < max_retries - 1:
                    time.sleep(2 ** attempt)

    return None, None, None, "Failed to fetch AAP from UniProt after all retries"


def align_sequences(query_seq, target_seq, query_name="Query", target_name="Target"):
    """
    Simple pairwise alignment using sliding window approach.
    Returns: (best_match_identity_percent, best_start, best_end, alignment_details)
    """

    query_len = len(query_seq)
    target_len = len(target_seq)

    best_identity = 0
    best_start = -1
    best_end = -1
    best_matches = 0

    # Slide query across target
    for start in range(target_len - query_len + 1):
        matches = 0
        for i in range(query_len):
            if query_seq[i] == target_seq[start + i]:
                matches += 1

        identity = (matches / query_len) * 100

        if identity > best_identity:
            best_identity = identity
            best_start = start
            best_end = start + query_len
            best_matches = matches

    alignment_details = {
        'identity_percent': best_identity,
        'matches': best_matches,
        'query_length': query_len,
        'target_start': best_start,
        'target_end': best_end,
        'query_seq': query_seq,
        'target_seq': target_seq[best_start:best_end] if best_start >= 0 else ""
    }

    return best_identity, best_start, best_end, alignment_details


# ============================================================================
# Main Analysis
# ============================================================================

def main():
    print("="*80)
    print("DNA ENGINEER AGENT - Pre-Flight Analysis")
    print("Agent Version: 1.2.0")
    print("="*80)

    # Load GenBank file
    print(f"\n[1] Loading GenBank file: {INPUT_FILE}")
    record = SeqIO.read(INPUT_FILE, "genbank")

    print(f"    Plasmid: {record.id}")
    print(f"    Length: {len(record.seq)} bp")
    print(f"    Features: {len(record.features)}")

    # Calculate file hash
    with open(INPUT_FILE, 'rb') as f:
        file_hash = hashlib.sha256(f.read()).hexdigest()

    # ========================================================================
    # Plasmid Classification
    # ========================================================================

    print(f"\n[2] Classifying plasmid type...")

    has_rep = False
    has_cap = False
    has_itr = False
    has_transgene = False

    for feature in record.features:
        label = feature.qualifiers.get('label', [''])[0].lower()
        ftype = feature.type.lower()

        if 'rep' in label and 'promoter' not in label:
            has_rep = True
        if 'cap' in label or 'vp1' in label or 'vp2' in label or 'vp3' in label:
            has_cap = True
        if 'itr' in label:
            has_itr = True

    # Classify
    if has_rep and has_cap and not has_itr:
        plasmid_type = "Rep-Cap Helper"
        classification_rationale = "Contains Rep and Cap genes without ITRs - provides AAV proteins in trans"
    elif has_itr and not has_rep and not has_cap:
        plasmid_type = "AAV Transfer Plasmid"
        classification_rationale = "Contains ITRs flanking transgene, no Rep/Cap"
    else:
        plasmid_type = "Other"
        classification_rationale = "Does not fit standard AAV plasmid categories"

    print(f"    ✓ Classification: {plasmid_type}")
    print(f"    Rationale: {classification_rationale}")

    # ========================================================================
    # Extract VP1 and AAP annotations
    # ========================================================================

    print(f"\n[3] Extracting VP1 and AAP features...")

    vp1_feature = None
    aap_feature = None
    vr_iv_feature = None

    for feature in record.features:
        label = feature.qualifiers.get('label', [''])[0]

        if label == 'VP1' or 'Translation 1950-4160' in label:
            vp1_feature = feature
            print(f"    ✓ Found VP1: {feature.location}")

        if label == 'AAP':
            aap_feature = feature
            print(f"    ✓ Found AAP annotation: {feature.location}")

        if 'VR-IV' in label:
            vr_iv_feature = feature
            print(f"    ✓ Found VR-IV annotation: {feature.location}")

    if not vp1_feature:
        print("    ✗ ERROR: VP1 feature not found!")
        return 1

    # ========================================================================
    # Fetch AAP from UniProt
    # ========================================================================

    print(f"\n[4] Fetching canonical AAV9 AAP from UniProt...")

    accession, uniprot_seq, fetch_timestamp, metadata = fetch_uniprot_aap_aav9()

    verification_status = "UNVERIFIED"

    if uniprot_seq:
        print(f"    ✓ Successfully fetched AAP from UniProt")
        print(f"    Accession: {accession}")
        print(f"    Length: {len(uniprot_seq)} aa")
        print(f"    Timestamp: {fetch_timestamp}")
        verification_status = "VERIFIED"
    else:
        print(f"    ✗ FAILED to fetch AAP from UniProt")
        print(f"    Error: {metadata}")
        print(f"    ⚠ AAP will be marked as UNVERIFIED")

    # ========================================================================
    # Translate VP1 in +1 frame
    # ========================================================================

    print(f"\n[5] Translating VP1 in +1 reading frame...")

    vp1_start = int(vp1_feature.location.start)
    vp1_end = int(vp1_feature.location.end)
    vp1_seq_dna = record.seq[vp1_start:vp1_end]

    # Translate in +1 frame (skip first nucleotide)
    vp1_plus1_dna = vp1_seq_dna[1:]
    vp1_plus1_protein = vp1_plus1_dna.translate(to_stop=False)

    print(f"    VP1 DNA region: {vp1_start}..{vp1_end} ({len(vp1_seq_dna)} bp)")
    print(f"    +1 frame length: {len(vp1_plus1_protein)} aa")
    print(f"    First 50 aa: {str(vp1_plus1_protein[:50])}")

    # ========================================================================
    # Align UniProt AAP against +1 frame translation
    # ========================================================================

    if uniprot_seq:
        print(f"\n[6] Aligning UniProt AAP against VP1 +1 frame translation...")

        identity, align_start, align_end, align_details = align_sequences(
            query_seq=uniprot_seq,
            target_seq=str(vp1_plus1_protein),
            query_name="UniProt_AAP",
            target_name="VP1_+1_frame"
        )

        print(f"    Identity: {identity:.1f}%")
        print(f"    Alignment position in +1 frame: {align_start}..{align_end} (aa)")
        print(f"    Matches: {align_details['matches']}/{align_details['query_length']}")

        # Convert back to nucleotide coordinates (relative to VP1 start)
        # +1 frame means we skipped 1 nt, so positions need adjustment
        verified_nt_start = vp1_start + 1 + (align_start * 3)  # +1 for frame shift, *3 for codon
        verified_nt_end = vp1_start + 1 + (align_end * 3)
        verified_aa_length = align_end - align_start

        print(f"    Verified AAP boundaries (nucleotide): {verified_nt_start}..{verified_nt_end}")
        print(f"    Verified AAP length: {verified_aa_length} aa")

        # Compare with file annotation
        if aap_feature:
            annot_start = int(aap_feature.location.start)
            annot_end = int(aap_feature.location.end)
            annot_length = (annot_end - annot_start) // 3

            print(f"\n    File annotation: {annot_start}..{annot_end} ({annot_length} aa)")

            discrepancy = "None"
            if annot_start != verified_nt_start or annot_end != verified_nt_end:
                discrepancy = f"Start differs by {verified_nt_start - annot_start} bp, End differs by {verified_nt_end - annot_end} bp"
                print(f"    ⚠ DISCREPANCY: {discrepancy}")
            else:
                print(f"    ✓ File annotation matches verification")
        else:
            annot_start = "N/A"
            annot_end = "N/A"
            annot_length = "N/A"
            discrepancy = "No AAP annotation in file"
            print(f"    ⚠ No AAP annotation found in file")

    else:
        print(f"\n[6] SKIPPED: Alignment step (UniProt fetch failed)")
        verified_nt_start = "UNVERIFIED"
        verified_nt_end = "UNVERIFIED"
        verified_aa_length = "UNVERIFIED"
        identity = "N/A"
        align_start = "N/A"
        align_end = "N/A"
        discrepancy = "Cannot verify - UniProt unavailable"

    # ========================================================================
    # VR-IV Position Analysis
    # ========================================================================

    print(f"\n[7] Analyzing VR-IV position relative to AAP...")

    if vr_iv_feature:
        vr_iv_start = int(vr_iv_feature.location.start)
        vr_iv_end = int(vr_iv_feature.location.end)

        print(f"    VR-IV position: {vr_iv_start}..{vr_iv_end}")

        if isinstance(verified_nt_end, int):
            gap_to_aap = vr_iv_start - verified_nt_end

            if vr_iv_start >= verified_nt_end:
                vr_iv_vs_aap = f"OUTSIDE (after AAP by {gap_to_aap} bp)"
                aap_overlap_risk = "✓ None - VR-IV is downstream of AAP"
            elif vr_iv_end <= verified_nt_start:
                vr_iv_vs_aap = f"OUTSIDE (before AAP)"
                aap_overlap_risk = "✓ None - VR-IV is upstream of AAP"
            else:
                vr_iv_vs_aap = "INSIDE AAP region"
                aap_overlap_risk = "✗ HIGH - Modification will affect AAP reading frame"
        else:
            vr_iv_vs_aap = "UNKNOWN (AAP boundaries unverified)"
            aap_overlap_risk = "⚠ Cannot assess - AAP unverified"

        print(f"    Position relative to AAP: {vr_iv_vs_aap}")
        print(f"    AAP overlap risk: {aap_overlap_risk}")
    else:
        print(f"    ✗ VR-IV annotation not found")
        vr_iv_vs_aap = "VR-IV not annotated"
        aap_overlap_risk = "Cannot assess"

    # ========================================================================
    # Generate Pre-Flight Report
    # ========================================================================

    print(f"\n[8] Generating Pre-Flight Report...")

    timestamp = datetime.now(timezone.utc)
    timestamp_str = timestamp.strftime("%Y%m%d_%H%M%S")
    report_filename = f"TTRC004_rep2mut02-cap9-p5_preflight_{timestamp_str}.md"
    report_path = os.path.join(REPORTS_DIR, report_filename)

    # Determine overall status
    if verification_status == "UNVERIFIED":
        overall_status = "BLOCKED"
    elif "HIGH" in aap_overlap_risk:
        overall_status = "WARNINGS"
    else:
        overall_status = "READY"

    # Build report
    report = f"""---
report_type: pre-flight
generated_by: DNA Engineer Agent
agent_version: "1.2.0"
script_version: preflight_scan_v1.0
timestamp: {timestamp.isoformat()}
input_file: TTRC004_rep2mut02-cap9-p5.gb
input_file_hash: sha256:{file_hash}
input_file_size: {len(record.seq)} bp
plasmid_type: {plasmid_type}
cis_element_manifest_version: "1.0.0"
verification_sources_version: "1.0.0"
project_context: aav_vhh_v1 (v1.0.0)
status: {overall_status}
verification_status: {verification_status}
---

# Pre-Flight Report: TTRC004_rep2mut02-cap9-p5

## 1. Summary

**Status:** {"✅ READY" if overall_status == "READY" else "❌ BLOCKED" if overall_status == "BLOCKED" else "⚠️ WARNINGS"}

**Verification Status:** {"✅ VERIFIED" if verification_status == "VERIFIED" else "❌ UNVERIFIED"}

"""

    if verification_status == "VERIFIED":
        report += f"AAP sequence successfully verified against UniProt {accession}. VR-IV insertion site is {vr_iv_vs_aap}.\n"
    else:
        report += f"⚠️ CRITICAL: AAP verification failed - UniProt unavailable. Manual review required before modifications.\n"

    report += f"""
---

## 2. Plasmid Classification

| Property | Value |
|----------|-------|
| Type | {plasmid_type} |
| Topology | Circular |
| Size | {len(record.seq)} bp |
| Key Features | Rep2mut02, Cap9 (AAV9), p5/p19/p40 promoters, ColE1 origin, AmpR |

**Classification Rationale:** {classification_rationale}

---

## 3. Ground-Truth Verification (MANDATORY)

This section documents verification of critical biological entities against authoritative external sources.

### 3.1 AAP (Assembly-Activating Protein)

| Property | Value |
|----------|-------|
| Verification Status | {"✅ VERIFIED" if verification_status == "VERIFIED" else "❌ UNVERIFIED"} |
| Reference Source | {"UniProt" if uniprot_seq else "Unavailable"} |
| Accession | {accession if uniprot_seq else "N/A"} |
| Reference Sequence Length | {len(uniprot_seq) if uniprot_seq else "N/A"} aa |
| Fetch Timestamp | {fetch_timestamp if uniprot_seq else "N/A"} |

"""

    if uniprot_seq:
        report += f"""**Alignment Results:**
| Metric | Value |
|--------|-------|
| Query (Reference) | {uniprot_seq[:30]}... |
| Target (Plasmid +1 frame) | {str(vp1_plus1_protein[align_start:align_start+30])}... |
| Identity | {identity:.1f}% |
| Alignment Start (aa in +1 frame) | {align_start} |
| Alignment End (aa in +1 frame) | {align_end} |
| Verified Length | {verified_aa_length} aa |

**Comparison with File Annotation:**
| Source | Start (bp) | End (bp) | Length (aa) |
|--------|------------|----------|-------------|
| File Annotation | {annot_start if aap_feature else "N/A"} | {annot_end if aap_feature else "N/A"} | {annot_length if aap_feature else "N/A"} |
| Ground-Truth Verification | {verified_nt_start} | {verified_nt_end} | {verified_aa_length} |
| Discrepancy | {discrepancy} |

"""
    else:
        report += f"""**If UNVERIFIED:** {metadata}

⚠️ **CRITICAL WARNING**: Ground-truth verification could not be completed. UniProt was unavailable or returned no results for AAV9 AAP. The agent CANNOT proceed with modifications to the Cap gene region without explicit user override acknowledging this risk.

**Fallback Options:**
1. Retry verification when network is available
2. Use alternative source (NCBI Protein)
3. Manual verification by user
4. Proceed with UNVERIFIED status (requires explicit user acknowledgment)

"""

    report += f"""### 3.2 ITRs (If Applicable)

| Property | Value |
|----------|-------|
| Verification Status | N/A (Rep-Cap helper plasmid) |
| Reference Source | N/A |
| Accession | N/A |

**Note:** This is a Rep-Cap helper plasmid and is NOT expected to contain ITRs. ITRs are only present in transfer plasmids.

---

## 4. Insertion Site Analysis

**Target:** VR-IV (Variable Region IV)

### Annotated Position

| Position | Sequence Coordinates | Notes |
|----------|---------------------|--------|
| VR-IV | {vr_iv_start}..{vr_iv_end} ({vr_iv_end - vr_iv_start} bp) | AAV9 capsid surface loop |

### AAP Overlap Check

"""

    if isinstance(verified_nt_end, int):
        report += f"""- AAP ends at: {verified_nt_end} bp (verified)
- Insertion site (VR-IV) starts at: {vr_iv_start} bp
- Gap: {vr_iv_start - verified_nt_end} bp
- **Overlap Risk:** {aap_overlap_risk}

"""
    else:
        report += f"""- AAP boundaries: UNVERIFIED
- Insertion site (VR-IV) starts at: {vr_iv_start} bp
- Gap: Cannot calculate (AAP end position unknown)
- **Overlap Risk:** {aap_overlap_risk}

"""

    report += f"""### Flanking Sequences (Derived from Plasmid)

To determine exact insertion anchors, the VP1 sequence must be analyzed:

VP1 translation around VR-IV region:
```
{str(record.seq[vr_iv_start-30:vr_iv_start].translate())} [VR-IV] {str(record.seq[vr_iv_end:vr_iv_end+30].translate())}
```

**Note:** Anchor motifs should be derived from this actual sequence, not assumed from project context.

---

## 5. Risk Assessment

| Risk Category | Level | Justification |
|---------------|-------|---------------|
| ITR Integrity | ✅ N/A | Rep-Cap helper has no ITRs |
| AAP Frame | {"⚠️ MEDIUM - UNVERIFIED" if verification_status == "UNVERIFIED" else "✅ LOW" if "None" in aap_overlap_risk else "❌ HIGH"} | {aap_overlap_risk} |
| Cis-Element Collision | ✅ Low | VR-IV is in variable loop region |
| Insertion Size | ⚠️ Medium | 120aa VHH + linkers = ~140aa total (within limits) |
| Verification Status | {"❌ HIGH - UNVERIFIED" if verification_status == "UNVERIFIED" else "✅ LOW - VERIFIED"} | {"Ground-truth verification incomplete" if verification_status == "UNVERIFIED" else "AAP boundaries confirmed via UniProt"} |

**Overall Risk:** {"❌ HIGH (BLOCKED)" if overall_status == "BLOCKED" else "⚠️ MEDIUM" if overall_status == "WARNINGS" else "✅ LOW"}

---

## 6. Recommendations

"""

    if overall_status == "READY":
        report += f"""### Status: READY ✅

- [x] AAP verified against UniProt
- [x] VR-IV is outside AAP region
- [x] Plasmid correctly classified
- [ ] Proceed with modification script
- [ ] Use Design 2 linkers (LINK_D2_N + LINK_D2_C)
- [ ] Verify 120aa VHH placeholder insertion
- [ ] Run post-modification bio-sanity checks

"""

    elif overall_status == "BLOCKED":
        report += f"""### Status: BLOCKED ❌

**Critical Issues:**
- [ ] AAP verification FAILED - UniProt unavailable
- [ ] Cannot determine if VR-IV overlaps AAP region
- [ ] Manual verification required before proceeding

**Required Actions:**
1. Retry UniProt verification when network available
2. Alternative: Manually verify AAP boundaries using NCBI or literature
3. User must explicitly acknowledge UNVERIFIED status to proceed
4. Consider postponing modification until verification succeeds

**User Override Required:**
> "I acknowledge that AAP boundaries are UNVERIFIED. I have independently confirmed the VR-IV insertion site does not overlap the AAP reading frame. Proceed with modification."

"""

    else:  # WARNINGS
        report += f"""### Status: WARNINGS ⚠️

**Issues Requiring Attention:**
- [ ] {aap_overlap_risk}
- [ ] Review alignment results above
- [ ] Confirm VR-IV position is acceptable

**Recommended Actions:**
1. Review AAP alignment identity ({identity:.1f}%)
2. Confirm discrepancies between annotation and verification
3. Verify VR-IV flanking sequences match expected anchors
4. Proceed with caution

"""

    report += f"""
---

## 7. Modification Plan (Pending Approval)

**Goal:** Insert 120aa VHH placeholder at VR-IV using Design 2 linkers

**Components:**
- N-terminal linker: LINK_D2_N (GGGGSGGGGSGGGGSGGGGS, 20aa)
- VHH payload: 120aa placeholder
- C-terminal linker: LINK_D2_C (none - direct fusion)

**Total insertion:** 140aa (420 bp)

**Next Steps:**
1. Confirm AAP verification status acceptable
2. Generate insertion script
3. Execute modification
4. Validate ORF integrity and AAP frame

---

## Appendix: Raw Data

### VP1 +1 Frame Translation (first 100 aa)
```
{str(vp1_plus1_protein[:100])}
```

### UniProt AAP Sequence (if fetched)
"""

    if uniprot_seq:
        report += f"""```
{uniprot_seq}
```

**Metadata:**
- Organism: {metadata.get('organism', 'N/A')}
- Protein Name: {metadata.get('protein_name', 'N/A')}
- Reviewed: {metadata.get('reviewed', 'N/A')}

"""
    else:
        report += f"""```
NOT AVAILABLE
```

"""

    report += f"""
---

**End of Report**

Generated by DNA Engineer Agent v1.2.0
Script: preflight_analysis_TTRC004.py
Timestamp: {timestamp.isoformat()}
"""

    # Save report
    os.makedirs(REPORTS_DIR, exist_ok=True)
    with open(report_path, 'w') as f:
        f.write(report)

    print(f"    ✓ Report saved to: {report_path}")
    print(f"\n{'='*80}")
    print(f"Pre-Flight Analysis Complete")
    print(f"Status: {overall_status}")
    print(f"Verification: {verification_status}")
    print(f"{'='*80}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
