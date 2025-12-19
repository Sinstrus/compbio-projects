#!/usr/bin/env python3
"""
Analyze GenBank backbone files and generate a catalog.

Extracts:
1. Total length (bp)
2. Restriction sites in MCS/polylinker region (especially HindIII and XbaI)
3. ITR positions and type (for AAV backbones)
4. Selection marker (Amp or Kan)
5. Type IIS sites that would conflict with Golden Gate
"""

import json
import re
from pathlib import Path

# Restriction enzyme recognition sequences
RESTRICTION_SITES = {
    # Standard MCS enzymes
    'EcoRI': 'GAATTC',
    'HindIII': 'AAGCTT',
    'XhoI': 'CTCGAG',
    'NotI': 'GCGGCCGC',
    'AscI': 'GGCGCGCC',
    'SalI': 'GTCGAC',
    'AgeI': 'ACCGGT',
    'SphI': 'GCATGC',
    'PstI': 'CTGCAG',
    'EcoRV': 'GATATC',
    'XbaI': 'TCTAGA',
    'KpnI': 'GGTACC',
    'SpeI': 'ACTAGT',
    'NheI': 'GCTAGC',
    'SacI': 'GAGCTC',
    # Type IIS enzymes
    'BsaI': 'GGTCTC',
    'BsmBI': 'CGTCTC',
    'BbsI': 'GAAGAC',
    'Esp3I': 'CGTCTC',  # Same as BsmBI
    'SapI': 'GCTCTTC',
    'PaqCI': 'CACCTGC',
}

# Type IIS enzymes that conflict with Golden Gate
TYPE_IIS_ENZYMES = ['BsaI', 'BsmBI', 'BbsI', 'Esp3I', 'PaqCI', 'SapI']

# Standard cloning sites for pGS AAV backbones
STANDARD_CLONING_SITES = ['HindIII', 'XbaI']

def parse_genbank(filepath):
    """Parse GenBank file to extract sequence and features."""
    with open(filepath, 'r') as f:
        content = f.read()

    # Extract LOCUS line for length
    locus_match = re.search(r'LOCUS\s+\S+\s+(\d+)\s+bp', content)
    length = int(locus_match.group(1)) if locus_match else 0

    # Extract definition
    definition_match = re.search(r'DEFINITION\s+(.+?)(?:\n\S|\nFEATURES)', content, re.DOTALL)
    definition = definition_match.group(1).strip().replace('\n', ' ') if definition_match else ''

    # Extract sequence from ORIGIN section
    origin_match = re.search(r'ORIGIN\s*\n(.*?)\/\/', content, re.DOTALL)
    if origin_match:
        seq_lines = origin_match.group(1)
        # Remove line numbers and spaces
        sequence = re.sub(r'[\s\d]+', '', seq_lines).upper()
    else:
        sequence = ''

    # Extract features
    features = []
    features_section = re.search(r'FEATURES\s+Location/Qualifiers\n(.*?)ORIGIN', content, re.DOTALL)
    if features_section:
        feature_text = features_section.group(1)
        # Parse features (simplified parser)
        current_feature = None
        for line in feature_text.split('\n'):
            if line and not line.startswith('     '):
                continue
            if re.match(r'^\s{5}\w+\s+', line):
                # New feature
                if current_feature:
                    features.append(current_feature)
                feature_match = re.match(r'^\s{5}(\w+)\s+([\w\(\).,]+)', line)
                if feature_match:
                    current_feature = {
                        'type': feature_match.group(1),
                        'location': feature_match.group(2),
                        'qualifiers': {}
                    }
            elif current_feature and '=' in line:
                # Qualifier
                qual_match = re.match(r'^\s+/(\w+)="?([^"]*)"?', line)
                if qual_match:
                    key, value = qual_match.groups()
                    current_feature['qualifiers'][key] = value.strip('"')
        if current_feature:
            features.append(current_feature)

    return {
        'length': length,
        'definition': definition,
        'sequence': sequence,
        'features': features,
    }

def find_restriction_sites(sequence, enzyme_name):
    """Find all occurrences of a restriction site in the sequence."""
    pattern = RESTRICTION_SITES.get(enzyme_name)
    if not pattern:
        return []

    positions = []
    pos = 0
    while True:
        pos = sequence.find(pattern, pos)
        if pos == -1:
            break
        positions.append(pos + 1)  # 1-indexed
        pos += 1

    return positions

def analyze_backbone(gb_file):
    """Analyze a single GenBank backbone file."""
    data = parse_genbank(gb_file)

    result = {
        'file': gb_file.name,
        'description': data['definition'],
        'length': data['length'],
        'selection_marker': None,
        'itrs': {},
        'mcs_sites': {},
        'type_iis_sites': {},
        'cloning_method': None,
    }

    # Determine selection marker
    if 'Amp' in gb_file.name:
        result['selection_marker'] = 'Ampicillin'
    elif 'Kan' in gb_file.name:
        result['selection_marker'] = 'Kanamycin'

    # Extract ITR information from features
    for feature in data['features']:
        label = feature['qualifiers'].get('label', '')
        if 'ITR' in label:
            # Parse location
            loc_str = feature['location']
            if '..' in loc_str:
                start, end = loc_str.split('..')
                result['itrs'][label] = {
                    'start': int(start),
                    'end': int(end),
                    'length': int(end) - int(start) + 1,
                }

    # Find MCS restriction sites
    for enzyme_name in RESTRICTION_SITES.keys():
        if enzyme_name in TYPE_IIS_ENZYMES:
            continue  # Handle Type IIS separately
        sites = find_restriction_sites(data['sequence'], enzyme_name)
        if sites:
            result['mcs_sites'][enzyme_name] = sites

    # Find Type IIS sites
    for enzyme_name in TYPE_IIS_ENZYMES:
        sites = find_restriction_sites(data['sequence'], enzyme_name)
        if sites:
            result['type_iis_sites'][enzyme_name] = sites

    # Determine cloning method
    if 'FLASH' in gb_file.name or 'pUC57' in gb_file.name:
        result['cloning_method'] = 'EcoRV_blunt'
    elif 'AAV' in gb_file.name:
        result['cloning_method'] = 'HindIII_XbaI'

    # Check for Key IIS Free claim
    if 'Key IIS Free' in gb_file.name:
        result['claimed_iis_free'] = True
        # Verify the claim
        key_iis = ['BsaI', 'BsmBI', 'BbsI']
        has_key_iis = any(enzyme in result['type_iis_sites'] and result['type_iis_sites'][enzyme]
                          for enzyme in key_iis)
        result['verified_iis_free'] = not has_key_iis

    return result

def main():
    """Analyze all backbone files and generate catalog."""
    backbone_dir = Path(__file__).parent.parent.parent / 'backbones' / 'genscript'

    catalog = {
        'aav_backbones': [],
        'flash_backbones': [],
        'metadata': {
            'standard_cloning_sites': STANDARD_CLONING_SITES,
            'type_iis_enzymes_checked': TYPE_IIS_ENZYMES,
            'restriction_sites': RESTRICTION_SITES,
        }
    }

    # Analyze AAV backbones
    aav_dir = backbone_dir / 'aav'
    for gb_file in sorted(aav_dir.glob('*.gb')):
        print(f"Analyzing {gb_file.name}...")
        result = analyze_backbone(gb_file)
        catalog['aav_backbones'].append(result)

    # Analyze FLASH backbones
    flash_dir = backbone_dir / 'flash'
    for gb_file in sorted(flash_dir.glob('*.gb')):
        print(f"Analyzing {gb_file.name}...")
        result = analyze_backbone(gb_file)
        catalog['flash_backbones'].append(result)

    # Write catalog
    output_file = backbone_dir / 'BACKBONE_CATALOG.json'
    with open(output_file, 'w') as f:
        json.dump(catalog, f, indent=2)

    print(f"\nCatalog written to {output_file}")

    # Print summary
    print(f"\nSummary:")
    print(f"  AAV backbones: {len(catalog['aav_backbones'])}")
    print(f"  FLASH backbones: {len(catalog['flash_backbones'])}")

    # Check for Type IIS sites
    print(f"\nType IIS Site Analysis:")
    for backbone in catalog['aav_backbones'] + catalog['flash_backbones']:
        if backbone['type_iis_sites']:
            print(f"  {backbone['file']}: {', '.join(backbone['type_iis_sites'].keys())}")
        else:
            print(f"  {backbone['file']}: CLEAN (no Type IIS sites)")

    # Print MCS summary for AAV backbones
    print(f"\nAAV Backbone MCS Sites (HindIII and XbaI):")
    for backbone in catalog['aav_backbones']:
        hindiii_sites = backbone['mcs_sites'].get('HindIII', [])
        xbai_sites = backbone['mcs_sites'].get('XbaI', [])
        print(f"  {backbone['file']}:")
        print(f"    HindIII: {hindiii_sites}")
        print(f"    XbaI: {xbai_sites}")

if __name__ == '__main__':
    main()
