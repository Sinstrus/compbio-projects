from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

def annotate_vr_loops(input_file, output_file):
    # Load the GenBank file
    record = SeqIO.read(input_file, "genbank")

    # Locate the VP1 CDS feature to establish the reading frame
    vp1_feature = None
    for feature in record.features:
        if feature.type == "CDS":
            # Check for 'capsid protein' or specific gene qualifiers
            # Casting to str converts the qualifier list to a string for easy searching
            if "VP1" in str(feature.qualifiers) or "major coat protein VP1" in str(feature.qualifiers):
                vp1_feature = feature
                break
    
    if not vp1_feature:
        print("Error: Could not locate VP1 CDS feature automatically. Please ensure the input file has a CDS labeled 'VP1'.")
        return

    # FIX: Cast directly to int. Biopython 1.78+ removed the .position attribute
    vp1_start = int(vp1_feature.location.start)
    print(f"VP1 CDS Start detected at nucleotide: {vp1_start + 1}")

    # AAV2 Variable Region Definitions
    # Residue numbering based on VP1
    vr_loops = {
        "VR-I":   (262, 269),
        "VR-II":  (326, 332),
        "VR-III": (380, 389),
        "VR-IV":  (450, 461),
        "VR-V":   (487, 505),
        "VR-VI":  (527, 540),
        "VR-VII": (546, 557),
        "VR-VIII":(581, 594),
        "VR-IX":  (705, 714)
    }

    # Create new features for each loop
    new_features = []
    for vr_name, (aa_start, aa_end) in vr_loops.items():
        # Convert Amino Acid position to Nucleotide position
        # Formula: NT_index = VP1_Start_Index + (AA_position - 1) * 3
        nt_start = vp1_start + (aa_start - 1) * 3
        nt_end = vp1_start + (aa_end * 3) 
        
        location = FeatureLocation(nt_start, nt_end)
        feature = SeqFeature(location=location, type="misc_feature", qualifiers={
            "label": [vr_name],
            "note": [f"Variable Region {vr_name}"]
        })
        new_features.append(feature)

    # Add new features to the record
    record.features.extend(new_features)
    
    # Sort features by start location for clean viewing
    record.features.sort(key=lambda x: int(x.location.start))

    # Write the annotated file
    SeqIO.write(record, output_file, "genbank")
    print(f"Success! Annotated file saved as: {output_file}")

if __name__ == "__main__":
    try:
        annotate_vr_loops("aav2_genome.gb", "aav2_genome_annotated.gb")
    except FileNotFoundError:
        print("File 'aav2_genome.gb' not found. Please ensure it is in the same directory.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")