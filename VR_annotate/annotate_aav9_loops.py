from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

def annotate_aav9_vr_loops(input_file, output_file):
    # Load the GenBank file
    try:
        record = SeqIO.read(input_file, "genbank")
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Locate the VP1 CDS feature to establish the reading frame
    vp1_feature = None
    for feature in record.features:
        if feature.type == "CDS":
            # Check for 'VP1' or 'capsid protein' identifiers
            # AAV9 GenBank files often label it "capsid protein VP1"
            if "VP1" in str(feature.qualifiers) or "major coat protein" in str(feature.qualifiers):
                vp1_feature = feature
                break
    
    if not vp1_feature:
        print("Error: Could not locate VP1 CDS feature automatically. Please ensure the input file has a CDS labeled 'VP1'.")
        return

    # FIX: Cast directly to int (resolves the AttributeError from previous versions)
    vp1_start = int(vp1_feature.location.start)
    print(f"AAV9 VP1 CDS Start detected at nucleotide: {vp1_start + 1}")

    # AAV9 Variable Region Definitions 
    # Coordinates taken from Table 1 of the provided Technical Report
    vr_loops = {
        "VR-I":   (262, 269), # Involved in 2/5-fold wall topography
        "VR-II":  (327, 332), # 5-fold pore constriction
        "VR-III": (382, 386), # Shoulder of protrusions
        "VR-IV":  (452, 460), # Major Engineering Hotspot (3-fold spike tip)
        "VR-V":   (488, 505), # Contains Galactose binding residues in AAV9
        "VR-VI":  (527, 539), # 3-fold base
        "VR-VII": (545, 558), # 2/5-fold wall
        "VR-VIII":(581, 593), # Critical functional domain; distinct from AAV2
        "VR-IX":  (704, 714)  # Near C-terminus
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
            "label": [f"AAV9_{vr_name}"],
            "note": [f"AAV9 Variable Region {vr_name} (AA {aa_start}-{aa_end})"]
        })
        new_features.append(feature)

    # Add new features to the record
    record.features.extend(new_features)
    
    # Sort features by start location for clean viewing
    record.features.sort(key=lambda x: int(x.location.start))

    # Write the annotated file
    SeqIO.write(record, output_file, "genbank")
    print(f"Success! Annotated AAV9 file saved as: {output_file}")

if __name__ == "__main__":
    # Ensure you have an AAV9 GenBank file named 'aav9_genome.gb'
    annotate_aav9_vr_loops("aav9_cap.gb", "aav9_cap_annotated.gb")