import csv
import gzip
import os
import sys

# --- Configuration ---
METADATA_FILE = "GSE221868_exon_skipping.filtered.psi.csv"
SEQUENCE_FILE = "output_sequences.tsv.gz"
OUTPUT_FILE = "output_sequences_full.tsv.gz"

# Increase CSV field size limit to handle massive intron sequences
csv.field_size_limit(sys.maxsize)

def get_delimiter(file_path):
    """Robustly detect if file is CSV or TSV."""
    with open(file_path, 'r', encoding='utf-8-sig') as f:
        line = f.readline()
        if '\t' in line: return '\t'
        if ',' in line: return ','
    return '\t' # Default

def main():
    if not os.path.exists(METADATA_FILE) or not os.path.exists(SEQUENCE_FILE):
        print("Error: Input files not found.")
        print(f"Ensure '{METADATA_FILE}' and '{SEQUENCE_FILE}' are in the folder.")
        return

    print("1. Loading metadata from original input...")
    metadata_map = {}
    
    # Target columns to extract
    target_cols = [
        "DMSO_PSI", "DMSO_PSI_q0.025", "DMSO_PSI_q0.975", 
        "Risdiplam_PSI", "Risdiplam_PSI_q0.025", "Risdiplam_PSI_q0.975"
    ]

    with open(METADATA_FILE, 'r', encoding='utf-8-sig') as f_in:
        delimiter = get_delimiter(METADATA_FILE)
        reader = csv.DictReader(f_in, delimiter=delimiter)
        
        # Strip whitespace from keys just in case
        reader.fieldnames = [name.strip() for name in reader.fieldnames]
        
        # Verify columns exist
        missing = [c for c in target_cols if c not in reader.fieldnames]
        if missing:
            print(f"Error: Original file is missing columns: {missing}")
            return

        for row in reader:
            eid = row['event_id']
            # Store the 6 columns as a list
            values = [row[c] for c in target_cols]
            metadata_map[eid] = values

    print(f"   Loaded metadata for {len(metadata_map)} events.")

    print("2. Merging with sequence file...")
    print("   (Using compression level 1 for maximum speed...)")
    
    processed_count = 0
    
    # CHANGE: Added compresslevel=1 to output file for 3x write speed
    with gzip.open(SEQUENCE_FILE, 'rt') as f_seq, gzip.open(OUTPUT_FILE, 'wt', compresslevel=1) as f_out:
        seq_reader = csv.reader(f_seq, delimiter='\t')
        
        # Read old header
        try:
            old_header = next(seq_reader)
        except StopIteration:
            print("Sequence file is empty.")
            return

        # Find where sequences start (Seq1_ExonA)
        try:
            seq_start_idx = old_header.index("Seq1_ExonA")
        except ValueError:
            print("Error: Could not find 'Seq1_ExonA' in sequence file header.")
            return
            
        # Create NEW header: event_id + 6 metadata cols + sequences
        new_header = ["event_id"] + target_cols + old_header[seq_start_idx:]
        f_out.write("\t".join(new_header) + "\n")
        
        for row in seq_reader:
            event_id = row[0]
            
            # Retrieve the full metadata for this ID
            if event_id in metadata_map:
                meta_values = metadata_map[event_id]
                
                # Get sequences from the existing file
                sequences = row[seq_start_idx:]
                
                # Construct new row
                new_row = [event_id] + meta_values + sequences
                f_out.write("\t".join(new_row) + "\n")
                processed_count += 1
            else:
                # If ID not found in metadata (rare, but possible if files mismatched), 
                # we skip it to ensure data integrity
                continue
            
            if processed_count % 1000 == 0:
                print(f"   Processed {processed_count} rows...", end='\r')
                
    print(f"\n3. Done! {processed_count} rows merged.")
    print(f"   Output saved to: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()