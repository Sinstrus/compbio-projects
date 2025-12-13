import pandas as pd
from pyfaidx import Fasta
import sys
import os
from tqdm import tqdm

# Files
PARSED_COORDS = "parsed_coordinates.tsv"
PSI_DATA = "GSE221868_exon_skipping.filtered.psi.csv"
GENOME_FILE = "hg38.fa"
OUTPUT_FILE = "risdiplam_full_dataset.csv.gz" # .gz automatically triggers compression

def get_seq(genome, location_str, strand):
    """
    Fetches sequence for a location string "chr:start-end".
    Handles strand: returns Reverse Complement if strand is '-'.
    """
    try:
        if pd.isna(location_str): return None
        
        chrom, span = location_str.split(':')
        start, end = map(int, span.split('-'))
        
        # Determine length
        length = end - start + 1
        
        # Safety check for massive introns (optional warning)
        # if length > 100000: return None 
        
        if strand == '+':
            seq = genome[chrom][start-1 : end].seq
        else:
            # -strand: Reverse Complement
            seq = genome[chrom][start-1 : end].reverse.complement.seq
            
        return seq.upper()
    except Exception as e:
        return None

def main():
    # Pre-flight cleanup
    if os.path.exists("risdiplam_full_dataset.csv"):
        print("Removing old uncompressed file to free space...")
        os.remove("risdiplam_full_dataset.csv")

    print("Loading Data...")
    if not pd.io.common.file_exists(PARSED_COORDS):
        print(f"Error: {PARSED_COORDS} not found. Run parse_coordinates.py first.")
        sys.exit(1)
        
    coords = pd.read_csv(PARSED_COORDS, sep='\t')
    psi = pd.read_csv(PSI_DATA)
    
    print("Merging Tables...")
    df = pd.merge(coords, psi, on='event_id')
    
    print(f"Loading Genome {GENOME_FILE}...")
    try:
        genome = Fasta(GENOME_FILE)
    except FileNotFoundError:
        print("ERROR: hg38.fa not found.")
        sys.exit(1)

    components = ['bio_ExonA', 'bio_Intron1', 'bio_PseudoExon', 'bio_Intron2', 'bio_ExonB']
    new_cols = {f"seq_{comp.replace('bio_', '')}": [] for comp in components}
    
    print("Fetching Sequences...")
    for idx, row in tqdm(df.iterrows(), total=len(df)):
        strand = row['strand']
        for comp in components:
            loc_str = row[comp]
            seq = get_seq(genome, loc_str, strand)
            col_name = f"seq_{comp.replace('bio_', '')}"
            new_cols[col_name].append(seq)

    for col_name, data in new_cols.items():
        df[col_name] = data

    df['dPSI'] = df['Risdiplam_PSI'] - df['DMSO_PSI']
    df = df.dropna(subset=new_cols.keys())
    
    # --- Statistics Report ---
    print("\n--- Component Statistics ---")
    for col in new_cols.keys():
        # Calculate mean length of strings in this column
        avg_len = df[col].astype(str).map(len).mean()
        print(f"  {col}: {avg_len:.0f} bp (Average)")

    psi_cols = ['DMSO_PSI', 'Risdiplam_PSI', 'dPSI']
    seq_cols = list(new_cols.keys())
    meta_cols = ['event_id', 'chrom', 'strand']
    final_cols = meta_cols + psi_cols + seq_cols
    
    out_df = df[final_cols]
    
    print(f"\nSaving to {OUTPUT_FILE} (Compressed)...")
    out_df.to_csv(OUTPUT_FILE, index=False, compression='gzip')
    print("Success!")

if __name__ == "__main__":
    main()