import pandas as pd
import re
from pathlib import Path

# Load the dataset
INPUT_FILE = "GSE221868_exon_skipping.filtered.psi.csv"
OUTPUT_FILE = "parsed_coordinates.tsv"

def parse_event_id(event_id):
    """
    Parses strings like: chr21:46635673-46635763,46636438-46636547,46636895-46636990+
    Returns a dictionary of the 5 genomic regions in BIOLOGICAL order (ExonA -> ExonB).
    """
    # Regex to grab Chrom, the 3 coordinate pairs, and the strand
    # Matches: chr:start-end,start-end,start-end[strand]
    pattern = r"([^:]+):(\d+)-(\d+),(\d+)-(\d+),(\d+)-(\d+)([\+\-])"
    match = re.match(pattern, event_id.strip('"')) # Strip quotes just in case
    
    if not match:
        return None

    chrom = match.group(1)
    strand = match.group(8)
    
    # Raw physical coordinates (Left to Right on Genome)
    # Block 1 (Leftmost)
    b1_start, b1_end = int(match.group(2)), int(match.group(3))
    # Block 2 (Middle)
    b2_start, b2_end = int(match.group(4)), int(match.group(5))
    # Block 3 (Rightmost)
    b3_start, b3_end = int(match.group(6)), int(match.group(7))

    # Calculate the physical Introns (Gaps)
    # Gap 1 is between Block 1 and Block 2
    g1_start, g1_end = b1_end + 1, b2_start - 1
    
    # Gap 2 is between Block 2 and Block 3
    g2_start, g2_end = b2_end + 1, b3_start - 1

    # --- Assign Biological Identity based on Strand ---
    
    if strand == '+':
        # 5' -> 3' is Left -> Right
        regions = {
            'chrom': chrom,
            'strand': strand,
            'bio_ExonA': f"{chrom}:{b1_start}-{b1_end}",
            'bio_Intron1': f"{chrom}:{g1_start}-{g1_end}",
            'bio_PseudoExon': f"{chrom}:{b2_start}-{b2_end}",
            'bio_Intron2': f"{chrom}:{g2_start}-{g2_end}",
            'bio_ExonB': f"{chrom}:{b3_start}-{b3_end}",
            # Metadata for fetching later
            'fetch_instructions': [
                (chrom, b1_start, b1_end, '+'), # Exon A
                (chrom, g1_start, g1_end, '+'), # Intron 1
                (chrom, b2_start, b2_end, '+'), # Pseudo
                (chrom, g2_start, g2_end, '+'), # Intron 2
                (chrom, b3_start, b3_end, '+')  # Exon B
            ]
        }
    else:
        # 5' -> 3' is Right -> Left
        # Exon A is the Rightmost block (Block 3)
        # Intron 1 is the Rightmost gap (Gap 2)
        regions = {
            'chrom': chrom,
            'strand': strand,
            'bio_ExonA': f"{chrom}:{b3_start}-{b3_end}",     # Block 3 (physically last, biologically first)
            'bio_Intron1': f"{chrom}:{g2_start}-{g2_end}",   # Gap 2
            'bio_PseudoExon': f"{chrom}:{b2_start}-{b2_end}", # Block 2
            'bio_Intron2': f"{chrom}:{g1_start}-{g1_end}",   # Gap 1
            'bio_ExonB': f"{chrom}:{b1_start}-{b1_end}",     # Block 1
             # Metadata for fetching later
            'fetch_instructions': [
                (chrom, b3_start, b3_end, '-'), # Exon A (RevComp required)
                (chrom, g2_start, g2_end, '-'), # Intron 1 (RevComp required)
                (chrom, b2_start, b2_end, '-'), # Pseudo (RevComp required)
                (chrom, g1_start, g1_end, '-'), # Intron 2 (RevComp required)
                (chrom, b1_start, b1_end, '-')  # Exon B (RevComp required)
            ]
        }
        
    return regions

def main():
    # Check if input file exists before attempting to read
    if not Path(INPUT_FILE).exists():
        print(f"Error: Input file '{INPUT_FILE}' not found.")
        print("Please ensure the CSV file is in the same directory.")
        return

    df = pd.read_csv(INPUT_FILE)
    print(f"Loaded {len(df)} rows.")

    parsed_data = []
    
    # Process all rows
    print("Parsing coordinates...")
    for idx, row in df.iterrows():
        try:
            res = parse_event_id(row['event_id'])
            if res:
                # We add the original ID to link back
                res['event_id'] = row['event_id']
                parsed_data.append(res)
        except Exception as e:
            print(f"Error parsing row {idx}: {e}")

    # Create new DataFrame
    result_df = pd.DataFrame(parsed_data)
    
    # Reorder columns for clarity
    cols = ['event_id', 'chrom', 'strand', 'bio_ExonA', 'bio_Intron1', 'bio_PseudoExon', 'bio_Intron2', 'bio_ExonB']
    result_df = result_df[cols]
    
    print("\n--- Verification (Positive Strand Example) ---")
    # Using your specific example from the chat
    ex_pos = "chr21:46635673-46635763,46636438-46636547,46636895-46636990+"
    if ex_pos in result_df['event_id'].values:
        print(result_df[result_df['event_id'] == ex_pos].iloc[0])
    else:
        print(f"Warning: Verification ID {ex_pos} not found in dataset.")

    print("\n--- Verification (Negative Strand Example) ---")
    # Using your specific example from the chat
    ex_neg = "chr21:46240810-46241017,46242801-46242931,46243464-46243722-"
    if ex_neg in result_df['event_id'].values:
        print(result_df[result_df['event_id'] == ex_neg].iloc[0])
    else:
         print(f"Warning: Verification ID {ex_neg} not found in dataset.")

    # Save
    result_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
    print(f"\nSaved parsed map to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()