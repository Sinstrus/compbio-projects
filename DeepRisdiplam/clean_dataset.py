import csv
import gzip
import sys

# --- Configuration ---
INPUT_FILE = "output_sequences_healed.tsv.gz"
OUTPUT_FILE = "dataset_X_sequences_strict.tsv"
OUTPUT_Y_FILE = "dataset_y_labels_strict.tsv"

# Constants for Stitching
EXON_A_LEN = 50
INTRON_1_PART_LEN = 250
INTRON_2_PART_LEN = 250
EXON_B_LEN = 50

# Filter Params
MAX_BASAL_PSI = 0.10     
MIN_ODDS_RATIO = 10.0    

csv.field_size_limit(sys.maxsize)

def calculate_odds_ratio(dmso_psi, ris_psi):
    epsilon = 1e-6
    d_p = max(min(dmso_psi, 1 - epsilon), epsilon)
    r_p = max(min(ris_psi, 1 - epsilon), epsilon)
    return (r_p / (1 - r_p)) / (d_p / (1 - d_p))

def stitch_minigene(seq_ea, seq_i1, seq_ps, seq_i2, seq_eb):
    part1 = seq_ea[-EXON_A_LEN:] if len(seq_ea) > EXON_A_LEN else seq_ea
    
    if len(seq_i1) <= (INTRON_1_PART_LEN * 2): part2 = seq_i1
    else: part2 = seq_i1[:INTRON_1_PART_LEN] + seq_i1[-INTRON_1_PART_LEN:]
        
    part3 = seq_ps
    
    if len(seq_i2) <= (INTRON_2_PART_LEN * 2): part4 = seq_i2
    else: part4 = seq_i2[:INTRON_2_PART_LEN] + seq_i2[-INTRON_2_PART_LEN:]
        
    part5 = seq_eb[:EXON_B_LEN]
    return part1 + part2 + part3 + part4 + part5

def main():
    print(f"--- Creating STRICT Dataset (GT-only) ---")
    
    count_total = 0
    count_kept = 0
    count_dropped = 0
    
    with gzip.open(INPUT_FILE, 'rt') as f_in, \
         open(OUTPUT_FILE, 'w') as f_x, \
         open(OUTPUT_Y_FILE, 'w') as f_y:
        
        reader = csv.DictReader(f_in, delimiter='\t')
        
        # Write Headers
        f_x.write("event_id\tminigene_sequence\n")
        f_y.write("event_id\tlabel\todds_ratio\tdmso_psi\tris_psi\n")

        for row in reader:
            count_total += 1
            
            # STRICT CHECK: Must start with GT
            seq_i2 = row['Seq4_Intron2'].upper()
            if not seq_i2.startswith("GT"):
                count_dropped += 1
                continue
                
            try:
                dmso_psi = float(row['DMSO_PSI'])
                ris_psi = float(row['Risdiplam_PSI'])
                
                if dmso_psi > MAX_BASAL_PSI: 
                    continue # Skip high basal (not a switch)
                
                odds_ratio = calculate_odds_ratio(dmso_psi, ris_psi)
                label = 1 if odds_ratio >= MIN_ODDS_RATIO else 0
                
                full_seq = stitch_minigene(
                    row['Seq1_ExonA'], row['Seq2_Intron1'], 
                    row['Seq3_Pseudoexon'], row['Seq4_Intron2'], 
                    row['Seq5_ExonB']
                )
                
                eid = row['event_id']
                f_x.write(f"{eid}\t{full_seq}\n")
                f_y.write(f"{eid}\t{label}\t{odds_ratio:.4f}\t{dmso_psi:.4f}\t{ris_psi:.4f}\n")
                
                count_kept += 1
                
            except ValueError:
                continue
                
            if count_total % 50000 == 0: 
                print(f"Processed {count_total}...", end='\r')

    print(f"\n\n--- Strict Filtering Complete ---")
    print(f"Total Processed: {count_total}")
    print(f"Dropped (Non-GT): {count_dropped} ({count_dropped/count_total*100:.2f}%)")
    print(f"Kept (Valid GT + Low Basal): {count_kept}")
    print(f"Saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()