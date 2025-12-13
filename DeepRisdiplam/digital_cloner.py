import csv
import gzip
import numpy as np
import sys

# --- Configuration ---
INPUT_FILE = "output_sequences_healed.tsv.gz"
OUTPUT_X_FILE = "dataset_X_sequences_healed.tsv"
OUTPUT_Y_FILE = "dataset_y_labels_healed.tsv"

# Filtering Thresholds
MAX_BASAL_PSI = 0.10     # Max allowed PSI in DMSO to be considered a candidate
MIN_ODDS_RATIO = 10.0    # Min fold change in Odds to be considered "Responsive"

# Sequence Limits
EXON_A_LEN = 50
INTRON_1_PART_LEN = 250
INTRON_2_PART_LEN = 250
EXON_B_LEN = 50

csv.field_size_limit(sys.maxsize)

def calculate_odds_ratio(dmso_psi, ris_psi):
    """
    Calculates Odds Ratio with epsilon for stability.
    OR = (Ris_PSI / (1-Ris_PSI)) / (DMSO_PSI / (1-DMSO_PSI))
    """
    epsilon = 1e-6
    
    # Clip probabilities to avoid exactly 0 or 1
    d_p = np.clip(dmso_psi, epsilon, 1 - epsilon)
    r_p = np.clip(ris_psi, epsilon, 1 - epsilon)
    
    odds_dmso = d_p / (1 - d_p)
    odds_ris = r_p / (1 - r_p)
    
    return odds_ris / odds_dmso

def stitch_minigene(seq_ea, seq_i1, seq_ps, seq_i2, seq_eb):
    """
    Implements the "Digital Minigene" logic.
    Stitches: ExonA(End) - Intron1(AC) - Pseudo - Intron2(AC) - ExonB(Start)
    AC = Anchor Context (Start + End)
    """
    # 1. Exon A (Last N bp)
    part1 = seq_ea[-EXON_A_LEN:] if len(seq_ea) > EXON_A_LEN else seq_ea
    
    # 2. Intron 1 (First N + Last N)
    if len(seq_i1) <= (INTRON_1_PART_LEN * 2):
        part2 = seq_i1
    else:
        part2 = seq_i1[:INTRON_1_PART_LEN] + seq_i1[-INTRON_1_PART_LEN:]
        
    # 3. Pseudoexon (Full)
    part3 = seq_ps
    
    # 4. Intron 2 (First N + Last N)
    if len(seq_i2) <= (INTRON_2_PART_LEN * 2):
        part4 = seq_i2
    else:
        part4 = seq_i2[:INTRON_2_PART_LEN] + seq_i2[-INTRON_2_PART_LEN:]
        
    # 5. Exon B (First N bp)
    part5 = seq_eb[:EXON_B_LEN]
    
    return part1 + part2 + part3 + part4 + part5

def main():
    print(f"--- Generating Minigene Dataset from {INPUT_FILE} ---")
    
    count_total = 0
    count_filtered_out = 0
    count_positives = 0
    count_negatives = 0
    
    with gzip.open(INPUT_FILE, 'rt') as f_in, \
         open(OUTPUT_X_FILE, 'w') as f_x, \
         open(OUTPUT_Y_FILE, 'w') as f_y:
        
        reader = csv.reader(f_in, delimiter='\t')
        
        try:
            header = next(reader)
        except StopIteration:
            print("Error: Empty file.")
            return

        # Map columns
        try:
            # Metadata
            idx_dmso = header.index("DMSO_PSI")
            idx_ris = header.index("Risdiplam_PSI")
            
            # Sequences
            idx_ea = header.index("Seq1_ExonA")
            idx_i1 = header.index("Seq2_Intron1")
            idx_ps = header.index("Seq3_Pseudoexon")
            idx_i2 = header.index("Seq4_Intron2")
            idx_eb = header.index("Seq5_ExonB")
            
        except ValueError as e:
            print(f"Error finding columns: {e}")
            return

        # Write Headers for Output
        f_x.write("event_id\tminigene_sequence\n")
        f_y.write("event_id\tlabel\todds_ratio\tdmso_psi\tris_psi\n")

        print("Processing rows...")
        
        for row in reader:
            count_total += 1
            if count_total % 10000 == 0:
                print(f"  Processed {count_total}...", end='\r')

            try:
                # 1. Get PSI values
                dmso_psi = float(row[idx_dmso])
                ris_psi = float(row[idx_ris])
                
                # 2. FILTER: Check Basal PSI (Must be low)
                if dmso_psi > MAX_BASAL_PSI:
                    count_filtered_out += 1
                    continue
                
                # 3. Calculate Odds Ratio and Label
                odds_ratio = calculate_odds_ratio(dmso_psi, ris_psi)
                
                if odds_ratio >= MIN_ODDS_RATIO:
                    label = 1
                    count_positives += 1
                else:
                    label = 0
                    count_negatives += 1
                
                # 4. Extract and Stitch Sequence
                full_seq = stitch_minigene(
                    row[idx_ea], 
                    row[idx_i1], 
                    row[idx_ps], 
                    row[idx_i2], 
                    row[idx_eb]
                )
                
                # 5. Save
                eid = row[0]
                f_x.write(f"{eid}\t{full_seq}\n")
                f_y.write(f"{eid}\t{label}\t{odds_ratio:.4f}\t{dmso_psi:.4f}\t{ris_psi:.4f}\n")
                
            except ValueError:
                # Handle missing/malformed PSI values
                count_filtered_out += 1
                continue

    print("\n\n--- Dataset Statistics ---")
    print(f"Total Rows Scanned: {count_total}")
    print(f"Excluded (High Basal PSI > {MAX_BASAL_PSI}): {count_filtered_out}")
    print(f"Candidates (Low Basal PSI): {count_positives + count_negatives}")
    print("-" * 30)
    print(f"POSITIVE Samples (Responder, OR >= {MIN_ODDS_RATIO}): {count_positives}")
    print(f"NEGATIVE Samples (Non-Responder, OR < {MIN_ODDS_RATIO}): {count_negatives}")
    
    if count_positives > 0:
        ratio = count_negatives / count_positives
        print(f"Class Imbalance Ratio: 1 Positive : {ratio:.1f} Negatives")
    else:
        print("Class Imbalance: No positives found!")

if __name__ == "__main__":
    main()