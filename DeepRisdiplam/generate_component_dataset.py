import csv
import gzip
import sys
import numpy as np

# --- Configuration ---
INPUT_FILE = "output_sequences_fully_healed.tsv.gz"  # NEW FILE
OUTPUT_FILE = "dataset_components_v8.tsv"

# Filter Params
MAX_BASAL_PSI = 0.10     
MIN_ODDS_RATIO = 10.0    

csv.field_size_limit(sys.maxsize)

def calculate_odds_ratio(dmso_psi, ris_psi):
    epsilon = 1e-6
    d_p = max(min(dmso_psi, 1 - epsilon), epsilon)
    r_p = max(min(ris_psi, 1 - epsilon), epsilon)
    return (r_p / (1 - r_p)) / (d_p / (1 - d_p))

def main():
    print(f"--- Generating Component Dataset V8 (Fully Healed) ---")
    
    count_pos = 0
    count_neg = 0
    
    with gzip.open(INPUT_FILE, 'rt') as f_in, open(OUTPUT_FILE, 'w') as f_out:
        reader = csv.DictReader(f_in, delimiter='\t')
        
        headers = ["event_id", "label", "ExonA", "Intron1", "Pseudoexon", "Intron2", "ExonB"]
        writer = csv.DictWriter(f_out, fieldnames=headers, delimiter='\t')
        writer.writeheader()

        for row in reader:
            try:
                dmso_psi = float(row['DMSO_PSI'])
                ris_psi = float(row['Risdiplam_PSI'])
                if dmso_psi > MAX_BASAL_PSI: continue
                
                odds_ratio = calculate_odds_ratio(dmso_psi, ris_psi)
                label = 1 if odds_ratio >= MIN_ODDS_RATIO else 0
                
                # DOUBLE STRICT CHECK: Both Introns must start with GT
                if not row['Seq2_Intron1'].upper().startswith("GT"): continue
                if not row['Seq4_Intron2'].upper().startswith("GT"): continue

                # Save Components
                ea = row['Seq1_ExonA'][-50:] if len(row['Seq1_ExonA']) > 50 else row['Seq1_ExonA']
                
                i1 = row['Seq2_Intron1']
                if len(i1) > 500: i1 = i1[:250] + i1[-250:]
                    
                ps = row['Seq3_Pseudoexon']
                
                i2 = row['Seq4_Intron2']
                if len(i2) > 500: i2 = i2[:250] + i2[-250:]
                    
                eb = row['Seq5_ExonB'][:50]
                
                writer.writerow({
                    "event_id": row['event_id'],
                    "label": label,
                    "ExonA": ea, "Intron1": i1, "Pseudoexon": ps, "Intron2": i2, "ExonB": eb
                })
                
                if label == 1: count_pos += 1
                else: count_neg += 1
                
            except ValueError: continue

    print(f"Saved {count_pos} Positives and {count_neg} Negatives.")

if __name__ == "__main__":
    main()