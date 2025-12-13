import pandas as pd
import gzip
import csv
import sys

# --- Configuration ---
INPUT_FILE = "output_sequences_full.tsv.gz"
INPUT_LABELS = "dataset_y_labels.tsv" 

csv.field_size_limit(sys.maxsize)

def main():
    print("--- Verifying Data Integrity (Splice Site Check) ---")
    
    # 1. Load Labels (to focus on the 'Positives' we care about)
    try:
        labels_df = pd.read_csv(INPUT_LABELS, sep='\t')
        positive_ids = set(labels_df[labels_df['label'] == 1]['event_id'])
        print(f"Checking {len(positive_ids)} positive samples...")
    except FileNotFoundError:
        print("Labels file not found. Checking ALL rows...")
        positive_ids = None

    # Stats Counters
    total = 0
    intron2_starts_GT = 0  # Correct
    intron2_starts_T = 0   # Off-by-one (Skipped G)
    intron2_starts_A = 0   # Random?
    intron2_starts_C = 0   # Random?
    
    # Example Storage
    examples = []

    try:
        with gzip.open(INPUT_FILE, 'rt') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            for row in reader:
                eid = row['event_id']
                
                # Filter for positives if labels exist
                if positive_ids and eid not in positive_ids:
                    continue
                
                # Check Intron 2 Start (The Donor)
                seq_i2 = row['Seq4_Intron2']
                if not seq_i2: continue
                
                # Check first 2 bases
                start_2 = seq_i2[:2].upper()
                
                if start_2 == "GT":
                    intron2_starts_GT += 1
                elif start_2.startswith("T"): # Starts with T... might be TG or TA
                    intron2_starts_T += 1
                elif start_2.startswith("A"):
                    intron2_starts_A += 1
                elif start_2.startswith("C"):
                    intron2_starts_C += 1
                
                total += 1
                
                # Save first 5 for visual inspection
                if len(examples) < 5:
                    examples.append({
                        "id": eid,
                        "Pseudo_End": row['Seq3_Pseudoexon'][-10:],
                        "Intron2_Start": seq_i2[:10]
                    })
                    
    except FileNotFoundError:
        print(f"Error: {INPUT_FILE} not found.")
        return

    # --- REPORT ---
    print(f"\nTotal Sequences Checked: {total}")
    print("-" * 30)
    print(f"Starts with GT (Canonical): {intron2_starts_GT} ({intron2_starts_GT/total*100:.1f}%)")
    print(f"Starts with T  (Off-by-1?): {intron2_starts_T} ({intron2_starts_T/total*100:.1f}%)")
    print(f"Starts with A             : {intron2_starts_A} ({intron2_starts_A/total*100:.1f}%)")
    print(f"Starts with C             : {intron2_starts_C} ({intron2_starts_C/total*100:.1f}%)")
    print("-" * 30)
    
    print("\n--- Visual Inspection (First 5 Examples) ---")
    print("Format: ...[Exon End] | [Intron Start]...")
    for ex in examples:
        print(f"\nID: {ex['id']}")
        print(f"Sequence: ...{ex['Pseudo_End']} | {ex['Intron2_Start']}...")
        
        # Diagnosis
        if ex['Intron2_Start'].startswith("GT"):
            print("Status:   ✅ OK")
        elif ex['Intron2_Start'].startswith("T"):
            print("Status:   ❌ Suspicious (Starts with T)")
        else:
            print(f"Status:   ❌ Non-Canonical ({ex['Intron2_Start'][:2]})")

if __name__ == "__main__":
    main()