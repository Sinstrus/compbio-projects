import pandas as pd
import sys

# --- Configuration ---
INPUT_X_FILE = "dataset_X_sequences_strict.tsv"
INPUT_Y_FILE = "dataset_y_labels_strict.tsv"
OFFSET_END_FROM_RIGHT = 550

def main():
    print("--- Verifying Training Positives ---")
    try:
        df_x = pd.read_csv(INPUT_X_FILE, sep='\t')
        df_y = pd.read_csv(INPUT_Y_FILE, sep='\t')
        df = pd.merge(df_x, df_y, on="event_id")
    except Exception as e:
        print(f"Error loading files: {e}")
        return

    positives = df[df['label'] == 1]['minigene_sequence'].values
    print(f"Loaded {len(positives)} positive sequences.")
    
    if len(positives) == 0:
        print("No positives found.")
        return

    # Check GT at offset
    correct = 0
    total = 0
    
    print(f"Checking for GT at offset {OFFSET_END_FROM_RIGHT} from right...")
    
    for seq in positives:
        # Calculate index: len - 550
        idx = len(seq) - OFFSET_END_FROM_RIGHT
        
        # Extract 2 bases
        dinuc = seq[idx : idx+2].upper()
        
        if dinuc == "GT":
            correct += 1
        else:
            # Print failure for inspection
            if total < 5:
                print(f"FAILURE: Found '{dinuc}' at index {idx} (Len: {len(seq)})")
        total += 1
        
    print("-" * 30)
    print(f"Total Checked: {total}")
    print(f"Correct GT: {correct}")
    print(f"Success Rate: {correct/total*100:.2f}%")
    
    if correct < total:
        print("\nCRITICAL ERROR: The training script assumes GT is at (len - 550).")
        print("But your data has GT somewhere else.")
        print("This breaks 'mutate_dead_donor', so Negatives still have GT.")
        print("Result: Model sees Positive (GT) vs Negative (GT). Learns nothing.")

if __name__ == "__main__":
    main()