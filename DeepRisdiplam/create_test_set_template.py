import pandas as pd
import sys

# --- Configuration ---
INPUT_FILE = "dataset_components.tsv"
OUTPUT_FILE = "test_set_template.csv"
NUM_CONTROLS = 10

def main():
    print(f"--- Creating Test Set Template ---")
    
    try:
        df = pd.read_csv(INPUT_FILE, sep='\t')
    except FileNotFoundError:
        print("Error: dataset_components.tsv not found.")
        return

    # 1. Select Positive Controls (High Odds Ratio)
    # We sort by odds_ratio (calculated implicitly during generation, but we only have label here)
    # We'll just take random ones since we don't have the OR column in components file easily accessible
    # Wait, we do not have OR in components file. We will just take the first N.
    
    positives = df[df['label'] == 1].head(NUM_CONTROLS).copy()
    positives['type'] = 'Positive Control'
    
    # 2. Select Negative Controls (Non-Responders)
    negatives = df[df['label'] == 0].head(NUM_CONTROLS).copy()
    negatives['type'] = 'Negative Control'
    
    # 3. Combine
    controls = pd.concat([positives, negatives])
    
    # 4. Prepare Output DataFrame
    # We want columns: id, type, ExonA, Intron1, Pseudoexon, Intron2, ExonB
    output_df = controls[['event_id', 'type', 'ExonA', 'Intron1', 'Pseudoexon', 'Intron2', 'ExonB']]
    
    # 5. Add Placeholder Rows for Your New Sequences
    new_rows = []
    for i in range(5):
        new_rows.append({
            'event_id': f'My_New_Seq_{i+1}',
            'type': 'Test Sequence',
            'ExonA': 'Insert_ExonA_Here',
            'Intron1': 'Insert_Intron1_Here',
            'Pseudoexon': 'Insert_Pseudoexon_Here',
            'Intron2': 'Insert_Intron2_Here', 
            'ExonB': 'Insert_ExonB_Here'
        })
        
    final_df = pd.concat([output_df, pd.DataFrame(new_rows)], ignore_index=True)
    
    # 6. Save
    final_df.to_csv(OUTPUT_FILE, index=False)
    print(f"Saved template to {OUTPUT_FILE}")
    print("Open this file, add your sequences, and save it.")
    print("NOTE: Ensure Intron2 starts with 'GT' or 'gt' if you want a fair test.")

if __name__ == "__main__":
    main()