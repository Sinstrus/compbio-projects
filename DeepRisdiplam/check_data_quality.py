import gzip
import csv
import sys

INPUT_FILE = "output_sequences_full.tsv.gz"

def peek_data(n=5):
    print(f"--- Checking {INPUT_FILE} ---\n")
    try:
        with gzip.open(INPUT_FILE, 'rt') as f:
            reader = csv.reader(f, delimiter='\t')
            
            # 1. Check Header
            try:
                header = next(reader)
            except StopIteration:
                print("Error: File is empty.")
                return

            print(f"HEADER COLUMNS: {header}")
            
            # Verify critical columns exist
            if "DMSO_PSI" in header and "Risdiplam_PSI" in header:
                print("✅ STATUS: 'DMSO_PSI' and 'Risdiplam_PSI' columns found.")
            else:
                print("❌ WARNING: Critical PSI columns missing.")
            
            print("-" * 50)

            # 2. Check Rows
            for i, row in enumerate(reader):
                if i >= n:
                    break
                
                print(f"\nRow {i+1} ID: {row[0]}")
                print(f"  DMSO PSI: {row[1]}")
                print(f"  DMSO PSI_q0.025: {row[2]}")
                print(f"  DMSO PSI_q0.975: {row[3]}")
                print(f"  Risdiplam PSI: {row[4]}")
                print(f"  Risdiplam PSI_q0.025: {row[5]}")
                print(f"  Risdiplam PSI_q0.975: {row[6]}")
                
                # Check sequence lengths (indices 3 to 7)
                if len(row) >= 8:
                    seq_names = ["ExonA", "Intron1", "Pseudo", "Intron2", "ExonB"]
                    seqs = row[3:]
                    print("  Sequence Lengths:")
                    for name, seq in zip(seq_names, seqs):
                        print(f"    - {name}: {len(seq)} bp {'(Empty!)' if len(seq)==0 else ''}")
                        # Optional: Print first 10 chars to visually verify 5' start
                        if len(seq) > 0:
                            print(f"      Start: {seq[:10]}...")

    except FileNotFoundError:
        print(f"Error: {INPUT_FILE} not found. Please run genomic_sequence_extractor.py first.")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    peek_data()