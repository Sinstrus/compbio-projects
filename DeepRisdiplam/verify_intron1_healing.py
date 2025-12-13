import gzip
import csv
import sys

INPUT_FILE = "output_sequences_fully_healed.tsv.gz" 

# Increase CSV limit
csv.field_size_limit(sys.maxsize)

def main():
    print(f"Checking Intron 1 Integrity in {INPUT_FILE}...")
    
    total = 0
    gt_starts = 0
    
    try:
        with gzip.open(INPUT_FILE, 'rt') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                i1 = row['Seq2_Intron1'].upper()
                if i1.startswith("GT"):
                    gt_starts += 1
                total += 1
                
                if total % 50000 == 0: print(f"Checked {total}...", end='\r')
                
    except FileNotFoundError:
        print("File not found. Run heal_dataset_full.py first.")
        return

    print(f"\nTotal Sequences: {total}")
    print(f"Intron 1 starts with GT: {gt_starts} ({gt_starts/total*100:.2f}%)")

if __name__ == "__main__":
    main()