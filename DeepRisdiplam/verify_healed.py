import gzip
import csv
import sys

INPUT_FILE = "output_sequences_healed.tsv.gz"
csv.field_size_limit(sys.maxsize)

def main():
    print(f"Checking {INPUT_FILE}...")
    
    gt_count = 0
    total = 0
    
    with gzip.open(INPUT_FILE, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            seq = row['Seq4_Intron2'].upper()
            if seq.startswith("GT"):
                gt_count += 1
            total += 1
            
            if total % 50000 == 0: print(f"Checked {total}...", end='\r')
            
    print(f"\nTotal: {total}")
    print(f"Starts with GT: {gt_count}")
    print(f"Success Rate: {gt_count/total*100:.2f}%")

if __name__ == "__main__":
    main()