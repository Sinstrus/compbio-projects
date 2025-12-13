import gzip
import sys
import io

# --- Configuration ---
INPUT_FILE = "output_sequences_full.tsv.gz" 
OUTPUT_FILE = "output_sequences_fully_healed.tsv.gz"

def heal_donor_junction_bytes(upstream_seq, intron_seq):
    """
    Optimized healer working directly on bytes.
    Goal: Ensure 'intron_seq' starts with b'GT'.
    """
    if len(upstream_seq) < 5 or len(intron_seq) < 5:
        return upstream_seq, intron_seq, "TOO_SHORT"
        
    # Case 0: Already Good (Starts with GT)
    start_2 = intron_seq[:2].upper()
    if start_2 == b"GT":
        return upstream_seq, intron_seq, "OK"
        
    # Case 1: Shifted Left (The 'G' is in the Upstream Exon)
    # Pattern: ...G|T...
    up_last = upstream_seq[-1:].upper()
    intron_first = intron_seq[:1].upper()
    
    if up_last == b'G' and intron_first == b'T':
        new_up = upstream_seq[:-1]
        new_intron = upstream_seq[-1:] + intron_seq
        return new_up, new_intron, "HEALED_LEFT"
        
    # Case 2: Shifted Right (The 'G' or 'GT' is deep in Intron)
    # Pattern: ...|GGT...
    intron_next_2 = intron_seq[1:3].upper()
    if intron_next_2 == b"GT":
        new_up = upstream_seq + intron_seq[:1]
        new_intron = intron_seq[1:]
        return new_up, new_intron, "HEALED_RIGHT"

    return upstream_seq, intron_seq, "UNFIXABLE"

def main():
    print("--- Ultra-Fast Auto-Healing ALL Junctions ---")
    
    stats_j1 = {"OK": 0, "HEALED_LEFT": 0, "HEALED_RIGHT": 0, "UNFIXABLE": 0, "TOO_SHORT": 0}
    stats_j2 = {"OK": 0, "HEALED_LEFT": 0, "HEALED_RIGHT": 0, "UNFIXABLE": 0, "TOO_SHORT": 0}
    
    count = 0
    
    # Use io.BufferedReader/Writer for speed
    with gzip.open(INPUT_FILE, 'rb') as f_in, gzip.open(OUTPUT_FILE, 'wb') as f_out:
        # Wrap in buffered readers/writers for performance
        reader = io.BufferedReader(f_in)
        writer = io.BufferedWriter(f_out)
        
        # Read header
        header = reader.readline()
        if not header:
            print("Error: Empty file.")
            return
            
        writer.write(header)
        
        # Parse indices from header to be robust
        header_parts = header.strip().split(b'\t')
        try:
            # Note: Headers might be bytes like b'Seq1_ExonA'
            idx_ea = -1
            idx_i1 = -1
            idx_ps = -1
            idx_i2 = -1
            
            for i, h in enumerate(header_parts):
                if h == b'Seq1_ExonA': idx_ea = i
                elif h == b'Seq2_Intron1': idx_i1 = i
                elif h == b'Seq3_Pseudoexon': idx_ps = i
                elif h == b'Seq4_Intron2': idx_i2 = i
                
            if -1 in [idx_ea, idx_i1, idx_ps, idx_i2]:
                print(f"Error: Columns not found. Header: {header_parts}")
                return
                
        except ValueError:
            print("Error: Columns not found in header.")
            return

        for line in reader:
            parts = line.strip().split(b'\t')
            # Safety check for short lines
            if len(parts) <= max(idx_ea, idx_i1, idx_ps, idx_i2): 
                continue
            
            # 1. Heal Junction 1 (Exon A -> Intron 1)
            ea = parts[idx_ea]
            i1 = parts[idx_i1]
            new_ea, new_i1, s1 = heal_donor_junction_bytes(ea, i1)
            parts[idx_ea] = new_ea
            parts[idx_i1] = new_i1
            stats_j1[s1] += 1
            
            # 2. Heal Junction 2 (Pseudo -> Intron 2)
            ps = parts[idx_ps]
            i2 = parts[idx_i2]
            new_ps, new_i2, s2 = heal_donor_junction_bytes(ps, i2)
            parts[idx_ps] = new_ps
            parts[idx_i2] = new_i2
            stats_j2[s2] += 1
            
            # Write back
            writer.write(b'\t'.join(parts) + b'\n')
            
            count += 1
            if count % 10000 == 0: 
                print(f"Processed {count}...", end='\r')
                
    print(f"\nDone! Processed {count} sequences.")
    print("-" * 40)
    print("Junction 1 Stats (Exon A -> Intron 1):")
    total = max(1, count)
    for k, v in stats_j1.items(): print(f"  {k}: {v} ({v/total*100:.1f}%)")
    print("-" * 40)
    print("Junction 2 Stats (Pseudo -> Intron 2):")
    for k, v in stats_j2.items(): print(f"  {k}: {v} ({v/total*100:.1f}%)")
    print("-" * 40)
    print(f"Saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()