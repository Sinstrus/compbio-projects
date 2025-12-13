import sys
import os
import csv
import gzip

# --- Configuration ---
FASTA_PATH = "hg38.fa"
FASTA_INDEX_PATH = "hg38.fa.fai"
INPUT_FILE = "GSE221868_exon_skipping.filtered.psi.csv"
OUTPUT_FILE = "output_sequences.tsv.gz"

class FastaReader:
    """
    A helper class to read random genomic regions from an indexed FASTA file
    without loading the entire genome into memory.
    """
    def __init__(self, fasta_path, index_path):
        self.fasta_path = fasta_path
        self.index = {}
        
        if not os.path.exists(index_path):
            raise FileNotFoundError(f"Index file not found: {index_path}. Please ensure .fai file exists.")

        print(f"Loading index from {index_path}...")
        with open(index_path, 'r') as f:
            for line in f:
                parts = line.split('\t')
                chrom = parts[0]
                length = int(parts[1])
                offset = int(parts[2])
                linebases = int(parts[3])
                linewidth = int(parts[4])
                self.index[chrom] = {
                    'length': length,
                    'offset': offset,
                    'linebases': linebases,
                    'linewidth': linewidth
                }
        self.file_handle = open(self.fasta_path, 'rb')

    def close(self):
        if self.file_handle:
            self.file_handle.close()

    def get_file_offset(self, chrom, coord, chrom_info):
        """Calculates the byte offset in the file for a given 0-based coordinate."""
        offset = chrom_info['offset']
        linebases = chrom_info['linebases']
        linewidth = chrom_info['linewidth']
        
        return offset + (coord // linebases) * linewidth + (coord % linebases)

    def fetch_seq(self, chrom, start, end):
        """
        Fetches sequence from chrom:start-end (1-based inclusive).
        """
        if chrom not in self.index:
            raise KeyError(f"Chromosome {chrom} not found in index.")
        
        chrom_info = self.index[chrom]
        
        # Convert 1-based inclusive to 0-based half-open for Python slicing logic
        start_0 = start - 1
        end_0 = end # because Python end is exclusive, which matches 1-based inclusive end in count
        
        if start_0 < 0 or end_0 > chrom_info['length']:
             # Basic bounds check, though often simple clipping is preferred in bioinformatics
             raise ValueError(f"Coordinates {chrom}:{start}-{end} out of bounds.")

        # Calculate byte offsets
        start_byte = self.get_file_offset(chrom, start_0, chrom_info)
        end_byte = self.get_file_offset(chrom, end_0, chrom_info)
        
        read_length = end_byte - start_byte
        
        self.file_handle.seek(start_byte)
        raw_bytes = self.file_handle.read(read_length)
        
        # Remove newlines/carriage returns to get pure sequence
        seq = raw_bytes.replace(b'\n', b'').replace(b'\r', b'')
        return seq.decode('utf-8')

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def parse_coordinate_string(coord_str):
    """
    Parses string like 'chr21:46635673-46635763,46636438-46636547,46636895-46636990+'
    Returns chrom, list of (start, end) tuples, and strand.
    """
    # Split chromosome
    chrom, rest = coord_str.split(':')
    
    # The strand is the last character
    strand = rest[-1]
    rest = rest[:-1]
    
    # Split the blocks
    blocks = rest.split(',')
    ranges = []
    for b in blocks:
        s, e = map(int, b.split('-'))
        ranges.append((s, e))
        
    return chrom, ranges, strand

def main():
    if not os.path.exists(INPUT_FILE):
        print(f"Error: {INPUT_FILE} not found. Please create it per instructions.")
        return

    try:
        fasta = FastaReader(FASTA_PATH, FASTA_INDEX_PATH)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return

    print(f"Processing {INPUT_FILE}...")
    
    # Open with utf-8-sig to handle Byte Order Marks (BOM) often found in Excel exports
    with open(INPUT_FILE, 'r', encoding='utf-8-sig', newline='') as f_in, gzip.open(OUTPUT_FILE, 'wt') as f_out:
        
        # Robust delimiter detection
        first_line = f_in.readline()
        f_in.seek(0) # Reset to start
        
        # Heuristic: if there are tabs, it's likely TSV. If no tabs but commas, likely CSV.
        delimiter = '\t'
        if '\t' not in first_line and ',' in first_line:
            delimiter = ','
            print(f"Detected delimiter: comma (',')")
        else:
            print(f"Detected delimiter: tab ('\\t') or whitespace")
            
        reader = csv.reader(f_in, delimiter=delimiter)
        
        try:
            header_parts = next(reader)
        except StopIteration:
            print("Error: Input file is empty.")
            return

        # Clean header parts (strip whitespace)
        header_parts = [h.strip() for h in header_parts]

        try:
            # Find indices dynamically
            id_idx = header_parts.index("event_id")
            
            # Simple matching for columns
            dmso_idx = -1
            ris_idx = -1
            
            for i, h in enumerate(header_parts):
                if "DMSO_PSI" == h:
                    dmso_idx = i
                elif "Risdiplam_PSI" == h:
                    ris_idx = i
            
            if dmso_idx == -1 or ris_idx == -1:
                print("Warning: Could not find exact 'DMSO_PSI' or 'Risdiplam_PSI' columns. Using columns 1 and 4 as fallback.")
                dmso_idx = 1
                ris_idx = 4
                
        except ValueError:
            print(f"Error: Could not find 'event_id' in header.")
            print(f"Header found: {header_parts}")
            print("Please check your input file format.")
            return

        # Write Output Header
        output_header = ["event_id", "DMSO_PSI", "Risdiplam_PSI", 
                         "Seq1_ExonA", "Seq2_Intron1", "Seq3_Pseudoexon", "Seq4_Intron2", "Seq5_ExonB"]
        f_out.write("\t".join(output_header) + "\n")

        count = 0
        skipped = 0
        missing_chroms = set()

        for parts in reader:
            if not parts: continue
            
            try:
                event_id = parts[id_idx]
                dmso_psi = parts[dmso_idx]
                ris_psi = parts[ris_idx]
            except IndexError:
                # Handle empty lines or malformed rows gracefully
                continue

            # Ensure event_id looks valid before processing
            if not event_id or ':' not in event_id:
                continue

            try:
                chrom, ranges, strand = parse_coordinate_string(event_id)
            except Exception as e:
                # print(f"Skipping malformed event_id '{event_id}': {e}")
                skipped += 1
                continue

            if len(ranges) != 3:
                # print(f"Skipping {event_id}: Expected 3 coordinate blocks, found {len(ranges)}.")
                skipped += 1
                continue

            # ranges[0] = ExonA_loc, ranges[1] = Pseudo_loc, ranges[2] = ExonB_loc (in genome coordinates)
            
            r1_start, r1_end = ranges[0]
            r2_start, r2_end = ranges[1]
            r3_start, r3_end = ranges[2]
            
            # Define the 5 regions in Genomic Forward (+) coordinates
            # 1. Exon 1
            loc1 = (r1_start, r1_end)
            # 2. Intron 1 (inferred)
            loc2 = (r1_end + 1, r2_start - 1)
            # 3. Pseudoexon
            loc3 = (r2_start, r2_end)
            # 4. Intron 2 (inferred)
            loc4 = (r2_end + 1, r3_start - 1)
            # 5. Exon 2
            loc5 = (r3_start, r3_end)

            locs = [loc1, loc2, loc3, loc4, loc5]
            
            # Extract sequences (Always 5' to 3' relative to the strand)
            sequences = []
            
            try:
                if strand == '+':
                    # Straightforward 1-2-3-4-5
                    for s, e in locs:
                        seq = fasta.fetch_seq(chrom, s, e)
                        sequences.append(seq)
                else:
                    # Negative Strand Logic
                    # The genomic coordinates read left-to-right are 3' to 5' for this gene.
                    # We want the output 1-2-3-4-5 to represent 5'-ExonA...ExonB-3'.
                    # On (-) strand, the "last" genomic block (loc5) is Exon A (5').
                    # The "first" genomic block (loc1) is Exon B (3').
                    
                    # Reverse the list of locations to go from High Genomic Coords (5') to Low Genomic Coords (3')
                    # locs becomes [loc5, loc4, loc3, loc2, loc1]
                    reversed_locs = locs[::-1]
                    
                    for s, e in reversed_locs:
                        # Fetch from genome (+)
                        seq_plus = fasta.fetch_seq(chrom, s, e)
                        # Reverse Complement to get (-) strand sequence
                        seq_minus = reverse_complement(seq_plus)
                        sequences.append(seq_minus)

            except KeyError as e:
                # This catches the "Chromosome not found" error
                missing_chroms.add(chrom)
                skipped += 1
                continue
            except ValueError as e:
                # This catches coordinates out of bounds
                skipped += 1
                continue
            
            # Write row
            row = [event_id, dmso_psi, ris_psi] + sequences
            f_out.write("\t".join(row) + "\n")
            count += 1
            
            if count % 100 == 0:
                print(f"Processed {count} records...", end='\r')

    print(f"\nDone! Extracted {count} sequences to {OUTPUT_FILE}")
    if skipped > 0:
        print(f"Skipped {skipped} entries (malformed or missing chromosomes).")
    if missing_chroms:
        print(f"Warning: The following chromosomes were not found in your FASTA file: {list(missing_chroms)[:5]}...")
    fasta.close()

if __name__ == "__main__":
    main()