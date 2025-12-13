import pandas as pd

FILE = "risdiplam_full_dataset.csv.gz"

print(f"Reading first 5 rows of {FILE}...")

# Pandas reads gzip automatically
df = pd.read_csv(FILE, nrows=5)

# Set pandas to display full sequence width for inspection
pd.set_option('display.max_colwidth', 50) 

# Display
print(df[['event_id', 'dPSI', 'seq_ExonA', 'seq_Intron1']].to_string())

print("\n(Sequences truncated for display, but full length in file)")