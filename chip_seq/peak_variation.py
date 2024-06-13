import pandas as pd

# Load data
df = pd.read_csv("peak_counts.csv")

# Check for variations in Peak_Counts for each combination
variations = df.groupby(['Project', 'Control', 'Peak_Type', 'Aligner', 'Peak_Caller', 'Deduplicator']).agg(
    Peak_Counts_Unique=('Peak_Counts', 'nunique'),
    Peak_Counts_Range=('Peak_Counts', lambda x: f"{x.min()}-{x.max()}"),
    Range_in_Called_Peaks=('Peak_Counts', lambda x: x.max() - x.min())
)

# Filter combinations with variation in Peak_Counts
combinations_with_variation = variations[variations['Peak_Counts_Unique'] > 1]

# Print the combinations with their peak counts range and actual range
print(combinations_with_variation)
