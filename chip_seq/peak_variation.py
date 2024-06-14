import pandas as pd

# Load data
df = pd.read_csv("peak_counts.csv")

# Check for variations in Peak_Counts for each combination
variations = df.groupby(['Project', 'Control', 'Peak_Type', 'Aligner', 'Peak_Caller', 'Deduplicator']).agg(
    Peak_Counts_Unique=('Peak_Counts', 'nunique'),
    Peak_Counts_Range=('Peak_Counts', lambda x: f"{x.min()}-{x.max()}"),
    Range_in_Called_Peaks=('Peak_Counts', lambda x: x.max() - x.min())
)

# Total number of combinations
total_combinations = variations.shape[0]

# Combinations with variation in Peak_Counts
combinations_with_variation = variations[variations['Peak_Counts_Unique'] > 1]

# Number of combinations with variability
num_combinations_with_variation = combinations_with_variation.shape[0]

# Number of combinations with 100% reproducibility
num_combinations_reproducible = total_combinations - num_combinations_with_variation

# Percentages
percent_with_variation = (num_combinations_with_variation / total_combinations) * 100
percent_reproducible = (num_combinations_reproducible / total_combinations) * 100

# Print the results
print(f"Total number of combinations: {total_combinations}")
print(f"Number of combinations with variability: {num_combinations_with_variation} ({percent_with_variation:.2f}%)")
print(f"Number of combinations with 100% reproducibility: {num_combinations_reproducible} ({percent_reproducible:.2f}%)")
print("\nCombinations with variation in Peak_Counts:")
print(combinations_with_variation)