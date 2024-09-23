#!/bin/bash python3

'''
python3 anova.py
'''

import pandas as pd
from scipy import stats

# Load data
df = pd.read_csv('new_peak_counts.csv')

# Define the columns representing software combinations
software_columns = ['Project', 'Control', 'Peak_Type', 'Aligner', 'Peak_Caller', 'Deduplicator']

# Pivot the data so that each software combination is a unique row
df_pivot = df.pivot_table(index=software_columns, columns='Test', values='Peak_Counts').reset_index()

# Prepare empty lists to store the ANOVA results
f_statistics = []
p_values = []

# Loop over each row in the pivoted DataFrame (each software combination)
for index, row in df_pivot.iterrows():
    # Extract the peak counts for test1 to test100
    peak_counts = row.filter(like='test').values
    
    # Ensure that we have more than one test (trial) to perform ANOVA
    if len(peak_counts) > 1 and not all(pd.isna(peak_counts)):  # Check for NaNs
        # Remove NaN values before performing ANOVA
        valid_peak_counts = peak_counts[~pd.isna(peak_counts)]
        
        # Check if the valid peak counts are constant
        if len(valid_peak_counts) > 1 and len(set(valid_peak_counts)) == 1:
            # All values are the same, cannot perform ANOVA
            f_statistics.append(None)
            p_values.append(None)
        else:
            if len(valid_peak_counts) > 1:  # Ensure we have more than one group
                anova_result = stats.f_oneway(*[valid_peak_counts for _ in range(len(valid_peak_counts))])
                # Store the F-statistic and p-value
                f_statistics.append(anova_result.statistic)
                p_values.append(anova_result.pvalue)
            else:
                # If we don't have enough valid peak counts, append None
                f_statistics.append(None)
                p_values.append(None)
    else:
        # If there are not enough trials, append None
        f_statistics.append(None)
        p_values.append(None)

# Add the F-statistic and p-value as new columns in the pivoted DataFrame
df_pivot['ANOVA_F-statistic'] = f_statistics
df_pivot['ANOVA_p-value'] = p_values

# Save the results to a new CSV file
df_pivot.to_csv('anova_results.csv', index=False)

print("ANOVA results saved to anova_results.csv")
