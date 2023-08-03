#!/usr/bin/env python3

'''
Note: I ran this script in the exp_vs_obs/ directory:
    python3 02_visualize_exp_vs_obs.py
'''

####################
## Import Modules ##
####################

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read in data
df = pd.read_csv('expected_vs_observed_peaks_master.csv')

#######################
## Overall Histogram ##
#######################

# Plot histogram 
plt.hist(df['Observed_Peaks'], bins=10, range=(0, 100), color='skyblue', edgecolor='black')

# Set x-axis tick positions and labels
plt.xticks(range(0, 101, 10)) 
plt.xlabel('Number of Observed Peaks')
plt.ylabel('Frequency')
plt.title('Total Distribution of Observed Peaks')
plt.show()

# Save figure
plt.savefig('tables_and_figures/total_distribution_of_observed_peaks.pdf')

##########################
## Table of Proportions ##
##########################

# Condense tests for each condition and calculate basic stats 
proportions_per_condition = df.groupby(["Endedness", "Peak_Type", "Aligner", "Peak_Caller", "Deduplicator", "Control"]).agg(
    Min_Observed_Peaks=("Observed_Peaks", "min"),
    Max_Observed_Peaks=("Observed_Peaks", "max"),
    Std_Observed_Peaks=("Observed_Peaks", "std"),
    Mean_Observed_Peaks=("Observed_Peaks", "mean")
).reset_index()

# Save to CSV
proportions_per_condition.to_csv("tables_and_figures/proportions_per_condition.csv", index=False)

##############
## Heatmaps ##
##############

'''
# Filter rows so peak caller is unique for peak type and endedness groups
heatmap_df = df.groupby(["Endedness", "Peak_Type"]).filter(lambda x: x["Peak_Caller"].nunique() == 1)

# Generate a heatmap for each unique combination of peak type, endedness, and peak caller
for (readtype, peaktype, peakcaller), data in heatmap_df.groupby(["Endedness", "Peak_Type", "Peak_Caller"]):
    
    # Plot size
    plt.figure(figsize=(10, 6))
    
    # Create a matrix for the heatmap
    pivot_df = data.pivot(index="Deduplicator", columns="Aligner", values="Observed_Peaks")
    sns.heatmap(pivot_df, cmap="coolwarm", annot=True, fmt=".2f", cbar_kws={'label': 'Observed Peaks'})
    
    # Format naming
    if peakcaller == "macs3":
        peakcaller_name = "MACS3"
    elif peakcaller == "genrich":
        peakcaller_name = "Genrich"
    elif peakcaller == "cisgenome":
        peakcaller_name = "CisGenome"
    elif peakcaller == "pepr":
        peakcaller_name = "PePr"
    
    # Configure titles and axis labels
    plt.title(f'Number of Peaks for {readtype.title()}-End Data with {peaktype.title()} Peaks using {peakcaller}')
    plt.xlabel("Aligner")
    plt.ylabel("Deduplicator")
    
    # Show figure
    plt.show()
    
    # Save figure
    plt.savefig(f'tables_and_figures/heatmap_{readtype}_{peaktype}_{peakcaller}.pdf')
'''