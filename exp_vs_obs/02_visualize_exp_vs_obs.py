#########################
## Table of Statistics ##
#########################

# Condense tests for each condition and calculate basic stats
proportions_per_condition = df.groupby(["Endedness", "Peak_Type", "Aligner", "Peak_Caller", "Deduplicator", "Control"]).agg(
    Min_Observed_Peaks = ("Observed_Peaks", "min"),
    Max_Observed_Peaks = ("Observed_Peaks", "max"),
    Std_Observed_Peaks = ("Observed_Peaks", "std"),
    Mean_Observed_Peaks = ("Observed_Peaks", "mean"),
    Percent_Accuracy = ("Observed_Peaks", lambda x: (x == 50).sum() / len(x) * 100),
).reset_index()

# Transpose "Test_Dataset" values into columns for new dataframe
pivot_df = df.pivot_table(index=["Endedness", "Peak_Type", "Aligner", "Peak_Caller", "Deduplicator", "Control"],
                          columns="Test_Dataset",
                          values="Observed_Peaks",
                          aggfunc='first')

# Rename the columns to add the prefix "Test_Data_"
pivot_df.columns = [f"Test_Data_{col}" for col in pivot_df.columns]

# Merge 'proportions_per_condition' with the transposed dataframe
merged_df = pd.merge(proportions_per_condition, pivot_df, on=["Endedness", "Peak_Type", "Aligner", "Peak_Caller", "Deduplicator", "Control"])

# Add column
merged_df["Expected_Peaks"] = 50

# Add true_positives, true_negatives, false_positives, and false_negatives
print(merged_df)


# Calculate sensitivity



# Specificity

# Precision

# Calculate read density









'''
# Delineate column order
desired_order = ["Endedness", "Peak_Type", "Aligner", "Peak_Caller", "Deduplicator", "Control", "Expected_Peaks"] + \
                [f"Test_Data_{i}" for i in range(1, 7)] + \
                ["True_Positives", "True_Negatives", "False_Positives", "False_Negatives", "Min_Observed_Peaks", "Max_Observed_Peaks", "Std_Observed_Peaks", "Mean_Observed_Peaks", "Sensitivity", "Specificity", "Precision", "Percent_Accuracy"]

# Reorder columns 
merged_df = merged_df.reindex(columns=desired_order)

# Sort by percent accuracy 
merged_df = merged_df.sort_values(by="Percent_Accuracy", ascending = False)

# Save to CSV
merged_df.to_csv("02_tables_and_figures/02_expected_vs_observed_results.csv", index = False)
'''