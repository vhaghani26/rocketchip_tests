#!/usr/bin/env python3

'''
Usage:
python3 peak_count_csv.py --in_dir peak_counts/ --out_dir .
'''

####################
## Import Modules ##
####################

import argparse
import os
import pandas as pd
import re

#####################
## Set Up Argparse ##
#####################

# Initialize argparse
parser = argparse.ArgumentParser(
    description='Count peaks using peak output directory')

parser.add_argument('--in_dir', required=True, type=str,
    metavar='<str>', help='Directory containing the files with the peak data in it')
    
parser.add_argument('--out_dir', required=True, type=str,
    metavar='<str>', help='Directory to output peak count CSV file')
    
# Finalization of argparse
arg = parser.parse_args()

#################
## Count Peaks ##
#################

# Initialize list for dataframe
data = []

# Iterate over each file in the directory
for filename in os.listdir(arg.in_dir):
    # Construct the full file path
    filepath = os.path.join(arg.in_dir, filename)
    # Check if it is a file (not a directory)
    if os.path.isfile(filepath):
        # Open the file and read its content
        with open(filepath, 'r') as file:
            value = file.read().strip()
            # Append the filename and value to the data list
            data.append((filename, value))

# Convert the data list to a dataframe
df = pd.DataFrame(data, columns=['Filename', 'Peak_Counts'])

# Sort alphabetically by file name
df = df.sort_values(by='Filename')

# Separate variables into columns
pattern = r'(?P<Project>[^_]+)_(?P<Control>no_control|with_control)_(?P<Peak_Type>[^_]+)_(?P<Aligner>bowtie2|STAR|bwa_mem)_(?P<Peak_Caller>[^_]+)_(?P<Deduplicator>no_deduplication|picard|samtools|sambamba)_(?P<Test>[^_]+)_(?P<Peaks>.+)'

df = df.join(df['Filename'].str.extract(pattern))

# Save the dataframe
df.to_csv(f'{arg.out_dir}/peak_counts.csv', index=False)