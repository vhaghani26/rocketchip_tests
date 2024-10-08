#!/bin/bash python3 

'''
python3 deep_dive_yaml.py
'''

import yaml

# Variables
controltypes = ["with_control", "no_control"]
projects = ["Rube"]
peaktypes = ["narrow", "broad"]
aligners = ["bwa_mem", "bowtie2", "STAR"]
peakcallers = ["macs3", "cisgenome", "genrich", "pepr"]
deduplicators = ["samtools", "no_deduplication", "sambamba", "picard"]
num_tests = 100

# Store sample names
samples = []

# Generate sample names, skipping invalid combinations
for control in controltypes:
    for project in projects:
        for peaktype in peaktypes:
            for aligner in aligners:
                for peakcaller in peakcallers:
                    for deduplicator in deduplicators:
                        for i in range(1, num_tests + 1):

                            # Skip illegal combinations
                            if control == "no_control" and (peakcaller == "cisgenome" or peakcaller == "pepr"):
                                continue

                            # Generate valid combinations
                            sample_name = f"{project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}"
                            samples.append(sample_name)

# Structure the data for YAML
data = {
    'samples': samples
}

# Write to snakefiles.yaml
with open('snakefiles.yaml', 'w') as file:
    yaml.dump(data, file, default_flow_style=False)

print("snakefiles.yaml has been created with valid combinations!")
