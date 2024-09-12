#!/usr/bin/env python3

'''
Note: I ran this script in the functional_test/ directory:
    python3 functional_test.py
'''

####################
## Import Modules ##
####################

import sys
from collections import defaultdict
import pandas as pd
import os
import textwrap
import subprocess 

###########################
## Set Working Variables ##
###########################

# User-specific variables
authors = 'Viktoria_Haghani'
working_dir = '/share/korflab/home/viki/rocketchip_tests/functional_test' # Do NOT end the directory name with / here
output_path = 'peaks_counts_before_updates.csv'

# Combinatorial testing variables
controltypes = ["with_control", "no_control"]  
readtypes = ["paired", "single"]
peaktypes = ["narrow", "broad"]
aligners = ["bwa_mem", "bowtie2", "STAR"]
peakcallers = ["macs3", "cisgenome", "genrich", "pepr"]
deduplicators = ["samtools", "no_deduplication", "sambamba", "picard"]
num_tests = 6

# Create DataFrame for peak counting 
df = pd.DataFrame(columns=["Endedness", "Peak_Type", "Aligner", "Peak_Caller", "Deduplicator", "Test_Dataset", "Control", "Observed_Peaks"])

#########################
## Delineate Functions ##
#########################

# Make project files 
def generate_project_files(authors, working_dir, controltypes, readtypes, peaktypes, aligners, peakcallers, deduplicators, num_tests):
    # Set up directory structure for project files if needed
    if not os.path.exists(f'{working_dir}/project_files/'):
        print(f'Directory {working_dir}/project_files/ not found. Creating {working_dir}/project_files/')
        os.system(f'mkdir {working_dir}/project_files')
    if not os.path.exists(f'{working_dir}/project_files/with_control/'):
        print(f'Directory {working_dir}/project_files/with_control/ not found. Creating {working_dir}/project_files/with_control/')
        os.system(f'mkdir {working_dir}/project_files/with_control')    
    if not os.path.exists(f'{working_dir}/project_files/no_control/'):
        print(f'Directory {working_dir}/project_files/no_control/ not found. Creating {working_dir}/project_files/no_control/')
        os.system(f'mkdir {working_dir}/project_files/no_control')  
    
    # Start combinatorial project file generation
    for control in controltypes:
        for readtype in readtypes:
            for peaktype in peaktypes:
                for aligner in aligners:
                    for peakcaller in peakcallers:
                        for deduplicator in deduplicators:
                            for i in range(1, num_tests + 1):
                                if control == "with_control":
                                    proj_file_info = textwrap.dedent(f"""
                                    Author: {authors}
                                    Project: functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}
                                    Genome:
                                        Name: genome
                                        Location: '{working_dir}/seq_data/{readtype}_{peaktype}/test_{i}/genome.fa'
                                    Reads:
                                        Samples:
                                            grp1: 
                                                - '{working_dir}/seq_data/{readtype}_{peaktype}/test_{i}/exp_a'
                                                - '{working_dir}/seq_data/{readtype}_{peaktype}/test_{i}/exp_b'
                                        Controls:
                                            ctl1: 
                                                - '{working_dir}/seq_data/{readtype}_{peaktype}/test_{i}/input'
                                    Readtype: {readtype}
                                    Peaktype: {peaktype}
                                    Aligner: {aligner}
                                    Deduplicator: {deduplicator}
                                    Peakcaller: {peakcaller}
                                    Threads: 1
                                    """)
                                elif control == "no_control":
                                    if peakcaller == "cisgenome" or peakcaller == "pepr": continue
                                    proj_file_info = textwrap.dedent(f"""
                                    Author: {authors}
                                    Project: functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}
                                    Genome:
                                        Name: genome
                                        Location: '{working_dir}/seq_data/{readtype}_{peaktype}/test_{i}/genome.fa'
                                    Reads:
                                        Samples:
                                            grp1: 
                                                - '{working_dir}/seq_data/{readtype}_{peaktype}/test_{i}/exp_a'
                                                - '{working_dir}/seq_data/{readtype}_{peaktype}/test_{i}/exp_b'
                                        Controls:
                                    Readtype: {readtype}
                                    Peaktype: {peaktype}
                                    Aligner: {aligner}
                                    Deduplicator: {deduplicator}
                                    Peakcaller: {peakcaller}
                                    Threads: 1
                                    """)
                                print(f'Generating {working_dir}/project_files/{control}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml...')
                                os.system(f'touch {working_dir}/project_files/{control}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml')
                                with open(f'{working_dir}/project_files/{control}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml', 'w') as f:
                                    f.write(f'{proj_file_info}')
                                os.system(f'sed -i \'1d\' {working_dir}/project_files/{control}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml')

# Make Snakefiles
def generate_snakefiles(working_dir, controltypes, readtypes, peaktypes, aligners, peakcallers, deduplicators, num_tests):
    # Set up directory structure for Snakefiles if needed
    if not os.path.exists(f'{working_dir}/snakefiles/'):
        print(f'Directory {working_dir}/snakefiles/ not found. Creating {working_dir}/snakefiles/')
        os.system(f'mkdir {working_dir}/snakefiles')
    if not os.path.exists(f'{working_dir}/snakefiles/with_control/'):
        print(f'Directory {working_dir}/snakefiles/with_control/ not found. Creating {working_dir}/snakefiles/with_control/')
        os.system(f'mkdir {working_dir}/snakefiles/with_control')    
    if not os.path.exists(f'{working_dir}/snakefiles/no_control/'):
        print(f'Directory {working_dir}/snakefiles/no_control/ not found. Creating {working_dir}/snakefiles/no_control/')
        os.system(f'mkdir {working_dir}/snakefiles/no_control') 

    # Start combinatorial Snakefile generation
    for control in controltypes:
        for readtype in readtypes:
            for peaktype in peaktypes:
                for aligner in aligners:
                    for peakcaller in peakcallers:
                        for deduplicator in deduplicators:
                            for i in range(1, num_tests + 1):
                                if (control == "no_control") and (peakcaller == "cisgenome" or peakcaller == "pepr"):
                                    continue
                                else:
                                    # Make the directory structure if it does not already exist
                                    snakefile_dir = f'{working_dir}/snakefiles/{control}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}'
                                    if not os.path.exists(snakefile_dir):
                                        os.makedirs(snakefile_dir)
                                    
                                    # Create the snakefiles using Rocketchip
                                    print(f"Generating {snakefile_dir}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}")
                                    os.system(f'rocketchip {working_dir}/project_files/{control}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml --data {working_dir}/seq_data/{readtype}_{peaktype}/test_{i} --output_file {snakefile_dir}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                    
                                    # MACS3 has an issue with small read files requiring the --nomodel flag, so I will manually add it for the single-end data that are having problems with peak-calling
                                    if readtype == "single" and peakcaller == "macs3":
                                        print(f'Adding --nomodel flag in MACS3 for {snakefile_dir}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                        file_to_open = f'{snakefile_dir}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}'
                                        # Read in the file
                                        with open(file_to_open, 'r') as file:
                                            filedata = file.read()
                                        # Replace the target string
                                        filedata = filedata.replace('macs3 callpeak ', 'macs3 callpeak --nomodel ')
                                        # Write the file out again
                                        with open(file_to_open, 'w') as file:
                                            file.write(filedata)

# Run Snakefiles
def run_snakefiles(working_dir, controltypes, readtypes, peaktypes, aligners, peakcallers, deduplicators, num_tests):
    # Run Snakemake for all combinations
    for control in controltypes:
        for readtype in readtypes:
            for peaktype in peaktypes:
                for aligner in aligners:
                    for peakcaller in peakcallers:
                        for deduplicator in deduplicators:
                            for i in range(1, num_tests + 1):
                                if (control == "no_control") and (peakcaller == "cisgenome" or peakcaller == "pepr"):
                                    continue
                                else:
                                    # Change into snakefile directory
                                    os.chdir(f'{working_dir}/snakefiles/{control}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                    os.system('pwd')
                                    # Run snakefile
                                    os.system(f'snakemake -j 4 -s functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                    # Go back to original directory
                                    os.chdir(f'../../../')

# Count peaks and compare called vs. true peaks
def count_peaks(working_dir, controltypes, readtypes, peaktypes, aligners, peakcallers, deduplicators, num_tests, output_path):
    for control in controltypes:
        for readtype in readtypes:
            for peaktype in peaktypes:
                for aligner in aligners:
                    for peakcaller in peakcallers:
                        for deduplicator in deduplicators:
                            for i in range(1, num_tests + 1):                            
                                # Handle illegal combinations 
                                if (control == "no_control") and (peakcaller == "cisgenome" or peakcaller == "pepr"):                                                          
                                    # Assign variables
                                    obs_peak_num = "NA"
                                    
                                    # Add null data for Cisgenome and Pepr to generate heatmap later
                                    df.loc[len(df)] = [readtype, peaktype, aligner, peakcaller, deduplicator, i, control, obs_peak_num]
                                
                                else:
                                    # Change into snakefile directory
                                    os.chdir(f'{working_dir}/snakefiles/{control}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                    print(f'Counting peaks for functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}...')
                                    
                                    # Count peaks and determine peak locations for MACS3 outputs
                                    if peakcaller == "macs3":
                                        os.chdir('06_macs3_peaks')
                                        if control == "with_control":
                                            if peaktype == "narrow":
                                                if os.path.isfile('grp1_ctl1_peaks.narrowPeak'):
                                                    # Determine peak locations
                                                    detected_peaks = []
                                                    with open('grp1_ctl1_peaks.narrowPeak', 'r') as file:
                                                        # Read the file line by line
                                                        for line in file:
                                                            # Split the line into columns
                                                            columns = line.strip().split('\t')
                                                            # Access the values in the second and third columns
                                                            peak_range = (int(columns[1]), int(columns[2]))
                                                            detected_peaks.append(peak_range)
                                                    # Count peaks
                                                    obs_peak_num = len(detected_peaks)
                                                else:
                                                    obs_peak_num = 0

                                            elif peaktype == "broad":
                                                if os.path.isfile('grp1_ctl1_peaks.broadPeak'):
                                                    # Determine peak locations
                                                    detected_peaks = []
                                                    with open('grp1_ctl1_peaks.broadPeak', 'r') as file:
                                                        # Read the file line by line
                                                        for line in file:
                                                            # Split the line into columns
                                                            columns = line.strip().split('\t')
                                                            # Access the values in the second and third columns
                                                            peak_range = (int(columns[1]), int(columns[2]))
                                                            detected_peaks.append(peak_range)
                                                    # Count peaks
                                                    obs_peak_num = len(detected_peaks)
                                                else:
                                                    obs_peak_num = 0
                                        elif control == "no_control":
                                            if peaktype == "narrow":
                                                if os.path.isfile('grp1_peaks.narrowPeak'):
                                                    # Determine peak locations
                                                    detected_peaks = []
                                                    with open('grp1_peaks.narrowPeak', 'r') as file:
                                                        # Read the file line by line
                                                        for line in file:
                                                            # Split the line into columns
                                                            columns = line.strip().split('\t')
                                                            # Access the values in the second and third columns
                                                            peak_range = (int(columns[1]), int(columns[2]))
                                                            detected_peaks.append(peak_range)
                                                    # Count peaks
                                                    obs_peak_num = len(detected_peaks)
                                                else:
                                                    obs_peak_num = 0
                                            elif peaktype == "broad":
                                                if os.path.isfile('grp1_peaks.broadPeak'):
                                                    # Determine peak locations
                                                    detected_peaks = []
                                                    with open('grp1_peaks.broadPeak', 'r') as file:
                                                        # Read the file line by line
                                                        for line in file:
                                                            # Split the line into columns
                                                            columns = line.strip().split('\t')
                                                            # Access the values in the second and third columns
                                                            peak_range = (int(columns[1]), int(columns[2]))
                                                            detected_peaks.append(peak_range)
                                                    # Count peaks
                                                    obs_peak_num = len(detected_peaks)
                                                else:
                                                    obs_peak_num = 0
                                        os.chdir('..')
                                    
                                    # Count peaks and determine peak locations for Cisgenome outputs
                                    elif peakcaller == "cisgenome":
                                        os.chdir('06_cisgenome_peaks')
                                        if os.path.isfile('grp1_ctl1_peak.cod'):
                                            # Determine peak locations
                                            detected_peaks = []
                                            # Open the file for reading
                                            with open("grp1_ctl1_peak.cod", "r") as file:
                                                # Skip the header line
                                                next(file)
                                                # Read the file line by line
                                                for line in file:
                                                    # Split the line into columns
                                                    columns = line.strip().split("\t")
                                                    # Access the values in the third and fourth columns
                                                    peak_range = (int(columns[2]), int(columns[3]))
                                                    detected_peaks.append(peak_range)
                                                    
                                            # Count peaks
                                            obs_peak_num = len(detected_peaks)
                                        else:
                                            obs_peak_num = 0
                                        os.chdir('..')

                                    # Count peaks and determine peak locations for Genrich outputs
                                    elif peakcaller == "genrich":
                                        os.chdir('06_genrich_peaks')
                                        if control == "with_control":
                                            if os.path.isfile('grp1_ctl1_peak.narrowPeak'):
                                                # Determine peak locations
                                                detected_peaks = []
                                                with open('grp1_ctl1_peak.narrowPeak', 'r') as file:
                                                    # Read the file line by line
                                                    for line in file:
                                                        # Split the line into columns
                                                        columns = line.strip().split('\t')
                                                        # Access the values in the second and third columns
                                                        peak_range = (int(columns[1]), int(columns[2]))
                                                        detected_peaks.append(peak_range)
                                                # Count peaks
                                                obs_peak_num = len(detected_peaks)
                                            else:
                                                obs_peak_num = 0
                                        elif control == "no_control":
                                            if os.path.isfile('grp1_peak.narrowPeak'):
                                                # Determine peak locations
                                                detected_peaks = []
                                                with open('grp1_peak.narrowPeak', 'r') as file:
                                                    # Read the file line by line
                                                    for line in file:
                                                        # Split the line into columns
                                                        columns = line.strip().split('\t')
                                                        # Access the values in the second and third columns
                                                        peak_range = (int(columns[1]), int(columns[2]))
                                                        detected_peaks.append(peak_range)
                                                # Count peaks
                                                obs_peak_num = len(detected_peaks)
                                            else:
                                                obs_peak_num = 0
                                        os.chdir('..')
                                    
                                    # Count peaks and determine peak locations for PePr outputs
                                    elif peakcaller == "pepr":
                                        os.chdir('06_pepr_peaks')
                                        if os.path.isfile('grp1_ctl1__PePr_peaks.bed'):
                                            # Determine peak locations
                                            detected_peaks = []
                                            with open('grp1_ctl1__PePr_peaks.bed', 'r') as file:
                                                # Read the file line by line
                                                for line in file:
                                                    # Split the line into columns
                                                    columns = line.strip().split('\t')
                                                    # Access the values in the second and third columns
                                                    peak_range = (int(columns[1]), int(columns[2]))
                                                    detected_peaks.append(peak_range)
                                            # Count peaks
                                            obs_peak_num = len(detected_peaks)
                                        else:
                                            obs_peak_num = 0
                                        os.chdir('..')
                                    
                                    # View peak number
                                    print(f'Observed peaks: {obs_peak_num}')
                                    
                                    # Go back to original directory
                                    os.chdir(f'../../../')
    
                                    # Add test to dataframe 
                                    df.loc[len(df)] = [readtype, peaktype, aligner, peakcaller, deduplicator, i, control, obs_peak_num]
                                
                                # Save to CSV
                                df.to_csv(output_path, index=False)
    
####################
## Run Everything ##
####################

# Generate project_files
generate_project_files(authors = authors, working_dir = working_dir, controltypes = controltypes, readtypes = readtypes, peaktypes = peaktypes, aligners = aligners, peakcallers = peakcallers, deduplicators = deduplicators, num_tests = num_tests)

# Generate Snakefiles
generate_snakefiles(working_dir = working_dir, controltypes = controltypes, readtypes = readtypes, peaktypes = peaktypes, aligners = aligners, peakcallers = peakcallers, deduplicators = deduplicators, num_tests = num_tests)

# Run Snakefiles
run_snakefiles(working_dir = working_dir, controltypes = controltypes, readtypes = readtypes, peaktypes = peaktypes, aligners = aligners, peakcallers = peakcallers, deduplicators = deduplicators, num_tests = num_tests)

# Count peaks and calculate statistics based on peak calling
count_peaks(working_dir = working_dir, controltypes = controltypes, readtypes = readtypes, peaktypes = peaktypes, aligners = aligners, peakcallers = peakcallers, deduplicators = deduplicators, num_tests = num_tests, output_path = output_path)