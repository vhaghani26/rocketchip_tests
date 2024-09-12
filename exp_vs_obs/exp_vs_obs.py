#!/usr/bin/env python3

'''
Note: I ran this script in the exp_vs_obs/ directory:
    python3 exp_vs_obs.py
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
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

###########################
## Set Working Variables ##
###########################

# User-specific variables
authors = 'Viktoria_Haghani_and_Aditi_Goyal_and_Alan_Zhang'
working_dir = '/share/korflab/home/viki/rocketchip_tests/exp_vs_obs' # Do NOT end the directory name with / here

# Combinatorial testing variables
controltypes = ["with_control", "no_control"]  
readtypes = ["paired", "single"]
peaktypes = ["narrow", "broad"]
aligners = ["bwa_mem", "bowtie2", "STAR"]
peakcallers = ["macs3", "cisgenome", "genrich", "pepr"]
deduplicators = ["samtools", "no_deduplication", "sambamba", "picard"]
num_tests = 6

# Create DataFrame for peak counting 
df = pd.DataFrame(columns=["Endedness", "Peak_Type", "Aligner", "Peak_Caller", "Deduplicator", "Test_Dataset", "Control", "Synthetic_Genome_Path", "Synthetic_Forward_Read_1_Path", "Synthetic_Reverse_Read_1_Path", "Synthetic_Forward_Read_2_Path", "Synthetic_Reverse_Read_2_Path", "Reads_per_Peak", "Padding", "Reads_STD_Dev", "Width", "Read_Length", "Paired", "Flank", "Expected_Peaks", "Observed_Peaks", "True_Positives", "True_Negatives", "False_Positives", "False_Negatives"])

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
                                    Project: exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}
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
                                    Project: exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}
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
                                print(f'Generating {working_dir}/project_files/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml...')
                                os.system(f'touch {working_dir}/project_files/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml')
                                with open(f'{working_dir}/project_files/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml', 'w') as f:
                                    f.write(f'{proj_file_info}')
                                os.system(f'sed -i \'1d\' {working_dir}/project_files/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml')

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
                                    snakefile_dir = f'{working_dir}/snakefiles/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}'
                                    if not os.path.exists(snakefile_dir):
                                        os.makedirs(snakefile_dir)
                                    
                                    # Create the snakefiles using Rocketchip
                                    print(f"Generating {snakefile_dir}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}")
                                    os.system(f'rocketchip {working_dir}/project_files/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml --data {working_dir}/seq_data/{readtype}_{peaktype}/test_{i} --output_file {snakefile_dir}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                    
                                    # MACS3 has an issue with small read files requiring the --nomodel flag, so I will manually add it for the single-end data that are having problems with peak-calling
                                    if readtype == "single" and peakcaller == "macs3":
                                        print(f'Adding --nomodel flag in MACS3 for {snakefile_dir}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                        file_to_open = f'{snakefile_dir}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}'
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
                                    os.chdir(f'{working_dir}/snakefiles/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                    os.system('pwd')
                                    # Run snakefile
                                    os.system(f'snakemake -j 4 -s exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                    # Go back to original directory
                                    os.chdir(f'../../../')

# Calculate true positive, true negative, false positive, and false negative peak regions
def get_stats(real_peaks, detected_peaks, genome_size):
	coors = defaultdict(lambda: defaultdict(str))

	for peak in real_peaks:
		coors[peak[0]]['real'] = 'start'
		coors[peak[1]]['real'] = 'end'

	for peak in detected_peaks:
		coors[peak[0]]['dtct'] = 'start'
		coors[peak[1]]['dtct'] = 'end'

	sorted(coors.items(), key=lambda x: x[0])

	poss = list(coors.keys())
	poss.sort()

	true_pos = 0
	true_neg = 0
	false_pos = 0
	false_neg = 0

	in_real_peak = False
	in_detected_peak = False

	previous = 1
	previous_start = True
	current_start = True
	reverse = False
	for i, pos in enumerate(poss):
		if len(coors[pos]) > 1:
			if coors[pos]['real'] == coors[pos]['dtct']:
				if coors[pos]['real'] == 'start': current_start = True
				else: current_start = False
			else:
				true_pos += 1
				current_start = True
				reverse = True
		elif 'real' in coors[pos]:
			if coors[pos]['real'] == 'start': current_start = True
			else: current_start = False
		else:
			if coors[pos]['dtct'] == 'start': current_start = True
			else: current_start = False
		
		if previous_start == current_start: buf = 0
		elif previous_start and not current_start: buf = 1
		else: buf = -1
		
		if not in_real_peak and not in_detected_peak:
			true_neg += pos - previous + buf
		if not in_real_peak and in_detected_peak:
			false_pos += pos - previous + buf
		if in_real_peak and not in_detected_peak:
			false_neg += pos - previous + buf
		if in_real_peak and in_detected_peak:
			true_pos += pos - previous + buf
		
		previous_start = current_start
		if reverse:
			previous_start = not previous_start
			reverse = False
			
		if len(coors[pos]) > 1:
			in_detected_peak = not in_detected_peak
			in_real_peak = not in_real_peak
		elif 'real' in coors[pos]:
			in_real_peak = not in_real_peak
		else:
			in_detected_peak = not in_detected_peak

		previous = pos 

	if poss[-1] < genome_size: true_neg += genome_size - poss[-1]
	
	return true_pos, true_neg, false_pos, false_neg

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
                                    genome_path = "NA"
                                    read_1_for_path = "NA"
                                    read_1_rev_path = "NA"
                                    read_2_for_path = "NA"
                                    read_2_rev_path = "NA"
                                    reads_per_peak = "NA"
                                    padding = "NA"
                                    reads_std_dev = "NA"
                                    width = "NA"
                                    length = "NA"
                                    paired = "NA"
                                    flank = "NA"
                                    expected_peaks = "NA"
                                    obs_peak_num = "NA"
                                    true_positives = "NA"
                                    true_negatives = "NA"
                                    false_positives = "NA"
                                    false_negatives = "NA"
                                    
                                    # Add null data for Cisgenome and Pepr to generate heatmap later
                                    df.loc[len(df)] = [readtype, peaktype, aligner, peakcaller, deduplicator, i, control, genome_path, read_1_for_path, read_1_rev_path, read_2_for_path, read_2_rev_path, reads_per_peak, padding, reads_std_dev, width, length, paired, flank, expected_peaks, obs_peak_num, true_positives, true_negatives, false_positives, false_negatives]
                                else:
                                    # Change into snakefile directory
                                    os.chdir(f'{working_dir}/snakefiles/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                    print(f'Counting peaks for exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}...')
                                    
                                    # Assign global data frame variables 
                                    genome_path = f'{working_dir}/seq_data/{readtype}_{peaktype}/test_{i}/genome.fa'
                                    
                                    # Determine the genome size
                                    with open(genome_path, "r") as fasta_file:
                                        # Read all lines in the file (excluding header lines)
                                        sequence = "".join(line.strip() for line in fasta_file if not line.startswith(">"))
                                    
                                    # Calculate the length of the sequence
                                    genome_size = len(sequence)
                                   
                                    # Assign global data frame variables 
                                    if readtype == "paired":
                                        read_1_for_path = f'{working_dir}/seq_data/{readtype}_{peaktype}/test_{i}/exp_a_1.fastq.gz'
                                        read_1_rev_path = f'{working_dir}/seq_data/{readtype}_{peaktype}/test_{i}/exp_a_2.fastq.gz'
                                        read_2_for_path = f'{working_dir}/seq_data/{readtype}_{peaktype}/test_{i}/exp_b_1.fastq.gz'
                                        read_2_rev_path = f'{working_dir}/seq_data/{readtype}_{peaktype}/test_{i}/exp_b_2.fastq.gz'
                                    elif readtype == "single":
                                        read_1_for_path = f'{working_dir}/seq_data/{readtype}_{peaktype}/test_{i}/exp_a.fastq.gz'
                                        read_1_rev_path = "NA"
                                        read_2_for_path = f'{working_dir}/seq_data/{readtype}_{peaktype}/test_{i}/exp_b.fastq.gz'
                                        read_2_rev_path = "NA"
                                    
                                    # Determine real peaks (note exp_a.peaks and exp_b.peaks are identical)
                                    real_peak_path = f'{working_dir}/seq_data/{readtype}_{peaktype}/test_{i}/exp_a.peaks'
                                    real_peaks = []
                                    with open(real_peak_path, "r") as file:
                                        # Read the file line by line
                                        for line in file:
                                            # Split the line into columns
                                            peak_summit = line.strip()
                                            peak_summit = int(peak_summit)
                                            # Access the values in the second and third columns
                                            if peaktype == "narrow":
                                                peak_length = 400
                                            elif peaktype == "broad":
                                                peak_length = 1500
                                            peak_range = (peak_summit - peak_length/2, peak_summit + peak_length/2)
                                            real_peaks.append(peak_range)
                                    
                                    # Assign sequence metadata
                                    reads_per_peak = 2**(i+1)
                                    expected_peaks = 1000
                                    reads_std_dev = 0.1
                                    length = 80
                                    flank = 2000 
                                    
                                    if peaktype == "narrow":
                                        padding = 2500
                                        width = 400
                                    elif peaktype == "broad":
                                        padding = 5000
                                        width = 1500
                                        
                                    if readtype == "paired":
                                        paired = 10
                                    elif readtype == "single":
                                        paired = "NA"
                                    
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
                                                    # Calculate true positives, true negatives, false positives, and false negatives
                                                    true_positives, true_negatives, false_positives, false_negatives = get_stats(real_peaks = real_peaks, detected_peaks = detected_peaks, genome_size = genome_size)
                                                else:
                                                    obs_peak_num = 0
                                                    true_positives = "NA"
                                                    true_negatives = "NA"
                                                    false_positives = "NA"
                                                    false_negatives = "NA"
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
                                                    # Calculate true positives, true negatives, false positives, and false negatives
                                                    true_positives, true_negatives, false_positives, false_negatives = get_stats(real_peaks = real_peaks, detected_peaks = detected_peaks, genome_size = genome_size)
                                                else:
                                                    obs_peak_num = 0
                                                    true_positives = "NA"
                                                    true_negatives = "NA"
                                                    false_positives = "NA"
                                                    false_negatives = "NA"
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
                                                    # Calculate true positives, true negatives, false positives, and false negatives
                                                    true_positives, true_negatives, false_positives, false_negatives = get_stats(real_peaks = real_peaks, detected_peaks = detected_peaks, genome_size = genome_size)
                                                else:
                                                    obs_peak_num = 0
                                                    true_positives = "NA"
                                                    true_negatives = "NA"
                                                    false_positives = "NA"
                                                    false_negatives = "NA"
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
                                                    # Calculate true positives, true negatives, false positives, and false negatives
                                                    true_positives, true_negatives, false_positives, false_negatives = get_stats(real_peaks = real_peaks, detected_peaks = detected_peaks, genome_size = genome_size)
                                                else:
                                                    obs_peak_num = 0
                                                    true_positives = "NA"
                                                    true_negatives = "NA"
                                                    false_positives = "NA"
                                                    false_negatives = "NA"
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
                                            # Calculate true positives, true negatives, false positives, and false negatives
                                            true_positives, true_negatives, false_positives, false_negatives = get_stats(real_peaks = real_peaks, detected_peaks = detected_peaks, genome_size = genome_size)
                                        else:
                                            obs_peak_num = 0
                                            true_positives = "NA"
                                            true_negatives = "NA"
                                            false_positives = "NA"
                                            false_negatives = "NA"
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
                                                # Calculate true positives, true negatives, false positives, and false negatives
                                                true_positives, true_negatives, false_positives, false_negatives = get_stats(real_peaks = real_peaks, detected_peaks = detected_peaks, genome_size = genome_size)
                                            else:
                                                obs_peak_num = 0
                                                true_positives = "NA"
                                                true_negatives = "NA"
                                                false_positives = "NA"
                                                false_negatives = "NA"
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
                                                # Calculate true positives, true negatives, false positives, and false negatives
                                                true_positives, true_negatives, false_positives, false_negatives = get_stats(real_peaks = real_peaks, detected_peaks = detected_peaks, genome_size = genome_size)
                                            else:
                                                obs_peak_num = 0
                                                true_positives = "NA"
                                                true_negatives = "NA"
                                                false_positives = "NA"
                                                false_negatives = "NA"
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
                                            # Calculate true positives, true negatives, false positives, and false negatives
                                            true_positives, true_negatives, false_positives, false_negatives = get_stats(real_peaks = real_peaks, detected_peaks = detected_peaks, genome_size = genome_size)
                                        else:
                                            obs_peak_num = 0
                                            true_positives = "NA"
                                            true_negatives = "NA"
                                            false_positives = "NA"
                                            false_negatives = "NA"
                                        os.chdir('..')
                                    
                                    # View peak number
                                    print(f'Observed peaks: {obs_peak_num}')
                                    
                                    # Go back to original directory
                                    os.chdir(f'../../../')
    
                                    # Add test to dataframe 
                                    df.loc[len(df)] = [readtype, peaktype, aligner, peakcaller, deduplicator, i, control, genome_path, read_1_for_path, read_1_rev_path, read_2_for_path, read_2_rev_path, reads_per_peak, padding, reads_std_dev, width, length, paired, flank, expected_peaks, obs_peak_num, true_positives, true_negatives, false_positives, false_negatives]
                                
                                # Save to CSV
                                df.to_csv(output_path, index=False)

# Compute sensitivity, precision, and F1 scores
def calculate_stats(dataframe):
    # Read in data 
    df = pd.read_csv(dataframe)
    
    # Calculate sensitivity [TP/(TP + FN)]
    df['Sensitivity'] = df['True_Positives'] / (df['True_Positives'] + df['False_Negatives'])
    
    # Calculate precision [TP/(TP + FP)]
    df['Precision'] = df['True_Positives'] / (df['True_Positives'] + df['False_Positives'])
    
    # Calculate F1 score [2TP/(2TP + FP + FN)]
    df['F1_Score'] = (2 * df['True_Positives']) / (2 * df['True_Positives'] + df['False_Positives'] + df['False_Negatives'])
    
    # Fill NA values with "NA"
    df = df.fillna("NA")
    
    # Save to CSV (overwrites previously saved file)
    df.to_csv(dataframe, index=False)

# Generate histogram
def observed_peaks_histogram(dataframe, output_file):
    # Read in data 
    df = pd.read_csv(dataframe)

    # Drop rows with NA values in the 'Observed_Peaks' column
    df = df.dropna(subset=['Observed_Peaks'])
    
    # Find the maximum value of observed peaks in order to adjust visualization
    max_observed_peaks = df['Observed_Peaks'].max()
    print(f'{max_observed_peaks} was the highest observed peak number.')
    
    # Round range up to the nearest multiple of 10
    max_range = ((max_observed_peaks + 9) // 10) * 10
    
    # Divide by number into every 10
    bin_num = int(max_range/10)
    
    # Plot histogram 
    plt.hist(df['Observed_Peaks'], bins = bin_num, range = (0, max_range), color = 'skyblue', edgecolor = 'black')
    
    # Set x-axis tick positions and labels
    plt.xlabel('Number of Observed Peaks')
    plt.ylabel('Frequency')
    plt.title('Total Distribution of Observed Peaks')
    plt.show()
    
    # Save figure
    plt.savefig(output_file)
    plt.close()

# Linear regression to compare reads per peak vs. observed number of peaks
def reads_per_peak_vs_obs_peaks(dataframe, output_file):
    # Read in data 
    df = pd.read_csv(dataframe)
    
    # Drop rows with NA values in the 'Observed_Peaks' column
    df = df.dropna(subset=['Observed_Peaks'])

    # Perform linear regression
    fit = np.polyfit(df['Reads_per_Peak'], df['Observed_Peaks'], 1)
    fit_fn = np.poly1d(fit)
    
    # Create a scatter plot
    plt.scatter(df['Reads_per_Peak'], df['Observed_Peaks'], color = 'blue', marker = 'o')
    
    # Add best-fit line
    plt.plot(df['Reads_per_Peak'], fit_fn(df['Reads_per_Peak']), color='red') 
    
    # Add labels
    plt.xlabel('Reads per Peak')
    plt.ylabel('Observed Peaks')
    plt.title('Read Coverage vs. Observed Peaks')
    plt.grid(False)
    
    # Save figure
    plt.savefig(output_file)
    plt.close()

# Linear regression to compare reads per peak vs. F1 scores
def reads_per_peak_vs_F1_score(dataframe, output_file):
    # Read in data 
    df = pd.read_csv(dataframe)
    
    # Drop rows with NA values in the 'Reads_per_Peak' column
    df = df.dropna(subset=['Reads_per_Peak'])

    # Perform linear regression
    fit = np.polyfit(df['Reads_per_Peak'], df['F1_Score'], 1)
    fit_fn = np.poly1d(fit)
    
    # Create a scatter plot
    plt.scatter(df['Reads_per_Peak'], df['F1_Score'], color = 'blue', marker = 'o')
    
    # Add best-fit line
    plt.plot(df['Reads_per_Peak'], fit_fn(df['Reads_per_Peak']), color='red') 
    
    # Add labels
    plt.xlabel('Reads per Peak')
    plt.ylabel('F1 Score')
    plt.title('Read Coverage vs. F1 Score')
    plt.grid(False)
    
    # Save figure
    plt.savefig(output_file)
    plt.close()
    
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
count_peaks(working_dir = working_dir, controltypes = controltypes, readtypes = readtypes, peaktypes = peaktypes, aligners = aligners, peakcallers = peakcallers, deduplicators = deduplicators, num_tests = num_tests, output_path = 'tables_and_figures/expected_vs_observed_peaks_master.csv')

# Compute sensitivity, precision, and F1 scores
calculate_stats(dataframe = 'tables_and_figures/expected_vs_observed_peaks_master.csv')

# Generate histogram
observed_peaks_histogram(dataframe = 'tables_and_figures/expected_vs_observed_peaks_master.csv', output_file = 'tables_and_figures/total_distribution_of_observed_peaks.pdf')

# Linear regression to compare reads per peak vs. observed number of peaks
reads_per_peak_vs_obs_peaks(dataframe = 'tables_and_figures/expected_vs_observed_peaks_master.csv', output_file = 'tables_and_figures/reads_per_peak_vs_observed_peaks.pdf')

# Linear regression to compare reads per peak vs. F1 scores
reads_per_peak_vs_F1_score(dataframe = 'tables_and_figures/expected_vs_observed_peaks_master.csv', output_file = 'tables_and_figures/reads_per_peak_vs_f1_score.pdf')