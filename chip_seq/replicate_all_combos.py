#!/usr/bin/env python3

'''
python3 replicate_all_combos.py
'''

####################
## Import Modules ##
####################

import sys
import pandas as pd
import os
import textwrap
import subprocess 

###########################
## Set Working Variables ##
###########################

# User-specific variables
authors = 'Viktoria_Haghani'
working_dir = '/share/korflab/home/viki/rocketchip_tests/chip_seq' # Do NOT end the directory name with / here

# Combinatorial testing variables
controltypes = ["with_control", "no_control"] 
projects = ["Rube", "Namani"] 
peaktypes = ["narrow", "broad"]
aligners = ["bwa_mem", "bowtie2", "STAR"]
peakcallers = ["macs3", "cisgenome", "genrich", "pepr"]
deduplicators = ["samtools", "no_deduplication", "sambamba", "picard"]
num_tests = 3

# Create DataFrame for peak counting 
df = pd.DataFrame(columns=["Endedness", "Project", "Peak_Type", "Aligner", "Peak_Caller", "Deduplicator", "Sample_ID", "Control_ID", "Trial_1_Peaks", "Trial_2_Peaks", "Trial_3_Peaks"])

#########################
## Delineate Functions ##
#########################

def download_genome(genome, out_dir):
    genome_links = {
        'dm6': 'https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz',
        'hg19': 'https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz',
        'hg38': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz',
        'mm9': 'https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz',
        'mm10': 'https://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz',
        'rn6': 'https://hgdownload.cse.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz',
        'ce11': 'https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/chromFa.tar.gz',
        'sacCer3': 'https://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz',
        'danRer11': 'https://hgdownload.cse.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.fa.gz'
    }
    
    link = genome_links[genome]
    output_directory = os.path.join(out_dir, genome)
    os.makedirs(output_directory, exist_ok=True)
    
    os.system(f'wget {link} -P {output_directory}')
    
    if link.endswith(f'{genome}/bigZips/{genome}.fa.gz'):
        output_file = os.path.join(output_directory, f'{genome}.fa.gz')
        os.system(f'gunzip -f {output_file}')
        
    if link.endswith(f'{genome}/bigZips/chromFa.tar.gz'):
        output_file = os.path.join(output_directory, 'chromFa.tar.gz')
        os.system(f'tar -xzf {output_file} -C {output_directory}')
        os.system(f'cat {output_directory}/*.fa > {output_directory}/{genome}.fa')
        os.system(f'rm -f {output_directory}/chr*.fa')
        
    if link.endswith('hg38.chromFa.tar.gz'):
        output_file = os.path.join(output_directory, 'chromFa.tar.gz')
        os.system(f'tar -xzf {output_file} -C {output_directory}')
        os.system(f'cat {output_directory}/chroms/*.fa > {output_directory}/{genome}.fa')
        os.system(f'rm -rf {output_directory}/chroms/')

# Make project files 
def generate_project_files(working_dir, controltypes, projects, peaktypes, aligners, peakcallers, deduplicators, num_tests):
    # Set up directory structure for project files if needed
    if not os.path.exists(f'{working_dir}/project_files/'):
        print(f'Directory {working_dir}/project_files/ not found. Creating {working_dir}/project_files/')
        os.system(f'mkdir {working_dir}/project_files')
    
    # Start combinatorial project file generation
    for control in controltypes:
        for project in projects:
            for peaktype in peaktypes:
                for aligner in aligners:
                    for peakcaller in peakcallers:
                        for deduplicator in deduplicators:
                            for i in range(1, num_tests + 1):
                            
                                if control == "with_control":
                                
                                    # Rube project files
                                    if project == "Rube":
                                        proj_file_info = textwrap.dedent(f"""
                                        Author: {authors}
                                        Project: Rube_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}
                                        Genome:
                                            Name: mm9
                                            Location: '{working_dir}/mm9/mm9.fa'
                                        Reads:
                                            Samples:
                                                grp1: 
                                                    - SRR2119601
                                                    - SRR2119602
                                            Controls:
                                                ctl1: 
                                                    - SRR2119603
                                                    - SRR2119604
                                        Readtype: single
                                        Peaktype: {peaktype}
                                        Aligner: {aligner}
                                        Deduplicator: {deduplicator}
                                        Peakcaller: {peakcaller}
                                        Threads: 6
                                        """)
                                    
                                    # Namani project files
                                    elif project == "Namani":
                                        proj_file_info = textwrap.dedent(f"""
                                        Author: {authors}
                                        Project: Namani_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}
                                        Genome:
                                            Name: hg38
                                            Location: '{working_dir}/hg38/hg38.fa'
                                        Reads:
                                            Samples:
                                                grp1: 
                                                    - SRR10588628
                                            Controls:
                                                ctl1: 
                                                    - SRR10588629
                                        Readtype: single
                                        Peaktype: {peaktype}
                                        Aligner: {aligner}
                                        Deduplicator: {deduplicator}
                                        Peakcaller: {peakcaller}
                                        Threads: 6
                                        """)
                                    
                                elif control == "no_control":
                                    if peakcaller == "cisgenome" or peakcaller == "pepr": continue
                                    
                                    # Rube project files
                                    if project == "Rube":
                                    
                                    proj_file_info = textwrap.dedent(f"""
                                    Author: {authors}
                                    Project: Rube_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}
                                    Genome:
                                        Name: mm9
                                        Location: '{working_dir}/mm9/mm9.fa'
                                    Reads:
                                        Samples:
                                            grp1: 
                                                - SRR2119601
                                                - SRR2119602
                                        Controls:
                                    Readtype: single
                                    Peaktype: {peaktype}
                                    Aligner: {aligner}
                                    Deduplicator: {deduplicator}
                                    Peakcaller: {peakcaller}
                                    Threads: 6
                                    """)
                                    
                                    # Namani project files
                                    elif project == "Namani":
                                        proj_file_info = textwrap.dedent(f"""
                                        Author: {authors}
                                        Project: Namani_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}
                                        Genome:
                                            Name: hg38
                                            Location: '{working_dir}/hg38/hg38.fa'
                                        Reads:
                                            Samples:
                                                grp1: 
                                                    - SRR10588628
                                        Readtype: single
                                        Peaktype: {peaktype}
                                        Aligner: {aligner}
                                        Deduplicator: {deduplicator}
                                        Peakcaller: {peakcaller}
                                        Threads: 6
                                        """)
                                
                                print(f'Generating {working_dir}/project_files/{project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}.yaml...')
                                os.system(f'touch {working_dir}/project_files/{project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}.yaml')
                                with open(f'{working_dir}/project_files/{project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}.yaml', 'w') as f:
                                    f.write(f'{proj_file_info}')
                                os.system(f'sed -i \'1d\' {working_dir}/project_files/{project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}.yaml')

# Make Snakefiles
def generate_snakefiles(working_dir, controltypes, projects, peaktypes, aligners, peakcallers, deduplicators, num_tests):
    # Set up directory structure for Snakefiles if needed
    if not os.path.exists(f'{working_dir}/snakefiles/'):
        print(f'Directory {working_dir}/snakefiles/ not found. Creating {working_dir}/snakefiles/')
        os.system(f'mkdir {working_dir}/snakefiles')
    
    # Start combinatorial project file generation
    for control in controltypes:
        for project in projects:
            for peaktype in peaktypes:
                for aligner in aligners:
                    for peakcaller in peakcallers:
                        for deduplicator in deduplicators:
                            for i in range(1, num_tests + 1):
                                if (control == "no_control") and (peakcaller == "cisgenome" or peakcaller == "pepr"):
                                    continue
                                    
                                else:
                                    # Make the directory structure if it does not already exist
                                    snakefile_dir = f'{working_dir}/snakefiles/{project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}'
                                    if not os.path.exists(snakefile_dir):
                                        os.makedirs(snakefile_dir)
                                    
                                    # Create the snakefiles using Rocketchip
                                    print(f"Generating {snakefile_dir}/{project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}")
                                    os.system(f'rocketchip {working_dir}/project_files/{project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}.yaml --data {snakefile_dir} --output_file {snakefile_dir}/{project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}')

# Run Snakefiles
def run_snakefiles(working_dir, controltypes, projects, peaktypes, aligners, peakcallers, deduplicators, num_tests):
    # Run Snakemake for all combinations
    for control in controltypes:
        for project in projects:
            for peaktype in peaktypes:
                for aligner in aligners:
                    for peakcaller in peakcallers:
                        for deduplicator in deduplicators:
                            for i in range(1, num_tests + 1):
                                if (control == "no_control") and (peakcaller == "cisgenome" or peakcaller == "pepr"):
                                    continue
                                else:
                                    # Change into snakefile directory
                                    os.chdir(f'{working_dir}/snakefiles/{project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}')
                                    os.system('pwd')
                                    # Run snakefile
                                    os.system(f'snakemake -j 4 -s {project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}')
                                    # Go back to original directory
                                    os.chdir(f'../../')

# Count peaks and compare called vs. true peaks
def count_peaks(working_dir, controltypes, projects, peaktypes, aligners, peakcallers, deduplicators, num_tests, output_path):
    for control in controltypes:
        for project in projects:
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
                                    df.loc[len(df)] = [readtype, project, aligner, peakcaller, deduplicator, i, control, genome_path, read_1_for_path, read_1_rev_path, read_2_for_path, read_2_rev_path, reads_per_peak, padding, reads_std_dev, width, length, paired, flank, expected_peaks, obs_peak_num, true_positives, true_negatives, false_positives, false_negatives]
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

####################
## Run Everything ##
####################

# Download genomes
download_genome(mm9, working_dir)
download_genome(hg38, working_dir)

# Generate project_files
generate_project_files(working_dir, controltypes, projects, peaktypes, aligners, peakcallers, deduplicators, num_tests)

# Generate Snakefiles
generate_snakefiles(working_dir, controltypes, projects, peaktypes, aligners, peakcallers, deduplicators, num_tests)

# Run Snakefiles
run_snakefiles(working_dir, controltypes, projects, peaktypes, aligners, peakcallers, deduplicators, num_tests)

'''
# Count peaks
count_peaks(working_dir = working_dir, controltypes = controltypes, readtypes = readtypes, peaktypes = peaktypes, aligners = aligners, peakcallers = peakcallers, deduplicators = deduplicators, num_tests = num_tests, output_path = 'peaks.csv')
'''