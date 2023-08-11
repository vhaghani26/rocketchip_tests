#!/usr/bin/env python3

'''
Note: I ran this script in the exp_vs_obs/ directory:
    python3 exp_vs_obs.py
'''

####################
## Import Modules ##
####################

import pandas as pd
import os
import textwrap
import subprocess 
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import seaborn as sns

# User-specific variables
working_dir = '/share/korflab/home/viki/rocketchip_tests'
authors = 'Viktoria_Haghani_and_Aditi_Goyal_and_Alan_Zhang'

#############################################
## Set Variables for Combinatorial Testing ##
#############################################

controltypes = ["no_control", "with_control"]  
readtypes = ["paired", "single"]
peaktypes = ["narrow", "broad"]
aligners = ["bwa_mem", "bowtie2", "STAR"]
peakcallers = ["macs3", "cisgenome", "genrich", "pepr"]
deduplicators = ["samtools", "no_deduplication", "sambamba", "picard"]
num_tests = 6

# Create DataFrame for peak counting 
df = pd.DataFrame(columns=["Endedness", "Peak_Type", "Aligner", "Peak_Caller", "Deduplicator", "Test_Dataset", "Control", "Synthetic_Genome_Path", "Synthetic_Forward_Read_1_Path", "Synthetic_Reverse_Read_1_Path", "Synthetic_Forward_Read_2_Path", "Synthetic_Reverse_Read_2_Path", "Reads_per_Peak", "Padding", "Reads_STD_Dev", "Width", "Read_Length", "Paired", "Flank", "Expected_Peaks", "Observed_Peaks"])

################################
## Set up Directory Structure ##
################################

# Project files
if not os.path.exists('project_files/'):
    print('Directory project_files/ not found. Creating project_files/')
    os.system(f'mkdir project_files')
if not os.path.exists('project_files/with_control/'):
    print('Directory project_files/with_control/ not found. Creating project_files/with_control/')
    os.system(f'mkdir project_files/with_control')    
if not os.path.exists('project_files/no_control/'):
    print('Directory project_files/no_control/ not found. Creating project_files/no_control/')
    os.system(f'mkdir project_files/no_control')  

# Snakefiles
if not os.path.exists('snakefiles/'):
    print('Directory snakefiles/ not found. Creating snakefiles/')
    os.system(f'mkdir snakefiles')
if not os.path.exists('snakefiles/with_control/'):
    print('Directory snakefiles/with_control/ not found. Creating snakefiles/with_control/')
    os.system(f'mkdir snakefiles/with_control')    
if not os.path.exists('snakefiles/no_control/'):
    print('Directory snakefiles/no_control/ not found. Creating snakefiles/no_control/')
    os.system(f'mkdir snakefiles/no_control')   

'''
########################
## Make Project Files ##
########################

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
                                    Location: '{working_dir}/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/genome.fa'
                                Reads:
                                    Samples:
                                        grp1: 
                                            - '{working_dir}/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/exp_a'
                                            - '{working_dir}/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/exp_b'
                                    Controls:
                                        ctl1: 
                                            - '{working_dir}/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/input'
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
                                    Location: '{working_dir}/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/genome.fa'
                                Reads:
                                    Samples:
                                        grp1: 
                                            - '{working_dir}/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/exp_a'
                                            - '{working_dir}/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/exp_b'
                                    Controls:
                                Readtype: {readtype}
                                Peaktype: {peaktype}
                                Aligner: {aligner}
                                Deduplicator: {deduplicator}
                                Peakcaller: {peakcaller}
                                Threads: 1
                                """)
                            print(f'Generating project_files/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml...')
                            os.system(f'touch project_files/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml')
                            with open(f'project_files/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml', 'w') as f:
                                f.write(f'{proj_file_info}')
                            os.system(f'sed -i \'1d\' project_files/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml')

#####################
## Make Snakefiles ##
#####################  

                            # Make the directory structure if it does not already exist
                            if not os.path.exists(f'snakefiles/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}'):
                                os.system(f'mkdir snakefiles/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                            
                            # Create the snakefiles using Rocketchip 
                            print(f"Generating snakefiles/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}")
                            os.system(f'python3 ../rocketchip {working_dir}/exp_vs_obs/project_files/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml --data {working_dir}/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i} --src {working_dir} --output_file {working_dir}/exp_vs_obs/snakefiles/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                            
                            # MACS3 has an issue with small read files requiring the --nomodel flag, so I will manually add it for the single-end data that are having problems with peak-calling
                            if readtype == "single" and peakcaller == "macs3":
                                print(f'Adding --nomodel flag in MACS3 for snakefile snakefiles/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                file_to_open = f'snakefiles/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}'
                                # Read in the file
                                with open(file_to_open, 'r') as file :
                                    filedata = file.read()
                                # Replace the target string
                                filedata = filedata.replace('macs3 callpeak ', 'macs3 callpeak --nomodel ')
                                # Write the file out again
                                with open(file_to_open, 'w') as file:
                                    file.write(filedata)
                       
####################
## Run Snakefiles ##
####################
                            if (control == "no_control") and (peakcaller == "cisgenome" or peakcaller == "pepr"):
                                continue
                            else:
                                # Change into snakefile directory
                                os.chdir(f'snakefiles/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                os.system('pwd')
                                # Run snakefile
                                os.system(f'snakemake -j 4 -s exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                # Go back to original directory
                                os.chdir(f'../../../')
                            
#################
## Count Peaks ##
#################

                            if (control == "no_control") and (peakcaller == "cisgenome" or peakcaller == "pepr"):
                                continue
                            else:
                                # Change into snakefile directory
                                os.chdir(f'snakefiles/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                print(f'Counting peaks for exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}...')
                                
                                # Assign global data frame variables 
                                genome_path = f'{working_dir}/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/genome.fa'
                                
                                if readtype == "paired":
                                    read_1_for_path = f'{working_dir}/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/exp_a_1.fastq.gz'
                                    read_1_rev_path = f'{working_dir}/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/exp_a_2.fastq.gz'
                                    read_2_for_path = f'{working_dir}/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/exp_b_1.fastq.gz'
                                    read_2_rev_path = f'{working_dir}/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/exp_b_2.fastq.gz'
                                elif readtype == "single":
                                    read_1_for_path = f'{working_dir}/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/exp_a.fastq.gz'
                                    read_1_rev_path = "NA"
                                    read_2_for_path = f'{working_dir}/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/exp_b.fastq.gz'
                                    read_2_rev_path = "NA"
                                
                                # Assign sequence metadata
                                reads_per_peak = 2**i
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
                                
                                # Count peaks and determine peak locations
                                if peakcaller == "macs3":
                                    os.chdir('06_macs3_peaks')
                                    if control == "with_control":
                                        if peaktype == "narrow":
                                            result = subprocess.run('less grp1_ctl1_peaks.narrowPeak | wc -l', shell = True, stdout = subprocess.PIPE, text = True)
                                        elif peaktype == "broad":
                                            result = subprocess.run('less grp1_ctl1_peaks.broadPeak | wc -l', shell = True, stdout = subprocess.PIPE, text = True)
                                    elif control == "no_control":
                                        if peaktype == "narrow":
                                            result = subprocess.run('less grp1_peaks.narrowPeak | wc -l', shell = True, stdout = subprocess.PIPE, text = True)
                                        elif peaktype == "broad":
                                            result = subprocess.run('less grp1_peaks.broadPeak | wc -l', shell = True, stdout = subprocess.PIPE, text = True)
                                    obs_peak_num = int(result.stdout.strip())
                                    os.chdir('..')
                                elif peakcaller == "cisgenome":
                                    os.chdir('06_cisgenome_peaks')
                                    result = subprocess.run('less grp1_ctl1_peak.cod | wc -l', shell = True, stdout = subprocess.PIPE, text = True)
                                    obs_peak_num = int(result.stdout.strip())
                                    obs_peak_num = obs_peak_num - 1
                                    os.chdir('..')
                                elif peakcaller == "genrich":
                                    os.chdir('06_genrich_peaks')
                                    if control == "with_control":
                                        result = subprocess.run('less grp1_ctl1_peak.narrowPeak | wc -l', shell = True, stdout = subprocess.PIPE, text = True)
                                        obs_peak_num = int(result.stdout.strip())
                                    elif control == "no_control":
                                        result = subprocess.run('less grp1_peak.narrowPeak | wc -l', shell = True, stdout = subprocess.PIPE, text = True)
                                        obs_peak_num = int(result.stdout.strip())
                                    os.chdir('..')
                                elif peakcaller == "pepr":
                                    os.chdir('06_pepr_peaks')
                                    if os.path.isfile('grp1_ctl1__PePr_peaks.bed'):
                                        result = subprocess.run('less grp1_ctl1__PePr_peaks.bed | wc -l', shell = True, stdout = subprocess.PIPE, text = True)
                                        obs_peak_num = int(result.stdout.strip())
                                    else:
                                        obs_peak_num = 0
                                    os.chdir('..')
                                print(f'Observed peaks: {obs_peak_num}')
                                
                                # Go back to original directory
                                os.chdir(f'../../../')

                                # Add test to dataframe 
                                df.loc[len(df)] = [readtype, peaktype, aligner, peakcaller, deduplicator, i, control, genome_path, read_1_for_path, read_1_rev_path, read_2_for_path, read_2_rev_path, reads_per_peak, padding, reads_std_dev, width, length, paired, flank, expected_peaks, obs_peak_num]
                                
# Save to CSV
df.to_csv("tables_and_figures/expected_vs_observed_peaks_master.csv", index=False)
'''

###############
## Visualize ##
###############

# Read in data
df = pd.read_csv('tables_and_figures/expected_vs_observed_peaks_master.csv')

'''
#######################
## Overall Histogram ##
#######################

# Find the maximum value of observed peaks in order to adjust visualization
max_observed_peaks = df['Observed_Peaks'].max()
print(f'{max_observed_peaks} was the highest observed peak number.')

# Plot histogram 
plt.hist(df['Observed_Peaks'], bins=200, range=(0, 2000), color='skyblue', edgecolor='black')

# Set x-axis tick positions and labels
#plt.xticks(range(0, 101, 10)) 
plt.xlabel('Number of Observed Peaks')
plt.ylabel('Frequency')
plt.title('Total Distribution of Observed Peaks')
plt.show()

# Save figure
plt.savefig('tables_and_figures/total_distribution_of_observed_peaks.pdf')
'''

###############################
## Visualize Data Parameters ##
###############################

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
plt.savefig('tables_and_figures/reads_per_peak_vs_observed_peaks.pdf')

'''
##############
## Heatmaps ##
##############

# Filter rows so peak caller is unique for peak type and endedness groups
heatmap_df = df.groupby(["Endedness", "Peak_Type"]).filter(lambda x: x["Peak_Caller"].nunique() == 1)

# Generate a heatmap for each unique combination of peak type, endedness, and peak caller
for (readtype, peaktype, peakcaller), data in heatmap_df.groupby(["Endedness", "Peak_Type", "Peak_Caller"]):
    
    # Plot size
    plt.figure(figsize=(10, 6))
    
    # Create a matrix for the heatmap
    pivot_df = data.pivot(index = "Deduplicator", columns = "Aligner", values = "Observed_Peaks")
    sns.heatmap(pivot_df, cmap = "coolwarm", annot = True, fmt = ".2f", cbar_kws = {'label': 'Observed Peaks'})
    
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
    plt.title(f'Number of Peaks for {readtype.title()}-End Data with {peaktype.title()} Peaks using {peakcaller_name}')
    plt.xlabel("Aligner")
    plt.ylabel("Deduplicator")
    
    # Show figure
    plt.show()
    
    # Save figure
    plt.savefig(f'02_tables_and_figures/heatmap_{readtype}_{peaktype}_{peakcaller}.pdf')
'''