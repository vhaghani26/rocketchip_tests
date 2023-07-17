#!/usr/bin/env python3

'''
Note: I ran this script in the exp_vs_obs/ directory:
    python3 exp_vs_obs.py
'''

####################
## Import Modules ##
####################

import textwrap
import os
import pandas as pd
import sys # only for sys exit

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
#df = pd.DataFrame(columns=["Endedness", "Peak_Type", "Aligner", "Peak_Caller", "Deduplicator", "Test_Dataset", "Control", "Synthetic_Genome_Path", "Synthetic_Forward_Read_1_Path", "Synthetic_Reverse_Read_1_Path", "Synthetic_Forward_Read_2_Path", "Synthetic_Reverse_Read_2_Path", "Expected_Peaks", "Observed_Peaks"])

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
                            if "control" == "with_control":
                                proj_file_info = textwrap.dedent(f"""
                                Author: Viktoria_Haghani_and_Aditi_Goyal_and_Alan_Zhang
                                Project: exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}
                                Genome:
                                Name: genome
                                Location: 'seq_data/{readtype}_{peaktype}/test_{i}/genome.fa'
                                Reads:
                                Samples:
                                    grp1: 
                                    - 'seq_data/{readtype}_{peaktype}/test_{i}/exp_a'
                                    - 'seq_data/{readtype}_{peaktype}/test_{i}/exp_b'
                                Controls:
                                    ctl1: 
                                    - 'seq_data/{readtype}_{peaktype}/test_{i}/input'
                                Readtype: {readtype}
                                Peaktype: {peaktype}
                                Aligner: {aligner}
                                Deduplicator: {deduplicator}
                                Peakcaller: {peakcaller}
                                Threads: 1
                                """)
                            if "control" == "no_control":
                                proj_file_info = textwrap.dedent(f"""
                                Author: Viktoria_Haghani_and_Aditi_Goyal_and_Alan_Zhang
                                Project: exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}
                                Genome:
                                Name: genome
                                Location: 'seq_data/{readtype}_{peaktype}/test_{i}/genome.fa'
                                Reads:
                                Samples:
                                    grp1: 
                                    - 'seq_data/{readtype}_{peaktype}/test_{i}/exp_a'
                                    - 'seq_data/{readtype}_{peaktype}/test_{i}/exp_b'
                                Controls:
                                    ctl1: 
                                    - 'seq_data/{readtype}_{peaktype}/test_{i}/input'
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
                            os.system(f'python3 ../rocketchip project_files/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml --data seq_data/{readtype}_{peaktype}/test_{i} --src . --output_file snakefiles/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                            
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

                            # Change into snakefile directory
                            os.chdir(f'exp_vs_obs/snakefiles/{control}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                            os.system('pwd')
                            os.system('ls')
                            # Run snakefile
                            os.system(f'snakemake -j 4 -s exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                            # Go back to original directory
                            os.chdir(f'../../../../')
                            
#################
## Count Peaks ##
#################

                            '''
                            # Change into snakefile directory
                            os.chdir(f'exp_vs_obs/snakefiles/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}')
                            os.system('pwd')
                            os.system('ls')
                            # Run snakefile
                            os.system(f'snakemake -j 4 -s exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}')
                            # Go back to original directory
                            os.chdir(f'../../../')

                            genome_path = 
                            read_1_for_path =
                            read_1_rev_path = 
                            read_2_for_path =
                            read_2_rev_path = 
                            obs_peak_num = 
                            df.loc[len(df)] = [readtype, peaktype, aligner, peakcaller, deduplicator, i, control, genome_path, read_1_for_path, read_1_rev_path, read_2_for_path, read_2_rev_path, 50, obs_peak_num]
                            '''