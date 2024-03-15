#!/usr/bin/env python3

'''
python3 count_peaks.py
'''

####################
## Import Modules ##
####################

import pandas as pd

# Run Snakefiles 
def run_snakefiles(working_dir, controltypes, projects, peaktypes, aligners, peakcallers, deduplicators, num_tests, output_path):
    # Run Snakemake for all combinations
    for control in controltypes:
        for project in projects:
            for peaktype in peaktypes:
                for aligner in aligners:
                    for peakcaller in peakcallers:
                        for deduplicator in deduplicators:
                            for i in range(1, num_tests + 1):
                            
                                # Create DataFrame for peak counting 
                                df = pd.DataFrame(columns=["Endedness", "Project", "Peak_Type", "Aligner", "Peak_Caller", "Deduplicator", "Sample_ID", "Control_ID", "Trial_1_Peaks", "Trial_2_Peaks", "Trial_3_Peaks"])
                                
                                if project == "Rube" or project == "Namani":
                                    readtype = "single"
                            
                                # Handle illegal combinations
                                if (control == "no_control") and (peakcaller == "cisgenome" or peakcaller == "pepr"):
                                    # Assign variables
                                    sample_id = "NA"
                                    control_id = "NA"
                                    trial_1_peaks = "NA"
                                    trial_2_peaks = "NA"
                                    trial_3_peaks = "NA"
                                    # Add data to dataframe
                                    df.loc[len(df)] = [readtype, project, aligner, peakcaller, deduplicator, i, control, genome_path, read_1_for_path, read_1_rev_path, read_2_for_path, read_2_rev_path, reads_per_peak, padding, reads_std_dev, width, length, paired, flank, expected_peaks, obs_peak_num, true_positives, true_negatives, false_positives, false_negatives]
                                
                                # Run snakefiles and count peaks
                                else:
                                    # Change into snakefile directory
                                    os.chdir(f'{working_dir}/snakefiles/{project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}')
                                    os.system('pwd')
                                    # Run snakefile
                                    os.system(f'snakemake -j 4 -s {project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}')
                                    
                                    # Count peaks



                                    
                                    
                                    
                                    # Go back to original directory
                                    os.chdir(f'../../../')








                                    
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