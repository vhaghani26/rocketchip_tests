#!/usr/bin/env python3

'''
Usage:
python3 count_peaks.py --in_dir {Snakefile directory} --out_dir {peak counts}
'''

####################
## Import Modules ##
####################

import os

#####################
## Set Up Argparse ##
#####################

# Initialize argparse
parser = argparse.ArgumentParser(
    description='Count peaks using peak output directory')

parser.add_argument('--in_dir', required=True, type=str,
    metavar='<str>', help='Absolute directory containing the subdirectory with the peak data in it (aka Snakefile directory)')
    
parser.add_argument('--out_dir', required=True, type=str,
    metavar='<str>', help='Absolute directory to output peak counts')
    
# Finalization of argparse
arg = parser.parse_args()

#################
## Count Peaks ##
#################

def count_peaks(in_dir, out_dir):                 
    # Count peaks for MACS3 outputs
    if "macs3" in in_dir:
        os.chdir('06_macs3_peaks')
        if "with_control" in in_dir:
            if "narrow" in in_dir:
                if os.path.isfile('grp1_ctl1_peaks.narrowPeak'):
                    with open('grp1_ctl1_peaks.narrowPeak', 'r') as file:
                        obs_peak_num = sum(1 for line in file)
                else:
                    obs_peak_num = 0
            elif "broad" in in_dir:
                if os.path.isfile('grp1_ctl1_peaks.broadPeak'):
                    with open('grp1_ctl1_peaks.broadPeak', 'r') as file:
                        obs_peak_num = sum(1 for line in file)
                else:
                    obs_peak_num = 0
        elif "no_control"  in in_dir:
            if "narrow" in in_dir:
                if os.path.isfile('grp1_peaks.narrowPeak'):
                    with open('grp1_peaks.narrowPeak', 'r') as file:
                        obs_peak_num = sum(1 for line in file)
                else:
                    obs_peak_num = 0
            elif "broad" in in_dir:
                if os.path.isfile('grp1_peaks.broadPeak'):
                    with open('grp1_peaks.broadPeak', 'r') as file:
                        obs_peak_num = sum(1 for line in file)
                else:
                    obs_peak_num = 0
        os.chdir('..')

    # Count peaks for Cisgenome outputs
    elif "cisgenome" in in_dir:
        os.chdir('06_cisgenome_peaks')
        if os.path.isfile('grp1_ctl1_peak.cod'):
            # Open the file for reading
            with open("grp1_ctl1_peak.cod", "r") as file:
                obs_peak_num = sum(1 for line in file)
                obs_peak_num = obs_peak_num - 1
        else:
            obs_peak_num = 0
        os.chdir('..')

    # Count peaks for Genrich outputs
    elif "genrich" in in_dir:
        os.chdir('06_genrich_peaks')
        if "with_control" in in_dir:
            if os.path.isfile('grp1_ctl1_peak.narrowPeak'):
                with open('grp1_ctl1_peak.narrowPeak', 'r') as file:
                    obs_peak_num = sum(1 for line in file)
            else:
                obs_peak_num = 0
        elif "no_control" in in_dir:
            if os.path.isfile('grp1_peak.narrowPeak'):
                with open('grp1_peak.narrowPeak', 'r') as file:
                    obs_peak_num = sum(1 for line in file)
            else:
                obs_peak_num = 0
        os.chdir('..')
    
    # Count peaks for PePr outputs
    elif "pepr"  in in_dir:
        os.chdir('06_pepr_peaks')
        if os.path.isfile('grp1_ctl1__PePr_peaks.bed'):
            with open('grp1_ctl1__PePr_peaks.bed', 'r') as file:
                obs_peak_num = sum(1 for line in file)
        else:
            obs_peak_num = 0
        os.chdir('..')
    
    # Get basename of input directory only
    input_dir = os.path.basename(in_dir)
    
    # Check if input_dir ends with '/' and remove it if so
    if input_dir.endswith('/'):
        input_dir = input_dir[:-1]
    
    # Set up output directory 
    output_file = os.join(out_dir, input_dir + "_peaks")
    
    # Save peak count
    os.system(f'touch {output_file}')
    with open(f'{output_file}', 'w') as f:
        f.write(f'{obs_peak_num}')

######################
## Execute Function ##
######################

count_peaks(arg.in_dir, arg.out_dir)