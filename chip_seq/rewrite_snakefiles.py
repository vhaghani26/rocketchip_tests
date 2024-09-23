#!/usr/bin/env python3

'''
Note: I ran this script in the chip_seq/ directory:
    python3 rewrite_snakefiles.py --indir snakefiles/
'''

####################
## Import Modules ##
####################

import argparse
import os

#####################
## Set Up Argparse ##
#####################

# Initialize argparse
parser = argparse.ArgumentParser(
    description='Replace unnecessary files for peak-calling')

parser.add_argument('--indir', required=True, type=str,
    metavar='<str>', help='Input directory containing snakefiles')
    
# Finalization of argparse
arg = parser.parse_args()

#####################
## Edit Snakefiles ##
#####################

# Walk through file system starting at input directory
for root, dirs, files in os.walk(arg.indir):
    for file_name in files:
        # Use only snakefiles with the following conditions
        if "_control_" in file_name and "_test" in file_name:
        
            # Get the full path of the file
            file_to_open = os.path.join(root, file_name)
            
            # Read in the file
            with open(file_to_open, 'r') as file:
                filedata = file.read()
                
            # Remove empty lines
            filedata = '\n'.join([line for line in filedata.split('\n') if line.strip()])
            
            # Make general changes (apply to all files)
            filedata = filedata.replace('expand("02_fastqc_analysis/', '#expand("02_fastqc_analysis/')
            filedata = filedata.replace('expand("05_bigwig_files', '#expand("05_bigwig_files')
            filedata = filedata.replace('output: "03_sam_files/{sample}.sam"', 'output: temp("03_sam_files/{sample}.sam")')
            filedata = filedata.replace('output: "04_bam_files/{sample}.bam"', 'output: temp("04_bam_files/{sample}.bam")')
            filedata = filedata.replace('output: "04_bam_files/{sample}.fixmate.bam"', 'output: temp("04_bam_files/{sample}.fixmate.bam")')
            filedata = filedata.replace('output: "04_bam_files/{sample}.sorted.fixmate.bam"', 'output: temp("04_bam_files/{sample}.sorted.fixmate.bam")')
            filedata = filedata.replace('output: "04_bam_files/{sample}.sorted.dedup.bam"', 'output: temp("04_bam_files/{sample}.sorted.dedup.bam")')
            filedata = filedata.replace('output: "04_bam_files/{sample}.sorted.dedup.bam.bai"', 'output: temp("04_bam_files/{sample}.sorted.dedup.bam.bai")')
            
            # MACS3 has an issue with small read files (this data is subsetted) requiring the --nomodel flag, so I will manually add it
            if "macs3 callpeak" in filedata:
                filedata = filedata.replace('macs3 callpeak', 'macs3 callpeak --nomodel')
                       
            # Write the file out again
            with open(file_to_open, 'w') as file:
                file.write(filedata)
                
            # Add verbosity
            print(f'Changes complete for {file_to_open}')