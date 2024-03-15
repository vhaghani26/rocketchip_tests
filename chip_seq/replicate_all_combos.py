#!/usr/bin/env python3

'''
python3 replicate_all_combos.py
'''

####################
## Import Modules ##
####################

import sys
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
                                    os.system(f'rocketchip {working_dir}/project_files/{project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}.yaml --data {working_dir} --output_file {snakefile_dir}/{project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}')

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

                                # Handle illegal combinations
                                if (control == "no_control") and (peakcaller == "cisgenome" or peakcaller == "pepr"):
                                    continue
                                    
                                # Run snakefiles and count peaks
                                else:
                                    # Change into snakefile directory
                                    os.chdir(f'{working_dir}/snakefiles/{project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}')
                                    os.system('pwd')
                                    # Run snakefile
                                    os.system(f'snakemake -j 4 -s {project}_{control}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test{i}')
                                    
                                    # Go back to original directory
                                    os.chdir(f'../../')

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