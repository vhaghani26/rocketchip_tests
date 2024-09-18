# Replicability 

In order to prove that Rocketchip is capable of replicating experimental results and running all software combinations, we ran ChIP-seq data from the study ["Sequence features accurately predict genome-wide MeCP2 binding in vivo"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4820824/) by Rube et al. 2016 and ["Genome-wide global identification of NRF2 binding sites in A549 non-small cell lung cancer cells by ChIP-Seq reveals NRF2 regulation of genes involved in focal adhesion pathways"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6949066/#SD1) by Namani et al. 2019 using Rocketchip.

## Instructions

Before you begin, make sure you have properly installed Rocketchip and that you are in the activated `rocketchip` Conda environment. I ran the following to conduct the analysis:

```
# Run all combinations of software (this needs a MASSIVE amount of time/storage)
python3 replicate_all_combos.py

# Count all the peaks and create a CSV file documenting the results
python3 peak_count_csv.py --in_dir peak_counts/ --out_dir .

# Print out the combinations that yielded variation in peak-calling
python3 peak_variation.py
```

In order to generate the heatmaps, I ran the code seen in `generate_heatmap.ipynb` using the Conda environment for figure generation detailed in the home directory of this repository.

## Deeper Dive

### Downsampling Data

The results revealed some variation in the Rube data set. In order to better understand what may be contributing to the non-deterministic results, we ran a subset of the data through Rocketchip 100 times. We downloaded the raw data from the Rube data set and subsampled it to 5000 reads using Seqtk (v1.4):

```
seqtk sample -s100 SRR2119601.fastq.gz 5000 > SRR2119601_5000.fastq
seqtk sample -s100 SRR2119602.fastq.gz 5000 > SRR2119602_5000.fastq
seqtk sample -s100 SRR2119603.fastq.gz 5000 > SRR2119603_5000.fastq
seqtk sample -s100 SRR2119604.fastq.gz 5000 > SRR2119604_5000.fastq
```

The `-s100` is an assigned random seed to ensure the results can be reproduced if needed. To format the data, we then gzipped the files and renamed them to be compatible with what Rocketchip is expecting:

```
# Gzip
gzip SRR2119601_5000.fastq
gzip SRR2119602_5000.fastq
gzip SRR2119603_5000.fastq
gzip SRR2119604_5000.fastq

# Rename
mv SRR2119601_5000.fastq.gz SRR2119601.fastq.gz
mv SRR2119602_5000.fastq.gz SRR2119602.fastq.gz
mv SRR2119603_5000.fastq.gz SRR2119603.fastq.gz
mv SRR2119604_5000.fastq.gz SRR2119604.fastq.gz
```

### Generate YAML File

We started by generating the yaml file (`snakefiles.yaml`) that has all the combinations we want to test. We generated it by running:

```
python3 deep_dive_yaml.py
```

This contains 14400 "samples" (these are each independent analysis combinations and trials) to test. Essentially, there are 144 combinations of software being tested with 100 trials each

### Generate Project and Snakefiles

We changed the combinatorial testing variables section of `replicate_all_combos.py` to:

```
# Combinatorial testing variables
controltypes = ["with_control", "no_control"] 
projects = ["Rube"]
peaktypes = ["narrow", "broad"]
aligners = ["bwa_mem", "bowtie2", "STAR"]
peakcallers = ["macs3", "cisgenome", "genrich", "pepr"]
deduplicators = ["samtools", "no_deduplication", "sambamba", "picard"]
num_tests = 100
```

We changed the executed functions at the bottom of the `rewrite_snakefiles.py` script to:

```
####################
## Run Everything ##
####################

# Create master dataframe
#create_csv(working_dir)

# Download genomes
download_genome('mm9', working_dir)
#download_genome('hg38', working_dir)

# Generate project_files
generate_project_files(working_dir, controltypes, projects, peaktypes, aligners, peakcallers, deduplicators, num_tests)

# Generate Snakefiles
generate_snakefiles(working_dir, controltypes, projects, peaktypes, aligners, peakcallers, deduplicators, num_tests)

# Run Snakefiles
#run_snakefiles(working_dir, controltypes, projects, peaktypes, aligners, peakcallers, deduplicators, num_tests)

# Make and run SLURM scripts 
#run_slurm_scripts(working_dir, controltypes, projects, peaktypes, aligners, peakcallers, deduplicators, num_tests)
```

Then we ran the script:

```
python3 replicate_all_combos.py
```

This generates the project files and snakefiles needed for Rocketchip. Because we are also not interested in storing the data or generating the bigwig files, as we are just focusing on peak counts in this analysis, we ran `rewrite_snakefiles.py` to make the intermediate outputs temporary up to when peaks are called and removed the quality control and bigwig steps. This was run using:

```
python3 rewrite_snakefiles.py --indir snakefiles/
```

### Running the Analysis

While Snakemake natively supports parallelization, we'd like to further parellelize the analyses since there are 14400 analyses to be conducted. As such, we ran the entire analysis by running `snakefile_runner`, a Snakefile that parallelizes the jobs per analysis while also parallelizing the analyses themselves. It stores the final peak count in a file, deleting all the intermediates and stored data along the way so as not to overload the disc storage space. By combining everything into one rule, it also forces the every analysis to run to completion. This was done because we previously experienced issues running multiple analyses at once, where several hundred analyses would be half completed, filling our disc space, but none had any results documented yet. This portion was run using:

```
snakemake -j 5 -s snakefile_runner
```
