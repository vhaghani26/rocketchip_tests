# Validating CUT&Tag and CUT&RUN

In order to validate Rocketchip's ability to analyze CUT&Tag and CUT&RUN data, we are running data from the study ["Identification of chromatin states during zebrafish gastrulation using CUT&RUN and CUT&Tag"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8976701/) by Akdogan-Ozdilek et al. through Rocketchip using Bowtie2, Picard, and MACS3.

Note: For reasons that elude me (presumably processing at different times), the CUT&RUN data set had a mix of single- and paired-end reads, which require different processing steps and were therefore split into analyses corresponding to their endedness.

## Instructions

Before you begin, make sure you have exported `ROCKETCHIP_SRC` and/or `ROCKETCHIP_DATA` or add them onto the ends of the commands as arguments (as seen below). Also ensure that you are in the activated `rocketchip` Conda environment. I ran the following to conduct the analysis for the CUT&Tag and CUT&RUN data.

CUT&Tag:

```
rocketchip cut_and_tag.yaml --output_file cut_and_tag --data .
snakemake -j 4 -s cut_and_tag
```

CUT&RUN for single-end reads:

```
rocketchip cut_and_run_se.yaml --output_file cut_and_run_se --data .
snakemake -j 4 -s cut_and_run_se
```

CUT&RUN for paired-end reads:

```
rocketchip cut_and_run_pe.yaml --output_file cut_and_run_pe --data .
snakemake -j 4 -s cut_and_run_pe
```

## Results

| SRA ID       | Origin   | Endedness | Metadata               | Raw Reads    | Aligned Reads | % Aligned (Akdogan-Ozdilek et. al) | % Aligned (Rocketchip)  |
| :----------: | :------: | :-------: | :--------------------: | :----------: | :-----------: | :--------------------------------: | :---------------------: |
| SRR14850825  | CUT&RUN  | Single    | 6hpf_H3K4me3_rep1      | 31,343,063   | 23,967,960    | 77.19%                             | 95.68%                  |
| SRR14850826  | CUT&RUN  | Single    | 6hpf_H3K4me3_rep2      | 70,960,918   | 53,529,863    | 76.16%                             | 95.41%                  |
| SRR14850827  | CUT&RUN  | Single    | 6hpf_H3K27me3_rep1     | 26,960,471   | 18,939,978    | 70.93%                             | 95.97%                  |
| SRR14850828  | CUT&RUN  | Single    | 6hpf_H3K27me3_rep2     | 27,851,277   | 19,356,459    | 70.16%                             | 95.34%                  |
| SRR14850829  | CUT&RUN  | Single    | 6hpf_H3K9me3_rep1      | 22,214,067   | 4,830,123     | 22.03%                             | 96.25%                  |
| SRR14850830  | CUT&RUN  | Single    | 6hpf_H3K9me3_rep2      | 43,462,568   | 9,922,627     | 23.13%                             | 96.19%                  |
| SRR14850831  | CUT&RUN  | Paired    | 6hpf_pol2_rep1         | 9,842,850    | 7,214,488     | 73.38%                             | 92.91%                  |
| SRR14850832  | CUT&RUN  | Paired    | 6hpf_pol2_rep2         | 5,707,930    | 4,137,608     | 72.58%                             | 92.45%                  |
| SRR14850833  | CUT&RUN  | Paired    | 6hpf_IgG               | 8,094,912    | 4,777,942     | 59.14%                             | 78.45%                  |
| SRR14870792  | CUT&Tag  | Paired    | 6hpf_H2AZ_Rep1         | 1,233,216    | 1,126,474     | 91.34%                             | 84.53%                  |
| SRR14870793  | CUT&Tag  | Paired    | 6hpf_H2AZ_Rep2         | 1,447,412    | 1,305,818     | 90.22%                             | 82.65%                  |
| SRR14870794  | CUT&Tag  | Paired    | 6hpf_H2AZ_Rep3         | 3,135,635    | 2,847,270     | 90.80%                             | 78.47%                  |
| SRR14870795  | CUT&Tag  | Paired    | 24hpf_H2AZ_Rep1        | 23,376,962   | 21,779,997    | 93.17%                             | 85.90%                  |
| SRR14870796  | CUT&Tag  | Paired    | 24hpf_H2AZ_Rep2        | 28,652,092   | 26,947,558    | 94.05%                             | 85.68%                  |
| SRR14870797  | CUT&Tag  | Paired    | 24hpf_H2AZ_Rep3        | 23,810,168   | 22,288,633    | 93.61%                             | 85.51%                  |

The original study describes different protocols for processing CUT&RUN vs. CUT&Tag data. It is important to note that the purpose of this experiment is not to reproduce identical results, but rather to show that Rocketchip can be applied to CUT&Tag and CUT&RUN data relatively successfully. As such, we have chosen to compare alignment rates for each sample, which were reported as Supplementary Table 2 of the original study. To carry out our analysis, we have tried to match software as appropriately as possible for each analysis. For CUT&RUN, we separated the single- and paired-end reads into separate analyses as required by Rocketchip, aligned the data with Bowtie2, deduplicated with Samtools, and ran broad peak-calling with MACS3. For the CUT&Tag data, we aligned the data with Bowtie2, deduplicated with Picard, and ran broad peak-calling with MACS3. One of the primary differences was that trimming was not performed with Rocketchip, whereas the original study conducted trimming. Furthermore, there were some differences in command line arguments, as Rocketchip was employed using the preset parameters. 

Rocketchip outputs logs containing the standard error of every step. For the alignment, these are named `{sample}_align_reads_err.log`. It is the only log output containing "align" in the file name. Therefore, to assess the alignment, I ran:

```
for file in *align*; do
    echo "File: $file"
    cat "$file"
    echo -e "\n"
done
```
