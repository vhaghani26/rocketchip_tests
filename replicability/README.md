# Replicability 
In order to prove that Rocketchip is capable of replicating experimental results, we are running ChIP-seq data from the study ["Sequence features accurately predict genome-wide MeCP2 binding in vivo"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4820824/) by Gong et al. through Rocketchip using Bowtie2, Samtools, and MACS3.

## Instructions

Before you begin, make sure you have exported `ROCKETCHIP_SRC` and/or `ROCKETCHIP_DATA` or add them onto the ends of the commands as arguments (as seen below). Also ensure that you are in the activated `rocketchip` Conda environment. I ran the following to conduct the analysis:

```
rocketchip replicability.yaml --output_file replicability --data .
snakemake -j 4 -s replicability
```

The `narrowPeak` file was inspected in order to count the number of peaks each analysis yielded.

## Results

| SRA ID    | Sample Metadata    |
| :-------: | :---------------:  |
|SRR2119601 | MeCP2_ChIP_WT_rep1 |
|SRR2119602 | MeCP2_ChIP_WT_rep2 |
|SRR2119603 | Input_WT_rep1      |
|SRR2119604 | Input_WT_rep2      |

As per the original study, the experimental replicates were combined, and peak calling was conducted using MACS2, which is an earlier version of MACS3. The input control replicates were also combined and utilized in the peak-calling analysis. While there was no explicit mention of the deduplication software employed, it was clarified that the data underwent a deduplication process. Bowtie2 was used for alignment in the original study. Consequently, the analysis was conducted using Bowtie2, Samtools, and MACS3 with the mm9 reference mouse genome. Each of the three runs of Rocketchip consistently yielded a total of 1,329,347 peaks. This outcome highlights Rocketchip's ability to replicate experimental results during the data analysis stage.