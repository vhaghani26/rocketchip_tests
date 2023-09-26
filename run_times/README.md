# Run Times

Before you begin, make sure you have exported `ROCKETCHIP_SRC` and/or `ROCKETCHIP_DATA` or add them onto the ends of the commands as arguments. Also ensure that you are in the activated `rocketchip` Conda environment.

| SRA ID      | Raw Read Count | Genome               | Run Time |
| :---------: | :------------: | :------------------: | :------: |
| DRR345782   |       NA       | Human (hg19)         |    NA    |
| SRR17514595 |       NA       | Human (hg19)         |    NA    | 
| SRR17409984 |       NA       | Human (hg19)         |    NA    |
| SRR17887304 |       NA       | Mouse (mm10)         |    NA    |
| SRR26041601 |       NA       | Mouse (mm10)         |    NA    |
| SRR14407118 |       NA       | Mouse (mm10)         |    NA    |
| ERR6356102  |       NA       | Rat (rn6)            |    NA    |
| ERR6356099  |       NA       | Rat (rn6)            |    NA    |
| ERR6356096  |       NA       | Rat (rn6)            |    NA    |
| SRR15509782 |       NA       | Zebrafish (danRer11) |    NA    |
| SRR15509781 |       NA       | Zebrafish (danRer11) |    NA    |
| SRR15046104 |       NA       | Zebrafish (danRer11) |    NA    |
| SRR16638475 |       NA       | Fruitfly (dm6)       |    NA    |
| SRR16638474 |       NA       | Fruitfly (dm6)       |    NA    |
| SRR15243009 |       NA       | Fruitfly (dm6)       |    NA    |
| SRR13125172 |       NA       | Worm (ce11)          |    NA    |
| SRR13125170 |       NA       | Worm (ce11)          |    NA    |
| SRR13125168 |       NA       | Worm (ce11)          |    NA    |
| SRR17329288 |       NA       | Yeast (sacCer3)      |    NA    |
| SRR17329314 |       NA       | Yeast (sacCer3)      |    NA    |
| SRR17329289 |       NA       | Yeast (sacCer3)      |    NA    |

These experiments were run on an HPC with 64 CPUs and 250 GB of memory available. However, jobs were run without being parallelized (i.e. one job at a time). Additionally, genome copies were deleted between runs using the same genome to ensure that run time accounts for the full workflow.

All `yaml` files are formatted for use, as data is downloaded from the NCBI SRA. The only thing that should be edited is the author name. The commands to be run for each test are documented below.

## Fruitfly Genome

## Human Genome

## Mouse Genome

## Rat Genome

## Worm Genome

## Yeast Genome 

### SRR17329288

```
python3 ../../rocketchip yeast_SRR17329288.yaml --output_file  yeast_SRR17329288 --data ..
```

This will create a file called Snakefile. To run it, simply run:

```
snakemake -j 1 -s yeast_SRR17329288
```

### SRR17329289

### SRR17329314

## Zebrafish Genome
