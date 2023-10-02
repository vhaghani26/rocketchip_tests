# Run Times

Before you begin, make sure you have exported `ROCKETCHIP_SRC` and/or `ROCKETCHIP_DATA` or add them onto the ends of the commands as arguments. Also ensure that you are in the activated `rocketchip` Conda environment.

| SRA ID      | Raw Read Count | Genome               | Real Time (s) | User Time (s) | System Time (s) |
| :---------: | :------------: | :------------------: | :-------: | :-------: | :---------: |
| DRR345782   |       NA       | Human (hg19)         |    NA     |    NA     |      NA     |
| SRR17514595 |       NA       | Human (hg19)         |    NA     |    NA     |      NA     | 
| SRR17409984 |       NA       | Human (hg19)         |    NA     |    NA     |      NA     |
| SRR17887304 |       NA       | Mouse (mm10)         |    NA     |    NA     |      NA     |
| SRR26041601 |       NA       | Mouse (mm10)         |    NA     |    NA     |      NA     |
| SRR14407118 |       NA       | Mouse (mm10)         |    NA     |    NA     |      NA     |
| ERR6356102  |       NA       | Rat (rn6)            |    NA     |    NA     |      NA     |
| ERR6356099  |       NA       | Rat (rn6)            |    NA     |    NA     |      NA     |
| ERR6356096  |       NA       | Rat (rn6)            |    NA     |    NA     |      NA     |
| SRR15509782 |       NA       | Zebrafish (danRer11) |    NA     |    NA     |      NA     |
| SRR15509781 |       NA       | Zebrafish (danRer11) |    NA     |    NA     |      NA     |
| SRR15046104 |       NA       | Zebrafish (danRer11) |    NA     |    NA     |      NA     |
| SRR16638475 |       NA       | Fruitfly (dm6)       |    NA     |    NA     |      NA     |
| SRR16638474 |       NA       | Fruitfly (dm6)       |    NA     |    NA     |      NA     |
| SRR15243009 |       NA       | Fruitfly (dm6)       |    NA     |    NA     |      NA     |
| SRR13125172 |       NA       | Worm (ce11)          |    NA     |    NA     |      NA     |
| SRR13125170 |       NA       | Worm (ce11)          |    NA     |    NA     |      NA     |
| SRR13125168 | 12569340 | Worm (ce11)          | 182m30.046s | 183m30.985s | 2m35.470s |
| SRR17329288 | 76503112 | Yeast (sacCer3)      | 164m58.914s | 169m55.851s | 3m8.378s |
| SRR17329314 |       NA       | Yeast (sacCer3)      |    NA     |    NA     |      NA     |
| SRR17329289 |       NA       | Yeast (sacCer3)      |    NA     |    NA     |      NA     |

These experiments were run on an HPC with 64 CPUs and 250 GB of memory available. However, jobs were run without being parallelized (i.e. one job at a time). Additionally, genome copies were deleted between runs using the same genome to ensure that run time accounts for the full workflow.

Raw read counts were calculated using:

```
zcat {SRA_ID}_1.fastq.gz | awk 'NR%4==1{count++} END{print count}'
```

This method of assessing read count was verified based on the "total sequences" section of the the FastQC report outputs.

All `yaml` files are formatted for use, as data is downloaded from the NCBI SRA. The only thing that should be edited is the author name. The commands to be run for each test are documented below.

## Fruitfly Genome

## Human Genome

## Mouse Genome

## Rat Genome

## Worm Genome

### SRR13125168

```
python3 ../../rocketchip worm_SRR13125168.yaml --output_file worm_SRR13125168 --data .
time snakemake -j 1 -s worm_SRR13125168
```

### SRR13125170

```
python3 ../../rocketchip worm_SRR13125170.yaml --output_file worm_SRR13125170 --data .
time snakemake -j 1 -s worm_SRR13125170
```

### SRR13125172

```
python3 ../../rocketchip worm_SRR13125172.yaml --output_file worm_SRR13125172 --data .
time snakemake -j 1 -s worm_SRR13125172
```

## Yeast Genome 

### SRR17329288

```
python3 ../../rocketchip yeast_SRR17329288.yaml --output_file yeast_SRR17329288 --data .
time snakemake -j 1 -s yeast_SRR17329288
```

### SRR17329289

```
python3 ../../rocketchip yeast_SRR17329289.yaml --output_file  yeast_SRR17329289 --data .
time snakemake -j 1 -s yeast_SRR17329289
```

### SRR17329314

```
python3 ../../rocketchip yeast_SRR17329314.yaml --output_file  yeast_SRR17329314 --data .
time snakemake -j 1 -s yeast_SRR17329314
```

## Zebrafish Genome
