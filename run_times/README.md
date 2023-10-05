# Run Times

Before you begin, make sure you have exported `ROCKETCHIP_SRC` and/or `ROCKETCHIP_DATA` or add them onto the ends of the commands as arguments. Also ensure that you are in the activated `rocketchip` Conda environment.

| SRA ID      | Raw Read Count | Genome               | Real Time (s) | User Time (s) | System Time (s) |
| :---------: | :------------: | :------------------: | :-------: | :-------: | :---------: |
| SRR15243009 | (running)       | Fruitfly (dm6)       |    NA     |    NA     |      NA     |
| SRR16638474 | NA       | Fruitfly (dm6)       |    NA     |    NA     |      NA     |
| SRR16638475 | NA       | Fruitfly (dm6)       |    NA     |    NA     |      NA     |
| DRR345782   | (running)  | Human (hg19)         |    NA     |    NA     |      NA     |
| SRR17409984 | NA       | Human (hg19)         |    NA     |    NA     |      NA     |
| SRR17514595 | NA       | Human (hg19)         |    NA     |    NA     |      NA     | 
| SRR14407118 | (running) | Mouse (mm10)         |    NA     |    NA     |      NA     |
| SRR17887304 | NA       | Mouse (mm10)         |    NA     |    NA     |      NA     |
| SRR26041601 | NA       | Mouse (mm10)         |    NA     |    NA     |      NA     |
| ERR6356096  | (running) | Rat (rn6)            |    NA     |    NA     |      NA     |
| ERR6356099  | NA       | Rat (rn6)            |    NA     |    NA     |      NA     |
| ERR6356102  | NA       | Rat (rn6)            |    NA     |    NA     |      NA     |
| SRR13125168 | 12569340 | Worm (ce11)          | 182m30.046s | 183m30.985s | 2m35.470s |
| SRR13125170 | 11144851 | Worm (ce11)          |    NA     |    NA     |      NA     |
| SRR13125172 | NA       | Worm (ce11)          |    NA     |    NA     |      NA     |
| SRR17329288 | 76503112 | Yeast (sacCer3)      | 164m58.914s | 169m55.851s | 3m8.378s |
| SRR17329289 |       NA       | Yeast (sacCer3)      |    NA     |    NA     |      NA     |
| SRR17329314 | 26073747 | Yeast (sacCer3)      | 117m8.917s | 120m59.169s | 1m55.124s |
| SRR15046104 | (running) | Zebrafish (danRer11) |    NA     |    NA     |      NA     |
| SRR15509781 | NA       | Zebrafish (danRer11) |    NA     |    NA     |      NA     |
| SRR15509782 | NA       | Zebrafish (danRer11) |    NA     |    NA     |      NA     |

These experiments were run on an HPC with 64 CPUs and 250 GB of memory available. However, jobs were run without being parallelized (i.e. one job at a time with one thread). Additionally, genome copies were deleted between runs using the same genome to ensure that run time accounts for the full workflow. All data selected was run for narrow-peak calling using paired-end data. The software used was BWA-MEM for alignment, Samtools for deduplication, and MACS3 for peak-calling. 

Raw read counts were calculated using:

```
zcat {SRA_ID}_1.fastq.gz | awk 'NR%4==1{count++} END{print count}'
```

This method of assessing read count was verified based on the "total sequences" section of the the FastQC report outputs.

All `yaml` files are formatted for use, as data is downloaded from the NCBI SRA. The only thing that should be edited is the author name. The commands to be run for each test are documented below each relevant section.

## Fruitfly Genome

### SRR15243009

```
python3 ../../rocketchip fruitfly_SRR15243009.yaml --output_file fruitfly_SRR15243009 --data .
time snakemake -j 1 -s fruitfly_SRR15243009
```

### SRR16638474

```
python3 ../../rocketchip fruitfly_SRR16638474.yaml --output_file fruitfly_SRR16638474 --data .
time snakemake -j 1 -s fruitfly_SRR16638474
```

### SRR16638475

```
python3 ../../rocketchip fruitfly_SRR16638475.yaml --output_file fruitfly_SRR16638475 --data .
time snakemake -j 1 -s fruitfly_SRR16638475
```

## Human Genome

### DRR345782

```
python3 ../../rocketchip human_DRR345782.yaml --output_file human_DRR345782 --data .
time snakemake -j 1 -s human_DRR345782
```

### SRR17409984

```
python3 ../../rocketchip human_SRR17409984.yaml --output_file human_SRR17409984 --data .
time snakemake -j 1 -s human_SRR17409984
```

### SRR17514595

```
python3 ../../rocketchip human_SRR17514595.yaml --output_file human_SRR17514595 --data .
time snakemake -j 1 -s human_SRR17514595
```

## Mouse Genome

### SRR14407118

```
python3 ../../rocketchip mouse_SRR14407118.yaml --output_file mouse_SRR14407118 --data .
time snakemake -j 1 -s mouse_SRR14407118
```

### SRR17887304

```
python3 ../../rocketchip mouse_SRR17887304.yaml --output_file mouse_SRR17887304 --data .
time snakemake -j 1 -s mouse_SRR17887304
```

### SRR26041601

```
python3 ../../rocketchip mouse_SRR26041601.yaml --output_file mouse_SRR26041601 --data .
time snakemake -j 1 -s mouse_SRR26041601
```

## Rat Genome

### ERR6356096

```
python3 ../../rocketchip rat_ERR6356096.yaml --output_file rat_ERR6356096 --data .
time snakemake -j 1 -s rat_ERR6356096
```

### ERR6356099

```
python3 ../../rocketchip rat_ERR6356099.yaml --output_file rat_ERR6356099 --data .
time snakemake -j 1 -s rat_ERR6356099
```

### ERR6356102

```
python3 ../../rocketchip rat_ERR6356102.yaml --output_file rat_ERR6356102 --data .
time snakemake -j 1 -s rat_ERR6356102
```

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

### SRR15046104

```
python3 ../../rocketchip zebrafish_SRR15046104.yaml --output_file  zebrafish_SRR15046104 --data .
time snakemake -j 1 -s zebrafish_SRR15046104
```

### SRR15509781

```
python3 ../../rocketchip zebrafish_SRR15509781.yaml --output_file  zebrafish_SRR15509781 --data .
time snakemake -j 1 -s zebrafish_SRR15509781
```

### SRR15509782

```
python3 ../../rocketchip zebrafish_SRR15509782.yaml --output_file  zebrafish_SRR15509782 --data .
time snakemake -j 1 -s zebrafish_SRR15509782
```