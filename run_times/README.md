# Run Times

Before you begin, make sure you have installed Rocketchip and that you are in the activated `rocketchip` Conda environment. 

| SRA ID      | Raw Read Count | Genome               | Real Time     | User Time     | System Time     |
| :---------: | :------------: | :------------------: | :-----------: | :-----------: | :-------------: |
| SRR15243009 | 13660923       | Fruitfly (dm6)       | 134m27.029s   | 87m27.505s    | 1m55.765s       |
| SRR16638474 | 3329820        | Fruitfly (dm6)       | 68m14.252s    | 30m12.775s    | 0m51.145s       |
| SRR16638475 | 2880763        | Fruitfly (dm6)       | 33m10.423s    | 26m14.124s    | 0m46.629s       |
| DRR345782   | 22649733       | Human (hg38)         | 428m14.225s   | 336m44.640s   | 9m48.430s       |
| SRR17409984 | 32055273       | Human (hg38)         | 554m36.029s   | 393m17.329s   | 13m27.333s      |
| SRR17514595 | 19332910       | Human (hg38)         | 782m17.619s   | 737m29.775s   | 17m35.879s      | 
| SRR14407118 | 102208451      | Mouse (mm10)         | 1232m38.595s  | 1089m12.640s  | 28m49.801s      |
| SRR17887304 | 20180276       | Mouse (mm10)         | 659m50.839s   | 545m20.644s   | 14m23.511s      |
| SRR26041601 | 3914984        | Mouse (mm10)         | 102m7.136s    | 85m24.658s    | 2m8.867s        |
| ERR6356096  | 25140500       | Rat (rn6)            | 629m4.017s    | 543m30.316s   | 15m31.332s      |
| ERR6356099  | 25802282       | Rat (rn6)            | 627m49.210s   | 558m25.806s   | 16m22.164s      |
| ERR6356102  | 35225234       | Rat (rn6)            | 924m12.801s   | 789m7.768s    | 27m24.284s      |
| SRR13125168 | 12569340       | Worm (ce11)          | 154m43.743s   | 98m54.013s    | 1m55.233s       |
| SRR13125170 | 11144851       | Worm (ce11)          | 150m22.033s   | 94m22.257s    | 1m41.009s       |
| SRR13125172 | 12265487       | Worm (ce11)          | 114m56.070s   | 99m46.430s    | 1m43.391s       |
| SRR17329288 | 76503112       | Yeast (sacCer3)      | 154m27.814s   | 78m15.281s    | 2m46.423s       |
| SRR17329289 | 74716042       | Yeast (sacCer3)      | 149m16.474s   | 73m52.756s    | 2m39.203s       |
| SRR17329314 | 26073747       | Yeast (sacCer3)      | 64m54.372s    | 52m8.798s     | 1m25.590s       |
| SRR15046104 | 19749132       | Zebrafish (danRer11) | 314m12.676s   | 261m21.603s   | 7m1.943s        |
| SRR15509781 | 27448764       | Zebrafish (danRer11) | 613m19.826s   | 546m51.961s   | 14m27.508s      |
| SRR15509782 | 22765779       | Zebrafish (danRer11) | 571m47.449s   | 446m46.142s   | 11m47.989s      |

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
rocketchip fruitfly_SRR15243009.yaml --output_file fruitfly_SRR15243009 --data .
time snakemake -j 1 -s fruitfly_SRR15243009
```

### SRR16638474

```
rocketchip fruitfly_SRR16638474.yaml --output_file fruitfly_SRR16638474 --data .
time snakemake -j 1 -s fruitfly_SRR16638474
```

### SRR16638475

```
rocketchip fruitfly_SRR16638475.yaml --output_file fruitfly_SRR16638475 --data .
time snakemake -j 1 -s fruitfly_SRR16638475
```

## Human Genome

### DRR345782

```
rocketchip human_DRR345782.yaml --output_file human_DRR345782 --data .
time snakemake -j 1 -s human_DRR345782
```

### SRR17409984

```
rocketchip human_SRR17409984.yaml --output_file human_SRR17409984 --data .
time snakemake -j 1 -s human_SRR17409984
```

### SRR17514595

```
rocketchip human_SRR17514595.yaml --output_file human_SRR17514595 --data .
time snakemake -j 1 -s human_SRR17514595
```

## Mouse Genome

### SRR14407118

```
rocketchip mouse_SRR14407118.yaml --output_file mouse_SRR14407118 --data .
time snakemake -j 1 -s mouse_SRR14407118
```

### SRR17887304

```
rocketchip mouse_SRR17887304.yaml --output_file mouse_SRR17887304 --data .
time snakemake -j 1 -s mouse_SRR17887304
```

### SRR26041601

```
rocketchip mouse_SRR26041601.yaml --output_file mouse_SRR26041601 --data .
time snakemake -j 1 -s mouse_SRR26041601
```

## Rat Genome

### ERR6356096

```
rocketchip rat_ERR6356096.yaml --output_file rat_ERR6356096 --data .
time snakemake -j 1 -s rat_ERR6356096
```

### ERR6356099

```
rocketchip rat_ERR6356099.yaml --output_file rat_ERR6356099 --data .
time snakemake -j 1 -s rat_ERR6356099
```

### ERR6356102

```
rocketchip rat_ERR6356102.yaml --output_file rat_ERR6356102 --data .
time snakemake -j 1 -s rat_ERR6356102
```

## Worm Genome

### SRR13125168

```
rocketchip worm_SRR13125168.yaml --output_file worm_SRR13125168 --data .
time snakemake -j 1 -s worm_SRR13125168
```

### SRR13125170

```
rocketchip worm_SRR13125170.yaml --output_file worm_SRR13125170 --data .
time snakemake -j 1 -s worm_SRR13125170
```

### SRR13125172

```
rocketchip worm_SRR13125172.yaml --output_file worm_SRR13125172 --data .
time snakemake -j 1 -s worm_SRR13125172
```

## Yeast Genome 

### SRR17329288

```
rocketchip yeast_SRR17329288.yaml --output_file yeast_SRR17329288 --data .
time snakemake -j 1 -s yeast_SRR17329288
```

### SRR17329289

```
rocketchip yeast_SRR17329289.yaml --output_file  yeast_SRR17329289 --data .
time snakemake -j 1 -s yeast_SRR17329289
```

### SRR17329314

```
rocketchip yeast_SRR17329314.yaml --output_file  yeast_SRR17329314 --data .
time snakemake -j 1 -s yeast_SRR17329314
```

## Zebrafish Genome

### SRR15046104

```
rocketchip zebrafish_SRR15046104.yaml --output_file  zebrafish_SRR15046104 --data .
time snakemake -j 1 -s zebrafish_SRR15046104
```

### SRR15509781

```
rocketchip zebrafish_SRR15509781.yaml --output_file  zebrafish_SRR15509781 --data .
time snakemake -j 1 -s zebrafish_SRR15509781
```

### SRR15509782

```
rocketchip zebrafish_SRR15509782.yaml --output_file  zebrafish_SRR15509782 --data .
time snakemake -j 1 -s zebrafish_SRR15509782
```

## Results

Results were visualized via [visualize_run_times.ipynb](https://github.com/vhaghani26/rocketchip_tests/blob/main/run_times/visualize_run_times.ipynb).
