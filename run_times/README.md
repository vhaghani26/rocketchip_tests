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

## Fruitfly Genome

## Human Genome

## Mouse Genome

## Rat Genome

## Worm Genome

## Yeast Genome 

### `yeast_SRR17329288.yaml`

Enter `yeast_SRR17329288/`. Because `yeast_SRR17329288.yaml` uses data from the NCBI SRA, the only thing that needs to be changed is the author name. Once you have edited the file, you can run the following:

```
rocketchip yeast_SRR17329288.yaml --output_file Snakefile
```

This will create a file called Snakefile. To run it, simply run:

```
snakemake -j 1
```

## Zebrafish Genome
