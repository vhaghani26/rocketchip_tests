# Yeast Genome Tests

Before you begin, make sure you have exported `ROCKETCHIP_SRC` and/or `ROCKETCHIP_DATA` or add them onto the ends of the commands as arguments. Also ensure that you are in the activated `rocketchip` Conda environment.

## `yeast_SRR17329288.yaml`

Enter `yeast_SRR17329288/`. Because `yeast_SRR17329288.yaml` uses data from the NCBI SRA, the only thing that needs to be changed is the author name. Once you have edited the file, you can run the following:

```
rocketchip yeast_SRR17329288.yaml --output_file Snakefile
```

This will create a file called Snakefile. To run it, simply run:

```
snakemake -j 1
```