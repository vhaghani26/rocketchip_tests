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