# Assessment of Synthetic Input

To ensure that the synthetic input control is capable of accurately modeling true input controls, we ran peak-calling for each input in two scenarios: (1) peak-calling for the input only (i.e. seeing the number of peaks in the input, but not running it with a control) and (2) peak-calling for the input against itself to ...(FIGURE OUT HOW TO PHRASE THIS GOAL)

https://academic.oup.com/bioinformatics/article/38/9/2617/6535233 for coverage 

| Sample ID              | Coverage         | Peak Number for Input Only | Peak Number for Input vs. Input | Study              | Genome    |
| :-------------:        | :--------------: | :------------------------: | :-----------------------------: | :----------------: | :-------: |
| paired_narrow_1        |                  | 0                          | 0                               | NA                 | Synthetic |
| paired_narrow_2        |                  | 0                          | 0                               | NA                 | Synthetic |
| paired_narrow_3        |                  | 0                          | 0                               | NA                 | Synthetic |
| paired_narrow_4        |                  | 0                          | 0                               | NA                 | Synthetic |
| paired_narrow_5        |                  | 0                          | 0                               | NA                 | Synthetic |
| paired_narrow_6        |                  | 0                          | 0                               | NA                 | Synthetic |
| paired_broad_1         |                  | 0                          | 0                               | NA                 | Synthetic |
| paired_broad_2         |                  | 0                          | 0                               | NA                 | Synthetic |
| paired_broad_3         |                  | 0                          | 0                               | NA                 | Synthetic |
| paired_broad_4         |                  | 0                          | 0                               | NA                 | Synthetic |
| paired_broad_5         |                  | 0                          | 0                               | NA                 | Synthetic |
| paired_broad_6         |                  | 0                          | 0                               | NA                 | Synthetic |
| single_narrow_1        |                  | 0                          | 0                               | NA                 | Synthetic |
| single_narrow_2        |                  | 0                          | 0                               | NA                 | Synthetic |
| single_narrow_3        |                  | 0                          | 0                               | NA                 | Synthetic |
| single_narrow_4        |                  | 0                          | 0                               | NA                 | Synthetic |
| single_narrow_5        |                  | 0                          | 0                               | NA                 | Synthetic |
| single_narrow_6        |                  | 0                          | 0                               | NA                 | Synthetic |
| single_broad_1         |                  | 0                          | 0                               | NA                 | Synthetic |
| single_broad_2         |                  | 0                          | 0                               | NA                 | Synthetic |
| single_broad_3         |                  | 0                          | 0                               | NA                 | Synthetic |
| single_broad_4         |                  | 0                          | 0                               | NA                 | Synthetic |
| single_broad_5         |                  | 0                          | 0                               | NA                 | Synthetic |
| single_broad_6         |                  | 0                          | 0                               | NA                 | Synthetic |
| SRR2119603, SRR2119604 | 2.01286, 1.70862 | 7349                       | 0                               | Rube et al. 2016   | mm10      |
| SRR9624455             | 7.96487          | 13338                      | 0                               | Rieder et al. 2019 | dm6       |
| SRR9624459             | 2.37311          | 4458                       | 0                               | Rieder et al. 2019 | dm6       |
| SRR9624464, SRR9624465 | 4.29152, 2.97774 | 11147                      | 0                               | Rieder et al. 2019 | dm6       |
| SRR9624471, SRR9624472 | 1.87005, 2.8268  | 8305                       | 0                               | Rieder et al. 2019 | dm6       |
| SRR9624480, SRR9624481 | 4.5775, 4.23696  | 34216                      | 0                               | Rieder et al. 2019 | dm6       |
| SRR9624487, SRR9624488 | 4.84733, 5.31379 | 15490                      | 0                               | Rieder et al. 2019 | dm6       |
| SRR9624496, SRR9624497 | 3.31872, 3.25337 | 17921                      | 0                               | Rieder et al. 2019 | dm6       |
| SRR9624503, SRR9624504 | 4.26115, 2.26914 | 11940                      | 0                               | Rieder et al. 2019 | dm6       |

## Drosophila Data

```
cd drosophila_inputs
rocketchip drosophila_inputs_ctrl.yaml --output_file drosophila_inputs_ctrl --data .
snakemake -j 10 -s drosophila_inputs_ctrl
rocketchip drosophila_inputs_no_ctrl.yaml --output_file drosophila_inputs_no_ctrl --data .
snakemake -j 10 -s drosophila_inputs_no_ctrl
cd 06_macs3_peaks
for file in *.narrowPeak; do filename=$(basename "$file"); line_count=$(wc -l < "$file"); echo "$filename: $line_count"; done

```

## Mouse Data

```
cd mouse_inputs
rocketchip mouse_inputs_ctrl.yaml --output_file mouse_inputs_ctrl --data .
snakemake -j 10 -s mouse_inputs_ctrl
rocketchip mouse_inputs_no_ctrl.yaml --output_file mouse_inputs_no_ctrl --data .
snakemake -j 10 -s mouse_inputs_no_ctrl
cd 06_macs3_peaks
for file in *.narrowPeak; do filename=$(basename "$file"); line_count=$(wc -l < "$file"); echo "$filename: $line_count"; done
```

## Synthetic Data

```
cd synthetic_inputs
python3 synthetic_inputs.py
```

## Calculating Coverage

```
for file in *.sorted.dedup.bam; do
    # Extract sample names
    sample=$(basename "$file" .sorted.dedup.bam)
    
    # Run samtools depth to calculate average depth
    avg_depth=$(samtools depth "${sample}.sorted.dedup.bam" | awk '{sum+=$3} END {print sum/NR}')

    # Print sample name and average depth
    echo "Sample: $sample, Average Depth: $avg_depth"
done
```
