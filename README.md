# Rocketchip Tests

## Overview

This repository contains the project yaml files and scripts required to run the different experiments outlined in the original Rocketchip manuscript. Each subdirectory contains explicit instructions and should be run in the directories as they are organized here. 

In order to replicate these tests, please install and set up [Rocketchip](https://github.com/vhaghani26/rocketchip) first. The figures were generated using a different Conda environment (`figure_env.yml`) than the actual Rocketchip environment.

Feel free to contact me at [vhaghani@ucdavis.edu](vhaghani@ucdavis.edu) if you have any questions, comments, or concerns, and I will get back to you as soon as possible.

## Extended Data

The [`Extended_Data/` subdirectory](https://github.com/vhaghani26/rocketchip_tests/tree/main/Extended_Data) contains the supplementary figures and tables associated with the Rocketchip manuscript.

## [`chip_seq/`](https://github.com/vhaghani26/rocketchip_tests/tree/main/chip_seq)

In order to demonstrate Rocketchip's utility in replicating experimental results and running all software combinations, we ran ChIP-seq data from the study ["Sequence features accurately predict genome-wide MeCP2 binding in vivo"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4820824/) by Rube et al. 2016 and ["Genome-wide global identification of NRF2 binding sites in A549 non-small cell lung cancer cells by ChIP-Seq reveals NRF2 regulation of genes involved in focal adhesion pathways"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6949066/#SD1) by Namani et al. 2019 using Rocketchip. This subdirectory contains all necessary analysis scripts and instructions associated with this experiment.

## [`cut_and_tag_run/`](https://github.com/vhaghani26/rocketchip_tests/tree/main/cut_and_tag_run)

In order to demonstrate that Rocketchip can be used to analyze CUT&Tag and CUT&RUN data, we ran data from the study ["Identification of chromatin states during zebrafish gastrulation using CUT&RUN and CUT&Tag"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8976701/) by Akdogan-Ozdilek et al. 2023 through Rocketchip. This subdirectory contains all necessary analysis scripts and instructions associated with this experiment.

## [`functional_test/`](https://github.com/vhaghani26/rocketchip_tests/tree/main/functional_test)

The purpose of this functional test is to ensure that all combinations of software are capable of running without compatibility issues. I intend to employ this functional test between updates to Rocketchip's software to identify possible bugs in the software or software intercompatibiity problems. This is not included in the manuscript; it is simply maintained separately from the source code repository.

## [`run_times/`](https://github.com/vhaghani26/rocketchip_tests/tree/main/run_times)

This experiment served dual purposes: (1) demonstrate Rocketchip's compatibility and use for various genomes and (2) give users a rough estimate for how long it takes Rocketchip to run based on different input data and genome sizes. This subdirectory contains all necessary analysis scripts and instructions associated with this experiment.

## Synthetic Data Sets ([`exp_vs_obs/`](https://github.com/vhaghani26/rocketchip_tests/tree/main/exp_vs_obs))

The synthetic sequence data sets, all analysis scripts, and usage instructions related to the synthetic datasets can be found in the [`exp_vs_obs/` subdirectory](https://github.com/vhaghani26/rocketchip_tests/tree/main/exp_vs_obs). This dataset was used for benchmarking algorithms in Rocketchip and evaluating its accuracy, precision, and sensitivity.

## Publicly Available Data Usage

All publicly available data we used in this study can be found using the following GEO Accession Numbers: GSE71126, GSE141497, GSE178343, GSE178559, GSE180812, GSE186853, GSE192929, GSE193317, GSE173740, GSE196170, GSE242974, GSE162063, GSE192592, and GSE182346 as well as the following BioProject Accession Numbers: PRJDB13012, PRJEB46586, and PRJNA657343.