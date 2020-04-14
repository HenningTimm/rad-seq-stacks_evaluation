# rad-seq-stacks_evaluation
Collection of evaluation workflows for the [rad-seq-stacks](https://github.com/snakemake-workflows/rad-seq-stacks) snakemake pipeline based on the [Stacks](http://catchenlab.life.illinois.edu/stacks/) software.


## Content
This repository contains five evaluations.
The first four use data simulated with our ddRADseq read simulator [ddRAGE](https://ddrage.readthedocs.io/en/latest/) and can be executed as is.
For the fifth evaluation workflow we used an unbuplished in-house dataset.
Hence, this workflow cannot be executed without access to the data set.
Once the data is published, we will add the link.

### 01_influence_of_pcr_deduplication
Illustrates the influence of PCR deduplication on the workflow.

### 02_influence_of_minimum_reads_per_locus
Simulates a low coverage dataset to analyze the impact of the minimum reads per locus parameter.

### 03_influence_of_locus_distance
Simulates a high diversity dataset to analyze the impact of locus distance parameters.

### 04_detection_of_snps
Simulates a data set with an increased number of mutations (all SNPs) to analyze the performance of different parameter sets for SNP detection.

### 05_analyzing_a_large_dataset
Analysis workflow for a dataset with 315 individuals of *Gammarus fossarum*.
This dataset is not included, since it is not yet published.
Once it is published, we will link to it here.
