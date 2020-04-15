# rad-seq-stacks_evaluation
Collection of evaluation workflows for the [rad-seq-stacks](https://github.com/snakemake-workflows/rad-seq-stacks) snakemake pipeline based on the [Stacks](http://catchenlab.life.illinois.edu/stacks/) software.

## Usage
Each folder contains a snakemake workflow that calls two instances of the rad-seq-stacks pipeline (each itself a snakemake workflow) on a simulated dataset.
To perform an evaluation, navigate into the corresponding folder and call
```bash
$ snakemake --use-conda --jobs 6
```
to run the pipeline.
All subworkflows are restricted to 3 cores using a parameter in the `config.yaml` file in the respective folder.
To increase the number of cores for each subworkflow, change the value
```yaml
cores_per_subworkflow: 3
```
and call the workflow with twice this amount, so that each subworkflow can use half of the total assigned cores.

Read data is simulated using [ddRAGE](https://ddrage.readthedocs.io/en/latest/).
To change the number of simulated loci (or other parameters of the simulation), change the respective values in the config file:
```yaml
ddrage:
  loci:
    10000
  individuals:
    # max number possible with standard BC set
    24
```
Note, that smaller instances (less loci, less individuals or less coverage) will execute faster, however some effects observed in the evaluation might not be visible with these parameters.


## Content
This repository contains five evaluations.
The first four use data simulated with our ddRADseq read simulator [ddRAGE](https://ddrage.readthedocs.io/en/latest/) and can be executed as is.
For the fifth evaluation workflow we used an unpublished in-house dataset.
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
