[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4420213.svg)](https://doi.org/10.5281/zenodo.4420213)

# rad-seq-stacks_evaluation
Collection of evaluation workflows for the [rad-seq-stacks](https://github.com/snakemake-workflows/rad-seq-stacks) snakemake workflow based on the [Stacks](http://catchenlab.life.illinois.edu/stacks/) software. Notice that this evaluation was performed for my PhD thesis. For the most recent version of this analysis workflow, please refer to the [rad-seq-stacks workflow](https://github.com/snakemake-workflows/rad-seq-stacks) that is now part of the Snakemake Workflows project.


## Requirements
The pipelines require a snakemake version of 5.10.0 or above and an installation of the [conda](https://docs.conda.io/en/latest/miniconda.html) package manager to run (we recommend using the miniconda Python 3.7 installer).
Lower versions result in a crash due to a `WorkflowError`.
Assuming conda is already installed (for installation instructions please refer to the [(mini)conda](https://docs.conda.io/en/latest/miniconda.html) website), the required version of snakemake can be installed into a conda environment as follows:
```bash
$ conda create -n snakemake-env "snakemake>=5.10.0" -c bioconda -c conda-forge
```
After installation, activate the environment using:
```bash
$ conda activate snakemake-env
(snakemake-env)$ snakemake --version
5.10.0
```

## Usage
Each folder contains a snakemake workflow that calls two instances of the rad-seq-stacks pipeline (each itself a snakemake workflow) on a simulated dataset.
To perform an evaluation, navigate into the corresponding folder and call
```bash
(snakemake-env)$ snakemake --use-conda --jobs 6
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
Analysis workflow for a dataset with 315 individuals of *Gammarus fossarum* with a total size of 103GB of gzipped FASTQ reads.
This dataset is not included, since it is not yet published by its owners.
Once it is published, we will provide a link to it here.
A list of individual names (id), barcoding information (p5 and p7 barcode), spacer lengths, and file names are contained in the respective [individuals.tsv](https://github.com/HenningTimm/rad-seq-stacks_evaluation/blob/master/05_analyzing_a_large_dataset/pipeline_with_dedup/individuals.tsv) and [units.tsv](https://github.com/HenningTimm/rad-seq-stacks_evaluation/blob/master/05_analyzing_a_large_dataset/pipeline_with_dedup/units.tsv) files. Used RADseq enzymes are documented in the [config.yaml](https://github.com/HenningTimm/rad-seq-stacks_evaluation/blob/master/05_analyzing_a_large_dataset/pipeline_with_dedup/config.yaml) file. The used DBR sequence is `NNNNNNMMGGACG`.
