"""Compare number and coverage of loci for different Stacks parameter sets.

Generate a violin plot of locus coverages and a point plot of locus number per individual.
Both plots compare runs with and without deduplication.
"""
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use("agg")
import seaborn as sns
import pandas as pd
import yaml
from collections import Counter


def get_simulated_counts():
    """Read in number of simulated loci from ddRAGE ground truth file.
    """
    locus_counts = Counter()
    with open(snakemake.input.gt, "r") as gt:
        try:
            # read all documents in the data
            inds, loci, *other = list(yaml.load_all(gt, Loader=yaml.FullLoader))
        except yaml.YAMLError as exc:
            print(exc)
    for loocus_id, locus in loci.items():
        for individual, alleles in locus["individuals"].items():
            if alleles:
                locus_counts[individual] += 1
    return locus_counts


def compare_pipeline_results():
    """Distribution of one or more
    """
    # thin lines
    mpl.rcParams['axes.linewidth'] = 0.5

    df_stack_sizes_no_dedup = pd.read_csv(snakemake.input.sizes_no_dedup)
    df_stack_sizes_with_dedup = pd.read_csv(snakemake.input.sizes_with_dedup)

    df_stack_counts_no_dedup = pd.read_csv(snakemake.input.counts_no_dedup)
    df_stack_counts_with_dedup = pd.read_csv(snakemake.input.counts_with_dedup)
    
    df_stack_sizes_no_dedup["dedup"] = "no"
    df_stack_sizes_with_dedup["dedup"] = "yes"

    df_stack_counts_no_dedup["dedup"] = "no"
    df_stack_counts_with_dedup["dedup"] = "yes"
    
    df_sizes = pd.concat([df_stack_sizes_no_dedup, df_stack_sizes_with_dedup])
    df_counts = pd.concat([df_stack_counts_no_dedup, df_stack_counts_with_dedup])
    
    # Set up plots to contain the histogram (left) and violin plots (right)
    # in approx 4:3 ratio
    scale = 1
    nr_param_sets = len(set(df_sizes["params"]))
    fig, ax_violin = plt.subplots(
        figsize=(12*scale, 6*nr_param_sets*scale)
    )

    thresholded_df = df_sizes[df_sizes["count"] < snakemake.params.threshold]

    sns.violinplot(data=thresholded_df, x="count", y="params", hue="dedup",
                   ax=ax_violin, split=True, b=0.05)

    ax_violin.xaxis.grid(True)
    
    plt.tight_layout()
    plt.savefig(snakemake.output.violin_path, format="pdf", dpi=300)

    plt.clf()
    fig, ax_scatter = plt.subplots(
        figsize=(4*scale, 8*scale)
    )

    sns.catplot(
        "dedup", col="params", y="nr_loci", hue="name", data=df_counts,
        kind="point", col_wrap=3
    )
    plt.savefig(snakemake.output.scatter_path, format="pdf", dpi=300)


compare_pipeline_results()
