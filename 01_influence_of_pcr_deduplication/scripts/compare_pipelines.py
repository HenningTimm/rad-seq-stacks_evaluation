"""
"""
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use("agg")
import seaborn as sns
import pandas as pd
import yaml
from collections import Counter
# def tuftefy(ax):
#     """Remove spines and tick position markers to reduce ink."""
#     ax.spines["top"].set_visible(False)
#     ax.spines["right"].set_visible(False)
#     ax.spines["left"].set_visible(False)
#     ax.spines["bottom"].set_visible(True)
#     ax.spines["bottom"].set_color('grey')
#     ax.grid(color="w", alpha=0.7)
#     ax.get_yaxis().grid(True)
#     ax.get_xaxis().grid(False)


# def dat_file_to_array(dat_file, threshold):
#     """Read dat files into a list. Split off everything over
#     threshold into a separate list.
#     """
#     data = []
#     outliers = []
#     for line in dat_file:
#         count = int(line.strip())
#         if count <= threshold:
#             data.append(int(count))
#         else:
#             outliers.append(int(count))
#     return data, outliers

def get_simulated_counts():
    # gt_df = pd.DataFrame(columns=["name", "simulated_loci"])
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

def compare_parameters():
    """Distribution of one or more
    """
    # thin lines
    mpl.rcParams['axes.linewidth'] = 0.5

    df_stack_sizes_no_dedup = pd.read_csv(snakemake.input.sizes_no_dedup)
    df_stack_sizes_with_dedup = pd.read_csv(snakemake.input.sizes_with_dedup)

    df_stack_counts_no_dedup = pd.read_csv(snakemake.input.counts_no_dedup)
    df_stack_counts_with_dedup = pd.read_csv(snakemake.input.counts_with_dedup)
    
    df_stack_sizes_no_dedup["Deduplication"] = "no"
    df_stack_sizes_with_dedup["Deduplication"] = "yes"

    df_stack_counts_no_dedup["Deduplication"] = "no"
    df_stack_counts_with_dedup["Deduplication"] = "yes"
    
    df_sizes = pd.concat([df_stack_sizes_no_dedup, df_stack_sizes_with_dedup])
    df_counts = pd.concat([df_stack_counts_no_dedup, df_stack_counts_with_dedup])

    
    
    # Set up plots to contain the histogram (left) and violin plots (right)
    # in approx 4:3 ratio
    scale = 3
    sns.set(font_scale=1.1, style="whitegrid")
    nr_param_sets = len(set(df_sizes["params"]))
    fig, ax_violin = plt.subplots(
        figsize=(4*scale, 1*nr_param_sets*scale)
    )

    thresholded_df = df_sizes[df_sizes["count"] < snakemake.params.threshold]

    sns.violinplot(data=thresholded_df, x="count", y="params", hue="Deduplication",
                   ax=ax_violin, split=True, b=0.05)

    # ax_violin.xaxis.grid(True)
    ax_violin.set(ylabel="Parameters", xlabel="Coverage")
    sns.despine(ax=ax_violin)
    plt.tight_layout()
    plt.savefig(snakemake.output.violin_path, format="pdf", dpi=300)

    plt.clf()
    scale=2
    fig, ax_scatter = plt.subplots(
        figsize=(4*scale, 10*scale)
    )

    df_counts = df_counts.rename(columns={
        "params": "Parameter set",
        "nr_loci": "Number of detected loci",
        })
    sns.set(font_scale=1.1, style="whitegrid")
    # NOTE: Also get total number of simulated loci
    # sim_counts = get_simulated_counts()
    # print(sim_counts)
    # sns.scatterplot(x="params", y="nr_loci", hue="name", style="Deduplication", data=df_counts)
    cat_plot = sns.catplot(
        "Deduplication", col="Parameter set", y="Number of detected loci", hue="name", data=df_counts,
        kind="point", col_wrap=5, aspect=0.55
    )
    cat_plot._legend.remove()
    plt.tight_layout()
    plt.savefig(snakemake.output.scatter_path, format="pdf", dpi=300)




compare_parameters()
