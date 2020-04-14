"""Plot the distributions of locus sizes and locus numbers between
different parameter sets.
"""
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


def compare_results_between_runs():
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

    ###########################################################################
    # Locus size violin plots
    ###########################################################################
    scale = 4
    sns.set(font_scale=1.5, style="whitegrid")
    nr_param_sets = len(set(df_sizes["params"]))
    fig, ax_violin = plt.subplots(
        figsize=(4*scale, 2*nr_param_sets*scale)
    )

    threshold = int(snakemake.wildcards.threshold)
    thresholded_df = df_sizes[df_sizes["count"] < threshold]
    thresholded_df["params"] = pd.Categorical(
        thresholded_df["params"],
        categories=[
            "n=2.M=2.m=3",
            "n=4.M=3.m=3",
            "n=5.M=5.m=3",
        ],
        ordered=True)
    
    g = sns.violinplot(
        data=thresholded_df, x="count", y="params",
        hue="Deduplication",
        ax=ax_violin, split=True,
        b=0.05,
    )
    ax_violin.xaxis.set_ticks([3] + list(np.arange(0, int(snakemake.wildcards.threshold) + 1, int(snakemake.wildcards.threshold) / 10)))
    
    ax_violin.set(ylabel="Parameters", xlabel="Coverage")
    sns.despine(ax=ax_violin)
    plt.tight_layout()
    plt.savefig(snakemake.output.violin_path, format="pdf", dpi=300)

    ###########################################################################
    # Locus number pointsplots
    ###########################################################################
    plt.clf()
    scale = 2
    fig, ax_scatter = plt.subplots(
        figsize=(10*scale, 8*scale)
    )

    df_counts["params"] = pd.Categorical(df_counts["params"],
                                         categories=[
                                             "n=2.M=2.m=3",
                                             "n=4.M=3.m=3",
                                             "n=5.M=5.m=3",
                                         ],
                                         ordered=True)

    df_counts = df_counts.rename(columns={
        "params": "Parameter set",
        "nr_loci": "Number of detected loci",
        })
    sns.set(font_scale=1.1, style="whitegrid")

    cat_plot = sns.catplot(
        "Deduplication", col="Parameter set", y="Number of detected loci",
        hue="name", data=df_counts,
        kind="point", col_wrap=3, height=10, aspect=0.6
    )
    for ax in cat_plot.axes:
        plt.setp(ax.lines, alpha=.3)
        ax.set_ylim(10000, 81000)
    cat_plot._legend.remove()
    plt.tight_layout()
    plt.savefig(snakemake.output.point_path, format="pdf", dpi=300)

    ###########################################################################
    # Locus number scatter plots
    ###########################################################################
    plt.clf()
    scale = 2
    fig, ax_scatter = plt.subplots(
        figsize=(10*scale, 8*scale)
    )

    sns.set(font_scale=1.0, style="whitegrid")

    pivoted_counts = pd.pivot_table(
        df_counts,
        index=["name", "Parameter set"],
        columns="Deduplication",
        values="Number of detected loci"
    ).reset_index()

    pivoted_counts = pivoted_counts.rename(columns={
        "yes": "Nr. of loci with deduplication",
        "no": "Nr. of loci without deduplication",
        })
    cat_plot = sns.relplot(
        x="Nr. of loci without deduplication", col="Parameter set",
        y="Nr. of loci with deduplication", data=pivoted_counts,
        kind="scatter", col_wrap=3, height=5, alpha=0.75
    )

    for ax in cat_plot.axes:
        ax.set_xlim(10000, 81000)
        ax.set_ylim(10000, 81000)
        ax.plot([10000, 81000], [10000, 81000], color="gray")
    plt.tight_layout()
    plt.savefig(snakemake.output.scatter_path, format="pdf", dpi=300)


compare_results_between_runs()
