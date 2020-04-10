import itertools
import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import yaml


def get_log_files():
    """Collect ustacks and cstacks log files for all individuals and
    parameter sets in both subworkflows.
    """
    # assert the individuals in both subworksflows are identical
    individuals_with = list(pd.read_csv(snakemake.input.individuals_with, sep="\t")["id"])
    individuals_no = list(pd.read_csv(snakemake.input.individuals_no, sep="\t")["id"])
    assert(individuals_with == individuals_no)
    individuals = individuals_with

    # assert both subworkflows used the same parameters
    stacks_params = None
    with open(snakemake.input.config_with, "r") as config_with_dedup, open(snakemake.input.config_no, "r") as config_no_dedup:
        params_with_dedup = yaml.safe_load(config_with_dedup)["params"]["stacks"]
        params_no_dedup = yaml.safe_load(config_no_dedup)["params"]["stacks"]
        assert(params_with_dedup == params_no_dedup)
        stacks_params = params_with_dedup

    # assemble all cstacks filenames without subfolder
    def cstacks_file(parameter_set):
        return f"n={parameter_set['max_locus_mm']}.M={parameter_set['max_individual_mm']}.m={parameter_set['min_reads']}.log"
    cstacks_filenames = [cstacks_file(param_set) for param_set in stacks_params]

    # assemble all ustacks filenames without subfolder
    def ustacks_file(parameter_set, individual):
        return f"M={parameter_set['max_individual_mm']}.m={parameter_set['min_reads']}/{individual}.log"
    ustacks_filenames = [ustacks_file(param_set, ind) for (param_set, ind) in itertools.product(stacks_params, individuals)]

    # assemble full filenames by adding subfolders
    paths_dict = {
        "logs_cstacks_with_dedup": [f"rad-seq-stacks_with_dedup/logs/cstacks/{cf}" for cf in cstacks_filenames],
        "logs_cstacks_no_dedup": [f"rad-seq-stacks_no_dedup/logs/cstacks/{cf}" for cf in cstacks_filenames],
        "logs_ustacks_with_dedup": [f"rad-seq-stacks_with_dedup/logs/ustacks/{uf}" for uf in ustacks_filenames],
        "logs_ustacks_no_dedup": [f"rad-seq-stacks_no_dedup/logs/ustacks/{uf}" for uf in ustacks_filenames],
        }
    return paths_dict



def parse_cstacks_file(text):
    """Parse a cstacks log file to extract the number of new loci added per sample.
    Return results as a dataframe.
    """
    segments = text.split("\n\n")
    newly_added_loci = []
    matched_to_known = []
    ind_nrs = []
    total = None

    # A segment is a group of lines in the log file that describe a certain action
    # this can be initialization, adding a sample to the catalog (processing),
    # or writing output.
    for segment in segments:
        
        if segment.startswith("Initializing"):
            nr = 1
            newly_added = int(re.search("(\d+) loci were newly added to the catalog", segment).groups()[0])
            ind_nr = re.search("/(Individual_\d+)\.", segment).groups()[0]
            ind_nrs.append(ind_nr)
            newly_added_loci.append(newly_added)
        elif segment.startswith("Processing"):
            nr, of = re.search("\[(\d+) of (\d+)\]", segment).groups()
            nr = int(nr)
            of = int(of)
            ind_nr = re.search("/(Individual_\d+)\.", segment).groups()[0]
            ind_nrs.append(ind_nr)
            # matched_to_known = int(re.search("(\d+) loci were matched to a catalog locus", segment).groups()[0])
            newly_added = int(re.search("(\d+) loci were newly added to the catalog", segment).groups()[0])
            newly_added_loci.append(newly_added)
        elif segment.startswith("Writing"):
            # total = int(re.search("Final catalog contains (\d+) loci.", segment).groups()[0])
            params = re.search("directory 'stacks/(.+)/'", segment).groups()[0]

    df = pd.DataFrame(
        {
            "newly_added_loci": newly_added_loci,
            "loci_in_catalog": np.cumsum(newly_added_loci),
            "Individual": ind_nrs,
        }
    )
    df["params"] = params

    return df


def parse_ustacks_file(text):
    """Parse a ustacks file to extract the number of blacklisted loci (HRLs)
    and the mean size of initial clusters.
    Return resulkts as a pandas data frame.
    """
    individual = re.search("(Individual_(\d+)).fq.gz", text).groups()[0]
    blacklisted = int(re.search("Blacklisted (\d+) stacks.", text).groups()[0])
    initial_mean_cov = float(re.search("Stack coverage: mean=(\d+.\d+);", text).groups()[0]) 

    df = pd.DataFrame(
        {
            "Individual": [individual],
            "Blacklisted Loci": [blacklisted],
            "Initial Mean Coverage": [initial_mean_cov],
        }
    )
    return df


def get_dataframe():
    """Parse all log files into dataframes and join them into one large dataframe.
    """
    log_files = get_log_files()
    cstacks_dataframes = []
    ustacks_dataframes = []
    for f in log_files["logs_cstacks_with_dedup"]:
        with open(f, "r") as log_file:
            df = parse_cstacks_file(log_file.read())
            df["dedup"] = "yes"
            cstacks_dataframes.append(df)

    for f in log_files["logs_cstacks_no_dedup"]:
        with open(f, "r") as log_file:
            df = parse_cstacks_file(log_file.read())
            df["dedup"] = "no"
            cstacks_dataframes.append(df)

    for f in log_files["logs_ustacks_with_dedup"]:
        with open(f, "r") as log_file:
            df = parse_ustacks_file(log_file.read())
            df["dedup"] = "yes"
            df["params"] = re.search("/(M=\d+\.m=\d+)/", f).groups()[0]
            ustacks_dataframes.append(df)

    for f in log_files["logs_ustacks_no_dedup"]:
        with open(f, "r") as log_file:
            df = parse_ustacks_file(log_file.read())
            df["dedup"] = "no"
            df["params"] = re.search("/(M=\d+\.m=\d+)/", f).groups()[0]
            ustacks_dataframes.append(df)

    return pd.concat(cstacks_dataframes), pd.concat(ustacks_dataframes)



def plot_accumulating_loci(df):
    """Plot the acuumulation of loci in the catalog with subsequent steps of cstacks.
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set(yscale="log")
    sns.lineplot(x="Individual", y="loci_in_catalog", style="dedup",
                 hue="params", markers=True, alpha=0.7, data=df, ax=ax)
    # NOTE: use number of individuals instead of 23?
    # use expected number of loci instead of 10000
    ax.plot([0, 23], [10000, 10000], linestyle=":", color="gray")
    ax.set_xticklabels(
        [str(i) for i, _ in enumerate(ax.get_xticklabels(), 1)]
    )
    ax.set_xlabel("Nr. of samples added to catalog (cumulative)")
    ax.set_ylabel("Nr. of loci in catalog")

    sns.despine()
    plt.tight_layout()
    plt.savefig(snakemake.output.accumulation, dpi=300)
    plt.clf()



def plot_endpoints(df):
    fig, (ax_no, ax_with) = plt.subplots(2,1, sharex=True, figsize=(8,12),)

    # ax_no.set(yscale="log")
    # ax_with.set(yscale="log")
    only_96 = df[df["Individual"] == "Individual_96"]

    sns.scatterplot(x="Individual", y="loci_in_catalog", style="dedup", hue="params", markers="X", alpha=0.7, data=only_96[only_96["dedup"] == "no"], ax=ax_no)
    sns.scatterplot(x="Individual", y="loci_in_catalog", style="dedup", hue="params", markers=True, alpha=0.7, data=only_96[only_96["dedup"] == "yes"], ax=ax_with)
    ax_with.plot([-0.1,0.1], [10000, 10000], linestyle=":", color="gray")
    # plt.xticks(rotation=45)
    # lp.set_xticklabels(rotations=30)
    sns.despine()
    plt.subplots_adjust(wspace=0, hspace=0.1)
    ax_with.set_xlabel("Catalog containing loci from all samples")
    ax_with.set_xticklabels([])
    ax_no.set_ylabel("Nr. of loci")
    ax_with.set_ylabel("Nr. of loci")
    # plt.show()
    plt.tight_layout()
    plt.savefig(snakemake.output.endpoints, dpi=300)

    plt.clf()


def plot_blacklisted_and_coverage(df):
    """Plot the number of HRl loci (blacklisted) and mean size as swarm plots.
    """
    # sort dedup calues so that they match the violin plot from before
    df = df.sort_values(by='dedup', ascending=True)

    df = df.rename(columns={
        "params": "Parameter set",
        "dedup": "Deduplication",
    })

    sns.set(style="whitegrid")

    fig, (ax_blacklisted, ax_mean) = plt.subplots(1,2, sharex=True, figsize=(10, 5),)
    p1 = sns.swarmplot(x="Parameter set", y="Blacklisted Loci", hue="Deduplication", data=df, ax=ax_blacklisted, size=4)
    p1.get_legend().set_visible(False)
    p2 = sns.swarmplot(x="Parameter set", y="Initial Mean Coverage", hue="Deduplication", data=df, ax=ax_mean, size=4)

    p2.get_legend().set_bbox_to_anchor((1.3, 1.0))
    ax_blacklisted.plot([-0.5, 4+0.5], [500, 500], color="gray", linestyle=":")

    sns.despine()
    plt.tight_layout()
    plt.savefig(snakemake.output.blacklisted, dpi=300)
    plt.clf()


# read in log files
cstacks_df, ustacks_df = get_dataframe()


plot_accumulating_loci(cstacks_df)
plot_endpoints(cstacks_df)
plot_blacklisted_and_coverage(ustacks_df)
