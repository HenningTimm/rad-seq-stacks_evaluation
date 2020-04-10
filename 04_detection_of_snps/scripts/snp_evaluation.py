import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import re


def get_data():
    # This data has been manually extracted from the logs of an
    # analysis run on compute (large)
    data = pd.DataFrame(
        columns=["Parameters", "Deduplication", "False Positive SNPs", "False Negative SNPs", "True Positive SNPs", "Total SNPs", "Split Loci", "Missed Loci"],
        )
    
    with open("../../../../pipelines/rss_snps/validation/final_eval.txt", "r") as results:
        for chunk in results.read().split("\n\n\n"):

            print("Chunk", chunk, "\n")
            if not chunk:
                continue
            dedup, n, M, m = re.search("Analyzed (with|no)_dedup\/n=(\d+)\.M=(\d+)\.m=(\d+).", chunk).groups()
            # n, M, m = (int(x) for x in (n, M, m))

            data_point = pd.DataFrame(
                {
                    "Parameters": f"n={n} M={M} m={m}",
                    "Deduplication": "yes" if dedup == "with" else "no",
                    "False Positive SNPs": int(re.search("Misidentified SNPs \(FP\): +(\d+)", chunk).groups()[0]),
                    "False Negative SNPs": int(re.search("Missed SNPs \(FN\): +(\d+)", chunk).groups()[0]),
                    "True Positive SNPs": int(re.search("Corectly called SNPs \(TP\): +(\d+)", chunk).groups()[0]),
                    "Total SNPs": int(re.search("Total SNPs: +(\d+)", chunk).groups()[0]),
                    "Split Loci": int(re.search("Split loci: +(\d+)", chunk).groups()[0]),
                    "Missed Loci": int(re.search("Missed loci: +(\d+)", chunk).groups()[0]),
                },
                index=[(1000 if dedup == "with" else 0) + int(n)]
            )
            # correctly_classified_loci = int(re.search("Correctly classified loci: +(\d+)", chunk).groups()[0])
            data = data.append(data_point)
    # print("Data:")
    # print(data)

    # data["Sensitivity"] = data["True Positive SNPs"] / (data["True Positive SNPs"] + data["False Negative SNPs"])
    data["Sensitivity"] = data["True Positive SNPs"] / data["Total SNPs"]

    data = data.sort_index()
    # print("Data:")
    # print(data)
    return data

def plot_snp_evaluation(out_path="snp_evaluation.pdf"):
    df = get_data()

    # fig, (ax_fp_fn, ax_tp, ax_loci) = plt.subplots(
    #     3,
    #     1,
    #     figsize=(12, 8),
    # )
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.6)
    # first plot: FP and FN
    g = sns.catplot(x="Parameters", y="False Positive SNPs", hue="Deduplication",
                    data=df,
                    height=4, aspect=3.5,
                    kind="point",
                    markers=["x", "x"],
                    legend_out=False
    )
    g.ax.set_xlabel('')
    sns.despine()
    plt.tight_layout()
    plt.savefig("FP_SNPs.pdf", dpi=3.500)
    plt.clf()

    # Second plot: TP
    g = sns.catplot(x="Parameters", y="True Positive SNPs", hue="Deduplication",
                    data=df,
                    height=4, aspect=3.5,
                    kind="point",
                    legend_out=False)
    
    mean = df["Total SNPs"].mean()
    g.ax.axhline(mean, ls="--", color="gray")
    # g.ax.legend(loc='upper right')
    g.ax.set_xlabel('')
    sns.despine()
    plt.tight_layout()
    plt.savefig("TP_SNPs.pdf", dpi=3.500)
    plt.clf()

    g = sns.catplot(x="Parameters", y="Sensitivity", hue="Deduplication",
                    data=df,
                    height=4, aspect=3.5,
                    kind="point",
                    legend_out=False
    )
    g.ax.set_xlabel('')
    sns.despine()
    plt.tight_layout()
    plt.savefig("Sensitivity.pdf", dpi=3.500)
    plt.clf()
    
    # Third plot: Split Loci, Missed Loci
    g = sns.catplot(x="Parameters", y="Split Loci", hue="Deduplication",
                    data=df,
                    height=4, aspect=3.5,
                    kind="point",
                    markers=["x", "x"],
                    legend_out=False
    )
    sns.despine()
    plt.tight_layout()
    plt.savefig("Split_Loci.pdf", dpi=300)
    plt.clf()



if __name__ == "__main__":
    plot_snp_evaluation()
