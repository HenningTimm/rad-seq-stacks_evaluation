def get_cargopath(w, input):
    """Check if called from within the .test folder to get the right path.
    """
    if ".test" in os.path.abspath(input.tsv):
        return "../scripts/stacks_analyzer/Cargo.toml"
    else:
        return "scripts/stacks_analyzer/Cargo.toml"


def get_plot_script(w, input):
    if isinstance(input.dats, (list, snakemake.io.Namedlist)):
        path = input.dats[0]
    else:
        path = input.dats
    prefix = "../" if ".test" in os.path.abspath(path) else ""
    return prefix + "scripts/plot_stack_sizes.py"


# compute the lengths of stacks loci by counting lines in the tags.tsv files
rule count:
    input:
        tsv="stacks/{parameter_set}/{individual}.tags.tsv.gz"
    output:
        dat="counts/{parameter_set}/{individual}.dat"
    params:
        cargo_path=get_cargopath
    conda:
        "../envs/rust.yaml"
    shell:
        """
        cargo run --release --manifest-path={params.cargo_path} -- count-locus-sizes {input.tsv} --dat {output.dat}
        """


# concatenate all plots for one parameter set into one pdf document
rule assemble_report:
    input:
        plots=expand("plots/{{parameter_set}}/{individual}.pdf",
                     individual=individuals.id,
    )
    output:
        pdf="plots/stacks_size_distribution_{parameter_set}.pdf"
    conda:
        "../envs/pdfunite.yaml"
    shell:
        "pdfunite {input.plots} {output.pdf}"


def get_all_parameter_sets():
    all_sets = []
    for parameters in config["params"]["stacks"]:
        all_sets.append(f"n={parameters['max_locus_mm']}.M={parameters['max_individual_mm']}.m={parameters['min_reads']}")
    return all_sets


rule plot_comparison:
    input:
        dats=expand("counts/{parameter_set}/{individual}.dat",
                    parameter_set=get_all_parameter_sets(),
                    individual=individuals.id
        ),
    output:
        violin_pdf="plots/distribution_comparison/stacks_size_distribution.pdf",
        scatter_pdf="plots/distribution_comparison/stacks_counts.pdf",
        sizes_dataframe="plots/distribution_comparison/sizes_dataframe.csv",
        counts_dataframe="plots/distribution_comparison/counts_dataframe.csv",
    conda:
        "../envs/plot_stacks_dist.yaml"
    params:
        plot_script=get_plot_script,
        threshold=100
    shell:
        """
        python {params.plot_script} compare_parameters {input.dats} --threshold {params.threshold}  --violin-path {output.violin_pdf} --scatter-path {output.scatter_pdf} --sizes-df-path {output.sizes_dataframe} --counts-df-path {output.counts_dataframe}
        """
