import pandas as pd

df = pd.read_csv(snakemake.input.bc_file, sep="\t", keep_default_na = False)

# join barcode and spacer
df["p5_stacks_bc"] = df["p5 bc"] + df["p5 spc"]
# replace spaces in individual names
df["id"] = df["#  Ind."].map(lambda x: x.replace(" ", "_"))
# compute length of p7 spacer
df["p7_spacer_len"] = df["p7 spc"].map(lambda x: len(x))
# set unit to 1, since only one p7 barcode is used here.
# for more p7 bcs this needs to be split up
df["unit"] = 1
# select all important columns for the individuals_file
df_ind = df.loc[:, ["id", "unit", "p5_stacks_bc"]]
# rename columns to match individuals.tsv template
df_ind.columns = ["id", "unit", "p5_barcode"]
df_ind.to_csv(snakemake.output.individuals, sep="\t", index=False)


# select important columns
df_units = df.loc[:, ["unit", "p7 bc", "p7_spacer_len"]]
# remove duplictaes. In this case only one columns should remain
# since only one p7 bc was used
df_units = df_units.drop_duplicates()
# insert input files from snakemake input
df_units["fq1"] = "../" + snakemake.input.fq1
df_units["fq2"] = "../" + snakemake.input.fq2
# rename columns to match units.tsv format
df_units.columns = ["id", "p7_barcode", "p7_spacer", "fq1", "fq2"]
df_units.to_csv(snakemake.output.units, sep="\t", index=False)


