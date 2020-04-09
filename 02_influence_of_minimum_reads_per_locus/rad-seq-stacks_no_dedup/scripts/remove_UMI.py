"""
"""
import sys
import dinopy
from shutil import copyfile

# redirect stderr to logfile
sys.stderr = open(snakemake.log[0], "w")

p7_reader = dinopy.FastqReader(snakemake.input.fq2)
umi_len = snakemake.params.umi["len"]

# copy p5 file without touching it
copyfile(snakemake.input.fq1, snakemake.output.fq1)

# trim the first UMI-len bases from the p7 read
with dinopy.FastqWriter(snakemake.output.fq2, force_overwrite=True) as writer:
    for seq, name, qual in p7_reader.reads(read_names=True, quality_values=True):
        writer.write(seq[umi_len:], name, qual[umi_len:])
