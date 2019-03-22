import pandas as pd
from snakemake.utils import validate, min_version
shell.prefix('source activate DPbase; source activate /hpcnfs/data/DP/SnakemakePipelines/RNAseq-snakemake/.snakemake/conda/cea34b3f; ')

##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"
# validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
# validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "lane"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
# validate(units, schema="schemas/units.schema.yaml")

SAMPLES = set(units["sample"])

rule all: 
	input:
		expand("04results/diffexp/{contrast}.diffexp.tsv", contrast = config["diffexp"]["contrasts"]),
		"04results/pca.pdf"

##### load rules #####

include: "trim.smk"
include: "align.smk"
include: "diffExp.smk"