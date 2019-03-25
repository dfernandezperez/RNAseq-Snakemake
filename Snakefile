import pandas as pd
from snakemake.utils import validate, min_version

singularity: "/hpcnfs/data/DP/Singularity/rnaseq-snakemake.simg"

##### set minimum snakemake version #####
min_version("5.4.3")


##### load config, cluster config and sample sheets #####

configfile: "config.yaml"
# validate(config, schema="schemas/config.schema.yaml")
CLUSTER = json.load(open(config['cluster_json']))

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
# validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "lane"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
# validate(units, schema="schemas/units.schema.yaml")

SAMPLES = set(units["sample"])

rule all: 
	input:
		expand("05results/diffexp/{contrast}.diffexp.tsv", contrast = config["diffexp"]["contrasts"]),
		"05results/pca.pdf",
		"01qc/multiqc_report.html"

##### load rules #####

include: "common.smk"
include: "trim.smk"
include: "align.smk"
include: "diffExp.smk"
include: "qc.smk"