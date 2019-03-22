rule deseq2:
    input:
        expand("02featureCounts/{sample}/{sample}.counts", sample = SAMPLES)
    output:
        "03deseq2/all.rds"
    params:
        samples=config["samples"]
    conda:
        "envs/deseq2.yaml"
    log:
        "00logs/deseq2/init.log"
    script:
        "scripts/DESeq2.R"

rule get_contrasts:
    input:
        "03deseq2/all.rds"
    output:
        table="results/diffexp/{contrast}.diffexp.tsv",
        ma_plot="results/diffexp/{contrast}.ma-plot.svg",
    conda:
        "envs/deseq2.yaml"
    params:
        contrast = lambda w: config["diffexp"]["contrasts"][w.contrast]
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    script:
        "scripts/get_DESeq2_contrasts.R"