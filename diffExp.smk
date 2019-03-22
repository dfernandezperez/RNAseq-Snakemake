rule deseq2:
    input:
        expand("02featureCounts/{sample}/{sample}.counts", sample = SAMPLES)
    output:
        "03deseq2/all.rds"
    params:
        samples=config["samples"]
    log:
        "00log/deseq2/init.log"
    script:
        "scripts/DESeq2.R"

rule get_contrasts:
    input:
        "03deseq2/all.rds"
    output:
        table="04results/diffexp/{contrast}.diffexp.tsv",
        ma_plot="04results/diffexp/{contrast}.ma-plot.pdf",
    params:
        contrast = lambda w: config["diffexp"]["contrasts"][w.contrast]
    log:
        "00log/deseq2/{contrast}.diffexp.log"
    script:
        "scripts/get_DESeq2_contrasts.R"

rule pca:
    input:
        "03deseq2/all.rds"
    output:
        "04results/pca.pdf"
    params:
        pca_labels=config["pca"]["labels"]
    log:
        "00log/pca.log"
    script:
        "scripts/plot-PCA.R"
