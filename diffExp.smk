rule deseq2:
    input:
        expand("03featureCounts/{sample}/{sample}.counts", sample = SAMPLES)
    output:
        "04deseq2/all.rds"
    params:
        samples=config["samples"]
    log:
        "00log/deseq2/init.log"
    script:
        "scripts/DESeq2.R"

rule get_contrasts:
    input:
        rules.deseq2.output
    output:
        table="05results/diffexp/{contrast}.diffexp.tsv",
        ma_plot="05results/diffexp/{contrast}.ma-plot.pdf",
    params:
        contrast = lambda w: config["diffexp"]["contrasts"][w.contrast]
    log:
        "00log/deseq2/{contrast}.diffexp.log"
    script:
        "scripts/get_DESeq2_contrasts.R"

rule pca:
    input:
        rules.deseq2.output
    output:
        "05results/pca.pdf"
    params:
        pca_labels=config["pca"]["labels"]
    log:
        "00log/pca.log"
    script:
        "scripts/plot-PCA.R"
