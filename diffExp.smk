rule deseq2:
    input:
        expand("03featureCounts/{sample}/{sample}.counts", sample = SAMPLES)
    output:
        rds         = "04deseq2/all.rds",
        norm_counts = "04deseq2/Normalized_counts.tsv"
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
        table   = "04deseq2/{contrast}/{contrast}_diffexp.tsv",
        ma_plot = "04deseq2/{contrast}/{contrast}_ma-plot.pdf",
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
        "04deseq2/pca.pdf"
    params:
        pca_labels = config["pca"]["labels"]
    log:
        "00log/pca.log"
    script:
        "scripts/plot-PCA.R"

# rule filter_deg:
#     input:
#         rules.get_contrasts.output.table
#     output:
#         "04deseq2/{contrast}/{log2fc}_{pvalue}/{contrast}_diffexp_{log2fc}_{pvalue}.tsv"
#     params:
#         pval   = config["diffexp"]["pvalue"],
#         log2fc = config["diffexp"]["log2fc"]
#     log:
#         "00log/deseq2/{contrast}.{log2fc}.{pvalue}.filter_deg.log"