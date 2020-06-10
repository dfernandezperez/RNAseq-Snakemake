rule create_tables:
    input:
        expand("results/03featureCounts/{sample}/{sample}.featureCounts", sample = SAMPLES)
    output:
        tpm         = "results/04deseq2/tpm.tsv",
        fpkm        = "results/04deseq2/fpkm.tsv",
        raw_counts  = "results/04deseq2/Raw_counts.tsv"
    params:
        exclude = config["diffexp"].get("exclude", None)
    log:
        "results/00log/deseq2/create_tables.log"
    script:
        "../scripts/createTables_count_rpkm_tpm.R"
 
rule deseq2:
    input:
        rules.create_tables.output.raw_counts
    output:
        rds         = "results/04deseq2/all.rds",
        norm_counts = "results/04deseq2/Normalized_counts.tsv",
    params:
        samples = config["samples"],
        exclude = config["diffexp"].get("exclude", None)
    log:
        "results/00log/deseq2/init.log"
    script:
        "../scripts/DESeq2.R"

rule get_contrasts:
    input:
        rds     = rules.deseq2.output.rds,
        fpkm    = rules.create_tables.output.fpkm
    output:
        table     = "results/04deseq2/{contrast}/{contrast}_diffexp.tsv",
        ma_plot   = "results/04deseq2/{contrast}/{contrast}_ma-plot.pdf",
        pval_hist = "results/04deseq2/{contrast}/{contrast}_pval-hist.pdf",
    params:
        contrast = lambda w: config["diffexp"]["contrasts"][w.contrast],
        samples   = config["samples"],
        exclude   = config["diffexp"].get("exclude", None)
    log:
        "results/00log/deseq2/{contrast}.diffexp.log"
    script:
        "../scripts/get_DESeq2_contrasts.R"

rule pca:
    input:
        rules.deseq2.output.rds
    output:
        "results/04deseq2/pca.pdf"
    params:
        pca_labels = config["pca"]["labels"]
    log:
        "results/00log/pca.log"
    script:
        "../scripts/plot-PCA.R"

rule filter_deg:
    input:
        diffExp = rules.get_contrasts.output.table,
    output:
        "results/04deseq2/{contrast}/log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}/{contrast}_diffexp_log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.tsv"
    params:
        pval      = lambda w: w.pvalue,
        log2fc    = lambda w: w.log2fc,
        fpkm_filt = lambda w: w.fpkm
    log:
        "results/00log/deseq2/{contrast}.{log2fc}.{pvalue}_fpkm{fpkm}.filter_deg.log"
    script:
        "../scripts/filter_deg.R"

rule deg_analysis:
    input:
        rules.filter_deg.output
    output:
        enrichments = "results/04deseq2/{contrast}/log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}/{contrast}_enrichments_log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.xls",
        volcano	    = "results/04deseq2/{contrast}/log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}/{contrast}_volcano_log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.pdf",
    params:
        pval     = lambda w: w.pvalue,
        log2fc   = lambda w: w.log2fc,
        contrast = lambda w: w.contrast,
        genome   = config["ref"]["genome"],
    log:
        "results/00log/deseq2/{contrast}.log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.deg_analysis.log"
    script:
        "../scripts/Volcano_GOenrichment.R"   
