rule create_tables_downsampled:
    input:
        expand("results/03featureCounts/{sample}/{sample}.featureCounts", sample = SAMPLES)
    output:
        tpm         = "results/04deseq2/downsampled/tpm.tsv",
        fpkm        = "results/04deseq2/downsampled/fpkm.tsv",
        raw_counts  = "results/04deseq2/downsampled/Raw_counts.tsv"
    params:
        exclude = config["diffexp"].get("exclude", None),
        seed    = config["seed"]
    log:
        "results/00log/deseq2/downsampled/create_tables.log"
    script:
        "../scripts/createTables_count_rpkm_tpm_downsampled.R"
 
rule deseq2_downsampled:
    input:
        rules.create_tables_downsampled.output.raw_counts
    output:
        rds         = "results/04deseq2/downsampled/all.rds",
        norm_counts = "results/04deseq2/downsampled/Normalized_counts.tsv",
    params:
        samples = config["samples"],
        exclude = config["diffexp"].get("exclude", None)
    log:
        "results/00log/deseq2/downsampled/init.log"
    script:
        "../scripts/DESeq2.R"

rule get_contrasts_downsampled:
    input:
        rds     = rules.deseq2_downsampled.output.rds,
        fpkm    = rules.create_tables_downsampled.output.fpkm
    output:
        table     = "results/04deseq2/downsampled/{contrast}/{contrast}_diffexp.tsv",
        ma_plot   = "results/04deseq2/downsampled/{contrast}/{contrast}_ma-plot.pdf",
        pval_hist = "results/04deseq2/downsampled/{contrast}/{contrast}_pval-hist.pdf",
    params:
        contrast        = lambda w: config["diffexp"]["contrasts"][w.contrast],
        lfcShrink       = config["lfcShrink"],
        samples         = config["samples"],
        exclude         = config["diffexp"].get("exclude", None),
        annot           = config["ref"]["geneInfo"].get("file", None),
        column_used     = config["ref"]["geneInfo"]["column_used"],
        column_toAdd    = config["ref"]["geneInfo"]["column_toAdd"],
        name_annotation = config["ref"]["geneInfo"]["name_annotation"],
    log:
        "results/00log/deseq2/downsampled/{contrast}.diffexp.log"
    script:
        "../scripts/get_DESeq2_contrasts.R"

rule pca_downsampled:
    input:
        rules.deseq2_downsampled.output.rds
    output:
        "results/04deseq2/downsampled/pca_elipse_names_top{ntop}.pdf",
        "results/04deseq2/downsampled/pca_top{ntop}.pdf",
        "results/04deseq2/downsampled/pca_names_top{ntop}.pdf"
    params:
        pca_labels = config["pca"]["labels"],
    log:
        "results/00log/downsampled/pca_{ntop}.log"
    script:
        "../scripts/plot-PCA.R"

rule filter_deg_downsampled:
    input:
        diffExp = rules.get_contrasts_downsampled.output.table,
    output:
        "results/04deseq2/downsampled/{contrast}/log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}/{contrast}_diffexp_log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.tsv"
    params:
        pval      = lambda w: w.pvalue,
        log2fc    = lambda w: w.log2fc,
        fpkm_filt = lambda w: w.fpkm
    log:
        "results/00log/deseq2/downsampled/{contrast}.{log2fc}.{pvalue}_fpkm{fpkm}.filter_deg.log"
    script:
        "../scripts/filter_deg.R"

rule enrichments_downsampled:
    input:
        rules.filter_deg_downsampled.output
    output:
        enrichments = "results/04deseq2/downsampled/{contrast}/log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}/{contrast}_enrichments_log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.xlsx",
    params:
        genome       = config["ref"]["genome"],
        pvalue       = config["enrichments"]["pval"],
        qvalue       = config["enrichments"]["qval"],
        set_universe = config["enrichments"]["set_universe"],
    log:
        "results/00log/deseq2/downsampled/{contrast}.log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.enrichments.log"
    script:
        "../scripts/enrichments.R"   


rule volcano_downsampled:
    input:
        rules.filter_deg_downsampled.output
    output:
        volcano_pdf	 = "results/04deseq2/downsampled/{contrast}/log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}/{contrast}_volcano_log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.pdf",
        volcano_png	 = "results/04deseq2/downsampled/{contrast}/log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}/{contrast}_volcano_log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.png",
    params:
        pval     = lambda w: w.pvalue,
        log2fc   = lambda w: w.log2fc,
        contrast = lambda w: w.contrast,
    log:
        "results/00log/deseq2/downsampled/{contrast}.log2fc{log2fc}_pval{pvalue}_fpkm{fpkm}.volcano.log"
    script:
        "../scripts/volcano.R"  
