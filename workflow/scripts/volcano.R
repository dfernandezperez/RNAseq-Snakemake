log <- file(snakemake@log[[1]], open= "wt")
sink(log)
sink(log, type = "message")

library(ggplot2)

source("workflow/scripts/custom_functions.R")

# Avoid stupid MASS package to mask dplyr functions
select <- dplyr::select
filter <- dplyr::filter


df <- read.delim(snakemake@input[[1]])

log2fc <- as.numeric(snakemake@params[["log2fc"]])
pval   <- as.numeric(snakemake@params[["pval"]])

#------------------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------------------

p <- VolcanoPlot(df, 
                    xlim      = c(-8,8),
                    ylim      = c(0,18),
                    main      = gsub("-", " ", snakemake@params[["contrast"]]),
                    labelSize = 8,
                    pval      = pval,
                    log2FC    = log2fc)

ggsave(snakemake@output[["volcano_pdf"]], p, width = 6, height = 6)
ggsave(snakemake@output[["volcano_png"]], p, width = 6, height = 6)