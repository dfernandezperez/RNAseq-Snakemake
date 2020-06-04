log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

source("scripts/custom_functions.R")

contrast <- read.delim(snakemake@input[[1]])

contrast.filt <- Annot_DE(contrast, log2FC = snakemake@params[["log2fc"]], padjust = snakemake@params[["pval"]])

write.table( contrast.filt, file=snakemake@output[[1]], sep = "\t", quote = F, row.names = F )