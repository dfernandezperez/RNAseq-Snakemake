log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("DESeq2")
library("dplyr")
library("ggplot2")
library("tibble")
source("scripts/custom_functions.R")
dds <- readRDS(snakemake@input[[1]])

contrast <- c("condition", snakemake@params[["contrast"]])
res <- results(dds, contrast = contrast)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast = contrast, res = res)
# sort by p-value
res <- res[order(res$padj),]
# TODO explore IHW usage

# Geneid to column, round to 3 deciamls all columns except pvalues
pval            <- res$pvalue
padjust         <- res$padj
res.filt        <- as.data.frame(res) %>% rownames_to_column(var = "Geneid") %>% round_df(2)
res.filt$pvalue <- pval
res.filt$padj   <- padjust

# store results
pdf(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-4,4))
dev.off()

pdf(snakemake@output[["pval_hist"]])
qplot(res.filt$pvalue, xlab = "p-value", ylab = "count")
dev.off()

write.table( res.filt, file = snakemake@output[["table"]], sep = "\t", quote = F, row.names = F ) 