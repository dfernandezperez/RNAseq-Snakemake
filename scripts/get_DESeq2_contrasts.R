log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("DESeq2")
library("IHW")
library("apeglm")
library("dplyr")
library("ggplot2")
library("tibble")
source("scripts/custom_functions.R")
dds <- readRDS(snakemake@input[[1]])

contrast <- c("condition", snakemake@params[["contrast"]])

# Apply IHW to weight p-values based on baseMean: https://bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html
res <- results(dds, contrast = contrast, filterFun = ihw, alpha = 0.05)

# shrink fold changes for lowly expressed genes, use the new ashr method: https://bioconductor.org/packages/release/bioc/vignettes/apeglm/inst/doc/apeglm.html#references
# Would be nice to use the apglm method also described in the paper, but can't be used with 'contrast', but with 'coef'	
# Finally I decided to keep the normal mode, we'll have to test more this
res <- lfcShrink(dds, contrast = contrast, res = res)

# sort by p-value
res <- res[order(res$padj),]

# Geneid to column, round to 2 deciamls all columns except pvalues
pval            <- res$pvalue
padjust         <- res$padj
res.filt        <- as.data.frame(res) %>% rownames_to_column(var = "Geneid") %>% round_df(2)
res.filt$pvalue <- pval
res.filt$padj   <- padjust

# Transform padj values that are 1 to NA (just to make nicer the volcano),
# this was the default behaviour of DESeq2, but IHW doesn't do it
res <- res %>% mutate(padj = ifelse(padj == 1, NA, padj))

# store results
pdf(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-4,4))
dev.off()

pdf(snakemake@output[["pval_hist"]])
qplot(res.filt$pvalue, xlab = "p-value", ylab = "count")
dev.off()

write.table( res.filt, file = snakemake@output[["table"]], sep = "\t", quote = F, row.names = F ) 