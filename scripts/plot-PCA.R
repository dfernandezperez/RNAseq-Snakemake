log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("DESeq2")
library("pcaExplorer")
# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

# obtain normalized counts
counts <- rlog(dds, blind = FALSE)

pdf(snakemake@output[[1]], width = 7, height = 7.5)
pcaplot(counts, intgroup = snakemake@params[["pca_labels"]], ellipse = FALSE, text_labels = TRUE)
pcaplot(counts, intgroup = snakemake@params[["pca_labels"]], ellipse = FALSE, text_labels = FALSE)
dev.off()