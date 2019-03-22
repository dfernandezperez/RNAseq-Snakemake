log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(DESeq2)
library(dplyr)
source("/hpcnfs/scratch/DP/dfernand/Scripts/DESeq2_customFunctions.R")


#------------------------- Preapre count table ------------------------------
read_files <- lapply(snakemake@input[[1]], function(x) read.delim(x, header = TRUE, skip = 1))
counts <- plyr::join_all(read_files, type='inner', by = "Geneid") %>% 
  setNames(c("GeneSymbol", sample_names)) %>%
  tibble::column_to_rownames("GeneSymbol")


#------------------------- DESeq2 workflow ------------------------------
colData <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample", check.names=FALSE)

# Counts and colData must have the same order
counts_ordered <- counts[colData$sample]
stopifnot(identical(colnames(counts_ordered), colData$sample))

dds <- DESeqDataSetFromMatrix(countData = counts_ordered,
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds)

saveRDS(dds, file=snakemake@output[[1]])