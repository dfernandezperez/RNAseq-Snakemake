log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


library(DESeq2)
library(dplyr)
library(tibble)


#------------------------------------------------------------------------------------------
# Read count table
#------------------------------------------------------------------------------------------
counts <- read.delim(snakemake@input[[1]], check.names = FALSE) %>%
    column_to_rownames("Geneid")


#------------------------------------------------------------------------------------------
# DESeq2 workflow
#------------------------------------------------------------------------------------------
colData <- read.table(snakemake@params[["samples"]], header=TRUE)
# Remove unwanted samples (outliers, for example)
if(!is.null(snakemake@params[["exclude"]])) {
    colData <- colData %>% filter( !sample %in% snakemake@params[["exclude"]] )
}

# Counts and colData must have the same order and samples
counts_ordered <- counts[as.character(colData$sample)]
stopifnot(identical(colnames(counts_ordered), as.character(colData$sample)))


dds <- DESeqDataSetFromMatrix(countData = counts_ordered,
                              colData   = colData,
                              design    = ~ condition)
dds <- DESeq(dds)

norm_counts <- counts(dds, normalized = T) %>% 
                data.frame %>% 
                round(3) %>% 
                rownames_to_column(var = "Geneid")


#------------------------------------------------------------------------------------------
# Save output
#------------------------------------------------------------------------------------------
write.table(norm_counts, snakemake@output[["norm_counts"]], sep = "\t", quote = F, row.names = FALSE)
saveRDS(dds, file=snakemake@output[["rds"]])