log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(DESeq2)
library(dplyr)

#------------------------- Preapre count table ------------------------------
read_files <- lapply(snakemake@input, function(x) read.delim(x, header = TRUE))

sample_names <- unlist(snakemake@input) %>%
				basename %>%
				gsub(pattern = ".counts", replacement = "")

counts <- plyr::join_all(read_files, type='inner', by = "Geneid") %>% 
  setNames(c("GeneSymbol", sample_names)) %>%
  tibble::column_to_rownames("GeneSymbol")


#------------------------- DESeq2 workflow ------------------------------
colData <- read.table("samples.tsv", header=TRUE)

# Counts and colData must have the same order
counts_ordered <- counts[as.character(colData$sample)]
stopifnot(identical(colnames(counts_ordered), as.character(colData$sample)))

dds <- DESeqDataSetFromMatrix(countData = counts_ordered,
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds)

write.table(counts(dds, normalized = T), snakemake@output[["norm_counts"]], sep = "\t", quote = F)
saveRDS(dds, file=snakemake@output[["rds"]])