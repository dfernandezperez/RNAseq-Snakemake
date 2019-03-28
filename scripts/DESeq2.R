log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


library(DESeq2)
library(dplyr)
library(tibble)

#------------------------- Preapre count table ------------------------------
read_files <- lapply(snakemake@input, function(x) read.delim(x, header = TRUE))

sample_names <- unlist(snakemake@input) %>%
				basename %>%
				gsub(pattern = ".counts", replacement = "")

counts <- plyr::join_all(read_files, type='inner', by = "Geneid") %>% 
  setNames(c("GeneSymbol", sample_names)) %>%
  column_to_rownames("GeneSymbol")


#------------------------- DESeq2 workflow ------------------------------
colData <- read.table(snakemake@params[["samples"]], header=TRUE)

# Counts and colData must have the same order
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

write.table(norm_counts, snakemake@output[["norm_counts"]], sep = "\t", quote = F, row.names = FALSE)

saveRDS(dds, file=snakemake@output[["rds"]])