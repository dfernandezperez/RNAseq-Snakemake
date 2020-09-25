log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


library(DESeq2)
library(dplyr)
library(tibble)
library(GenomicFeatures)
library(tximport)

#------------------------------------------------------------------------------------------
# Read abundances from salmon and create tx2gene
#------------------------------------------------------------------------------------------
files        <- unlist(snakemake@input)

names(files) <- as.character(snakemake@params[["sample_names"]])

txdb    <- makeTxDbFromGFF(file=snakemake@params[["tx2gene"]])
k       <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")


#------------------------------------------------------------------------------------------
# DESeq2 workflow
#------------------------------------------------------------------------------------------
colData <- read.table(snakemake@params[["samples"]], header=TRUE)
# Remove unwanted samples (outliers, for example)
if(!is.null(snakemake@params[["exclude"]])) {
    colData <- colData %>% dplyr::filter( !sample %in% snakemake@params[["exclude"]] )
}

# Filter samples that must be excluded from the list of quant.sf files.
# Also force the same order of samples in colData and quant.sf
files <- files[as.character(colData$sample)]
stopifnot(identical(names(files), as.character(colData$sample)))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene) # Load salmon quant files
dds <- DESeqDataSetFromTximport(txi       = txi,
                                colData   = colData,
                                design    = ~ condition)

dds <- DESeq(dds)


Geneid <- as.character(snakemake@params[["annot_col"]])
norm_counts <- counts(dds, normalized = T) %>% 
                data.frame %>% 
                round(3) %>% 
                rownames_to_column(var = Geneid)


#------------------------------------------------------------------------------------------
# Save output
#------------------------------------------------------------------------------------------
write.table(norm_counts, snakemake@output[["norm_counts"]], sep = "\t", quote = F, row.names = FALSE)
saveRDS(dds, file=snakemake@output[["rds"]])