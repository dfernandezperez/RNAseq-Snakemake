log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

source("workflow/scripts/custom_functions.R")
library(dplyr)


#-----------------------------------------------------------------------------
# Annot DE analysis and add fpkm information
#-----------------------------------------------------------------------------
contrast <- read.delim(snakemake@input[["diffExp"]])

log2fc    <- as.numeric(snakemake@params[["log2fc"]])
pval      <- as.numeric(snakemake@params[["pval"]])
fpkm_filt <- as.numeric(snakemake@params[["fpkm_filt"]])

contrast_annot <- contrast %>% 
  Annot_DE(log2FC = log2fc, padjust = pval, fpkm = fpkm_filt) %>%
  dplyr::select(-baseMean, -lfcSE, -pvalue)

write.table( contrast_annot, file=snakemake@output[[1]], sep = "\t", quote = F, row.names = F )