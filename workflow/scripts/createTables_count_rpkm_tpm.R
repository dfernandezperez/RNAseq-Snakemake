log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dplyr)
library(tibble)

source("workflow/scripts/custom_functions.R")

#------------------------------------------------------------------------------------------
# Preapre count table
#------------------------------------------------------------------------------------------
read_files <- snakemake@input %>%
  purrr::map(read.delim, header = TRUE) %>%
  purrr::map(setNames, c("Geneid", "Length", "Effective_Length", "TPM", "Counts")) %>%
  purrr::map(~ mutate(.x, fpkm = do_tpm(Counts,Length)))
  
sample_names <- as.character(snakemake@params[["sample_names"]])

# Create a df with raw counts, fpkm and tpm
counts <- read_files %>%
  purrr::map(select, Geneid, Counts) %>%
  plyr::join_all(type='inner', by = "Geneid") %>% 
  setNames(c("Geneid", sample_names))

fpkm <- read_files %>%
  purrr::map(select, Geneid, fpkm) %>%
  plyr::join_all(type='inner', by = "Geneid") %>% 
  setNames(c("Geneid", sample_names))

tpm <- read_files %>%
  purrr::map(select, Geneid, TPM) %>%
  plyr::join_all(type='inner', by = "Geneid") %>% 
  setNames(c("Geneid", sample_names)) 


#------------------------------------------------------------------------------------------
# Remove excluded samples in case they are defined in config file
#------------------------------------------------------------------------------------------
if(!is.null(snakemake@params[["exclude"]])) {
    counts <- counts %>% select(-one_of(snakemake@params[["exclude"]]))
    fpkm   <- fpkm %>% select(-one_of(snakemake@params[["exclude"]]))
    tpm    <- tpm %>% select(-one_of(snakemake@params[["exclude"]]))
}


#------------------------------------------------------------------------------------------
# Write output
#------------------------------------------------------------------------------------------
write.table(counts, snakemake@output[["raw_counts"]], sep = "\t", quote = F, row.names = FALSE)
write.table(fpkm, snakemake@output[["fpkm"]], sep = "\t", quote = F, row.names = FALSE)
write.table(tpm, snakemake@output[["tpm"]], sep = "\t", quote = F, row.names = FALSE)