log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dplyr)
library(tibble)

#------------------------- Functions to tpm and rpkm ------------------------------
do_fpkm = function (counts, effective_lengths) {
  exp(log(counts) - log(effective_lengths) - log(sum(counts)) + log(1E9))
}

do_tpm = function (counts, effective_lengths) {
  rate = log(counts) - log(effective_lengths)
  exp(rate - log(sum(exp(rate))) + log(1E6))
}


#------------------------- Preapre count table ------------------------------
read_files <- snakemake@input %>%
  purrr::map(read.delim, header = TRUE, skip = 1) %>%
  purrr::map(select, 1,6,7) %>% # Select geneid, length and counts
  purrr::map(setNames, c("Geneid", "Length", "Counts")) %>%
  purrr::map(~ mutate(.x, fpkm = do_fpkm(Counts,Length))) %>%
  purrr::map(~ mutate(.x, tpm  = do_tpm(Counts,Length)))
  
sample_names <- unlist(snakemake@input) %>%
  basename %>%
  gsub(pattern = ".featureCounts", replacement = "")

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
  purrr::map(select, Geneid, tpm) %>%
  plyr::join_all(type='inner', by = "Geneid") %>% 
  setNames(c("Geneid", sample_names)) 


#------------------------- Write output ------------------------------
write.table(counts, snakemake@output[["raw_counts"]], sep = "\t", quote = F, row.names = FALSE)
write.table(fpkm, snakemake@output[["fpkm"]], sep = "\t", quote = F, row.names = FALSE)
write.table(tpm, snakemake@output[["tpm"]], sep = "\t", quote = F, row.names = FALSE)