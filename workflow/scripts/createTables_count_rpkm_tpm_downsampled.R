# This script is a little bit different from the other createTables because
# we need to have a matrix with all samples to perform the downsampling.
# In the other script fpkm and tpm are first calculated for every sample
# and later the different data.frames are merged. Here first we have to merge
# all the counts and then perform the downsampling and fpkm/tpm.


log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dplyr)
library(tibble)

source("workflow/scripts/custom_functions.R")

# Seed for downsampling
set.seed(snakemake@params[["seed"]])


#------------------------------------------------------------------------------------------
# Preapre count table
#------------------------------------------------------------------------------------------
sample_names <- unlist(snakemake@input) %>%
  basename %>%
  gsub(pattern = ".featureCounts", replacement = "")

# Merge all samples keeping geneid and length columns to calculate fpkm/tpm
counts <- snakemake@input %>%
  purrr::map(read.delim, header = TRUE, skip = 1) %>% # Ignore first line of featureCounts output
  purrr::map(select, 1,6,7) %>% # Select geneid, length and counts
  purrr::map(setNames, c("Geneid", "Length", "Counts"))  %>%
  plyr::join_all(type='inner', by = c("Geneid", "Length")) %>%
  setNames(c("Geneid", "Length", sample_names))


#------------------------------------------------------------------------------------------
# Remove excluded samples in case they are defined in config file prior to downsample
#------------------------------------------------------------------------------------------
if(!is.null(snakemake@params[["exclude"]])) {
    counts <- counts %>% select(-one_of(snakemake@params[["exclude"]]))
}


#------------------------------------------------------------------------------------------
# Downsample counts and calculate fpkm / tpm
#------------------------------------------------------------------------------------------
# Downsample count matrix to the sample with less reads
counts_downsampled <- cbind(
  select(counts, Geneid, Length),
  Down_Sample_Matrix( select(counts, -Geneid, -Length) )
)

fpkm_downsampled <- counts_downsampled %>%
  mutate_at( vars(sample_names), list(~ do_fpkm(., Length)) ) %>%
  select(-Length)


tpm_downsampled <- counts_downsampled %>%
  mutate_at( vars(sample_names), list(~ do_tpm(., Length)) ) %>%
  select(-Length)

counts_downsampled <- select(counts_downsampled, - Length)


#------------------------------------------------------------------------------------------
# Write output
#------------------------------------------------------------------------------------------
write.table(counts_downsampled, snakemake@output[["raw_counts"]], sep = "\t", quote = F, row.names = FALSE)
write.table(fpkm_downsampled, snakemake@output[["fpkm"]], sep = "\t", quote = F, row.names = FALSE)
write.table(tpm_downsampled, snakemake@output[["tpm"]], sep = "\t", quote = F, row.names = FALSE)