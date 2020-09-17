log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("DESeq2")
library("IHW")
library("apeglm")
library("dplyr")
library("ggplot2")
library("tibble")
library("tidyr")
source("workflow/scripts/custom_functions.R")


#------------------------------------------------------------------------------------------
# Read data and set contrast
#------------------------------------------------------------------------------------------
dds      <- readRDS(snakemake@input[["rds"]])
contrast <- paste("condition_", snakemake@params[["contrast"]][1], "_vs_", snakemake@params[["contrast"]][2], sep = "")

# Set as reference for the GLM model the sample that is going to act as a control in the contrast
dds$condition <- relevel(dds$condition, snakemake@params[["contrast"]][2])
dds           <- nbinomWaldTest(dds)


#------------------------------------------------------------------------------------------
# Apply IHW to weight p-values based on baseMean
#------------------------------------------------------------------------------------------
# https://bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html
res <- results(dds, name = contrast, filterFun = ihw, alpha = 0.05)


#------------------------------------------------------------------------------------------
# Shrink fold changes with apeglm method of DESeq2 if it was set in config file.
#------------------------------------------------------------------------------------------
do_lfcShrink <- as.logical(snakemake@params[["lfcShrink"]])

if (do_lfcShrink == TRUE) {
  # https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty895/5159452
  res <- lfcShrink(dds = dds, coef = contrast, res = res, type = "apeglm")
}

#-----------------------------------------------------------------------------
# Read fpkm table and calculate the average of replicates
#-----------------------------------------------------------------------------
# Get samples from the conditions of the contrast
colData <- read.table(snakemake@params[["samples"]], header=TRUE) %>%
  filter( condition %in% snakemake@params[["contrast"]] ) 

# Remove unwanted samples (outliers, for example)
if(!is.null(snakemake@params[["exclude"]])) {
    colData <- colData %>% filter( !sample %in% snakemake@params[["exclude"]] )
}

fpkm <- read.delim(snakemake@input[["fpkm"]], header = TRUE, check.names = FALSE)

fpkm_table <- fpkm %>%
  pivot_longer(cols =  -"Geneid") %>%
  filter(name %in% colData$sample) %>%
  left_join(colData, by = c("name" = "sample")) %>%
  group_by(condition, Geneid) %>%
  summarise(FPKM = mean(value)) %>%
  ungroup %>%
  mutate(condition = paste(condition, "_FPKM", sep = "")) %>%
  pivot_wider(names_from = condition,
              values_from = FPKM)


#------------------------------------------------------------------------------------------
# Tidy output: Merge FPKM data, sort by p-value, Geneid to column, round to 2 deciamls
#------------------------------------------------------------------------------------------
res.tidy <- as.data.frame(res) %>%
  rownames_to_column(var = "Geneid") %>%
  left_join(fpkm_table) %>%
  arrange(padj) %>%
  mutate_at(vars(-Geneid, -pvalue, -padj), list(~ round(., 2))) %>%
# Transform padj values that are 1 to NA (just to make nicer the volcano),
# this was the default behaviour of DESeq2, but IHW doesn't do it
  mutate(padj = ifelse(padj == 1, NA, padj)) 


#-----------------------------------------------------------------------------
# Add other annotations (ENSEMBL, ENTREZ...)
#-----------------------------------------------------------------------------
if(!is.null(snakemake@params[["annot"]])) {
  file     <- as.character(snakemake@params[["annot"]])
  col_used <- as.numeric(snakemake@params[["column_used"]])
  col_add  <- as.numeric(snakemake@params[["column_toAdd"]])
  name_add <- as.character(snakemake@params[["name_annotation"]])

  res.tidy <- read.delim(file, skip = 1) %>%
                select(!!col_used, !!col_add) %>%
                setNames(c("Geneid", name_add)) %>%
                right_join(res.tidy) %>%
                # Remove any potential duplicates due to
                # the addition of the new annotation
                distinct(Geneid, .keep_all = TRUE) %>%
                arrange(padj) # Sort by pvalue because this changes this order
                
}

#-----------------------------------------------------------------------------
# Save data and plot log2fc with plotMA from deseq2. Save pvalue distribution
#-----------------------------------------------------------------------------
# store results
pdf(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-4,4))
dev.off()

pdf(snakemake@output[["pval_hist"]])
qplot(res.tidy$pvalue, xlab = "p-value", ylab = "count")
dev.off()

write.table( res.tidy, file = snakemake@output[["table"]], sep = "\t", quote = F, row.names = F ) 