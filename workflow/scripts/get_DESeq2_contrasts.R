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

# Apply IHW to weight p-values based on baseMean: https://bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html
res <- results(dds, name = contrast, filterFun = ihw, alpha = 0.05)


#------------------------------------------------------------------------------------------
# Get DEG genes, lfcshrink and volcano
#------------------------------------------------------------------------------------------
# shrink fold changes for lowly expressed genes, use the new apglm method: https://bioconductor.org/packages/release/bioc/vignettes/apeglm/inst/doc/apeglm.html#references
# https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty895/5159452
res <- lfcShrink(dds = dds, coef = contrast, res = res, type = "apeglm")

# Put nice the DEG table: sort by p-value, Geneid to column, round to 2 deciamls all columns except pvalues
res             <- res[order(res$padj),]
pval            <- res$pvalue
padjust         <- res$padj
res.filt        <- as.data.frame(res) %>% rownames_to_column(var = "Geneid") %>% round_df(2)
res.filt$pvalue <- pval
res.filt$padj   <- padjust

# Transform padj values that are 1 to NA (just to make nicer the volcano),
# this was the default behaviour of DESeq2, but IHW doesn't do it
res.filt <- res.filt %>% mutate(padj = ifelse(padj == 1, NA, padj))


#-----------------------------------------------------------------------------
# Read fpkm table and calculate the average of replicates
#-----------------------------------------------------------------------------
colData <- read.table(snakemake@params[["samples"]], header=TRUE)

# Remove unwanted samples (outliers, for example)
if(!is.null(snakemake@params[["exclude"]])) {
    colData <- colData %>% filter( !sample %in% snakemake@params[["exclude"]] )
}


fpkm         <- read.delim(snakemake@input[["fpkm"]], header = TRUE, check.names = FALSE)
colData_filt <- colData %>% 
  dplyr::filter(condition %in% snakemake@params[["contrast"]] ) 


fpkm_table <- fpkm %>%
  tidyr::pivot_longer(cols =  -"Geneid") %>%
  dplyr::filter(name %in% colData_filt$sample) %>%
  left_join(colData_filt, by = c("name" = "sample")) %>%
  group_by(condition, Geneid) %>%
  summarise(FPKM = mean(value)) %>%
  ungroup %>%
  mutate(condition = paste(condition, "_FPKM", sep = "")) %>%
  pivot_wider(names_from = condition,
              values_from = FPKM)

#-----------------------------------------------------------------------------
# Merge deseq2 results with FPKM data
#-----------------------------------------------------------------------------
res.filt <- res.filt %>%
    left_join(fpkm_table)


#-----------------------------------------------------------------------------
# Save data and plot log2fc with plotMA from deseq2. Save pvalue distribution
#-----------------------------------------------------------------------------
# store results
pdf(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-4,4))
dev.off()

pdf(snakemake@output[["pval_hist"]])
qplot(res.filt$pvalue, xlab = "p-value", ylab = "count")
dev.off()

write.table( res.filt, file = snakemake@output[["table"]], sep = "\t", quote = F, row.names = F ) 