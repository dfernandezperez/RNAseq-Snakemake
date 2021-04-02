log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("DESeq2")
library("pcaExplorer")
library("ggplot2")
library("ggforce")


#------------------------------------------------------------------------------------------
# Read parameters
#------------------------------------------------------------------------------------------
# Define number of genes to perform PCA
ntop <- as.numeric(snakemake@wildcards[["ntop"]])

# group variable set in config.yaml
grouping <- snakemake@params[["pca_labels"]]


#------------------------------------------------------------------------------------------
# Load data and define custom theme
#------------------------------------------------------------------------------------------
# Define theme for plot
theme_custom <-   theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"),
        legend.key.size  = unit(0.4, units = "cm"),
  ) +
  theme(plot.title = element_text(hjust = 0.5)
  ) +
  theme(
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  theme(
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(fill="white")
  )

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

# obtain normalized counts
counts <- rlog(dds, blind = FALSE)


#------------------------------------------------------------------------------------------
# Plotting
#------------------------------------------------------------------------------------------
p1 <- plotPCA(counts, intgroup = grouping, ntop = ntop) +
  theme_custom + theme(legend.position = "none") +
  geom_mark_ellipse(aes(fill = group, label = group), con.cap = 0, expand = unit(3, "mm"), label.fontsize = 14)

p2 <- plotPCA(counts, intgroup = grouping, ntop = ntop) +
  theme_custom

# PCA with sample labels. We use this other function from pcaExplorer for simplicity
p3 <- pcaplot(counts, intgroup = grouping, ellipse = FALSE, text_labels = TRUE, ntop = ntop)


#------------------------------------------------------------------------------------------
# Save to file
#------------------------------------------------------------------------------------------
ggsave(filename = snakemake@output[[1]], p1, width = 6, height = 6)
ggsave(filename = snakemake@output[[2]], p2, width = 6, height = 6)
ggsave(filename = snakemake@output[[3]], p3)
