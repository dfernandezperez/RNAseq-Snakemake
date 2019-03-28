##-------------------- Filtering and annotation of DE tables --------------------##
#Function to round all numeric colums of dataframe
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}


# Annot DE table. Useful to filter the excel file or to do a volcano plot coloring the DE and non DE genes
Annot_DE <- function(df, log2FC = 2, padjust = 0.05) {
  require(dplyr)
  require(tibble)
  # df <- data.frame(df)
  df$DEG <- "NotDE"
  df$DEG[which(df$log2FoldChange >= log2FC & df$padj <= padjust)] <- "Upregulated" #Annotate UPregulated genes
  df$DEG[which(df$log2FoldChange <= -log2FC & df$padj <= padjust)] <- "Downregulated" #Annotate DOWNregulated genes
  
  # df <- round_df(df,3) %>%
  # rownames_to_column(var = "Geneid")
  return(df)
}


##-------------------- Enrichment analyses functions --------------------##
goEnrichment <- function(df, ont = "BP", db = org.Mm.eg.db) {
  require(clusterProfiler)
  require(db)
  ego <- enrichGO(gene          = df$ENTREZID,
                  OrgDb         = db,
                  keyType       = 'ENTREZID',
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.1,
                  readable = TRUE)
  return(ego)
}

KEGGenrichment <- function(df, org = "mmu", db = org.Mm.eg.db) {
  require(clusterProfiler)
  require(db)
  ekgg <- enrichKEGG(gene          = df$ENTREZID,
                     organism      = org,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.1)
  ekgg <- DOSE::setReadable(ekgg, OrgDb = db, keytype="ENTREZID")
  return(ekgg)
}


PAenrichment <- function(df, org = "mouse") {
  require(clusterProfiler)
  require(ReactomePA)
  ePA <- enrichPathway(gene         = df$ENTREZID,
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.1,
                       organism     = org,
                       readable     = TRUE)
  return(ePA)
}

GSEA_enrichment <- function(df, pathways.gmt) {
  require(clusterProfiler)
  pathways <- fgsea::gmtPathways(pathways.gmt)
  
  # Create a list containing a named vector (with genenames) of log2fc of each PcgfsKO vs WT
  geneList <- df$log2FoldChange
  names(geneList) <- toupper(df$Geneid)
  
  # Run GSEA algorithm
  fgseaRes <- fgsea::fgsea(pathways = pathways, 
                           stats   = geneList,
                           minSize = 15,
                           maxSize = 2000,
                           nperm   = 10000)
  return(fgseaRes)
}



##-------------------- Plots --------------------##
VolcanoPlot <- function(df, xlim=NULL, ylim=NULL, main = NULL, labelSize = 7, pval = 0.1, log2FC = 1) {
  require(ggplot2)
  
  p <-  ggplot(data = df, aes(x=log2FoldChange, y=-log10(padj), colour=DEG) ) +
    geom_point(alpha=0.7, size=2) +
    
    annotate("text", label = sum(df$DEG == "Upregulated"), color = "red", y = 0, x = xlim[2], 
             vjust="inward",hjust="inward", size = labelSize) +
    annotate("text", label = sum(df$DEG == "Downregulated"), color = "darkgreen", y = 0, x = xlim[1],
             vjust="inward",hjust="inward", size = labelSize) +
    
    theme_classic() +
    theme(legend.title = element_blank()) +
    theme(legend.position = "top") +
    
    
    ggtitle(main) +
    theme(plot.title = element_text(lineheight=.8, face="bold", hjust = .5)) +
    
    xlim(xlim) + ylim(ylim) + 
    
    geom_hline(yintercept = -log10(pval), linetype = 2) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype = 2) +
    
    xlab("log2 fold change") + ylab("-log10 p-value") +
    scale_colour_manual(values=c("darkgreen", "gray", "red") )
  
  p
}