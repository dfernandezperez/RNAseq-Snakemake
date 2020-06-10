##-------------------- Filtering and annotation of DE tables --------------------##
#Function to round all numeric colums of dataframe
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}


# Annot DE table. Useful to filter the excel file or to do a volcano plot coloring the DE and non DE genes
Annot_DE <- function(df, log2FC = 2, padjust = 0.05, fpkm = 0) {
  require(dplyr)
  require(purrr)
  
  df <- df %>%
    
    mutate(max_fpkm = reduce(select(., contains("_FPKM")), pmax)) %>% # Add column with max fpkm value between conditions
    
    mutate(DEG = ifelse(log2FoldChange >= log2FC & padj <= padjust & max_fpkm >= fpkm, "Upregulated",
                        ifelse(log2FoldChange <= -log2FC & padj <= padjust & max_fpkm >= fpkm, "Downregulated", "NS"))) %>%
    
    select(-max_fpkm) # Remove temporary max fpkm column
  
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
  ekgg <- DOSE::setReadable(ekgg, OrgDb = db, keyType="ENTREZID")
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
  require(ggrastr)

  df <- mutate(df, shape = "circle")

  df <- mutate(df, shape = ifelse(-log10(padj) >  ylim[2], "triangle", shape))
  df <- mutate(df, padj = ifelse(-log10(padj) >  ylim[2], 10^-ylim[2], padj))

  df <- mutate(df, shape = ifelse(log2FoldChange > xlim[2], "triangle", shape))
  df <- mutate(df, shape = ifelse(log2FoldChange < -xlim[2], "triangle", shape))

  df <- mutate(df, log2FoldChange = ifelse(log2FoldChange > xlim[2], xlim[2], log2FoldChange))
  df <- mutate(df, log2FoldChange = ifelse(log2FoldChange < -xlim[2], -xlim[2], log2FoldChange))


  p <-  ggplot(data = na.omit(df), aes(x=log2FoldChange, y=-log10(padj), colour=DEG, shape=shape) ) +

    # geom_point_rast(alpha=0.7, size=1.7, raster.height = 5.15, raster.width = 6, raster.dpi = 400) +
    geom_point(alpha=0.7, size=1.7) +

    annotate("text", label = sum(df$DEG == "Upregulated"), color = "red", y = 0, x = xlim[2],
             vjust="inward",hjust="inward", size = labelSize) +
    annotate("text", label = sum(df$DEG == "Downregulated"), color = "darkgreen", y = 0, x = xlim[1],
             vjust="inward",hjust="inward", size = labelSize) +

    theme_classic(base_size = 20) +
    theme(legend.title = element_blank()) +
    theme(legend.position = "top") +


    ggtitle(main) +
    theme(plot.title = element_text(lineheight=.8, face="bold", hjust = .5)) +

    xlim(xlim) + ylim(ylim) +

    geom_hline(yintercept = -log10(pval), linetype = 2) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype = 2) +

    xlab(bquote(~Log[2]~ "fold change")) + ylab(bquote(~-Log[10]~italic(P))) +

    scale_colour_manual(values=c("Downregulated" = "darkgreen", "NS" = "gray", "Upregulated" = "red"),
                        labels = c("Downregulated" = "Downregulated", "NotDE" = "NS", "Upregulated" = "Upregulated"),
                        drop = FALSE) + #Force legend to show always

    guides(shape=FALSE) # Remove legend for shapes

  return(p)
}



# VolcanoPlot <- function(df, xlim=NULL, ylim=NULL, main = NULL, labelSize = 7, pval = 0.1, log2FC = 1) {
#   require(ggplot2)
#   require(ggrastr)

#   df <- mutate(df, shape = "circle")
#   df <- mutate(df, shape = ifelse(-log10(padj) >  ylim[2], "triangle_up", shape))
#   df <- mutate(df, padj = ifelse(-log10(padj) >  ylim[2], 10^-ylim[2], padj))

#   df <- mutate(df, shape = ifelse(log2FoldChange > xlim[2], "triangle_right", shape))
#   df <- mutate(df, shape = ifelse(log2FoldChange < -xlim[2], "triangle_left", shape))

#   df <- mutate(df, log2FoldChange = ifelse(log2FoldChange > xlim[2], xlim[2], log2FoldChange))
#   df <- mutate(df, log2FoldChange = ifelse(log2FoldChange < -xlim[2], -xlim[2], log2FoldChange))


#   p <-  ggplot(data = na.omit(df), aes(x=log2FoldChange, y=-log10(padj), colour=DEG, shape=shape) ) +

#     # geom_point_rast(alpha=0.7, size=1.7, raster.height = 5.15, raster.width = 6, raster.dpi = 400) +
#     geom_point(alpha=0.7, size=1.7) +
    
#     annotate("text", label = sum(df$DEG == "Upregulated"), color = "red", y = 0, x = xlim[2],
#              vjust="inward",hjust="inward", size = labelSize) +
#     annotate("text", label = sum(df$DEG == "Downregulated"), color = "darkgreen", y = 0, x = xlim[1],
#              vjust="inward",hjust="inward", size = labelSize) +

#     theme_classic(base_size = 20) +
#     theme(legend.title = element_blank()) +
#     theme(legend.position = "top") +


#     ggtitle(main) +
#     theme(plot.title = element_text(lineheight=.8, face="bold", hjust = .5)) +

#     xlim(xlim) + ylim(ylim) +

#     geom_hline(yintercept = -log10(pval), linetype = 2) +
#     geom_vline(xintercept = c(-log2FC, log2FC), linetype = 2) +

#     xlab(bquote(~Log[2]~ "fold change")) + ylab(bquote(~-Log[10]~italic(P))) +

#     scale_colour_manual(values=c("Downregulated" = "darkgreen", "NotDE" = "gray", "Upregulated" = "red"),
#                         labels = c("Downregulated" = "Downregulated", "NotDE" = "NS", "Upregulated" = "Upregulated"),
#                         drop = FALSE) + #Force legend to show always

#     scale_shape_manual(values=c("triangle_up" = "\u25B2", "triangle_right" = "\u25BA", "triangle_left" = "\u25C4", "circle" = "\u25CF",
#                                 drop = FALSE)) +

#     guides(shape=FALSE) # Remove legend for shapes

#   return(p)
# }