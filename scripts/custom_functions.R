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
  df <- data.frame(df)
  df$DEG <- "NotDE"
  df$DEG[which(df$log2FoldChange >= log2FC & df$padj <= padjust)] <- "Upregulated" #Annotate UPregulated genes
  df$DEG[which(df$log2FoldChange <= -log2FC & df$padj <= padjust)] <- "Downregulated" #Annotate DOWNregulated genes
  
  df <- round_df(df,3) %>%
  rownames_to_column(var = "Geneid")
  return(df)
}


##-------------------- Enrichment analyses functions --------------------##
goEnrichment <- function(df, ont = "BP") {
  require(clusterProfiler)
  require(org.Mm.eg.db)
  ego <- enrichGO(gene         = df$ENTREZID,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = 'ENTREZID',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.1,
                  readable = TRUE)
  return(ego)
}

KEGGenrichment <- function(df, org = "mmu") {
  require(clusterProfiler)
  require(org.Mm.eg.db)
  ekgg <- enrichKEGG(gene         = df$ENTREZID,
                     organism         = "mmu",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.1)
  ekgg <- DOSE::setReadable(ekgg, OrgDb = org.Mm.eg.db, keytype="ENTREZID")
  return(ekgg)
}


PAenrichment <- function(df, org = "mouse") {
  require(clusterProfiler)
  require(ReactomePA)
  ePA <- enrichPathway(gene = df$ENTREZID,
                       pvalueCutoff=0.05,
                       qvalueCutoff  = 0.1,
                       organism = "mouse",
                       readable=TRUE)
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
                           stats = geneList,
                           minSize=15,
                           maxSize=2000,
                           nperm=10000)
  return(fgseaRes)
}