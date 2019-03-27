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