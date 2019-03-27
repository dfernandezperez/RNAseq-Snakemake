setwd("/hpcnfs/scratch/DP/dfernand/Spivetti/PCGFS_PROJECT/RNAseq/counts")
library(dplyr)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
require(openxlsx)

DEA.annot <- read.delim(snakemake@input[[1]])

# Get UP and DOWN-regulated 
UP <- DEA.annot %>% 
  dplyr::filter(DEG == "Upregulated") %>% 
  dplyr::select(Geneid) %>% 
  pull %>%
  bitr(fromType = "SYMBOL",toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)

DWN <- DEA.annot %>% 
  dplyr::filter(DEG == "Downregulated") %>% 
  dplyr::select(Geneid) %>% 
  pull %>%
  bitr(fromType = "SYMBOL",toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)

UP.go <- goEnrichment(UP, ont = "BP")
UP.kegg <- KEGGenrichment(UP, org = "mmu")
UP.pa <- PAenrichment(UP, org = "mouse")
DWN.go <- goEnrichment(DWN, ont = "BP")
DWN.kegg <- KEGGenrichment(DWN, org = "mmu")
DWN.pa <- PAenrichment(DWN, org = "mouse")
GSEA <- GSEA_enrichment(DEA.annot, gmt)

# ------ save outputs to xlsx ---------
list_of_datasets <- list( "GO Upregulated"        = UP.go@result %>% filter(p.adjust < 0.05) %>% add_row(), 
                          "GO Downregulated"      = DWN.go@result %>% filter(p.adjust < 0.05) %>% add_row(),
                          "KEGG Upregulated"      = UP.kegg@result %>% filter(p.adjust < 0.05) %>% add_row(),
                          "KEGG Downregulated"    = DWN.kegg@result %>% filter(p.adjust < 0.05) %>% add_row(),
                          "Reactome Upregulated"  = UP.pa@result %>% filter(p.adjust < 0.05) %>% add_row(),
                          "Reactome Downregulaed" = DWN.pa@result %>% filter(p.adjust < 0.05) %>% add_row(),
                          "GSEA Hallmarks"        = GSEA.hall %>% filter(padj < 0.25) %>% dplyr::select(-c(leadingEdge,nMoreExtreme)) %>% arrange(pval) %>% add_row(),
                          "GSEA c2all"            = GSEA.g2all %>% filter(padj < 0.25) %>% dplyr::select(-c(leadingEdge,nMoreExtreme)) %>% arrange(pval) %>% add_row())

write.xlsx(list_of_datasets, file = snakemake@output[[1]])