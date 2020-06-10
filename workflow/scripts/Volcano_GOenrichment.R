log <- file(snakemake@log[[1]], open= "wt")
sink(log)
sink(log, type = "message")


library(dplyr)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
library(openxlsx)
source("workflow/scripts/custom_functions.R")

DEA.annot <- read.delim(snakemake@input[[1]])

log2fc <- as.numeric(snakemake@params[["log2fc"]])
pval   <- as.numeric(snakemake@params[["pval"]])

##############################################
##                 Volcano                  ##
##############################################

p <- VolcanoPlot(DEA.annot, 
                    xlim      = c(-8,8),
                    ylim      = c(0,18),
                    main      = gsub("-", " ", snakemake@params[["contrast"]]),
                    labelSize = 8,
                    pval      = pval,
                    log2FC    = log2fc)

pdf(snakemake@output[["volcano"]], width = 6, height = 6)
p
dev.off()



##############################################
##        Gene enrichment analysis          ##
##############################################

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


# Set mouse or human references for the databases
if (snakemake@params[["genome"]] == "mouse") { 
  kegg.genome <- "mmu"
  pa.genome   <- "mouse"
  db          <- "org.Mm.eg.db"
} else if (snakemake@params[["genome"]] == "human") {
    kegg.genome <- "hsa"
    pa.genome   <- "human"
    db          <- "org.Hs.eg.db"
}


UP.go      <- goEnrichment(UP, ont = "BP", db = db)
UP.kegg    <- KEGGenrichment(UP, org = kegg.genome, db = db)
UP.pa      <- PAenrichment(UP, org = pa.genome)
DWN.go     <- goEnrichment(DWN, ont = "BP", db = db)
DWN.kegg   <- KEGGenrichment(DWN, org = kegg.genome, db = db)
DWN.pa     <- PAenrichment(DWN, org = pa.genome)
GSEA.hall  <- GSEA_enrichment(DEA.annot, "resources/h.all.v6.2.symbols.gmt")
GSEA.c2all <- GSEA_enrichment(DEA.annot, "resources/c2.all.v6.2.symbols.gmt")
GSEA.c3tft <- GSEA_enrichment(DEA.annot, "resources/c3.tft.v6.2.symbols.gmt")


# ------ save outputs to xlsx ---------
list_of_datasets <- list( "GO Upregulated"        = UP.go@result %>% filter(p.adjust < 0.05) %>% add_row(), 
                          "GO Downregulated"      = DWN.go@result %>% filter(p.adjust < 0.05) %>% add_row(),
                          "KEGG Upregulated"      = UP.kegg@result %>% filter(p.adjust < 0.05) %>% add_row(),
                          "KEGG Downregulated"    = DWN.kegg@result %>% filter(p.adjust < 0.05) %>% add_row(),
                          "Reactome Upregulated"  = UP.pa@result %>% filter(p.adjust < 0.05) %>% add_row(),
                          "Reactome Downregulaed" = DWN.pa@result %>% filter(p.adjust < 0.05) %>% add_row(),
                          "GSEA Hallmarks"        = GSEA.hall %>% filter(padj < 0.25) %>% dplyr::select(-c(leadingEdge,nMoreExtreme)) %>% arrange(pval) %>% add_row(),
                          "GSEA c2all"            = GSEA.c2all %>% filter(padj < 0.25) %>% dplyr::select(-c(leadingEdge,nMoreExtreme)) %>% arrange(pval) %>% add_row(),
                          "GSEA c3tft"            = GSEA.c3tft %>% filter(padj < 0.25) %>% dplyr::select(-c(leadingEdge,nMoreExtreme)) %>% arrange(pval) %>% add_row())

write.xlsx(list_of_datasets, file = snakemake@output[["enrichments"]])