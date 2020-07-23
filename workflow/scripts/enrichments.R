log <- file(snakemake@log[[1]], open= "wt")
sink(log)
sink(log, type = "message")

library(tidyverse)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
library(openxlsx)

source("workflow/scripts/custom_functions.R")

#------------------------------------------------------------------------------------------
# Read files and parameters
#------------------------------------------------------------------------------------------
DEA.annot <- read.delim(snakemake@input[[1]]) 

genome       <- snakemake@params[["genome"]]
pvalue       <- as.numeric(snakemake@params[["pvalue"]])
qvalue       <- as.numeric(snakemake@params[["qvalue"]])
fpkm         <- as.numeric(snakemake@wildcards[["fpkm"]])
set_universe <- as.logical(snakemake@params[["set_universe"]])
Geneid       <- as.character(snakemake@params[["annot_col"]])

#------------------------------------------------------------------------------------------
# Prepare data
#------------------------------------------------------------------------------------------
# If gene code is ensembl remove the decimals after the end of ensemblid (e.g: ENSMUSG00000000028.15)
if (Geneid == "ENSEMBL") {
  DEA.annot <- DEA.annot %>%
    mutate(!!Geneid = gsub(x = !!Geneid, pattern = "\\.\\d+$", replacement = ""))
}

# Get UP and DOWN-regulated 
UP <- DEA.annot %>% 
  dplyr::filter(DEG == "Upregulated") %>% 
  dplyr::select(!!Geneid) %>% 
  pull %>%
  bitr(fromType = Geneid, toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)

DWN <- DEA.annot %>% 
  dplyr::filter(DEG == "Downregulated") %>% 
  dplyr::select(!!Geneid) %>% 
  pull %>%
  bitr(fromType = Geneid, toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)

# Get universe of genes. Genes that have been considered for differential expression.
if (set_universe == TRUE) {
  universe <- DEA.annot %>%
    mutate(max_fpkm = purrr::reduce(dplyr::select(., contains("_FPKM")), pmax)) %>%
    dplyr::filter(max_fpkm > fpkm) %>%
    dplyr::select(!!Geneid) %>% 
    pull %>%
    bitr(fromType = Geneid, toType = c("ENTREZID"), OrgDb = org.Mm.eg.db) %>%
    pull(ENTREZID)
} else {
  universe <- DEA.annot %>%
    dplyr::select(!!Geneid) %>% 
    pull %>%
    bitr(fromType = Geneid, toType = c("ENTREZID"), OrgDb = org.Mm.eg.db) %>%
    pull(ENTREZID)
}

# Set mouse or human references for the databases
if (genome == "mouse") { 
  kegg.genome <- "mmu"
  pa.genome   <- "mouse"
  db          <- "org.Mm.eg.db"
} else if (genome == "human") {
  kegg.genome <- "hsa"
  pa.genome   <- "human"
  db          <- "org.Hs.eg.db"
}


#------------------------------------------------------------------------------------------
# Perform enrichments
#------------------------------------------------------------------------------------------
UP.go      <- goEnrichment(UP, ont = "BP", db = db, pvalue = pvalue, qvalue = qvalue, universe = universe)
DWN.go     <- goEnrichment(DWN, ont = "BP", db = db, pvalue = pvalue, qvalue = qvalue, universe = universe)
UP.kegg    <- KEGGenrichment(UP, org = kegg.genome, db = db, pvalue = pvalue, qvalue = qvalue, universe = universe)
DWN.kegg   <- KEGGenrichment(DWN, org = kegg.genome, db = db, pvalue = pvalue, qvalue = qvalue, universe = universe)
UP.pa      <- PAenrichment(UP, org = pa.genome, pvalue = pvalue, qvalue = qvalue, universe = universe)
DWN.pa     <- PAenrichment(DWN, org = pa.genome, pvalue = pvalue, qvalue = qvalue, universe = universe)
GSEA.hall  <- GSEA_enrichment(DEA.annot, "resources/h.all.v6.2.symbols.gmt")
GSEA.c2all <- GSEA_enrichment(DEA.annot, "resources/c2.all.v6.2.symbols.gmt")
GSEA.c3tft <- GSEA_enrichment(DEA.annot, "resources/c3.tft.v6.2.symbols.gmt")


#------------------------------------------------------------------------------------------
# save outputs to xls
#------------------------------------------------------------------------------------------
list_of_datasets <- list( "GO Upregulated"        = UP.go@result %>% dplyr::filter(p.adjust < !!pvalue) %>% add_row(), 
                          "GO Downregulated"      = DWN.go@result %>% dplyr::filter(p.adjust < !!pvalue) %>% add_row(),
                          "KEGG Upregulated"      = UP.kegg@result %>% dplyr::filter(p.adjust < !!pvalue) %>% add_row(),
                          "KEGG Downregulated"    = DWN.kegg@result %>% dplyr::filter(p.adjust < !!pvalue) %>% add_row(),
                          "Reactome Upregulated"  = UP.pa@result %>% dplyr::filter(p.adjust < !!pvalue) %>% add_row(),
                          "Reactome Downregulaed" = DWN.pa@result %>% dplyr::filter(p.adjust < !!pvalue) %>% add_row(),
                          "GSEA Hallmarks"        = GSEA.hall %>% dplyr::filter(padj < 0.25) %>% dplyr::select(-c(leadingEdge,nMoreExtreme)) %>% arrange(pval) %>% add_row(),
                          "GSEA c2all"            = GSEA.c2all %>% dplyr::filter(padj < 0.25) %>% dplyr::select(-c(leadingEdge,nMoreExtreme)) %>% arrange(pval) %>% add_row(),
                          "GSEA c3tft"            = GSEA.c3tft %>% dplyr::filter(padj < 0.25) %>% dplyr::select(-c(leadingEdge,nMoreExtreme)) %>% arrange(pval) %>% add_row())

write.xlsx(list_of_datasets, file = snakemake@output[["enrichments"]])