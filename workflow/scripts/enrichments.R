log <- file(snakemake@log[[1]], open= "wt")
sink(log)
sink(log, type = "message")

library(tidyverse)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(openxlsx)
library(msigdbr)

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


# Set mouse or human references for the databases
if (genome == "mouse") { 
  kegg.genome <- "mmu"
  pa.genome   <- "mouse"
  db          <- "org.Mm.eg.db"
  m_df.h      <- msigdbr(species = "Mus musculus", category = "H") %>%
                   dplyr::select(gs_name, entrez_gene) %>%
                   as.data.frame()
  m_df.c2     <- msigdbr(species = "Mus musculus", category = "C2") %>%
                   dplyr::select(gs_name, entrez_gene) %>%
                   as.data.frame()
} else if (genome == "human") {
  kegg.genome <- "hsa"
  pa.genome   <- "human"
  db          <- "org.Hs.eg.db"
  m_df.h      <- msigdbr(species = "Homo sapiens", category = "H") %>%
                   dplyr::select(gs_name, entrez_gene) %>%
                   as.data.frame()
  m_df.c2     <- msigdbr(species = "Homo sapiens", category = "C2") %>%
                   dplyr::select(gs_name, entrez_gene) %>%
                   as.data.frame()
}

#------------------------------------------------------------------------------------------
# Prepare data
#------------------------------------------------------------------------------------------
# Get UP and DOWN-regulated 
UP <- DEA.annot %>% 
  dplyr::filter(DEG == "Upregulated") %>% 
  dplyr::select(Geneid) %>% 
  pull %>%
  bitr(fromType = "SYMBOL",toType = c("ENTREZID"), OrgDb = db) %>%
  pull(ENTREZID)

DWN <- DEA.annot %>% 
  dplyr::filter(DEG == "Downregulated") %>% 
  dplyr::select(Geneid) %>% 
  pull %>%
  bitr(fromType = "SYMBOL",toType = c("ENTREZID"), OrgDb = db) %>%
  pull(ENTREZID)

# Get universe of genes. Genes that have been considered for differential expression.
if (set_universe == TRUE) {
  universe <- DEA.annot %>%
    mutate(max_fpkm = purrr::reduce(dplyr::select(., contains("_FPKM")), pmax)) %>%
    dplyr::filter(max_fpkm > fpkm) %>%
    dplyr::select(Geneid) %>% 
    pull %>%
    bitr(fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = db) %>%
    pull(ENTREZID)
} else {
  universe <- DEA.annot %>%
    dplyr::select(Geneid) %>% 
    pull %>%
    bitr(fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = db) %>%
    pull(ENTREZID)
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
UP.msig_h  <- msig_db_enrichment(UP, db = db, pvalue = pvalue, qvalue = qvalue, universe = universe, term2gene = m_df.h)
DWN.msig_h <- msig_db_enrichment(DWN, db = db, pvalue = pvalue, qvalue = qvalue, universe = universe, term2gene = m_df.h)
UP.msig_c2 <- msig_db_enrichment(UP, db = db, pvalue = pvalue, qvalue = qvalue, universe = universe, term2gene = m_df.c2)
DWN.msig_c2<- msig_db_enrichment(DWN, db = db, pvalue = pvalue, qvalue = qvalue, universe = universe, term2gene = m_df.c2)
GSEA.hall  <- GSEA_enrichment(DEA.annot, "resources/h.all.v6.2.symbols.gmt")
GSEA.c2all <- GSEA_enrichment(DEA.annot, "resources/c2.all.v6.2.symbols.gmt")
GSEA.c3tft <- GSEA_enrichment(DEA.annot, "resources/c3.tft.v6.2.symbols.gmt")


#------------------------------------------------------------------------------------------
# save outputs to xls
#------------------------------------------------------------------------------------------
head(UP)
head(universe)
head(UP.go@result)
list_of_datasets <- list( "GO Upregulated"                 = UP.go@result %>% dplyr::filter(p.adjust < !!pvalue) %>% add_row(), 
                          "GO Downregulated"               = DWN.go@result %>% dplyr::filter(p.adjust < !!pvalue) %>% add_row(),
                          "KEGG Upregulated"               = UP.kegg@result %>% dplyr::filter(p.adjust < !!pvalue) %>% add_row(),
                          "KEGG Downregulated"             = DWN.kegg@result %>% dplyr::filter(p.adjust < !!pvalue) %>% add_row(),
                          "Reactome Upregulated"           = UP.pa@result %>% dplyr::filter(p.adjust < !!pvalue) %>% add_row(),
                          "Reactome Downregulaed"          = DWN.pa@result %>% dplyr::filter(p.adjust < !!pvalue) %>% add_row(),
                          "MSigDB_Hallmarks Upregulated"   = UP.msig_h@result %>% dplyr::filter(p.adjust < !!pvalue) %>% add_row(),
                          "MSigDB_Hallmarks Downregulated" = DWN.msig_h@result %>% dplyr::filter(p.adjust < !!pvalue) %>% add_row(),
                          "MSigDB_C2all Upregulated"       = UP.msig_c2@result %>% dplyr::filter(p.adjust < !!pvalue) %>% add_row(),
                          "MSigDB_C2all Downregulated"     = DWN.msig_c2@result %>% dplyr::filter(p.adjust < !!pvalue) %>% add_row(),                          
                          "GSEA Hallmarks"                 = GSEA.hall %>% dplyr::filter(padj < 0.25) %>% dplyr::select(-c(leadingEdge,nMoreExtreme)) %>% arrange(pval) %>% add_row(),
                          "GSEA c2all"                     = GSEA.c2all %>% dplyr::filter(padj < 0.25) %>% dplyr::select(-c(leadingEdge,nMoreExtreme)) %>% arrange(pval) %>% add_row(),
                          "GSEA c3tft"                     = GSEA.c3tft %>% dplyr::filter(padj < 0.25) %>% dplyr::select(-c(leadingEdge,nMoreExtreme)) %>% arrange(pval) %>% add_row()
                          )

write.xlsx(list_of_datasets, file = snakemake@output[["enrichments"]])
