rm(list=ls())
library(CREAM)
###########
GeneMatch <- function(Regions, RefTSS, TSS_distance = 1e3){
  
  chr_vec_ref <- as.character(RefTSS$chr)
  tss_ref <- as.numeric(RefTSS$tss)
  gene_names <- as.character(RefTSS$gene)
  #####
  chr_vec <- as.character(Regions[,1])
  start_vec <- as.numeric(Regions[,2])
  end_vec <- as.numeric(Regions[,3])
  
  GeneMatch <- c()
  for(reg_iter in 1:nrow(Regions)){
    matched_ind <- which(chr_vec[reg_iter] == chr_vec_ref & 
                           tss_ref-TSS_distance < end_vec[reg_iter] & 
                           tss_ref+TSS_distance > start_vec[reg_iter])
    GeneMatch <- c(GeneMatch, gene_names[matched_ind])
  }
  ##
  GeneMatch <- unique(GeneMatch)
  return(GeneMatch)
}
#######
TSS_Dist <- 100

RefTSS <- read.table("~/OneDrive - University of Toronto/CREAM/refGene_hg38_TSS.txt")
input_path <- "~/OneDrive - University of Toronto/CREAM/TCGA_ATAC/COREs/TCGA-06-A5U0-01A_SigCut1.65_.bed"

COREs <- CREAM(in_path = input_path, optimize = FALSE)
COREs <- COREs[,1:3]

Genes_TargetFile <- GeneMatch(Regions = COREs, RefTSS = RefTSS, TSS_distance = TSS_Dist)

