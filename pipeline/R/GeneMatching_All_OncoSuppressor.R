## Clear the workspace and import the data
rm(list = ls())
library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(ChIPpeakAnno) 
########################
GeneMatch <- function(file_path, RefTSS, TSS_distance = 1e3){
  
  chr_vec_ref <- as.character(RefTSS$chr)
  tss_ref <- as.numeric(RefTSS$tss)
  gene_names <- as.character(RefTSS$gene)
  #####
  Table <- read.table(file_path)
  Table <- as.matrix(Table)
  # Table <- Table[,c(1:3)]
  # colnames(Table) <- c("chr", "start", "end")
  chr_vec <- as.character(Table[,1])
  start_vec <- as.numeric(Table[,2])
  end_vec <- as.numeric(Table[,3])
  ######
  # Table_Clust  <- as.data.frame(cbind(Table[,1], Table[,2], Table[,3]))
  # colnames(Table_Clust) <- c("seqnames", "start", "end")
  # GR_TargetFeat <- toGRanges(Table_Clust, format="BED", header=FALSE) 
  # peakAnno <- annotatePeak(GR_TargetFeat, tssRegion=c(-TSS_distance, TSS_distance),
  #                          TxDb=txdb, annoDb="org.Hs.eg.db")
  # GeneMatch <- unique(as.character(peakAnno@anno$SYMBOL))
  
  GeneMatch <- c()
  for(reg_iter in 1:nrow(Table)){
    matched_ind <- which(chr_vec[reg_iter] == chr_vec_ref & 
                           tss_ref-TSS_distance < end_vec[reg_iter] & 
                           tss_ref+TSS_distance > start_vec[reg_iter])
    GeneMatch <- c(GeneMatch, gene_names[matched_ind])
  }
  ##
  GeneMatch <- unique(GeneMatch)
  
  return(GeneMatch)
}


Genes_FindingInRef <- function(target_Genes, ref_genes){
  
  
}
########################
TSS_Dist <- 100
RefTSS <- read.table("~/OneDrive - University of Toronto/CREAM/refGene_hg38_TSS.txt")
ClustDir <- "~/OneDrive - University of Toronto/CREAM/TCGA_ATAC/COREs/"
OutDir <- "~/OneDrive - University of Toronto/Arvind/AaronS/COREs/IncludedPeaks/"
ClusFiles <- list.files(path = ClustDir,pattern = "_SigCut1.65_.bed")

Genes_TargetFile <- GeneMatch(file_path = paste(ClustDir, ClusFiles[1],
                                                sep = "", collapse = ""),
                              RefTSS = RefTSS, TSS_distance = TSS_Dist)

########################
oncogene_info <- read.csv("~/Desktop/Obel/Data/Onco_HUMAN.csv")
Oncogenes <- unique(unlist(lapply(oncogene_info$Gene.names,
                                  function(gene_names){
                                    strsplit(as.character(gene_names),
                                             " ")[[1]]})))

suppresor_info <- read.csv("~/Desktop/Obel/Data/TumorSuppressor_HUMAN.csv")
Suppresor_genes <- unique(unlist(lapply(suppresor_info$Gene.names,
                                        function(gene_names){
                                          strsplit(as.character(gene_names),
                                                   " ")[[1]]})))

Oncogenes_mathced <- intersect(Genes_TargetFile, Oncogenes)
Suppresor_genes_mathced <- intersect(Genes_TargetFile, Suppresor_genes)

matched_oncogenes_info <- oncogene_info[unlist(lapply(Oncogenes_mathced, function(target_gene){
  which(grepl(target_gene, oncogene_info$Gene.names))})),]

matched_suppresorgenes_info <- suppresor_info[unlist(lapply(Suppresor_genes_mathced, function(target_gene){
  which(grepl(target_gene, suppresor_info$Gene.names))})),]

############################################################
# matched_oncogenes_info and matched_suppresorgenes_info
# should be shown to the user and Genes_TargetFile
# should be provided as a downloadable file
############################################################
View(matched_oncogenes_info)
