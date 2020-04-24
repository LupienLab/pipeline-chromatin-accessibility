rm(list = ls())
####Library######################################
library(optparse)
library(reticulate)
library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(ChIPpeakAnno)
########################
####User Input######################################
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset bed file", metavar="character"),
  make_option(c("-t", "--TSSDist"), type="integer", default=NULL,
              help="TSS Dist", metavar="number"),
  make_option(c("-o", "--out"), type="character",           default="Matching_genes_proximity",
              help="output file [default= %default]", metavar="character"),
  make_option(c("-d", "--dir"), type="character", default="R_analysis",
              help="output directory [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("No arguments given.", call.=FALSE)
}

#####################################################

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
  GeneMatch <- c()
  for(reg_iter in 1:nrow(Table)){
    matched_ind <- which(chr_vec[reg_iter] == chr_vec_ref & 
                           tss_ref-TSS_distance < end_vec[reg_iter] & 
                           tss_ref+TSS_distance > start_vec[reg_iter])
    GeneMatch <- c(GeneMatch, gene_names[matched_ind])
  }

  GeneMatch <- unique(GeneMatch)
  
  return(GeneMatch)
}

########################
source_python("../pipeline/python/checkBedfileQuality.py")
checkBedFile(opt$file)
TSS_Dist <- opt$TSSDist
RefTSS <- read.table("../static/refGene_hg38_TSS.txt")
input_path <- opt$file

Genes_TargetFile <- GeneMatch(file_path = input_path,
                              RefTSS = RefTSS, TSS_distance = TSS_Dist)

########################
oncogene_info <- read.csv("../static/Onco_HUMAN.csv")
Oncogenes <- unique(unlist(lapply(oncogene_info$Gene.names,
                                  function(gene_names){
                                    strsplit(as.character(gene_names),
                                             " ")[[1]]})))

suppresor_info <- read.csv("../static/TumorSuppressor_HUMAN.csv")
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


matched_genes_list <- list(matched_oncogenes_info, matched_suppresorgenes_info)

outputpath=paste(opt$dir,"/",opt$out,sep="")
write.csv(matched_oncogenes_info, paste(outputpath,"_oncogene.txt",sep = ""),quote=FALSE)
write.csv(matched_suppresorgenes_info,paste(outputpath,"_tumor_suppressor.txt",sep = ""),quote=FALSE)
