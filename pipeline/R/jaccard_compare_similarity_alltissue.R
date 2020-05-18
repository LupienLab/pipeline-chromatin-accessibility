rm(list = ls())
####Library######################################
library(optparse)
library(reticulate)
library(GenomicRanges)
library(data.table)
library(parallel)
library(survcomp)
library(survival)
########################
####User Input######################################
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset bed file", metavar="character"),
  make_option(c("-t", "--top_SamNum"), type="integer", default=NULL,
              help="Top SamNum", metavar="number"),
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="sample name", metavar="character"),
  make_option(c("-c", "--tissue"), type="character", default=NULL,
              help="tissue to match similarity", metavar="character"),
  make_option(c("-d", "--dir"), type="character", default=".",
              help="output directory", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("No arguments given.", call.=FALSE)
}

#####################################################
###################
jaccard_sim <- function(gr1, gr2){
  #####
  CommonUnique <- length(which(countOverlaps(gr1, gr2) > 0))
  Jacard <- CommonUnique/(length(gr1)+length(gr1)-CommonUnique)
  #####
  return(Jacard)
}
####################################
jaccard_simtoref <- function(query_file, reference_Granges, reference_names){
  
  query_mat1 <- fread(query_file,select = c(1:3), data.table=FALSE)
  ########
  gr1 <- GRanges(seqnames=Rle(query_mat1$V1),
                 ranges = IRanges(query_mat1$V2,
                                  end=query_mat1$V3))
  ########
  n_core <- min(detectCores(all.tests = FALSE, logical = TRUE), 3)
  jaccard <- unlist(mclapply(1:length(reference_Granges),
                             function(SamIter){
                               jaccard_sim(gr1, reference_Granges[[SamIter]])},
                             mc.cores = n_core))
  
  # jaccard <- rep(0, length(reference_Granges))
  # for(SamIter in 1:length(reference_Granges)){
  #   ##########
  #   gr2 <- reference_Granges[[SamIter]]
  #   ##########
  #   jaccard[SamIter] <- jaccard_sim(gr1, gr2)
  # }
  ###
  names(jaccard) <- reference_names
  ###
  return(jaccard)
}
############################
############################
source_python("../pipeline/python/checkBedfileQuality.py")
checkBedFile(opt$file)
query_file <- opt$file
reference_path <- "../pipeline/static/TCGA_SigCut1.65_AllPeaks_Granges.RDS"


Grange_TCGA <- readRDS(reference_path)
reference_names <- names(Grange_TCGA)
jaccard_vec <- jaccard_simtoref(query_file=query_file,
                                reference_Granges = Grange_TCGA,
                                reference_names = reference_names)

quantile(jaccard_vec)
sort(jaccard_vec, decreasing = T)[1:5]
####################################
#Sample_Name <- gsub(".bed", "", query_file)
Sample_Name <-opt$sample
Output_Dir <- opt$dir
PhenoVec <- readRDS("../pipeline/static/TCGA_ATACSamples_Phenotype.rds")


Sorted_Sim <- sort(jaccard_vec, decreasing = T)
Sorted_Sim <- Sorted_Sim#/max(Sorted_Sim)
PhenoMat <- c()
for(SamIter in names(Sorted_Sim)){
  
  #######
  Matched_Ind <- which(Survival$sample == SamIter)
  #######
  if(length(Matched_Ind) == 1){
    PhenoMat <- rbind(PhenoMat, c(SamIter, PhenoVec[which(names(PhenoVec) == SamIter)],
                                  Sorted_Sim[SamIter]))
  }
  
}

Top_SamNum <- opt$top_SamNum

colnames(PhenoMat) <- c("sample","tissue","survival status","time (day)","similarity")

##########################################
setwd(Output_Dir)
write.csv(PhenoMat[1:Top_SamNum,], paste("Phenotypic_info_Top", Top_SamNum, "samples.csv", sep= ""),quote=FALSE,row.names=FALSE)

TopTissue_Vec <- unlist(PhenoMat[1:Top_SamNum,"tissue"])
write.csv(names(table(TopTissue_Vec)[
  which(table(TopTissue_Vec) ==max(table(TopTissue_Vec)))]),
  paste("Mostsimilar_tissue_Top", Top_SamNum, "samples.csv", sep= ""),quote=FALSE,row.names=FALSE)
##########################################
# plotting similarity and survival
##########################################
Target_Tissues <- opt$tissue
TargetTissue_Ind <- which(PhenoMat[,"tissue"] %in% Target_Tissues)
PhenoMat <- PhenoMat[TargetTissue_Ind,]
PhenoMat[,"similarity"] <- as.numeric(PhenoMat[,"similarity"])/max(as.numeric(PhenoMat[,"similarity"]))

Top_SamNum <- min(Top_SamNum, length(TargetTissue_Ind))
##########################################

pdf(paste(Sample_Name, ".pdf",
          sep = "", collapse = ""), width=10,height=5)
par(mar=c(5,10,2,3)+0.1,mgp=c(3,1,0))
par(mfrow = c(1,2))
barplot(rev(as.numeric(PhenoMat[1:Top_SamNum,"similarity"])),
        horiz = T, beside = T, las = 2,
        xlim = c(0,1),
        xlab = "similarity",names.arg = rev(PhenoMat[
          1:Top_SamNum,"sample"]), 
        col = "darkred",cex.axis = 1.5, cex.lab = 1.5)
#########