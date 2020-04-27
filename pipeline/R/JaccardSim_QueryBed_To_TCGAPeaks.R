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
Survival <- read.table("../pipeline/static/GDC-PANCAN.survival.tsv",
                       stringsAsFactors = F, check.names = F, header = T)



Sorted_Sim <- sort(jaccard_vec, decreasing = T)
Sorted_Sim <- Sorted_Sim#/max(Sorted_Sim)
PhenoMat <- c()
for(SamIter in names(Sorted_Sim)){
  
  #######
  Matched_Ind <- which(Survival$sample == SamIter)
  #######
  if(length(Matched_Ind) == 1){
    PhenoMat <- rbind(PhenoMat, c(SamIter, PhenoVec[which(names(PhenoVec) == SamIter)],
                                  Survival[Matched_Ind,c("OS","OS.time")],
                                  Sorted_Sim[SamIter]))
  }
  
}

colnames(PhenoMat) <- c("sample","tissue","survival status","time (day)","similarity")

##########################################
# plotting similarity and survival
##########################################
Target_Tissues <- c("GBM")
TargetTissue_Ind <- which(PhenoMat[,"tissue"] %in% Target_Tissues)
PhenoMat <- PhenoMat[TargetTissue_Ind,]
PhenoMat[,"similarity"] <- as.numeric(PhenoMat[,"similarity"])/max(as.numeric(PhenoMat[,"similarity"]))
Top_SamNum <- opt$top_SamNum
Top_SamNum <- min(Top_SamNum, length(TargetTissue_Ind))
##########################################
outputpath=paste(opt$dir,"/",opt$out,sep="")
write.csv(PhenoMat[1:Top_SamNum,], paste(outputpath,"Phenotypic_info_Top", Top_SamNum, "samples.csv", sep= ""),quote=FALSE)

TopTissue_Vec <- unlist(PhenoMat[1:Top_SamNum,"tissue"])
write.csv(names(table(TopTissue_Vec)[
  which(table(TopTissue_Vec) ==max(table(TopTissue_Vec)))]),
  paste(outputpath,"Mostsimilar_tissue_Top", Top_SamNum, "samples.csv", sep= ""),quote=FALSE)
##########################################
setwd(Output_Dir)
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
TimeVec <- as.numeric(PhenoMat[,"time (day)"])/365
EventVec <- as.character(PhenoMat[,"survival status"])
Class <- rep(0,length(TimeVec))
Class[1:Top_SamNum] <- 1
EventVec <- as.numeric(EventVec)
surv.obj <- survfit(Surv(TimeVec, EventVec) ~ Class)

Cox_Fit <- summary(coxph(Surv(time = TimeVec, event = EventVec) ~ Class))
Hazard <- Cox_Fit$coefficients[1]
pvalue <- min(Cox_Fit$coefficients[5], (1-Cox_Fit$coefficients[5]))
#######

plot(surv.obj,col =c("blue", "darkred"),
     lty = 1,lwd = 3,  cex.lab = 1.5, cex.axis = 1.5,frame.plot=F,
     xlab = "Time (years)",xlim = c(0,0.1*floor(max(TimeVec)*10)),
     ylab = "Probability of Overall Survival",
     main = paste("Hazard = ", round(Hazard,2),
                  "(p-value = ", formatC(pvalue, format = "e", digits = 2), ")",
                  sep = "", collapse = ""), xaxt = "n")
axis(side = 1, at = c(0, 0.05*floor(max(TimeVec)*10), 0.1*floor(max(TimeVec)*10)),
     labels = c(0, 0.05*floor(max(TimeVec)*10), 0.1*floor(max(TimeVec)*10)), cex.axis = 1.5)
dev.off()

