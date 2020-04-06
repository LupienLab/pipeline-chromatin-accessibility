rm(list = ls())
library(GenomicRanges)
library(data.table)
library(parallel)
library(survcomp)
library(survival)
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
setwd("~/OneDrive - University of Toronto/CREAM/GBM_Calgary/pGBM_peakData_for_Ali/")
query_file <- 'pat3_relapse.bed'
reference_path <- "~/Desktop/Obel/Data/TCGA_SigCut1.65_AllPeaks_Granges.RDS"


Grange_TCGA <- readRDS(reference_path)
reference_names <- names(Grange_TCGA)
jaccard_vec <- jaccard_simtoref(query_file=query_file,
                                reference_Granges = Grange_TCGA,
                                reference_names = reference_names)

quantile(jaccard_vec)
sort(jaccard_vec, decreasing = T)[1:5]
####################################
setwd("~/OneDrive - University of Toronto/CREAM/TCGA_ATAC/")

Sample_Name <- gsub(".bed", "", query_file)
Figure_Dir <- "~/OneDrive - University of Toronto/CREAM/GBM_Calgary/pGBM_peakData_for_Ali/Figures/"
PhenoVec <- readRDS("TCGA_ATACSamples_Phenotype.rds")
Survival <- read.table("GDC-PANCAN.survival.tsv",
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

colnames(PhenoMat) <- c("sample", "tissue", "survival status", "time (day)", "similarity")
##########################################
# plotting similarity and survival
##########################################
Target_Tissues <- c("GBM")
TargetTissue_Ind <- which(PhenoMat[,"tissue"] %in% Target_Tissues)
PhenoMat <- PhenoMat[TargetTissue_Ind,]
PhenoMat[,"similarity"] <- as.numeric(PhenoMat[,"similarity"])/max(as.numeric(PhenoMat[,"similarity"]))
Top_SamNum <- 5#floor(length(TargetTissue_Ind)/2)
Top_SamNum <- min(Top_SamNum, length(TargetTissue_Ind))

##########################################
setwd(Figure_Dir)
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

