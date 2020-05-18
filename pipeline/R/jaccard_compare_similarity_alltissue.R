rm(list = ls())
####Library######################################
library(reticulate)
library(GenomicRanges)
library(data.table)
library(parallel)
library(survcomp)
library(survival)
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
  
  ###
  names(jaccard) <- reference_names
  ###
  return(jaccard)
}
############################
############################
top_samnum <- 5
query_file <- "~/Desktop/Obel/Data/Roadmap_H3K27ac_peaks/E001-H3K27ac.imputed.narrowPeak.bed.nPk"
reference_path <- "~/Desktop/Obel/Data/H3K27ac_Roadmap_AllPeaks_Granges.RDS"
out_dir <- "~/Desktop/Obel/"
analysis_name <- "similarity_roadmap_H3K27ac"
PhenoMat <- readRDS("~/Desktop/Obel/Data/Roadmap_SamplePhenotypes.RDS")


Grange_reference <- readRDS(reference_path)
Grange_reference <- Grange_reference[unlist(lapply(PhenoMat[,"sample_name"],
                                                   function(x){which(grepl(
                                                     x, names(Grange_reference)))}))]

reference_names <- names(Grange_reference)
jaccard_vec <- sort(jaccard_simtoref(query_file=query_file,
                                     reference_Granges = Grange_reference,
                                     reference_names = reference_names),
                    decreasing = T)


write.csv(jaccard_vec, paste(out_dir, analysis_name,".csv",
                             sep = "", collapse = ""))
#######
jaccard_vec <- sort(jaccard_vec, decreasing = T)[1:top_samnum]


matched_samples <- lapply(PhenoMat[,"sample_name"],
                          function(x){which(grepl(
                            x, names(jaccard_vec)))})
names(jaccard_vec) <- paste(PhenoMat[which(unlist(lapply(matched_samples, function(X){length(X) != 0}))),"sample_name"],
                            PhenoMat[which(unlist(lapply(matched_samples, function(X){length(X) != 0}))),"phenotype"], sep = "_")


barplot(rev(jaccard_vec),
        horiz = T, beside = T, las = 2,
        xlim = c(0,1), xlab = "similarity",
        col = "darkred",cex.axis = 1.5, cex.lab = 1.5)

