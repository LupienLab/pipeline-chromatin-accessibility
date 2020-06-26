
# load dependencies
#----------------------------------------------------------------
library(chromVAR)
library(Matrix)
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BiocParallel)
#register(SerialParam())
library(BSgenome.Hsapiens.UCSC.hg38)

library(data.table)
library(tidyverse)


run_chromvar=function(narrowPeak,
                      repeatrds,
                      out_dir){
  
  #----------------------------------------------------------------
  ## Load Countdata
  #----------------------------------------------------------------
  set.seed(2019)
  
  print("read in count data")
  my_counts_matrix <- fread(narrowPeak, header = F)[,1:3]
  rownames(my_counts_matrix) <- paste(my_counts_matrix$seqname, my_counts_matrix$start, my_counts_matrix$end, sep="_")
  
  
  ## Remove random chromosomes
  
  toMatch <- c("random", "alt", "Un", "chrM", "EBV")
  my_counts_matrix <- subset(my_counts_matrix, !(grepl(paste(toMatch, collapse="|"), my_counts_matrix$seqname)))
  print("counts matrix")  
  print(head(my_counts_matrix))
  fragment_counts <- makeSummarizedExperimentFromDataFrame(my_counts_matrix)
  assayNames(fragment_counts) <- "counts"
  
  #----------------------------------------------------------------
  ## add gc content
  #----------------------------------------------------------------
  print("correcting for GC content")
  
  fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
  counts_filtered <- filterPeaks(fragment_counts,min_fragments_per_peak = 1, non_overlapping = TRUE)
  
  print("saving counts_filtered")
  saveRDS(counts_filtered, paste0(out_dir,"/counts_filtered.rds"))
  rm(fragment_counts)
  
  # #----------------------------------------------------------------
  # ## Get motifs and what peaks contain motifs --
  # #----------------------------------------------------------------
  
  rm(my_counts_matrix)
  
  print("generate annotation dataset")
  
  my_annotation_df <- as.data.frame(readRDS(repeatrds))
  rownames(my_annotation_df) <- paste(my_annotation_df$seqnames, my_annotation_df$start, my_annotation_df$end, sep="_")
  my_annotation_df <- my_annotation_df[rownames(counts_filtered),]
  
  print("generate annotation dataset")
  anno_ix <- getAnnotations(as.matrix(my_annotation_df[,4:ncol(my_annotation_df)]), rowRanges = rowRanges(counts_filtered))
  print("saving annotation")
  saveRDS(anno_ix, paste0(out_dir,"/anno_ix.rds"))
  
  #----------------------------------------------------------------
  ## compute deviation
  #----------------------------------------------------------------
  set.seed(2019)
  
  print("computing expectations")
  expectations <- computeExpectations(counts_filtered)
  print("saving expectations")
  saveRDS(expectations, paste0(out_dir,'/expectations.rds'))
  
  print("computing background peaks")
  background_peaks <- getBackgroundPeaks(counts_filtered)
  print("saving background peaks")
  saveRDS(background_peaks, paste0(out_dir,'/background_peaks.rds'))
  
  print("computing deviation")
  dev <- computeDeviations(object = counts_filtered, annotations = anno_ix, background_peaks = background_peaks)
  
  print("saving deviation results")
  saveRDS(dev, paste0(out_dir,"/devobj.rds"))
  
  z.scores = deviationScores(dev) ## deviation Z-score
  dev.scores = deviations(dev) ## bias corrected deviations
  
  write.table(z.scores, file=paste0(out_dir,"/zscore.tsv"), col.names=T, row.names=T, sep="\t")
  write.table(dev.scores, file=paste0(out_dir,"/deviations.tsv"), col.names=T, row.names=T, sep="\t")
  
  #----------------------------------------------------------------
  ## compute variablity
  #----------------------------------------------------------------
  print("Computing Variability")
  variability <- computeVariability(dev)
  write.table(variability, file=paste0(out_dir,"/Variability.tsv"), col.names=T, row.names=F, sep="\t")
  
}


out_dir <- 'out_dir_repats'
dir.create(out_dir)
run_chromvar(narrowPeak="/mnt/work1/users/lupiengroup/People/ankita/pipeline-chromatin-accessibility/data/Peaks/test-FNA-Notch-06-PM-pos-12-Bulk_L002_peaks.filtered.narrowPeak",
             repeatrds="bm_for_repeatElements.rds",
             out_dir = out_dir)

