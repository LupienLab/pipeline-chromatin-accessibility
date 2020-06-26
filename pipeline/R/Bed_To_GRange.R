rm(list = ls())
library(data.table)
library(GenomicRanges)
###
reference_name <- "TEonly_repeats_hg38"
file_name_pattern <- "bed"
reference_path <- '../static/TEonly_repeats_hg38/'
out_dir <- './'

reference_files <- list.files(reference_path, pattern = file_name_pattern)
###
Grange_reference_set <- list()
for(SamIter in 1:length(reference_files)){ 
  #print(SamIter)
  query_mat2 <- fread(paste(reference_path, reference_files[SamIter],
                            sep = "", collapse = ""),
                      select = c(1:3), data.table=FALSE)
  ##########
  gr <- GRanges(seqnames=Rle(query_mat2$V1),
                ranges = IRanges(query_mat2$V2, end=query_mat2$V3))
  ##########
  Grange_reference_set[[SamIter]] <- gr
}

names(Grange_reference_set) <- reference_files

saveRDS(Grange_reference_set,
        file = paste(out_dir, reference_name,
                     "_AllPeaks_Granges.RDS",
                     sep = "", collapse = ""))