rm(list = ls())

reference_path <- '~/OneDrive - University of Toronto/CREAM/TCGA_ATAC/Peaks/'
reference_files <- list.files(reference_path, pattern = "")

Grange_TCGA_Sig165 <- list()
for(SamIter in 1:length(reference_files)){ 
  print(SamIter)
  query_mat2 <- fread(paste(reference_path, reference_files[SamIter],
                            sep = "", collapse = ""),
                      select = c(1:3), data.table=FALSE)
  ##########
  gr <- GRanges(seqnames=Rle(query_mat2$V1),
                ranges = IRanges(query_mat2$V2, end=query_mat2$V3))
  ##########
  Grange_TCGA_Sig165[[SamIter]] <- gr
}


saveRDS(Grange_TCGA_Sig165, file = "H3K27ac_Roadmap_AllPeaks_Granges.RDS")
