rm(list = ls())

####Library######################################
library(optparse)
library(data.table)
library(GenomicRanges)

####User Input######################################
option_list = list(
  make_option(c("-f", "--narrowfile"), type="character", default=NULL,
              help="Narrow peaks file", metavar="character"),
  make_option(c("-r", "--repeats"), type="character", default=NULL,
              help="Repeats directory", metavar="character")
  
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$narrowfile) || is.null(opt$repeats)){
  print_help(opt_parser)
  stop("No arguments given.", call.=FALSE)
}

#####################################################

print(1)

#------------------------------------------------
narrowPeaksFile <- opt$narrowfile
rep_folder_dir <- opt$repeats

peaks <- fread(narrowPeaksFile, header = F)[,1:3]
colnames(peaks) <- c('seqnames', 'start', 'end')

rep_paths <- list.files(rep_folder_dir, full.names = T)

query_ranges <- GRanges(seqnames = peaks$seqnames, ranges = IRanges(start = peaks$start, end = peaks$end))

print(2)

for(this_rep_path in rep_paths){
  this_rep <- fread(this_rep_path)
  colnames(this_rep) <- c('seqnames', 'start', 'end')
  
  #remove rows with bad seqnames
  this_rep <- this_rep[which(this_rep$seqnames %in% paste0('chr',c(1:22,'X','Y')))]
  
  this_rep_name <- gsub('.bed.sorted', '', rev(strsplit(this_rep_path,'\\/')[[1]])[1])
  this_rep_ranges <- GRanges(seqnames = this_rep$seqnames, ranges = IRanges(start = this_rep$start, end = this_rep$end))
  
  ovlps <- as.data.table(findOverlaps(query = query_ranges, subject = this_rep_ranges))
  
  peaks[ , this_rep_name] <- 0
  peaks[unique(ovlps$queryHits), this_rep_name] <- 1
}

print(3)

saveRDS(peaks, 'bm_for_repeatElements.rds')
