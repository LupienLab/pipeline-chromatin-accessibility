rm(list=ls())
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(RColorBrewer)
####User Input######################################
library(optparse)
library(reticulate)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset bed file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="bedfile_statistics.pdf",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-d", "--dir"), type="character", default=".",
              help="output directory", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#####################################################


#####################################################
############ Human chromosome information ###########
#####################################################
chr_names = c('chr1', 'chr2', 'chr3', 'chr4',
              'chr5', 'chr6', 'chr7', 'chr8',
              'chr9', 'chr10', 'chr11', 'chr12',
              'chr13', 'chr14', 'chr15', 'chr16',
              'chr17', 'chr18', 'chr19', 'chr20',
              'chr21', 'chr22', 'chrX', 'chrY')
chrsize_list = list(hg19 = c(249250621, 243199373, 198022430, 191154276,
                             180915260, 171115067, 159138663, 146364022,
                             141213431, 135534747, 135006516, 133851895,
                             115169878, 107349540, 102531392, 90354753,
                             81195210, 78077248, 59128983, 63025520,
                             48129895, 51304566, 155270560, 59373566)/1e6,
                    hg38 = c(248956422, 242193529, 198295559, 190214555,
                             181538259, 170805979, 159345973, 145138636,
                             138394717, 133797422, 135086622, 133275309,
                             114364328, 107043718, 101991189, 90338345,
                             83257441, 80373285, 58617616, 64444167,
                             46709983, 50818468, 156040895, 57227415)/1e6)
#####################################################
############ Human chromosome information ###########
#####################################################
genomic_dist_perchr <- function(bed_frame = bedfile,
                                chr_names = chr_names,
                                genome='hg38',
                                loci_type = c("Promoters", "immediateDownstream", 
                                              "fiveUTRs", "threeUTRs", 
                                              "Exons", "Introns"),
                                promoter_proximity = 1000,
                                downstream_dist = 1000){
  bed <- bed_frame[,1:3]
  colnames(bed) <- c("seqnames", "start", "end")
  gr_allchr <- toGRanges(bed, format="BED", header=FALSE) 
  ########
  genomci_dist_list <- c()
  for(chr_iter in chr_names){
    print(chr_iter)
    
    matched_ind <- seqnames(gr_allchr) == chr_iter
    print(length(matched_ind))
    print(length(which(matched_ind)))
    # gr_targetchr <- toGRanges(bed[which(bed$chr == chr_iter),], format="BED", header=FALSE)
    gr_targetchr <- gr_allchr[matched_ind]
    if(length(which(matched_ind)) > 1){
      genomci_dist_list <- cbind(genomci_dist_list, genomic_dist(gr = gr_targetchr,
                                                                 genome = genome,
                                                                 loci_type = c("Promoters", "immediateDownstream", 
                                                                               "fiveUTRs", "threeUTRs", 
                                                                               "Exons", "Introns"),
                                                                 promoter_proximity = 1000,
                                                                 downstream_dist = 1000))
    }else{
      genomci_dist_list<- cbind(genomci_dist_list, rep(0, length(loci_type)))
    }
    
  }
  #######
  # names(genomci_dist_list) <- chr_names
  #######
  return(genomci_dist_list)
}
#####################################################
############ Extracting statistics of ############### 
############ regions in a bed file    ###############
#####################################################
bed_stat <- function(bed_frame, chr_names){
  
  
  bed_stat_frame <- data.frame(matrix(ncol = 3, nrow = 0))
  for(chriter in chr_names){
    matched_ind <- which(bed_frame$chr == chriter)
    if(length(matched_ind) > 0){
      size_dist <- as.numeric(bed_frame$end)[matched_ind]-as.numeric(bed_frame$start)[matched_ind]
      stat_vec <- c(length(matched_ind),
                    median(size_dist)/1e3,
                    max(size_dist)/1e3)
      
      bed_stat_frame <- rbind(bed_stat_frame, stat_vec)
    }else{
      bed_stat_frame <- rbind(bed_stat_frame, c(0, 0, 0))
    }
    
  }
  
  colnames(bed_stat_frame) <- c("region_numbers",
                                "size_median",
                                "size_max")
  ########
  return(bed_stat_frame)
  
}
#####################################################
############ Extracting statistics of ############### 
############ regions in a bed file    ###############
#####################################################
bed_stat_plot <- function(bed_frame, bed_stat_frame, chrsize_vec,
                          chr_names, genome = 'hg38',
                          loci_type = c("Promoters", "immediateDownstream", 
                                        "fiveUTRs", "threeUTRs", 
                                        "Exons", "Introns"),
                          promoter_proximity = 1000,
                          downstream_dist = 1000,
                          output_path){
  
  pdf(output_path, width = 20, height = 10)
  
  par(mfrow = c(1, 5))
  par(mar=c(5,5,2,3)+0.1,mgp=c(4,1.5,0))
  barplot(rev(bed_stat_frame$region_numbers),
          names.arg = rev(chr_names),
          horiz = T, las=2, col = 'forestgreen',
          cex.axis = 1.5, cex.names = 1.5, cex.lab=1.5,
          xlab = 'number of regions', ylab = '')
  ###
  barplot(rev(bed_stat_frame$size_median),
          names.arg = rev(chr_names),
          horiz = T, las=2, col = 'forestgreen',
          cex.axis = 1.5, cex.names = 1.5, cex.lab=1.5,
          xlab = 'average regions size (kbp)', ylab = '')
  ###
  barplot(rev(bed_stat_frame$size_max),
          names.arg = rev(chr_names),
          horiz = T, las=2, col = 'forestgreen',
          cex.axis = 1.5, cex.names = 1.5, cex.lab=1.5,
          xlab = 'maximum regions size (kbp)', ylab = '')
  ###
  barplot(rev(chrsize_vec),
          names.arg = rev(chr_names),
          horiz = T, las=2, col = 'forestgreen',
          cex.axis = 1.5, cex.names = 1.5, cex.lab=1.5,
          xlab = 'chromosome size (mbp)', ylab = '')
  ###
  barplot(genomic_dist_perchr(bed_frame = bed_frame,
                              chr_names = rev(chr_names),
                              genome = genome,
                              loci_type = c("Promoters", "immediateDownstream", 
                                            "fiveUTRs", "threeUTRs", 
                                            "Exons", "Introns"),
                              promoter_proximity = 1000,
                              downstream_dist = 1000),
          names.arg = rev(chr_names),
          horiz = T, las=2, col = col_vector[1:length(loci_type)],
          cex.axis = 1.5, cex.names = 1.5, cex.lab=1.5,
          xlab = 'chromosome size (mbp)', ylab = '')
  dev.off()
}


######################################################

source('../pipeline/R/genomic_dist.R')
source_python("../pipeline/python/checkBedfileQuality.py")
checkBedFile(opt$file)

input_file <- opt$file
bedfile <- read.table(input_file,stringsAsFactors = F, check.names = F)
bedfile <- bedfile[,c(1:3)]
colnames(bedfile) <- c('chr', 'start', 'end')
genome_assembly <- 'hg38'
loci_typevec <- c("Promoters", "immediateDownstream", 
                  "fiveUTRs", "threeUTRs", 
                  "Exons", "Introns")
promoter_proximity <- 1000
downstream_dist <- 1000
#########
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#########

stat_frame = bed_stat(bed_frame = bedfile[1:10,], chr_names = chr_names)


bed_stat_plot(bed_frame = bedfile[1:10,],
              bed_stat_frame = stat_frame,
              chrsize_vec = chrsize_list$hg38,
              chr_names = chr_names,
              genome = genome_assembly,
              loci_type = loci_typevec,
              promoter_proximity = promoter_proximity,
              downstream_dist = downstream_dist,
              output_path=paste(opt$dir,"/",opt$out,sep=""))
