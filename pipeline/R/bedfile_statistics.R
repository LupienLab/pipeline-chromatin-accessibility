####################################################
rm(list=ls())
####User Input######################################
library(optparse)
library(reticulate)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="sample.pdf", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

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
bed_stat_plot <- function(bed_stat_frame, chrsize_vec, chr_names, output_path){
  
  pdf(output_path, width = 15, height = 10)
  
  par(mfrow = c(1, 4))
  par(mar=c(5,5,2,3)+0.1,mgp=c(4,1.5,0))
  barplot(rev(bed_stat_frame$region_numbers),
          names.arg = rev(chr_names),
          horiz = T, las=2, col = 'forestgreen',
          cex.axis = 1.5, cex.names = 1.5, cex.lab=1.5,
          xlab = 'number of regions', ylab = '')
  barplot(rev(bed_stat_frame$size_median),
          names.arg = rev(chr_names),
          horiz = T, las=2, col = 'forestgreen',
          cex.axis = 1.5, cex.names = 1.5, cex.lab=1.5,
          xlab = 'average regions size (kbp)', ylab = '')
  barplot(rev(bed_stat_frame$size_max),
          names.arg = rev(chr_names),
          horiz = T, las=2, col = 'forestgreen',
          cex.axis = 1.5, cex.names = 1.5, cex.lab=1.5,
          xlab = 'maximum regions size (kbp)', ylab = '')
  barplot(rev(chrsize_vec),
          names.arg = rev(chr_names),
          horiz = T, las=2, col = 'forestgreen',
          cex.axis = 1.5, cex.names = 1.5, cex.lab=1.5,
          xlab = 'chromosome size (mbp)', ylab = '')
  dev.off()
}


######################################################
source_python("../checkBedfileQuality.py")
checkBedFile(opt$file)
bedfile <- read.table(opt$file, stringsAsFactors = F, check.names = F )
bedfile <- bedfile[,c(1:3)]
colnames(bedfile) <- c('chr', 'start', 'end')

stat_frame = bed_stat(bed_frame = bedfile, chr_names = chr_names)

bed_stat_plot(bed_stat_frame = stat_frame,
              chrsize_vec = chrsize_list$hg38,
              chr_names = chr_names,
              output_path = opt$out)
