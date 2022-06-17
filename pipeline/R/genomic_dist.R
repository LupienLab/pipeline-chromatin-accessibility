genomic_dist <- function(gr, genome="hg38", loci_type = c("Promoters", "immediateDownstream", 
                                                           "fiveUTRs", "threeUTRs", 
                                                           "Exons", "Introns"),
                         promoter_proximity = 1000, downstream_dist = 1000){
  
  if(genome=="hg19"){
    ref_genome <- TxDb.Hsapiens.UCSC.hg19.knownGene
  }else if(genome == "hg38"){
    ref_genome <- TxDb.Hsapiens.UCSC.hg38.knownGene
  }else{
    stop("Genome was not found!")
  }
  aCR <- assignChromosomeRegion(gr, nucleotideLevel = FALSE, 
                                precedence = loci_type,
                                TxDb = ref_genome,
                                proximal.promoter.cutoff = promoter_proximity,
                                immediate.downstream.cutoff = downstream_dist)
  return(aCR$percentage)
}
