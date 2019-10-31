## Code for getting 1000 genomes data using tabix from within R
library(GenomicRanges)
run_tabix <- function(region, add=10000) {
  
  get_1000_genomes_file <- function(chromosome, phase=3) {
    paste("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL", chromosome, "phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", sep=".")
  }
  
  getpositions <- function() {
    vv <- as.character(region)
    if (substring(vv,1,3)=="chr")
      vv <- substring(vv, 4)
    
    if (as.character(strand(region))!="*")
      vv <- substring(vv, 1, nchar(vv)-2)
    
    vv
  }
  
  chrom <- as.character(seqnames(region))
  
  ofile <- tempfile()
  position <- substring( as.character(region+add), first=4)
  cat("tabix -fh", get_1000_genomes_file(chrom), getpositions(),"\n")
  system2("tabix", args=c("-fh",
                          get_1000_genomes_file(chrom), 
                          getpositions()),
          stdout=ofile)
  return(ofile)
}


