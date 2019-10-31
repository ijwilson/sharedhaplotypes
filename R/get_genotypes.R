## Get genptypes from the genotype data file
## from an illumina file

library(data.table)
## need a map file first
if (!file.exists("SNPID.rda")) {
  mapfile <- fread("/users/nijw/lustre/BrainBank/A2271_Reports/A2271_Batch_01_SNP_Table.txt")
  SNPID <- mapfile[,.(Index, Name, Chr, Position, SNP)]
  SNPID[, A:= substr(SNP,2,2)]
  SNPID[, B:= substr(SNP,4,4)]
  head(SNPID)
  save(SNPID, file="SNPID.rda")
  rm(mapfile)
}


get_genotypes <- function(chrom, start, end, datafile= "/users/nijw/lustre/BrainBank/Genotypes/all.tsv.gz") {
  pos_string <- paste(chrom, ":", start,"-", end, sep="")
  ofile <- tempfile() 
  system2("tabix", c("-fh", datafile, pos_string), stderr="", stdout = ofile)
  res <- fread(ofile)
  unlink(ofile)
  info <- res[,.(name=V1, chrom=V2, position=V3)]
  res[,c("V1", "V2", "V3") := NULL] 
  m <- as.matrix(res)
  m[m=="AA"] <- "0"
  m[m=="AB"] <- "1"
  m[m=="BB"] <- "2"
  m <- array(as.numeric(m), dim=dim(m))
  list(info, g=t(m))
}

fn <- get_genotypes("1", "10000000", "11000000")
i
