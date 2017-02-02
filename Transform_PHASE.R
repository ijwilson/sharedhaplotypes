#!/usr/bin/Rscript
##################################################################
Splitmarkers <- function(str,SNP)
  {
    str <- paste(str)
     if (str=="") {
       if (SNP) return (c("?","?"))
       else return(c("-1","-1"))
     }
    if (grepl('/',str)==TRUE) return(unlist(strsplit(str, "/")))
    else return(unlist(strsplit(str, "")))
  }
##################################################################
## test for SNPs by determining whether the first character of the line is a base
locusSNP <- function(col) {
  issnp <- function(yy) {
    if (substr(yy,1,1)  %in% c("A","C","T","G")) return (TRUE)
    return (FALSE)
  }
  return(sapply(col,issnp)) 
}
####################################################################

transform_phase <- function(input_stem, target_data, ped, output_stem="phase1K", target_pops=c("CEU", "GBR")) {
  is_unrelated <- ped[,3]=="0" & ped[,4]=="0"
  ped = ped[is_unrelated & ped$Population %in% target_pops,]
  all_sample_names <- scan(paste(input_stem,"individuals", sep="."), what=character(0))
  keep <- all_sample_names %in% ped$Individual.ID
  used_sample_names <- all_sample_names[keep]
  
  ## open output files
  f <- file(paste(output_stem, "inp",sep="."), "w")
  known <- file(paste(output_stem,"known", sep="."), "w")        ## phase known
  
  n <- ncol(target_data) -2
  markers <- nrow(target_data)
  SNPs <- locusSNP(target_data[, 3])       ## the most complete data
  
  cat(n+length(used_sample_names), markers, sep="\n", file=f)
  cat("P ",target_data[,1],"\n",file=f)
  loc = rep("M", markers)
  loc[SNPs] <- "S"
  cat(loc,"\n",file=f)
  for (ind in 1:n) { 
    col=ind+2
    cat(colnames(target_data)[col],"\n", file=f)
    vals <- mapply(Splitmarkers, hmerf[, col], SNPs)
    cat(vals[1,],"\n",file=f)
    cat(vals[2,],"\n",file=f)
    cat(rep("*",ncol(vals)),"\n",file=known)
  }
  
  raw <- t(as.matrix(read.table(input_stem)))
  ## now add the mutation - first add the column
  raw <- raw[,c(1:14,1,15:ncol(raw))]
  raw[,15] <- "A"
  PhasedLine = paste(rep("0",ncol(raw)),collapse=" ")
  
  raw <- apply(raw,1,paste,collapse=" ")
  raw <- matrix(raw,ncol=2,byrow=T)
  raw <- raw[keep,]
  
  for (i in 1:length(used_sample_names)) {
    cat(used_sample_names[i], "\n", file=f)
    cat(raw[i,],sep="\n", file=f)
    cat(PhasedLine,"\n", file=known)
  }
  
  close(f)
  close(known)
  
  
}

ped3 <- read.csv("data/integrated_call_samples_v2.20130502.ALL.ped", sep="\t", header=TRUE, quote="", stringsAsFactors = FALSE)
hmerf <- read.csv("tmpfiles/HMERFmatch1K_phase3.csv", stringsAsFactors = FALSE)
transform_phase("tmpfiles/match1K_phase3", hmerf, ped3, "tmpfiles/phase1k_p3")


ped <- read.csv("data/phase1_samples_integrated_20101123.ped", sep="\t", header=TRUE, quote="", stringsAsFactors = FALSE)
hmerf <- read.csv("tmpfiles/HMERFmatch1K_phase1.csv", stringsAsFactors = FALSE)

transform_phase("tmpfiles/match1K_phase1", hmerf, ped, "tmpfiles/phase1k_p1")

