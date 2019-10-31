source("read_phase_functions.R")
# Read the pairs output.  We have columns
pairs <- readpairs("phase1K.out_pairs")
colnames(pairs) <- c("name", "hap1", "hap2", "pprob")
## Split the phased output 
splitter(pairs,c(4,8,11))
## The first row that begins HG 
w = min(grep("HG",pairs[,1]))
## Get the raw haplotype frequencies
## These are the frequencies from the 
raw <- as.matrix(read.table("match1K"))
raw <- t(raw)
b = apply(raw,1,paste,collapse="")
##
bcase <- mapply(function(x,y) x[1:(y[1]-1),], b, w, SIMPLIFY=FALSE)
casefreq(bcase[[1]],3,nsplit=4)
casefreq(bcase[[1]],2,nsplit=4)

bcontrol <- mapply(function(x,y){return(x[y[1]:y[2],])},b, w, SIMPLIFY=FALSE)
tb1 = table(c(paste(bcontrol[[1]][,7]),paste(bcontrol[[1]][,11])))
table(c(paste(bcontrol[[1]][,7]),paste(bcontrol[[1]][,12])))

b[[1]][1:10,]
table(c(paste(bcontrol[[1]][,6]),paste(bcontrol[[1]][,14])))

casefreq <- function(bbb,cl=1,nsplit=3) {
  fam <- paste(bbb[,1])
  hap <- paste(bbb[,5+cl])
  hap[bbb[,5]=="A"] <- paste(bbb[bbb[,5]=="A",6+nsplit+cl])
  tapply(bbb[,4],data.frame(fam=fam,hap=hap),sum)
}

raw <- as.matrix(read.table("match1K"))
raw <- t(raw)

b = apply(raw,1,paste,collapse="")
u <- substr(paste(b),15,23)=="ATGTGTCAC"

print(sum(u))
