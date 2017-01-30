

## script to read PHASE pairs output

readpairs <- function(filename)
{
  conn = file(filename,"r")
  res <- readLines(conn)
  inds = grep("IND",res)
  indnames <- substring(res[inds],6)
  n <- diff(c(inds,length(res)+1))-1
  res <- gsub("\\s","",res[-inds])
  namecol <- rep(indnames,n)
  data <- matrix(unlist(strsplit(res,",")),ncol=3,byrow=T)
  r = data.frame(cbind(namecol,data))
  r[,4] <- as.numeric(paste(r[,4]))
  close(conn)
  r
}

################### From Gerry Directory

splitter <- function(pairdata,wh=1) {
  sp <- function(x) {
    c(substr(x,15,15),substr(x,1,7),substr(x,8,8),substr(x,8,11),substr(x,12,13)
      ,paste(substr(x,14,14),substr(x,16,19),sep="")
      ,substr(x,20,22),substr(x,23,24))
  }
  t(sapply(paste(pairdata[,1+wh]),sp))
}

# function to split pairs data by a single positions
splitter <- function(pairdata,splits=c(1)) {
  nsplit <- length(splits)+1
  l <- nchar(paste(pairdata[1,3]))
  splitvecstart <- c(1,splits)
  splitvecend <- c(splits,l)
  pd <- rep(paste(pairdata[,2]),rep(nsplit,nrow(pairdata)))
  s1 <- substring(pd,splitvecstart,splitvecend)
  pd1 <- rep(paste(pairdata[,3]),rep(nsplit,nrow(pairdata)))
  s1 <- substring(pd,splitvecstart,splitvecend)
  data.frame(pairdata,hap1 = matrix(s1,ncol=nsplit,byrow=T),hap2 = matrix(s1,ncol=nsplit,byrow=T))
}