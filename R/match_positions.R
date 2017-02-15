
extract_positions <- function(target_file, dbstem, out, mutation_position=179410829, mutation_wild_type="A") {
  #  target_file= "data/HMERF_Haplotype.csv" 
  #  dbstem = "tmpfiles/impute1K1.impute" 
  #  out= "match1K_phase1" 
  #  mutation_position=179410829
  # mutation_wild_type="A"
  data_hmerf <- fread(target_file, showProgress = FALSE)
  setkey(data_hmerf, position)
  
  data_1K <- fread(paste(dbstem, "legend", sep="."), showProgress = FALSE)
  data_1K[, rownum := 1:nrow(data_1K)]
  setkey(data_1K, pos)
  matching_positions <- data_hmerf[data_1K, nomatch=0]  
  missing_positions <- data_hmerf[!(position %in% data_1K$pos)]
  
  code <- list()
  for (i in 1:nrow(matching_positions)) 
    code[[i]] <- c(matching_positions$allele0[i], matching_positions$allele1[i])
  
  sample_names <-  scan(paste(dbstem,"hap.indv",sep="."), what=character(0))
  cat(sample_names, "\n",sep="\n", file=paste("tmpfiles/", out,".individuals", sep="" ))
  fread(paste(dbstem, "hap.indv", sep="."), showProgress = FALSE)
  
  all_1K <- fread(paste(dbstem, "hap", sep="."), showProgress = FALSE)
  all_1K <- all_1K[matching_positions$rownum,]
  
  genotypes <- matrix("", nrow=nrow(matching_positions)+1, ncol=ncol(all_1K))
  for (index in 1:nrow(matching_positions)) {
    row <- as.integer(all_1K[index, ])+1
    genotypes[index,] <- code[[index]][row]
  }
  genotypes[nrow(matching_positions)+1, ] <- mutation_wild_type   ## add mutation
  
  mp <- copy(matching_positions)                         ## need to make a copy so that I can remove columns
  add <- data_hmerf[position==mutation_position]
  allele <- strsplit(as.character(add[,3]),'')[[1]]
  mutation_variant <- ifelse(mutation_wild_type==allele[1], allele[2], allele[1])
  add$allele0 <- mutation_wild_type
  add$allele1 <- mutation_variant 
  add$ID <- "targetmutation"
  
  mp <- rbind(mp, add, fill=TRUE )
  o <- order(mp$position)
  genotypes <- genotypes[o,]
  mp <- mp[o,]
  
  
  fwrite(
    mp[,.(ID, pos=position, allele0, allele1)],
    file=paste("tmpfiles/", out, ".legend", sep=""), showProgress = FALSE
  )
  
  mp[, c("ID","allele0","allele1","rownum") := NULL]     ## And put the mutation back
  fwrite(mp, file = paste("tmpfiles/HMERF", out, ".csv", sep=""), showProgress = FALSE)
  fwrite(data.table(genotypes), sep=" ", col.names = FALSE, file=paste("tmpfiles", out, sep="/"))
  list(missing=missing_positions, matching=matching_positions)
}
