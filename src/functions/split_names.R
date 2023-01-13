# Subsetting just the genename from the gff3 files

split.names <- function(x,split){
  split.genename <- unlist(strsplit(x, split = ';', fixed = TRUE))[1]
  split.genename2 <- unlist(strsplit(split.genename, split = ":", fixed = TRUE))[2]
  return(split.genename2)
}

split.names_m <- function(x,split){
  split.genename <- unlist(strsplit(x, split = ';', fixed = TRUE))[1]
  split.genename2 <- unlist(strsplit(split.genename, split = "=", fixed = TRUE))[2]
  return(split.genename2)
}
