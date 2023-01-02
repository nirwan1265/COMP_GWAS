pvalue.combine <- function(path, filename, chr, pca, organism){
  preprocess_data <- preprocess(path, filename, chr,  organism)
  geno <- geno <- geno_ld(chr)
  pca <- pca_ld()
  
  
  
}