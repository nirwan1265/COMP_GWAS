# Function for putting everything in a proper list
data_wrangle <- function(path, phenoname, chr, organism){
  preprocess_data <- preprocess(path, phenoname, chr,  organism)
  geno <- geno <- geno_ld(organism, chr)
  pca <- pca_ld(organism)
  return(list(preprocess = preprocess_data,genotype = geno,PCA = pca))
}


#Not use
path = "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data"
phenoname <- "tot"
organism <- "Zea"
chr <- 1
trial <- data_wrangle(path, phenoname, chr, organism)



