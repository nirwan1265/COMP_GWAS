# Function for putting everything in a proper list
data_wrangle <- function(path, phenoname, chr, organism){
  preprocess_data <- preprocess(path, phenoname, chr,  organism)
  geno <- geno <- geno_ld(organism, chr)
  pca <- pca_ld(organism)
  return(list(preprocess = preprocess_data,genotype = geno,PCA = pca))
}


#Not use
#path = "/Users/nirwan/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data"
##phenoname <- "tot"
#organism <- "Sorghum bicolor"
#chr <- 2
#trial <- data_wrangle(path, phenoname, chr, organism)



