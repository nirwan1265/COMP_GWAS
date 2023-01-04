pvalue.combine <- function(path, phenoname, chr, organism){
  preprocess_data <- preprocess(path, phenoname, chr,  organism)
  geno <- geno <- geno_ld(chr)
  pca <- pca_ld()
  return(list(preprocess = preprocess_data,genotype = geno,PCA = pca))
}



path = "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data"
phenoname <- "tot"
organism <- "Sorghum bicolor"
chr <- 2
trial <- pvalue.combine(path, phenoname, chr, organism)
