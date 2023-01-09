# Load PCA
pca_ld <- function(pca){
  pca <- vroom("PCA/sorghum/sorghum_PCA.txt") 
  pca <- pca[,-1]
}

#pca <- pca_ld()

