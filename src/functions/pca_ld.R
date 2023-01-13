# Load PCA
pca_ld <- function(organism){
  if (organism == "Sorghum"){
    pca <- vroom("PCA/sorghum/sorghum_PCA.txt") 
    pca <- pca[,-1]  
  }
  else{
    pca <- vroom("PCA/maize/maize_PCA.txt") 
    pca <- pca[,-1]  
  }
  return(pca)
}

#pca <- pca_ld()

