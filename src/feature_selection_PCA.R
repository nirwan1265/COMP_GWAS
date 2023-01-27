library(pcaL1)
library(vroom)
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
pheno <- read.csv("Maize_allphospho.csv")
help(pcaL1)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/geno/RomeroNavarro2017/MAF/")
snp_data <- vroom("c1_MAF_v5_sorghum.txt")
snp_data[1:4,1:4]
snp_data <- snp_data[,sample(colnames(snp_data), 500)]

##for 100x10 data matrix X, 
## lying (mostly) in the subspace defined by the first 2 unit vectors, 
## projects data into 1 dimension.
#X <- matrix(c(runif(100*2, -10, 10), rep(0,100*8)),nrow=100) + 
#  matrix(c(rep(0,100*2),rnorm(100*8,0,0.1)),ncol=10)
myl1pca <- l1pca(snp_data)

pca_scores <- myl1pca$scores

##projects data into 2 dimensions.
myl1pca <- l1pca(snp_data, projDim=2, center=FALSE, 
                 tolerance=0.00001, iterations=20)

## plot first two scores
plot(myl1pca$scores)




# First, install and load the pls package
install.packages("pls")
#install.packages("pls", dependencies = TRUE, repos = "http://cran.us.r-project.org")
library(pls)
library(caret)


# Next, load your phenotype data into R. In this example, the data is stored in a data frame called "pheno_data"
dim(pca_scores)
dim(pheno)
# Fit a PLS model using the PCs as the independent variables and the phenotypes as the dependent variables
pls_model <- plsr(pca_scores~as.matrix(pheno[,2:29]), ncomp = 1)

summary(pls_model)



# Print the variable importance in the projection (VIP) scores for the phenotypes
imp <- as.data.frame(varImp(pls_model))

imp_df <- as.data.frame(imp)
imp_df$pheno <- rownames(imp)
imp_sorted <- imp_df[order(-imp_df$Overall),]
top_n_pheno <- imp_sorted[1:28,]
imp_sorted <- imp[order(-imp$Overall),]
imp_sorted <- cbind(imp_sorted, rownames(imp))
imp_sorted <- imp_sorted[1:28,]

