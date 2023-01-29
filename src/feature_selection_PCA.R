library(pls)
library(caret)
library(pcaL1)
library(vroom)
library(glmnet)
library(modelStudio)
library(DALEX)
library(tidyverse)
library(tidymodels)
library(xgboost)



setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
system("ls")
pheno <- read.csv("Maize_allphospho.csv")


#setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/geno/RomeroNavarro2017/MAF/")
#snp_data <- vroom("c1_MAF_v5_sorghum.txt")
#snp_data[1:4,1:4]
#snp_data <- snp_data[,sample(colnames(snp_data), 500)]

##for 100x10 data matrix X, 
## lying (mostly) in the subspace defined by the first 2 unit vectors, 
## projects data into 1 dimension.
#X <- matrix(c(runif(100*2, -10, 10), rep(0,100*8)),nrow=100) + 
#  matrix(c(rep(0,100*2),rnorm(100*8,0,0.1)),ncol=10)
#myl1pca <- l1pca(snp_data)


#pca_scores <- myl1pca$scores

##projects data into 2 dimensions.
#myl1pca <- l1pca(snp_data, projDim=2, center=FALSE, 
#                 tolerance=0.00001, iterations=20)

## plot first two scores
#plot(myl1pca$scores)




# First, install and load the pls package
#install.packages("pls")
#install.packages("pls", dependencies = TRUE, repos = "http://cran.us.r-project.org")


# Next, load your phenotype data into R. In this example, the data is stored in a data frame called "pheno_data"
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data/PCA")
pca <- read.table("maize/maize_PCA.txt",header = T)
dim(pca)
#pca <- pca[,c(-1,-3:-11)]
# Giving weights to the Eigen Vector:
W = matrix(c(3,2,1,1,1,1,1,1,1,1), ncol = 1)

dim(pca)
dim(W)


#Combining Eigen Vector with weight 
composite_pca <- pca
composite_pca <- as.data.frame(rowSums(as.matrix(pca) %*% W))


dim(pca)
dim(pheno)

# Fit a PLS model using the PCs as the independent variables and the phenotypes as the dependent variables
pls_model <- plsr(as.matrix(composite_pca)~as.matrix(pheno[,2:29]), ncomp = 1)

summary(pls_model)

# Print the variable importance in the projection (VIP) scores for the phenotypes
imp <- as.data.frame(varImp(pls_model))
imp_df <- as.data.frame(imp)
imp_df$pheno <- rownames(imp)
imp_sorted <- imp_df[order(-imp_df$Overall),]
rownames(imp_sorted) <- NULL
imp_sorted$pheno <- sapply(strsplit(imp_sorted$pheno, ")"), "[",2)



##############
#XGBOOST
##############
data_tbl <- as_tibble(cbind(as.data.frame(pca[,2]), pheno[,2:29]))
names(data_tbl)[1] <- "pca"


#Model
fit_xgboost <- boost_tree(learn_rate = 0.3) %>%
  set_mode("regression") %>%
  set_engine("xgboost") %>%
  fit(pca ~. , data = data_tbl)
boost_
fit_xgboost


# Explainer
explainer <- DALEX::explain(
  model = fit_xgboost,
  data = data_tbl,
  y = data_tbl$pca,
  label = "XGBoost"
)


# MODEL STUDIO

modelStudio::modelStudio(explainer)
