######################
# Required packages
######################
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


######################
# Phenotype data
# Separated in columns
######################
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
pheno <- read.csv("Maize_allphospho.csv")


######################
# SNPs PCA eigen values
# Separated in columns
######################
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
pca <- read.csv("pca.csv")


######################
# Combining pheno and PCA
# Separated in columns
######################
data_tbl <- as_tibble(cbind(as.data.frame(pca[,2]), pheno[,2:29])) # Took the first PCA
names(data_tbl)[1] <- "pca"

# for weighted
# Giving weights to the Eigen Vector:
#W = matrix(c(3,2,1,1,1,1,1,1,1,1), ncol = 1)
#composite_pca <- as.data.frame(rowSums(as.matrix(composite_pca) %*% W))


###################
#XGBOOST
###################
#Model
fit_xgboost <- boost_tree(learn_rate = 0.3) %>%
  set_mode("regression") %>%
  set_engine("xgboost") %>%
  fit(pca ~. , data = data_tbl)

fit_xgboost


# Explainer
explainer <- DALEX::explain(
  model = fit_xgboost,
  data = data_tbl,
  y = data_tbl$pca,
  label = "XGBoost"
)

###################
# MODEL STUDIO
# Visualization 
###################
modelStudio::modelStudio(explainer)

