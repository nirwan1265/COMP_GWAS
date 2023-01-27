# ACAT
library(ACAT)
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination")
combined_omni_sorghum <- read.csv("combined_sorghum_omni.csv")


# Subset some rows
subs_omni <- as.data.frame(combined_omni_sorghum[1,c(1,3,4)])
omni_res <- as.data.frame(combined_omni_sorghum[1,c(2)])
names(omni_res) <- "pvalue"
pheno_name <- as.data.frame(combined_omni_sorghum[1,c(5)])
names(pheno_name) <- "phenotype"

subs_omni$no_of_SNPs <- lapply(strsplit(subs_omni$no_of_SNPs, ","), as.numeric)
subs_omni$list_of_pvalue <- lapply(strsplit(subs_omni$list_of_pvalue , ","), as.numeric)
omni_res$pvalue <- lapply(strsplit(omni_res$pvalue , ","), as.numeric)
pheno_name$phenotype <- lapply(strsplit(pheno_name$phenotype , ","), as.character)



#Creating a new dataframe
results_acat <- data.frame(start=integer(), end=integer(), acat_result=numeric())
hello darkness my old friend


# looping through rows and applying acat function for each value in col1
for(i in 1:nrow(subs_omni)){
  start <- 1
  for(j in 1:length(subs_omni$no_of_SNPs[[i]])){
    end <- start + subs_omni$no_of_SNPs[[i]][j] - 1
    acat_result <- ACAT(subs_omni$list_of_pvalue[[i]][start:end], weights = NULL, is.check = TRUE)
    results_acat <- rbind(results_acat, data.frame(start, end, acat_result))
    start <- end + 1
  }
}



results_acat$results_omni <- NA
# loop through rows of values_df and assign values to new_col in results_df
counter <- 1
for(i in 1:nrow(omni_res)){
  for(j in 1:length(omni_res$pvalue[[i]])){
    results_acat$results_omni[counter] <- omni_res$pvalue[[i]][j]
    counter <- counter + 1
  }
}


results_acat$phenotype <- NA
# loop through rows of values_df and assign values to new_col in results_df
counter <- 1
for(i in 1:nrow(pheno_name)){
  for(j in 1:length(pheno_name$phenotype[[i]])){
    results_acat$phenotype[counter] <- pheno_name$phenotype[[i]][j]
    counter <- counter + 1
  }
}


results_acat <- results_acat[,-c(1,2)]

