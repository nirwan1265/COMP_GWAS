# FINDING COUNT AND COMMON GENES BETWEEN OMNI AND MAGMA
# SORGHUM
library(dplyr)
setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination")
sorghum_omni <- read.csv("OMNI_sorghum.csv")
sorghum_omni <- sorghum_omni[which(sorghum_omni$pvalue < 0.05), ]
sorghum_omni <- sorghum_omni[,-4]

sorghum_omni <- aggregate(. ~ GeneName, data = sorghum_omni, paste, collapse = ",")

# Sorting by the number of genes in a particular phenotype
sorghum_omni <- sorghum_omni %>% 
  arrange(desc(count(phenotype)))
sorghum_omni$phenotype_count <- sapply(strsplit(sorghum_omni$phenotype,","), length)

sorghum_omni <- sorghum_omni %>% 
  arrange(desc(phenotype_count))

df <- sorghum_omni[,c(1,5,6)]
