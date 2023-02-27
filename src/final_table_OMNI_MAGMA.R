# FINDING COUNT AND COMMON GENES BETWEEN OMNI AND MAGMA
# SORGHUM
# OMNI
library(dplyr)
setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination")
#setwd("/Users/nirwan/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/pvalues.combination")
sorghum_omni <- read.csv("OMNI_sorghum.csv")
sorghum_omni <- sorghum_omni[which(sorghum_omni$pvalue < 0.05), ]
sorghum_omni <- sorghum_omni[,-4]

sorghum_omni <- aggregate(. ~ GeneName, data = sorghum_omni, paste, collapse = ",")

# Sorting by the number of genes in a particular phenotype
library(dplyr)
# sorghum_omni <- sorghum_omni %>% 
#   arrange(desc(count(phenotype)))
sorghum_omni$phenotype_count <- sapply(strsplit(sorghum_omni$phenotype,","), length)

sorghum_omni <- sorghum_omni %>% 
  arrange(desc(phenotype_count))


#MAGMA
library(dplyr)
system("ls")
# From combining_data_omni_magma
sorghum_MAGMA <- total_magma
#sorghum_MAGMA <- read.csv("MAGMA_sorghum.csv")
sorghum_MAGMA <- sorghum_MAGMA[which(sorghum_MAGMA$P_MULTI < 0.05), ]
sorghum_MAGMA <- sorghum_MAGMA[,c(-4,-5)]

sorghum_MAGMA <- aggregate(. ~ GENE, data = sorghum_MAGMA, paste, collapse = ",")

# Sorting by the number of genes in a particular phenotype
sorghum_MAGMA$phenotype_count <- sapply(strsplit(sorghum_MAGMA$phenotype,","), length)

sorghum_MAGMA <- sorghum_MAGMA %>% 
  arrange(desc(phenotype_count))

write.csv(sorghum_omni,"combined_sorghum_omni.csv",row.names = F)
write.csv(sorghum_MAGMA,"Combined_MAGMA_res_sorghum.csv",row.names = F)
### MAIZE ####
# FINDING COUNT AND COMMON GENES BETWEEN OMNI AND MAGMA
# maize
# OMNI
maize_omni <- read.csv("OMNI_maize.csv")
maize_omni <- maize_omni[which(maize_omni$pvalue < 0.05), ]
maize_omni <- maize_omni[,-4]

maize_omni <- aggregate(. ~ GeneName, data = maize_omni, paste, collapse = ",")

# Sorting by the number of genes in a particular phenotype
maize_omni <- maize_omni %>% 
  arrange(desc(count(phenotype)))
maize_omni$phenotype_count <- sapply(strsplit(maize_omni$phenotype,","), length)

maize_omni <- maize_omni %>% 
  arrange(desc(phenotype_count))


#MAGMA
maize_MAGMA <- read.csv("MAGMA_maize.csv")
maize_MAGMA <- maize_MAGMA[which(maize_MAGMA$P_MULTI < 0.05), ]
maize_MAGMA <- maize_MAGMA[,c(-4,-5)]

maize_MAGMA <- aggregate(. ~ GENE, data = maize_MAGMA, paste, collapse = ",")

# Sorting by the number of genes in a particular phenotype
maize_MAGMA$phenotype_count <- sapply(strsplit(maize_MAGMA$phenotype,","), length)

maize_MAGMA <- maize_MAGMA %>% 
  arrange(desc(phenotype_count))


