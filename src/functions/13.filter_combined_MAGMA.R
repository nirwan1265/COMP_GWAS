#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################
# SELECTING COMMON TOP GENES
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



################################################################################
# SORGHUM
################################################################################
# Loading data
dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/Top_MAGMA_hits/sorghum/"
combined_magma_sorghum <- read.csv(paste0(dir,"Combined_MAGMA_res_sorghum.csv"))

# Selecting the highest number of phenotypes 
max_pheno_count <- max(combined_magma_sorghum$phenotype_count)

# Selecting all genes with a cut off of half of max_pheno_count
combined_magma_sorghum <- combined_magma_sorghum[which(combined_magma_sorghum$phenotype_count >= 8), ]


################################################################################
# MAIZE
################################################################################
# Loading data
dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/Top_MAGMA_hits/maize/"
combined_magma_maize <- read.csv(paste0(dir,"Combined_MAGMA_res_maize.csv"))

# Selecting the highest number of phenotypes 
max_pheno_count <- max(combined_magma_maize$phenotype_count)

# Selecting all genes with a cut off of half of max_pheno_count
combined_magma_maize <- combined_magma_maize[which(combined_magma_maize$phenotype_count >= 8), ]





