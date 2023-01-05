library(tictoc)
#Register nodes
cluster <- makeCluster(parallel::detectCores() -1)
registerDoParallel(cluster)

#Run trial files
pvalue.combine <- function(path, phenoname, chr, organism){
  preprocess_data <- preprocess(path, phenoname, chr,  organism)
  geno <- geno <- geno_ld(chr)
  pca <- pca_ld()
  return(list(preprocess = preprocess_data,genotype = geno,PCA = pca))
}



path = "/Users/nirwan/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data"
phenoname <- "tot"
organism <- "Sorghum bicolor"
chr <- 2
trial <- pvalue.combine(path, phenoname, chr, organism)





#GBJ function
gbj_test <- function(gwas.zstat, gwas.marker, gwas.pvalue, geno, tab.pc){
  x <- as.data.frame(matrix(0, nrow = 1, ncol = 1))
  y <- vector()
  z <- vector()
  combined.test.statistics <- as.data.frame(matrix(NA, nrow = 1, ncol = 1))
  ref_genotype <- as.data.frame(matrix(NA, nrow = 1, ncol = 1))
  ref_genotype_skat <- as.data.frame(matrix(NA, nrow = 1, ncol = 1))
  for (i in 1:chr) {
    zstat_df <- as.data.frame(trial$preprocess$Zstat[[i]]) # i
    zstat_df <- zstat_df[1:5,1:5]
    pvalue_df <- as.data.frame(trial$preprocess$pvalue[[i]])
    marker_df <- as.data.frame(trial$preprocess$Marker[[i]])
    genontype_df <- as.data.frame(trial$genotype[[i]])
    
    #ref_genotype <- as.data.frame(genontype_df[,colnames(genontype_df)] %in% marker_df)
  }
  return(as.data.frame(zstat_df))
}

for(i in chr)

zstat_df <- as.data.frame(trial$preprocess$Zstat[[1]]) # i
zstat_df <- zstat_df[1:100,1:100]
pvalue_df <- as.data.frame(trial$preprocess$pvalue[[1]])
pvalue_df <- pvalue_df[1:100,1:100]
marker_df <- as.data.frame(trial$preprocess$Marker[[1]])
marker_df <- marker_df[1:100,1:100]
genotype_df <- as.data.frame(trial$genotype[[1]])



## Function for Subsetting columns with more than one element 
subset_element <- function (x) length(na.omit(x)) > 1

# Applying the function for each df
subset_zstat <- sapply(zstat_df, subset_element)
subset_pvalue <- sapply(pvalue_df, subset_element)
subset_marker <- sapply(marker_df, subset_element)

# Subsetting the data frame with more than 1 elements
zstat_df <- zstat_df[, subset_zstat]
pvalue_df <- pvalue_df[, subset_pvalue]
marker_df <- marker_df[, subset_marker]


# Subsetting reference genotype
sub_refgeno <- function(y,z){
  y <- y[!is.na(y)]
  ref_genotype <- as.data.frame(z[y])
  return(ref_genotype)
}

# Use foreach to apply the function to each column in parallel to subset ref_genotype for each column
ref_genotype = list()
tic()
ref_genotype <- foreach(i = 1:ncol(marker_df), .combine = c) %dopar% {
  list(sub_refgeno(marker_df[,i], genotype_df))
}
toc()



# PCA 
pca <- as.matrix(trial$PCA)


# Use foreach to apply the function to each column in parallel for calculating correlation matrix
corr_mat = list()
tic()
corr_mat <- foreach(i = 1:20, .combine = c) %dopar% {
  #corr_gbj(ref_genotype[i])
  list(estimate_ss_cor(ref_pcs = pca, ref_genotypes = as.data.frame(ref_genotype[i]), link_function = 'linear'))
}
toc()



# Use foreach to apply the function to each column in parallel for calculating correlation matrix
gbj_analysis = list()
tic()
gbj_analysis <- foreach(i = 1:20, .combine = c) %dopar% {
  #corr_gbj(ref_genotype[i])
  list(GBJ(test_stats = as.vector(unlist(na.omit(zstat_df[i]))), cor_mat=corr_mat[[i]])$GBJ_pvalue)
}
toc()

# Use foreach to apply the function to each column in parallel for calculating correlation matrix
omni_analysis = list()
tic()
omni_analysis <- foreach(i = 1:20, .combine = c) %dopar% {
  list(GBJ::OMNI_ss(test_stats = as.vector(unlist(na.omit(zstat_df[i]))), cor_mat=corr_mat[[i]], num_boots = 100)$OMNI_pvalue)
}
toc()


# Stop the parallel cluster
stopCluster(cluster)
