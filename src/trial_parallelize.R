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



path = "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data"
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

trial_add <- function(x){
  return(x+100)
}

trial_ref <- function(y,z){
  y <- y[!is.na(y)]
  #z <- z[!is.na(z)]
  ref_genotype <- as.data.frame(z[y])
  return(ref_genotype)
}


pca <- trial$PCA
trial_gbj <- function(ref_genotype){
  cor_mat <- estimate_ss_cor(ref_pcs = pca, ref_genotypes = ref_genotype, link_function = 'linear')
  return(cor_mat)
}

## NEED TO run ref_genotype individually first 
## Use that to run run GBJ


# Use foreach to apply the function to each column in parallel
#results <- foreach(i = 1:ncol(as.data.frame(trial$preprocess$Zstat[1])), .combine = c) %dopar% {
results = list()
tic()
results <- foreach(i = 1, .combine = c) %dopar% {
  list(trial_ref(marker_df[,i], genotype_df),ref_genotype = trial_add(zstat_df[,i]), trial_gbj(ref_genotype))
}
toc()
results

# Stop the parallel cluster
stopCluster(cluster)
