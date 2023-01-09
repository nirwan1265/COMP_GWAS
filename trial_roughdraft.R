#Required arguments
path = "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data"
phenoname <- "tot"
organism <- "Sorghum bicolor"
chr <- 2
res <- NULL
#trial <- data_wrangle(path, phenoname, chr, organism)
results = list()
colnames_df = NULL
x <- NULL
no_of_snps <- NULL
for(j in 1:2){
  
  # Subset the required data
  
  zstat_df <- as.data.frame(trial$preprocess$Zstat[[j]]) # i
  pvalue_df <- as.data.frame(trial$preprocess$pvalue[[j]])
  marker_df <- as.data.frame(trial$preprocess$Marker[[j]])
  
  ############################################################################
  
  # Subset the genotype
  
  genotype_df <- as.data.frame(trial$genotype[[j]])
  
  ############################################################################
  
  # Subset genes with only one column
  subset_one <- function (x) length(na.omit(x)) == 1
  
  # Applying the function for pvalue
  subset_pvalue <- sapply(pvalue_df, subset_one)
  
  # Subsetting the data frame with more than 1 elements
  single_snp <- NULL
  single_snp <- as.data.frame(pvalue_df[, subset_pvalue])
  single_snp <- as.data.frame(t(rbind(colnames(single_snp),single_snp[1,])))
  names(single_snp) <- c("GeneName","pvalue")
  single_snp$no_of_SNPs <- 1
  single_snp$pval_combination_GBJ_minP_GHC_SKAT <- "No"
  row.names(single_snp) <- NULL
  
  #############################################################################
  
  ## Function for Sub-setting columns with more than one element 
  
  subset_element <- function (x) length(na.omit(x)) > 1
  
  # Applying the function for each df
  
  subset_zstat <- sapply(zstat_df, subset_element)
  subset_pvalue <- sapply(pvalue_df, subset_element)
  subset_marker <- sapply(marker_df, subset_element)
  
  # Subsetting the data frame with more than 1 elements
  
  zstat_df <- zstat_df[, subset_zstat]
  pvalue_df <- pvalue_df[, subset_pvalue]
  marker_df <- marker_df[, subset_marker]
  
  ############################################################################
  
  # Subsetting reference genotype
  sub_refgeno <- function(y,z){
    y <- y[!is.na(y)]
    ref_genotype <- as.data.frame(z[y])
    return(ref_genotype)
  }
  
  # Parallizing each column to subset ref_genotype 
  
  ref_genotype = list()
  ref_genotype <- foreach(i = 1:ncol(marker_df), .combine = c) %dopar% {
    list(sub_refgeno(marker_df[,i], genotype_df))
  }
  
  ############################################################################
  
  # PCA 
  
  pca <- as.data.frame(trial$PCA)
  
  ############################################################################
  
  # Paralellizing for calculating correlation matrix
  
  corr_mat = list()
  corr_mat <- foreach(i = 1:ncol(marker_df), .combine = c) %dopar% {
    list(GBJ::estimate_ss_cor(ref_pcs = pca, 
                              ref_genotypes = as.data.frame(ref_genotype[i]), 
                              link_function = 'linear'))
  }
  
  ############################################################################
  
  # Parallelizing GBJ analysis
  
  gbj_analysis = list()
  gbj_analysis <- foreach(i = 1:ncol(marker_df), .combine = c) %dopar% {
    list(GBJ::GBJ(test_stats = as.vector(unlist(na.omit(zstat_df[i]))), 
                  cor_mat=corr_mat[[i]])$GBJ_pvalue)
  }
  
  results[[j]] <- gbj_analysis
  names(results)[[j]] <- paste0("chr",j)
  
  ############################################################################
  
  x <- as.data.frame(t(as.data.frame(results[[j]], header = FALSE)))
  res <- rbind(res,x)
  colnames_df <- c(colnames_df,colnames(zstat_df[1:ncol(zstat_df)]))
  for(i in 1:ncol(zstat_df)){
    no_of_snps <- c(no_of_snps,length(na.omit(zstat_df[,i])))
  }
}

res$names <- colnames_df
res <- res %>% dplyr::select("names","V1") 
names(res) <- c("GeneName","pvalue")



# omni_analysis = list()
# omni_analysis <- foreach(i = 1:20, .combine = c) %dopar% {
#   list(GBJ::OMNI_ss(test_stats = as.vector(unlist(na.omit(zstat_df[i]))), cor_mat=corr_mat[[i]], num_boots = 100)$OMNI_pvalue)
# }
# results[[j]] <- omni_analysis
# names(results)[[j]] <- paste0("chr",j)

# for(j in 1:length(zstat_df)){
#  names(results[[i]][[j]]) <- colnames(zstat_df)[j]
# }
