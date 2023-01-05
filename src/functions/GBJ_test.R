gbj_test <- function(chr){
  #Empty results list to store values for each chromosome
  results <- list()
  
  for(j in 1:chr){
    
    # Subset the required data
    zstat_df <- as.data.frame(trial$preprocess$Zstat[[j]]) # i
    zstat_df <- zstat_df[1:100,1:100]
    pvalue_df <- as.data.frame(trial$preprocess$pvalue[[j]])
    pvalue_df <- pvalue_df[1:100,1:100]
    marker_df <- as.data.frame(trial$preprocess$Marker[[j]])
    marker_df <- marker_df[1:100,1:100]
    genotype_df <- as.data.frame(trial$genotype[[j]])
    
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
    ref_genotype <- foreach(i = 1:ncol(marker_df), .combine = c) %dopar% {
      list(sub_refgeno(marker_df[,i], genotype_df))
    }
    
    # PCA 
    pca <- as.data.frame(trial$PCA)
    
    # Use foreach to apply the function to each column in parallel for calculating correlation matrix
    corr_mat = list()
    corr_mat <- foreach(i = 1:20, .combine = c) %dopar% {
      list(GBJ::estimate_ss_cor(ref_pcs = pca, ref_genotypes = as.data.frame(ref_genotype[i]), link_function = 'linear'))
    }
    
    # gbj_analysis = list()
    # gbj_analysis <- foreach(i = 1:20, .combine = c) %dopar% {
    #   list(GBJ::GBJ(test_stats = as.vector(unlist(na.omit(zstat_df[i]))), cor_mat=corr_mat[[i]])$GBJ_pvalue)
    # }
    # results[[j]] <- gbj_analysis
    # names(results)[[j]] <- paste0("chr",j)
    
    omni_analysis = list()
    omni_analysis <- foreach(i = 1:20, .combine = c) %dopar% {
      list(GBJ::OMNI_ss(test_stats = as.vector(unlist(na.omit(zstat_df[i]))), cor_mat=corr_mat[[i]], num_boots = 100)$OMNI_pvalue)
    }
    results[[j]] <- omni_analysis
    names(results)[[j]] <- paste0("chr",j)
  }
  return(results)
}