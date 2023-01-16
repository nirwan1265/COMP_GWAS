gbj_test <- function(path, phenoname, chr, organism){
  #Register nodes
  #cluster <- makeCluster(parallel::detectCores() - 1)
  #registerDoParallel(cluster)
  
  #Empty helper variables
  results <- list()
  res <- NULL
  colnames_df <- NULL
  no_of_snps <- NULL
  no_of_markers <- NULL
  no_of_pvalue <- NULL
  pvalue_df_list <- NULL
  marker_df_list <- NULL
  tot_single_snps <- NULL
  
  # Loops indices
  i <- 1
  k <- 1
  
  # Preprocess data
  trial <- data_wrangle(path, phenoname, chr, organism)
  
  # Looping through the chromosomes
  for(j in 1:chr){
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
    single_snp$list_of_pvalue <- single_snp$pvalue
    row.names(single_snp) <- NULL
    tot_single_snps <- rbind(tot_single_snps,single_snp)
    
    ############################################################################s
    
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
    # gbj_analysis = list()
    # gbj_analysis <- foreach(i = 1:ncol(marker_df), .combine = c) %dopar% {
    #   list(GBJ::GBJ(test_stats = as.vector(unlist(na.omit(zstat_df[i]))),
    #                 cor_mat=corr_mat[[i]])$GBJ_pvalue)
    # }
    # results[[j]] <- gbj_analysis
    # names(results)[[j]] <- paste0("chr",j)
    
    ############################################################################
    
    #Parallelizing OMNI analysis
    omni_analysis = list()
    omni_analysis <- foreach(i = 1:ncol(marker_df), .combine = c) %dopar% {
      list(GBJ::OMNI_ss(test_stats = as.vector(unlist(na.omit(zstat_df[i]))), cor_mat=corr_mat[[i]], num_boots = 100)$OMNI_pvalue)
    }
    results[[j]] <- omni_analysis
    names(results)[[j]] <- paste0("chr",j)
    
    
    ############################################################################
    
    # Combining the results
    x <- as.data.frame(t(as.data.frame(results[[j]], header = FALSE)))
    res <- rbind(res,x)
    
    ############################################################################
    
    # Getting the number ot SNPs and list of pvalues
    colnames_df <- c(colnames_df,colnames(zstat_df[1:ncol(zstat_df)]))
    for(i in 1:ncol(zstat_df)){
      no_of_snps <- c(no_of_snps,length(na.omit(zstat_df[,i])))
    }
    no_of_pvalue <- c(no_of_pvalue, apply(t(as.data.frame(pvalue_df)), 1, function(row) paste(na.omit(row), collapse = ",")))
    
    ############################################################################
    
  }
  ##############################################################################
  
  # Data wrangling
  res$names <- colnames_df
  res <- res %>% dplyr::select("names","V1")
  names(res) <- c("GeneName","pvalue")

  # Adding the total number of SNPs
  res$no_of_SNPs <- no_of_snps
  
  # Combination or not
  res$pval_combination_GBJ_minP_GHC_SKAT <- "Yes"
  
  # Adding number of pvalue
  res$list_of_pvalue <- no_of_pvalue
  
  #Empyting row name because ewww
  row.names(res) <- NULL
  
  # Combining tables with 1 and many SNPs- FINAL result
  res <- rbind(res, tot_single_snps) %>% dplyr::mutate(pvalue = as.double(pvalue))  %>% dplyr::arrange(pvalue) #%>% dplyr::mutate(pvalue = base::format(pvalue, scientific = TRUE, digits = 5))
  
  ##############################################################################
  
  # Returning the final results
  return(res)
  # Stop the parallel cluster
  #stopCluster(cluster)
  
}


# #Required arguments
# path = "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data"
# phenoname <- "tot"
# organism <- "Sorghum bicolor"
# chr <- 2
# 
# #Register nodes
# cluster <- makeCluster(parallel::detectCores() - 1)
# registerDoParallel(cluster)
# 
# #Getting the results with time
# tic()
# final_results <- gbj_test(path, phenoname, chr, organism)
# toc()
# 
# # Stop the parallel cluster
# stopCluster(cluster)



################################################################################
# for multiple phenotypes
################################################################################

# #Required arguments

path = "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data"
#Sorghum
#phenoname <- c("NPlim","occ","org","PBR1","PBR2","PHO1","PHO2","PMEH1","PMEH2","PNZ1","PNZ2","POL1","POL2","sec","solHi_","solLo_","solMo_","solmod_","solVL_","stp1", "stp2", "stp3","stp10_","tot","TP1","TP2")
#Maize
phenoname <- c("NPlim","occ","org","PBR1","PBR2","PHO1","PHO2","PMEH1","PMEH2","PNZ1","PNZ2","POL1","POL2","sec","sol_Hi","sol_Lo","sol_Mo","sol","sol_VL","stp10", "stp20", "stp30","stp100","tot","TP1","TP2")
phenoname <- c("NPlim","occ","org","PBR1","PBR2","PHO1","PHO2","PMEH1","PMEH2","PNZ1","PNZ2","POL1","POL2","sec","sol_Hi","sol_Lo","sol_VL","sol")#,"stp10","stp20","stp30","stp100", "tot","TP1","TP2")
phenoname <- c("tot","TP1","TP2","NPlim","occ","org","PBR1","PBR2","PHO1","PHO2","PMEH1","PMEH2","PNZ1","PNZ2","POL1","POL2","sec","sol_Hi","sol_Lo","sol_VL","sol")
#phenoname <- c("NPlim","occ","org","PBR1","PBR2","PHO1","PHO2","PMEH1","PMEH2","PNZ1","PNZ2","POL1","POL2","sec","sol_Hi","sol_Lo",)
#phenoname <- "tot"
organism <- "Zea"
chr <- 10

#Register nodes
cluster <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cluster)

#Getting the results with time
tic()
for (m in phenoname){
  print(m)
  assign(paste0(m,"_omni_maize"), gbj_test(path,m,chr,organism))
}
toc()

system("pwd")
# Save as CSV
for(i in phenoname){
  write.csv(get(paste0(i,"_omni_maize")),(paste0(i,"_omni_maize.csv")), row.names = FALSE)
}

#Save as RDS
for(i in phenoname){
  saveRDS(get(paste0(i,"_omni_maize")),(paste0(i,"omni_maize.RDS")))
}

# Stop the parallel cluster
stopCluster(cluster)
system("ls")
