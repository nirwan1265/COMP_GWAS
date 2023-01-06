
res <- NULL
for(i in 1:chr){
  x <- as.data.frame(t(as.data.frame(final_results[[i]], header = FALSE)))
  res <- rbind(res,x)
  #res <- res %>% dplyr::select("names","V1") %>% arrange(V1)# %>% mutate(no_of_snps = lofsnps("GeneName"))
  x <- NULL
}
res$names <- colnames(zstat_df[1:40])
res <- res %>% dplyr::select("names","V1") #%>% arrange(V1) %>% rename(GeneName = names, pvalue = V1) %>% mutate(no_of_snps = lofsnps("GeneName"))
names(res) <- c("GeneName","pvalue")

#Adding number of SNPs
no_of_snps <- NULL
for(i in 1:ncol(zstat_df)){
  no_of_snps <- c(no_of_snps,length(na.omit(zstat_df[,i])))
}
res$no_of_SNPs <- no_of_snps[1:40]

#Adding whether the had pvalue combination
res$pval_combination_GBJ_minP_GHC_SKAT <- "Yes"

#Empyting row name because ewww
row.names(res) <- NULL


#Combining tables with 1 and many SNPs
res <- rbind(res, single_snp) %>% arrange(desc(pvalue))


# Subset the required data
pvalue_df <- as.data.frame(trial$preprocess$pvalue[[1]])
#pvalue_df <- pvalue_df[1:100,1:100]
marker_df <- as.data.frame(trial$preprocess$Marker[[1]])


#Subsetting genes with only one column
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

typeof


