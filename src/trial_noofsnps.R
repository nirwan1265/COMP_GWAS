for(i in 1:chr){
  
  
}
res <- NULL
for(i in 1:chr){
  assign("res",)
}
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
res$pval_combination_GBJ_minP_GHC_SKAT <- "Yes"
res$no_of_snps <- length(na.omit(colnames(zstat_df[1])))


row.names(res) <- res[,1]
res <- res %>% mutate(no_of_snps = lofsnps(GeneName))



row.names(res) <- NULL
row.names(res)["SORBI_3001G007000"]



lofsnps <- function(x){
  len <- length(na.omit(zstat_df[x]))
  return(len)
}
zstat_df["SORBI_3001G000800"]


x <- res %>% mutate(no_of_snps = lofsnps(GeneName))

length(na.omit(zstat_df$SORBI_3001G000800))

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
single_snp$pval_combination_GBJ_minP_GHC_SKAT <- "No"
single_snp$no_of_SNPs <- 1
row.names(single_snp) <- NULL

x <- rbind(res,single_snp) %>% arrange(pvalue)
