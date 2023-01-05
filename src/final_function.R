#Required arguments
path = "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data"
phenoname <- "tot"
organism <- "Sorghum bicolor"
chr <- 3

#Register nodes
cluster <- makeCluster(parallel::detectCores() -1)
registerDoParallel(cluster)


#Final Function
final_function <- function(path, phenoname, chr, organism) {
  
  
  results <- gbj_test(chr)
  
  return(results)
  
}


#Getting the results with time
tic()
final_results <- final_function(path, phenoname, chr, organism)
toc()

# Stop the parallel cluster
stopCluster(cluster)
