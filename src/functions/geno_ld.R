#Load genotype

#Sorghum
# geno_ld <- function(chr){
#   geno <- list()
#   # it should be in the same folder, but the file is too large to push to git hub 
#   #geno <- vroom("Genotype/sorghum/allchrom_africa_filtered.MAF.txt")
#   
#   #use the above
#   for(i in 1:chr){
#     assign(paste0("geno",i), vroom(paste0("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/raw/africa.filtered/c",i,"_MAF_sorghum.txt")))
#     geno <- c(geno, list(get(paste0("geno",i))))
#     names(geno)[i] <- paste0("geno",i)
#   }
#   return(geno)
# }

#Maize
geno_ld <- function(chr){
  geno <- list()
  # it should be in the same folder, but the file is too large to push to git hub 
  #geno <- vroom("Genotype/sorghum/allchrom_africa_filtered.MAF.txt")
  
  #use the above
  for(i in 1:chr){
    assign(paste0("geno",i), vroom(paste0("/Users/nirwantandukar/Documents/RubenLab/Maize/perchr/numerical_MAF_pheno_filtered/geno",i,".txt")))
    geno <- c(geno, list(get(paste0("geno",i))))
    names(geno)[i] <- paste0("geno",i)
  }
  return(geno)
}

tic()
geno <- geno_ld(1)
toc()

