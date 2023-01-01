#Load genotype
geno_ld <- function(geno){
  # it should be in the same folder, but the file is too large to push to git hub 
  #geno <- vroom("Genotype/sorghum/allchrom_africa_filtered.MAF.txt")
  
  #use the above
  geno <- vroom("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/raw/africa.filtered/allchrom_africa_filtered.MAF.txt")
  return(geno)
}

geno <- geno_ld()
