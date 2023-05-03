# Combining data

################################################################################
# SORGHUM
################################################################################

# OMNIBUS
omni_path <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GBJ/sorghum/GBJ_sorghum_2kb_v3_RDS"
omni_files <- list.files(omni_path, pattern = "RDS")
omni_files <- c("PNZ1_omni_sorghum_2kbext.RDS","POL1_omni_sorghum_2kbext.RDS","PBR1_omni_sorghum_2kbext.RDS","stp10_omni_sorghum_2kbext.RDS")

for(file in omni_files){
  file_path <- file.path(omni_path, file)
  #assign(gsub(".RDS","",paste0(file)), readRDS(file_path))
  assign(paste0(file), readRDS(file_path))
}

omni_total <- ls()
total_omni <- data.frame()
for(all_omni in omni_files){
  omni <- get(all_omni)
  if(!is.data.frame(omni)){
    omni <- as.data.frame(omni)
  }
  omni$phenotype <- all_omni
  assign(all_omni,omni)
  total_omni <- rbind(total_omni, omni)
  #total_omni <- omni_total[-1,]
}
#total_omni[,1:ncol(total_omni)] <- apply(total_omni,2,function(x) gsub("2kbext.RDS","",x))
total_omni[,1:ncol(total_omni)] <- apply(total_omni,2,function(x) gsub("_","",x))
total_omni[,1:ncol(total_omni)] <- apply(total_omni,2,function(x) gsub("omnisorghum","",x))
total_omni[,1:ncol(total_omni)] <- apply(total_omni,2,function(x) gsub("kbext.RDS","",x))

total_omni <- total_omni %>% mutate(pvalue = as.numeric(as.character(pvalue))) %>% arrange(pvalue)

write.csv(total_omni,"Top_GBJ_hits_sorghum_4Pheno.csv",row.names = FALSE)

if("total_omni" %in% ls()) {
  # remove all variables except for "total_omni"
  rm(list = setdiff(ls(), "total_omni"))
  
}



# MAGMA
magma_path <- ("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/MAGMA/sorghum/MAGMA_sorghum_2kb_v3_csv")
magma_files <- list.files(magma_path, pattern = ".genes.out")
magma_files <- c("PBR1_magma_multi_snpwise.genes.out","PNZ1_magma_multi_snpwise.genes.out","stp10_magma_multi_snpwise.genes.out","POL1_magma_multi_snpwise.genes.out")

for(file in magma_files){
  file_path <- file.path(magma_path, file)
  #assign(gsub(".RDS","",paste0(file)), readRDS(file_path))
  assign(paste0(file), read.table(file_path))
}

magma_total <- ls()
total_magma <- data.frame()
for(all_magma in magma_files){
  magma <- get(all_magma)
  if(!is.data.frame(magma)){
    magma <- as.data.frame(magma)
  }
  magma$phenotype <- all_magma
  assign(all_magma,magma)
  total_magma <- rbind(total_magma, magma)
}
names(total_magma) <- total_magma[1,]
total_magma <- total_magma[-1,]
total_magma <- total_magma[,c(1,5,9:12)]
names(total_magma)[6] <- "phenotype"
total_magma[,1:ncol(total_magma)] <- apply(total_magma,2,function(x) gsub("_magma_multi_snpwise.genes.out","",x))
total_magma <- total_magma %>% arrange(P_MULTI)

write.csv(total_magma,"Top_MAGMA_hits_sorghum_4pheno.csv",row.names = FALSE)


################################################################################
# MAIZE
################################################################################

# OMNIBUS
omni_path <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GBJ/maize/GBJ_maize_2kb_v5_RDS"
omni_files <- list.files(omni_path, pattern = "RDS")
omni_files <- c("POL1_gbj_maize.RDS","PBR1_gbj_maize.RDS","PNZ1_gbj_maize.RDS","stp10_gbj_maize.RDS")

for(file in omni_files){
  file_path <- file.path(omni_path, file)
  #assign(gsub(".RDS","",paste0(file)), readRDS(file_path))
  assign(paste0(file), readRDS(file_path))
}


omni_total <- ls()
total_omni <- data.frame()
for(all_omni in omni_files){
  omni <- get(all_omni)
  if(!is.data.frame(omni)){
    omni <- as.data.frame(omni)
  }
  omni$phenotype <- all_omni
  assign(all_omni,omni)
  total_omni <- rbind(total_omni, omni)
  #total_omni <- omni_total[-1,]
}
#total_omni[,1:ncol(total_omni)] <- apply(total_omni,2,function(x) gsub("_omni_maize.RDS","",x))
total_omni[,1:ncol(total_omni)] <- apply(total_omni,2,function(x) gsub("_gbj_maize.RDS","",x))

total_omni <- total_omni %>% mutate(pvalue = as.numeric(as.character(pvalue))) %>% arrange(pvalue)

write.csv(total_omni,"Top_GBJ_hits_maize_4Pheno.csv",row.names = FALSE)



# MAGMA
magma_path <- ("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/MAGMA/maize_2kb")
magma_files <- list.files(magma_path, pattern = ".genes.out")
magma_files <- c("PBR1_maize_magma_multi_snpwise.genes.out","PNZ1_maize_magma_multi_snpwise.genes.out","stp10_maize_magma_multi_snpwise.genes.out","POL1_maize_magma_multi_snpwise.genes.out")

for(file in magma_files){
  file_path <- file.path(magma_path, file)
  #assign(gsub(".RDS","",paste0(file)), readRDS(file_path))
  assign(paste0(file), read.table(file_path))
}

magma_total <- ls()
total_magma <- data.frame()
for(all_magma in magma_files){
  magma <- get(all_magma)
  if(!is.data.frame(magma)){
    magma <- as.data.frame(magma)
  }
  magma$phenotype <- all_magma
  assign(all_magma,magma)
  total_magma <- rbind(total_magma, magma)
}
names(total_magma) <- total_magma[1,]
total_magma <- total_magma[-1,]
total_magma <- total_magma[,c(1,5,9:12)]
names(total_magma)[6] <- "phenotype"
total_magma[,1:ncol(total_magma)] <- apply(total_magma,2,function(x) gsub("_maize_magma_multi_snpwise.genes.out","",x))
total_magma <- total_magma %>% arrange(P_MULTI)

write.csv(total_magma,"Top_MAGMA_hits_maize_4Pheno.csv",row.names = FALSE)






