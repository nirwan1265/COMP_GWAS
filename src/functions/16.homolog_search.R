#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################
# HOMOLOGS SEARCH
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

################################################################################
# HOMOLOGS DATABASE
################################################################################
dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/sorghum_maize_homolog/"

homolog <- read.delim(paste0(dir,"Sorghum_Maize_Homolog_gene.txt"))
homolog <- homolog[,c("ensembl_gene_id1","ensembl_gene_id2")]


################################################################################
# FINDING ORTHOLOG between MAIZE and SORGHUM
################################################################################
dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/Top_MAGMA_hits/"

sorghum_top <- read.csv(paste0(dir,"sorghum/Top_MAGMA_hits_sorghum.csv")) %>% 
  filter(P_MULTI <= 0.05) %>% 
  dplyr::select(GENE) %>% 
  dplyr::rename(ensembl_gene_id1 = GENE) %>%
  left_join(homolog)

maize_top <- read.csv(paste0(dir,"maize/Top_MAGMA_hits_maize.csv")) %>% 
  filter(P_MULTI <= 0.05) %>% 
  dplyr::select(GENE) %>%
  dplyr::rename(ensembl_gene_id2 = GENE) %>%
  left_join(homolog)

# Common genes in sorghum from maize
maize_homologs_in_sorghum_gwas <- as.data.frame(intersect(sorghum_top$ensembl_gene_id2,maize_top$ensembl_gene_id2))
sorghum_homologs_in_maize_gwas <- as.data.frame(intersect(sorghum_top$ensembl_gene_id1,maize_top$ensembl_gene_id1))

write.csv(maize_homologs_in_sorghum_gwas, "MaizeOrthoForSorghumGWAS_from_MaizeGWAS.csv", row.names = F)
write.csv(sorghum_homologs_in_maize_gwas, "SorghumOrthoForMaizeGWAS_from_SorghumGWAS.csv", row.names = F)
system("pwd")

