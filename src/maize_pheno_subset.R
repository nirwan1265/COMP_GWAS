setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
maize_pheno <- read.csv("taxa_geoloc_pheno.csv")
maize_pheno <- maize_pheno[which(maize_pheno$sp == "Zea mays"), ]
maize_pheno <- maize_pheno[which(maize_pheno$GEO3 == "Meso-America" | maize_pheno$GEO3 == "South America"| maize_pheno$GEO3 == "Caribbean"), ]

maize_pheno$GEO3 <- ifelse(maize_pheno$GEO3 == "Meso-America", 1,
                      ifelse(maize_pheno$GEO3 == "South America", 2,
                             ifelse(maize_pheno$GEO3 == "Caribbean", 3, maize_pheno$GEO3)))
maize_pheno <- maize_pheno[,-c(1,3,5,6,7)]
maize_pheno <- maize_pheno[,-31]

maize_allpheno <- maize_pheno[,-2]
write.csv(maize_allpheno, "Maize_allphospho.csv", row.names = FALSE)

system("ls")
maize_cov <-as.data.frame(maize_pheno[,2])
names(maize_cov) <- "covariates"

write.csv(maize_cov, "maize_cov.csv", row.names = FALSE)
