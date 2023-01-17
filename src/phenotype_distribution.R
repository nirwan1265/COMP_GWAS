# Phenotype distribution
setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
pheno <- read.csv("Sorghum_allphospho_africa.csv")
pheno <- pheno[,-1]
pheno_t <- pheno[,c(1:21)]
library(ggplot2)
plots <- lapply(pheno_t,function(x){
  ggplot(data.frame(x),aes(x)) +
    geom_density() +
    ggtitle(as.character(substitute(x)))
})
library(cowplot)
#par(mfrow = c(2,5))
plot_grid(plotlist = plots, ncol = 7)
quartz()

