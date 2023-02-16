###PCA analysis required for GBJ
#Need a vcf file format of the hapmap which is converted to GDS format
#TASSEL or PLINK is used for converting hapmap to VCF file format
#Need a directory to  create the gds file. If working on the server, we might need to define this before starting
#Sorghum
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/raw/africa.filtered")

##Reading the vcf files
for(i in sprintf("%02d", 1)){
  assign(paste0("vcf.fn",i),paste0("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/raw/africa.filtered/allchrom.recode.vcf"))
}


#Maize
#setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/geno/RomeroNavarro2017/vcf")

##Reading the vcf files
# for(i in sprintf("%02d", 1)){
#   assign(paste0("vcf.fn",i),paste0("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/geno/RomeroNavarro2017/vcf/allchrom_maize_MAF0.01_pheno.vcf"))
# }

##Converting vcf to gds
#A bit time consuming
#Note: for some reason, you cannot run the next step twice if you make an error. you need to delete all this converted gds files, remove all your env variables and do it again.
j <- 1
for(i in paste0("vcf.fn",sprintf("%02d", 1))){
  d = get(i)
  snpgdsVCF2GDS(d, paste0("chr",sprintf("%02d",j),".gds"), method = "copy.num.of.ref")
  assign(i,d)
  j = j + 1
}


##Get the GDS file data
for(i in sprintf("%02d", 1)){
  assign(paste0("gdsfile",i), snpgdsOpen(paste0("chr",i,".gds")))
}

##LD-based SNP pruning
set.seed(1000)
# Try different LD thresholds for sensitivity analysis but read in a paper somewhere that 0.2 was used for GBJ
j <- 1
for(i in paste0("gdsfile",sprintf("%02d", 1))){
  d = get(i)
  assign(paste0("snpset",sprintf("%02d",j)), snpgdsLDpruning(d,ld.threshold = 0.2))
  assign(i,d)
  j = j + 1
}


## Get all selected snp id
j <- 1
for(i in paste0("snpset",sprintf("%02d",1))){
  d = get(i)
  assign(paste0("snpset.id",sprintf("%02d",j)), unlist(unname(d)))
  assign(i,d)
  j = j + 1
}

## Run PCA
j <- 1
for(i in paste0("snpset.id",sprintf("%02d",1))){
  d = get(i)
  assign(paste0("pca",sprintf("%02d",j)), snpgdsPCA(get(paste0("gdsfile",sprintf("%02d",j))), snp.id = d, num.thread = 10))
  assign(i,d)
  j = j+1
}


#In case there are population information
#https://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html
#In the case of no prior population information,
#First two shows max variance
#Make a table of eigen values
j <- 1
for(i in paste0("pca",sprintf("%02d",1))){
  d = get(i)
  assign(paste0("tab",sprintf("%02d",j)), data.frame(sample.id = d$sample.id,
                                                     EV1 = d$eigenvect[,1],
                                                     EV2 = d$eigenvect[,2],
                                                     EV3 = d$eigenvect[,3],
                                                     EV4 = d$eigenvect[,4],
                                                     EV5 = d$eigenvect[,5],
                                                     EV6 = d$eigenvect[,6],
                                                     EV7 = d$eigenvect[,7],
                                                     EV8 = d$eigenvect[,8],
                                                     EV9 = d$eigenvect[,9],
                                                     EV10 = d$eigenvect[,10],
                                                     EV11 = d$eigenvect[,11],
                                                     EV12 = d$eigenvect[,12],
                                                     EV13 = d$eigenvect[,13],
                                                     EV14 = d$eigenvect[,14],
                                                     EV15 = d$eigenvect[,15],
                                                     EV16 = d$eigenvect[,16],
                                                     EV17 = d$eigenvect[,17],
                                                     EV18 = d$eigenvect[,18],
                                                     EV19 = d$eigenvect[,19],
                                                     EV20 = d$eigenvect[,20],
                                                     EV21 = d$eigenvect[,21],
                                                     EV22 = d$eigenvect[,22],
                                                     EV23 = d$eigenvect[,23],
                                                     EV24 = d$eigenvect[,24],
                                                     EV25 = d$eigenvect[,25],
                                                     stringsAsFactors = FALSE))
  assign(i,d)
  j = j + 1
}

system("pwd")
saveRDS(pca01,"pca_sorghum.RDS")

# Variance Proportion (%)
pc.percent <- pca01$varprop*100
round(pc.percent)
plot(tab01$EV1,tab01$EV2)

#Population information 
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
pop_code <- read.csv("covariate_sorghum_region.csv")
head(pop_code)
pop_code <- as.data.frame(cbind(tab01$sample.id, pop_code$region))
names(pop_code) <- c("ID","Region")

plot(tab01$EV1,tab01$EV2, col = as.integer(pop_code$Region), xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomright", legend=levels(as.factor(pop_code$Region)), pch="o", col=1:nlevels(as.factor(pop_code$Region)))


lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca01$eigenvect[,1:4], col=pop_code$Region, labels=lbls)


########################
snp_selected <- pca01$snp.id
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/raw/africa.filtered")
maf_map <- vroom("allchrom_africa_filtered.MAF.txt")
maf_map[1:4,1:4]
maf_map2 <- maf_map[,snp_selected]
library(pcaL1)
mypcal1 <- pcal1(as.matrix(maf_map2), projDim=2, center=FALSE, projections="l2pca", initialize = "l2pca")

myawl1pca <- awl1pca(as.matrix(maf_map2))





########################################################################################


