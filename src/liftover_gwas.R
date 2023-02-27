setwd("~/Documents/RubenLab/Data for sorghum/sorghum/harvester")

#setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Data for sorghum/sorghum/harvester")
#BiocManager::install("liftOver")
library(liftOver)
library(dplyr)
library(GenomicRanges)
#BiocManager::install("regioneR")
library(regioneR)
library(rtracklayer)

# Loading "colm" dataframe with PCAdapt pvalues
load("pcadapt_corrected.Rimage", verbose = TRUE)
PCAdapt <- colm
colm <- NULL

# SEEDs data  (Romero Navarro 2017) is in AGPv2
# This version of pcadapt is also in AGPv2
# because PCAdapt data loads into a GRanges object with warnings
# about out of range SNPS if AGPv3 is used
#
# B73_RefGen_v2_Chr.bed
# file with each chromosome start and end coordinates
# this is useful for checking the genome version of the markers
#
# REMEMBER: BED file format specifications say that bed files are indexed
# at 0 (zero indexed coordinates)!!!!!!!!!!!
#
# head -n 10 /rsstu/users/r/rrellan/sara/B73_RefGen_v2.fa.fai | awk 'BEGIN {FS="\t"}; {gsub("Chr", "", $1); gsub(":.*", "", $1); print $1 FS "0" FS $2 FS "ws"}' | sort -k1,1n > B73_RefGen_v2_Chr.bed 
#
# double check SNP coordinates!
#
# I think we must uplift all coordinates to AGPv4
#
# Indexed with samtools faidx from
# https://ftp.maizegdb.org/MaizeGDB/FTP/B73_RefGen_v2/B73_RefGen_v2.fa.gz

AGPv2 <- toGRanges("B73_RefGen_v2_Chr.bed")

seqlengths(AGPv2) <- width(AGPv2)

genome(AGPv2) <- "AGPv2"

# Loading chainfile for liftover
# source:
# http://ftp.gramene.org/CURRENT_RELEASE/assembly_chain/zea_mays/
chain_file <- "AGPv2_to_B73_RefGen_v4.chain"
str(trial)

ch <-import.chain(chain_file)
PCAdapt_v2 <- makeGRangesFromDataFrame(PCAdapt,
                                       keep.extra.columns=TRUE,
                                       ignore.strand=TRUE,
                                       seqinfo = seqinfo(AGPv2),  # Sanity check
                                       seqnames.field="CHR",
                                       start.field="BP",
                                       end.field="BP") 

PCAdapt_v4 <- liftOver(PCAdapt_v2,ch) %>% unlist()

#saveRDS(PCAdapt_v4, file= "PCAdapt_v4.RDS")

#Converting to V5

chain_file2 <- "B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain"
ch2 <- import.chain(chain_file2)
PCAdapt_v5 <-  liftOver(PCAdapt_v4,ch2) %>% unlist()

#Extracting data
y <- as.data.frame(as.character(names(PCAdapt_v5)))
x <- data.frame(PCAdapt_v5)
nrow(y)
nrow(x)
pcadaptv5 <- cbind(y,x)
colnames(pcadaptv5)[1] <- "Marker"


#Importing gwas data
gwas.subpop <- read.delim("GLM_Q28_site_subpopulation.txt", header=T, sep="")
gwas.nosubpop <- read.delim("GLM_site", header=T, sep="")

#Extracting Marker and pvalues
pcadapt.pm <- pcadaptv5[,c(1,8)]
colnames(pcadapt.pm)[2] <- "p.pcadapt"
gwas.subpop.pm <- gwas.subpop[,c(2,6)]
colnames(gwas.subpop.pm)[2] <- "p.gwas.subpop"
gwas.nosubpop.pm <- gwas.nosubpop[,c(2,6)]
colnames(gwas.nosubpop.pm)[2] <- "p.gwas.nosubpop"

#Finding Common Makers
gwas.common.markers <- inner_join(gwas.subpop.pm,gwas.nosubpop.pm) 
common.all <- inner_join(pcadapt.pm,gwas.common.markers)
ch10.pvalues <- common.all[328430:352995,]

#Making genotype file for 
head(ch10.pvalues)
head(pcadaptv5)
geno.pca <- pcadaptv5[,c(1:8)] 
head(geno.pca)
