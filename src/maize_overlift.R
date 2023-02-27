# Required Packages
library(liftOver)
library(dplyr)
library(GenomicRanges)
library(regioneR)
library(rtracklayer)
library(vroom)
library(vcfR)
library(bedr)


# Chain files directory
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data/Overlift/maize")


# GRanges for V2
AGPv2 <- toGRanges("B73_RefGen_v2_Chr.bed")
seqlengths(AGPv2) <- width(AGPv2)
genome(AGPv2) <- "AGPv2"

# Loading chainfile for liftover
# source:
# http://ftp.gramene.org/CURRENT_RELEASE/assembly_chain/zea_mays/
# Chain V2 to V4
chain_file <- "AGPv2_to_B73_RefGen_v4.chain"
ch <-import.chain(chain_file)

# Chain V4 to V5
chain_file2 <- "B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain"
ch2 <- import.chain(chain_file2)



#Loading genotype:
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/geno/RomeroNavarro2017/MAF")
geno_ld <- function(organism,chr){
  geno <- list()
  if(organism == "Sorghum"){
    for(i in 1:chr){
      assign(paste0("geno",i), vroom(paste0("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/raw/africa.filtered/c",i,"_MAF_sorghum.txt")))
      geno <- c(geno, list(get(paste0("geno",i))))
      names(geno)[i] <- paste0("geno",i)
    }
  }
  else {
    for(i in 1:chr){
      assign(paste0("geno",i), vroom(paste0("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/geno/RomeroNavarro2017/hapmap/ch",i,".v5.hmp.txt")))
      geno <- c(geno, list(get(paste0("geno",i))))
      names(geno)[i] <- paste0("geno",i)
    }
  }
  return(geno)
}
system("ls")
geno <- geno_ld("Zea",1)


###############################################################################
# Overlifting Genome from v2 to v5

geno1 <- as.data.frame(geno[["geno1"]])
#geno1 <- as.data.frame(geno1[,c(1,3,4)])
geno1 <- as.data.frame(t(geno1))
geno1[1:4,1:4]

pheno_v2 <- makeGRangesFromDataFrame(geno1,
                                      keep.extra.columns=TRUE,
                                      ignore.strand=TRUE,
                                      seqinfo = seqinfo(AGPv2),  # Sanity check
                                      seqnames.field="chrom",
                                      start.field="pos",
                                      end.field="pos")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S1_",pheno_v5$start)


geno1 <- as.data.frame(geno[["geno1"]])
geno1 <- geno1[geno1$`rs#` %in% pheno_v5$`rs.`, ]
geno1$`rs#` <- pheno_v5$newSNP
geno1$pos <- pheno_v5$start
  
write.table(geno1,"geno1.hmp.txt", quote = F, row.names = F, sep = "\t")

###############################################################################
geno2 <- as.data.frame(geno[["geno2"]])
geno2 <- as.data.frame(geno2[,c(1,3,4)])


pheno_v2 <- makeGRangesFromDataFrame(geno2,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chrom",
                                     start.field="pos",
                                     end.field="pos")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S2_",pheno_v5$start)


geno2 <- as.data.frame(geno[["geno2"]])
geno2 <- geno2[geno2$`rs#` %in% pheno_v5$`rs.`, ]
geno2$`rs#` <- pheno_v5$newSNP
geno2$pos <- pheno_v5$start

write.table(geno2,"geno2.hmp.txt", quote = F, row.names = F, sep = "\t")

###############################################################################

geno3 <- as.data.frame(geno[["geno3"]])
geno3 <- as.data.frame(geno3[,c(1,3,4)])


pheno_v2 <- makeGRangesFromDataFrame(geno3,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chrom",
                                     start.field="pos",
                                     end.field="pos")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S3_",pheno_v5$start)


geno3 <- as.data.frame(geno[["geno3"]])
geno3 <- geno3[geno3$`rs#` %in% pheno_v5$`rs.`, ]
geno3$`rs#` <- pheno_v5$newSNP
geno3$pos <- pheno_v5$start

write.table(geno3,"geno3.hmp.txt", quote = F, row.names = F, sep = "\t")

###############################################################################

geno4 <- as.data.frame(geno[["geno4"]])
geno4 <- as.data.frame(geno4[,c(1,3,4)])


pheno_v2 <- makeGRangesFromDataFrame(geno4,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chrom",
                                     start.field="pos",
                                     end.field="pos")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S4_",pheno_v5$start)


geno4 <- as.data.frame(geno[["geno4"]])
geno4 <- geno4[geno4$`rs#` %in% pheno_v5$`rs.`, ]
geno4$`rs#` <- pheno_v5$newSNP
geno4$pos <- pheno_v5$start

write.table(geno4,"geno4.hmp.txt", quote = F, row.names = F, sep = "\t")

###############################################################################

geno5 <- as.data.frame(geno[["geno5"]])
geno5 <- as.data.frame(geno5[,c(1,3,4)])


pheno_v2 <- makeGRangesFromDataFrame(geno5,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chrom",
                                     start.field="pos",
                                     end.field="pos")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S5_",pheno_v5$start)


geno5 <- as.data.frame(geno[["geno5"]])
geno5 <- geno5[geno5$`rs#` %in% pheno_v5$`rs.`, ]
geno5$`rs#` <- pheno_v5$newSNP
geno5$pos <- pheno_v5$start

write.table(geno5,"geno5.hmp.txt", quote = F, row.names = F, sep = "\t")

###############################################################################

geno6 <- as.data.frame(geno[["geno6"]])
geno6 <- as.data.frame(geno6[,c(1,3,4)])


pheno_v2 <- makeGRangesFromDataFrame(geno6,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chrom",
                                     start.field="pos",
                                     end.field="pos")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S6_",pheno_v5$start)


geno6 <- as.data.frame(geno[["geno6"]])
geno6 <- geno6[geno6$`rs#` %in% pheno_v5$`rs.`, ]
geno6$`rs#` <- pheno_v5$newSNP
geno6$pos <- pheno_v5$start

write.table(geno6,"geno6.hmp.txt", quote = F, row.names = F, sep = "\t")

###############################################################################

geno7 <- as.data.frame(geno[["geno7"]])
geno7 <- as.data.frame(geno7[,c(1,3,4)])


pheno_v2 <- makeGRangesFromDataFrame(geno7,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chrom",
                                     start.field="pos",
                                     end.field="pos")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S7_",pheno_v5$start)


geno7 <- as.data.frame(geno[["geno7"]])
geno7 <- geno7[geno7$`rs#` %in% pheno_v5$`rs.`, ]
geno7$`rs#` <- pheno_v5$newSNP
geno7$pos <- pheno_v5$start

write.table(geno7,"geno7.hmp.txt", quote = F, row.names = F, sep = "\t")

###############################################################################

geno8 <- as.data.frame(geno[["geno8"]])
geno8 <- as.data.frame(geno8[,c(1,3,4)])


pheno_v2 <- makeGRangesFromDataFrame(geno8,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chrom",
                                     start.field="pos",
                                     end.field="pos")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S8_",pheno_v5$start)


geno8 <- as.data.frame(geno[["geno8"]])
geno8 <- geno8[geno8$`rs#` %in% pheno_v5$`rs.`, ]
geno8$`rs#` <- pheno_v5$newSNP
geno8$pos <- pheno_v5$start

write.table(geno8,"geno8.hmp.txt", quote = F, row.names = F, sep = "\t")

###############################################################################
geno9 <- as.data.frame(geno[["geno9"]])
geno9 <- as.data.frame(geno9[,c(1,3,4)])


pheno_v2 <- makeGRangesFromDataFrame(geno9,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chrom",
                                     start.field="pos",
                                     end.field="pos")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S9_",pheno_v5$start)


geno9 <- as.data.frame(geno[["geno9"]])
geno9 <- geno9[geno9$`rs#` %in% pheno_v5$`rs.`, ]
geno9$`rs#` <- pheno_v5$newSNP
geno9$pos <- pheno_v5$start

write.table(geno9,"geno9.hmp.txt", quote = F, row.names = F, sep = "\t")

###############################################################################

geno10 <- as.data.frame(geno[["geno10"]])
geno10 <- as.data.frame(geno10[,c(1,3,4)])


pheno_v2 <- makeGRangesFromDataFrame(geno10,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chrom",
                                     start.field="pos",
                                     end.field="pos")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S10_",pheno_v5$start)


geno10 <- as.data.frame(geno[["geno10"]])
geno10 <- geno10[geno10$`rs#` %in% pheno_v5$`rs.`, ]
geno10$`rs#` <- pheno_v5$newSNP
geno10$pos <- pheno_v5$start

write.table(geno10,"geno10.hmp.txt", quote = F, row.names = F, sep = "\t")

###############################################################################




# Loading GWAS files
folder <- "~/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data/GWAS/maize"
txt_files <- list.files(folder, pattern = "*.txt")

for (i in 1){#:length(txt_files)) {
  file <- paste0(folder, "/", txt_files[i])
  assign(gsub(".filtered.txt","",txt_files[i]),as.data.frame(vroom(file)))
}

#pheno <- read.table("apa.filtered.txt", header = T)
class(stp10)

#Making GRanges for the phenotype
for (i in 1:length(txt_files)) {
  assign(gsub(".filtered.txt","",txt_files[i]),makeGRangesFromDataFrame(get(paste0(gsub(".filtered.txt","",txt_files[i]))),
                                                                              keep.extra.columns=TRUE,
                                                                              ignore.strand=TRUE,
                                                                              seqinfo = seqinfo(AGPv2),  # Sanity check
                                                                              seqnames.field="chr",
                                                                              start.field="ps",
                                                                              end.field="ps") )
}

#pheno_v2 <- makeGRangesFromDataFrame(pheno,
#                                       keep.extra.columns=TRUE,
#                                       ignore.strand=TRUE,
#                                       seqinfo = seqinfo(AGPv2),  # Sanity check
#                                       seqnames.field="chr",
#                                       start.field="ps",
#                                       end.field="ps") 

# Liftover to V4
#pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
for (i in 1:length(txt_files)) {
  assign(gsub(".filtered.txt","",txt_files[i]), liftOver(get(paste0(gsub(".filtered.txt","",txt_files[i]))),ch) %>% unlist())
}


# Liftover to V5
#pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
for (i in 1:length(txt_files)) {
  assign(gsub(".filtered.txt","",txt_files[i]), liftOver(get(paste0(gsub(".filtered.txt","",txt_files[i]))),ch2) %>% unlist())
}


#Extracting data
for (i in 1:length(txt_files)) {
  assign(gsub(".filtered.txt","",txt_files[i]), data.frame(get(paste0(gsub(".filtered.txt","",txt_files[i])))) %>% dplyr::select(rs, seqnames, start, p_wald))
}
#pheno_v5 <- data.frame(pheno_v5)
#pheno_v5 <- phenov5 %>% dplyr::select(rs, seqnames, start, p_wald)
#colnames(phenov5)[1] <- "Marker"


#Saving data
for (i in 1:length(txt_files)) {
  write.table(get(paste0(gsub(".filtered.txt","",txt_files[i]))),paste0(gsub(".filtered.txt","",txt_files[i]),".txt"), row.names = F, quote = FALSE)
}






################################################################################
# Sorghum
################################################################################

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data/Overlift/sorghum")

# Loading chainfile for liftover
# source:
# http://ftp.gramene.org/CURRENT_RELEASE/assembly_chain/sorghum_bicolor/
# two versions of v1
# This paper: https://www.nature.com/articles/nature07723 referenced on Lasky says BTx642 v1.1

# Chain V1 to V3
chain_file <- "Sorbi1_to_Sorghum_bicolor_v2_tabs.chain"
# Before importing chain file:
# Have to convert spaces to tabs first. Do this in bash
# sed -r 's/^([0-9]+) ([0-9]+) ([0-9]+)$/\1\t\2\t\3/' Sorbi1_to_Sorghum_bicolor_NCBIv3.chain > Sorbi1_to_Sorghum_bicolor_
ch <-import.chain(chain_file)


#Granges for V1 refrence genome
NCBIv1 <- toGRanges("SbicolorBTx642_564_v1.0.bed")
seqlengths(NCBIv1) <- width(NCBIv1)
genome(NCBIv1) <- "NCBIv1"



# Granges for Genotypes
setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/raw/africa.filtered/")

# Reading vcf file
# The vcf file with all the chrom
x <- read.vcf("allchrom_v1.vcf")

# Converting to bed file
y <- vcf2bed(x)

# Creating the required Dataframe
geno_v1 <- data.frame(SNP = paste0("S",y$chr,"_",as.numeric(y$start)+1), CHR = as.numeric(y$chr), BP = (as.numeric(y$start)+1))
geno_v1 <- data.frame(SNP = paste0("S",y$chr,"_",as.numeric(y$start)+1), CHR = paste0("chromosome_",as.numeric(y$chr)), BP = (as.numeric(y$start)+1))
str(geno_v1)

# Making genomic ranges from dataframe
geno_v1 <- makeGRangesFromDataFrame(geno_v1,
                                       keep.extra.columns=TRUE,
                                       ignore.strand=TRUE,
                                       #seqinfo = seqinfo(NCBIv1),  # Sanity check
                                       seqnames.field="CHR",
                                       start.field="BP",
                                       end.field="BP") 

unique(seqnames(NCBIv1))
unique(seqnames(geno_v1))


# Uplift from V1 (geno) to V3 using ch file 
geno_v3 <- liftOver(geno_v1,ch) %>% unlist()
geno_v3 <- as.data.frame(geno_v3) %>% arrange(seqnames)




geno_v3_1 <- geno_v3[grepl("^S1_", geno_v3$SNP), ]
geno_v3_2 <- geno_v3[grepl("^S2_", geno_v3$SNP), ]
geno_v3_3 <- geno_v3[grepl("^S3_", geno_v3$SNP), ]
geno_v3_4 <- geno_v3[grepl("^S4_", geno_v3$SNP), ]
geno_v3_5 <- geno_v3[grepl("^S5_", geno_v3$SNP), ]
geno_v3_6 <- geno_v3[grepl("^S6_", geno_v3$SNP), ]
geno_v3_7 <- geno_v3[grepl("^S7_", geno_v3$SNP), ]
geno_v3_8 <- geno_v3[grepl("^S8_", geno_v3$SNP), ]
geno_v3_9 <- geno_v3[grepl("^S9_", geno_v3$SNP), ]
geno_v3_10 <- geno_v3[grepl("^S10_", geno_v3$SNP), ]

names(geno_v3)
unique(geno_v3$seqnames)

geno_v3_1[1:4,1:6]

geno_v3_1$V3 <- paste0("S1_",geno_v3_1$start)
geno_v3_2$V3 <- paste0("S2_",geno_v3_2$start)
geno_v3_3$V3 <- paste0("S3_",geno_v3_3$start)
geno_v3_4$V3 <- paste0("S4_",geno_v3_4$start)
geno_v3_5$V3 <- paste0("S5_",geno_v3_5$start)
geno_v3_6$V3 <- paste0("S6_",geno_v3_6$start)
geno_v3_7$V3 <- paste0("S7_",geno_v3_7$start)
geno_v3_8$V3 <- paste0("S8_",geno_v3_8$start)
geno_v3_9$V3 <- paste0("S9_",geno_v3_9$start)
geno_v3_10$V3 <- paste0("S10_",geno_v3_10$start)


# Splitting based on chromosome
#geno_v3_list <- split(geno_v3, geno_v3$seqnames)
#str(geno_v3_list)
# for (i in 1:length(geno_v3_list)) {
#   assign(paste0("geno_v3_", i), geno_v3_list[[i]])
# }



# Adding V1 column
# a = 1
# for(i in paste0("geno_v3_",1:10)){
#   d =get(i)
#   d$V3 = paste0("S",a,"_",d$start)
#   a = a+1
#   assign(i,d)
# }


# Load all the genotypes
setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/raw/africa.filtered/v1/")
for (i in 1:10){
  assign(paste0("geno", i), vroom(paste0("c",i,"_MAF_sorghum.txt")))
}

geno1[1:4,1:4]
head(geno_v3_1) 

# Indexing and replacing the column names
index1 <- geno_v3_1$SNP
geno1 <- geno1[,index1]
names(geno) <- geno_v3_1$V3



# for (i in 1:10) {
#   # get the index for this iteration
#   index <- get(paste0("geno_v3_", i))$SNP
#   
#   # subset the corresponding geno data frame
#   geno <- get(paste0("geno", i))[, index]
#   
#   # set the column names
#   names(geno) <- get(paste0("geno_v3_", i))$V3
#   
#   # assign the updated geno data frame back to the original variable name
#   assign(paste0("geno", i), geno)
# }
# #IGNORE THE ERROR
# 
# 
# # Saving the MAF files
# for(i in 1:10){
#   write.table(get(paste0("geno",i)), paste0("geno",i,".txt"), row.names = F, quote = F)
# }
# 

###
# Load the hapmap file
setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/raw/africa.filtered/v1/hapmap")
for(i in 1:10){
  assign(paste0("hmp.ch",i),vroom(paste0("ch",i,".hmp.txt")))
}

#Indexing
index1 <- geno_v3_1$SNP
index2 <- geno_v3_2$SNP
index3 <- geno_v3_3$SNP
index4 <- geno_v3_4$SNP
index5 <- geno_v3_5$SNP
index6 <- geno_v3_6$SNP
index7 <- geno_v3_7$SNP
index8 <- geno_v3_8$SNP
index9 <- geno_v3_9$SNP
index10 <- geno_v3_10$SNP

# Granges for Genotypes
hmp.ch3[1:4,1:4]
index3[1:4]
names(hmp.ch1)[1:4]
geno_v3_1$V3[1:4]
geno_v3_1[1:4,]
names(geno_v3_1)

hmp.ch1 <- hmp.ch1[hmp.ch1$`rs#` %in% index1, ]
hmp.ch1$`rs#` <- geno_v3_1$V3
hmp.ch1$pos <- geno_v3_1$start

hmp.ch2 <- hmp.ch2[hmp.ch2$`rs#` %in% index2, ]
hmp.ch2$`rs#` <- geno_v3_2$V3
hmp.ch2$pos <- geno_v3_2$start

hmp.ch3 <- hmp.ch3[hmp.ch3$`rs#` %in% index3, ]
hmp.ch3$`rs#` <- geno_v3_3$V3
hmp.ch3$pos <- geno_v3_3$start

hmp.ch4 <- hmp.ch4[hmp.ch4$`rs#` %in% index4, ]
hmp.ch4$`rs#` <- geno_v3_4$V3
hmp.ch4$pos <- geno_v3_4$start

hmp.ch5 <- hmp.ch5[hmp.ch5$`rs#` %in% index5, ]
hmp.ch5$`rs#` <- geno_v3_5$V3
hmp.ch5$pos <- geno_v3_5$start

hmp.ch6 <- hmp.ch6[hmp.ch6$`rs#` %in% index6, ]
hmp.ch6$`rs#` <- geno_v3_6$V3
hmp.ch6$pos <- geno_v3_6$start

hmp.ch7 <- hmp.ch7[hmp.ch7$`rs#` %in% index7, ]
hmp.ch7$`rs#` <- geno_v3_7$V3
hmp.ch7$pos <- geno_v3_7$start

hmp.ch8 <- hmp.ch8[hmp.ch8$`rs#` %in% index8, ]
hmp.ch8$`rs#` <- geno_v3_8$V3
hmp.ch8$pos <- geno_v3_8$start


hmp.ch9 <- hmp.ch9[hmp.ch9$`rs#` %in% index9, ]
hmp.ch9$`rs#` <- geno_v3_9$V3
hmp.ch9$pos <- geno_v3_9$start

hmp.ch10 <- hmp.ch10[hmp.ch10$`rs#` %in% index10, ]
hmp.ch10$`rs#` <- geno_v3_10$V3
hmp.ch10$pos <- geno_v3_10$start


# Save hmp file
for(i in 1:10){
  write.table(get(paste0("hmp.ch",i)), paste0("chr",i,".hmp.txt"), row.names = F, quote = F, sep = "\t")
}
