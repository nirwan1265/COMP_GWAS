# Chain files directory
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data/Overlift")

# Required Packages
library(liftOver)
library(dplyr)
library(GenomicRanges)
library(regioneR)
library(rtracklayer)
library(vroom)

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
      assign(paste0("geno",i), vroom(paste0("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/geno/RomeroNavarro2017/MAF/chr",i,"_MAF_MAF0.01_pheno.txt")))
      geno <- c(geno, list(get(paste0("geno",i))))
      names(geno)[i] <- paste0("geno",i)
    }
  }
  return(geno)
}
system("ls")
geno <- geno_ld("Zea",10)


x <- as.numeric(gsub("S1_","",names(geno[["geno1"]])))
x <- data.frame(names(geno[["geno1"]]),as.numeric(gsub("S1_","",names(geno[["geno1"]]))))
names(x) <- c("SNP","ps")
x$chr <- 1



pheno_v2 <- makeGRangesFromDataFrame(x,
                                      keep.extra.columns=TRUE,
                                      ignore.strand=TRUE,
                                      seqinfo = seqinfo(AGPv2),  # Sanity check
                                      seqnames.field="chr",
                                      start.field="ps",
                                      end.field="ps")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S1_",pheno_v5$start)

geno1 <- geno[["geno1"]]
geno1 <- geno1[pheno_v5$SNP]
names(geno1) <- pheno_v5$newSNP
geno1[1:4,1:4]
write.table(geno1,"geno1.txt", quote = F, row.names = F)


###################################################################################################
x <- as.numeric(gsub("S2_","",names(geno[["geno2"]])))
x <- data.frame(names(geno[["geno2"]]),as.numeric(gsub("S2_","",names(geno[["geno2"]]))))
names(x) <- c("SNP","ps")
x$chr <- 1



pheno_v2 <- makeGRangesFromDataFrame(x,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chr",
                                     start.field="ps",
                                     end.field="ps")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S2_",pheno_v5$start)

geno2 <- geno[["geno2"]]
geno2 <- geno2[pheno_v5$SNP]
names(geno2) <- pheno_v5$newSNP
geno2[1:4,1:4]
write.table(geno2,"geno2.txt", quote = F, row.names = F)

############
x <- as.numeric(gsub("S3_","",names(geno[["geno3"]])))
x <- data.frame(names(geno[["geno3"]]),as.numeric(gsub("S3_","",names(geno[["geno3"]]))))
names(x) <- c("SNP","ps")
x$chr <- 1



pheno_v2 <- makeGRangesFromDataFrame(x,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chr",
                                     start.field="ps",
                                     end.field="ps")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S3_",pheno_v5$start)

geno3 <- geno[["geno3"]]
geno3 <- geno3[pheno_v5$SNP]
names(geno3) <- pheno_v5$newSNP
geno3[1:4,1:4]
write.table(geno3,"geno3.txt", quote = F, row.names = F)
######3
x <- as.numeric(gsub("S4_","",names(geno[["geno4"]])))
x <- data.frame(names(geno[["geno4"]]),as.numeric(gsub("S4_","",names(geno[["geno4"]]))))
names(x) <- c("SNP","ps")
x$chr <- 1



pheno_v2 <- makeGRangesFromDataFrame(x,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chr",
                                     start.field="ps",
                                     end.field="ps")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S4_",pheno_v5$start)

geno4 <- geno[["geno4"]]
geno4 <- geno4[pheno_v5$SNP]
names(geno4) <- pheno_v5$newSNP
write.table(geno4,"geno4.txt", quote = F, row.names = F)

#####
x <- as.numeric(gsub("S5_","",names(geno[["geno5"]])))
x <- data.frame(names(geno[["geno5"]]),as.numeric(gsub("S5_","",names(geno[["geno5"]]))))
names(x) <- c("SNP","ps")
x$chr <- 1



pheno_v2 <- makeGRangesFromDataFrame(x,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chr",
                                     start.field="ps",
                                     end.field="ps")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S5_",pheno_v5$start)

geno5 <- geno[["geno5"]]
geno5 <- geno5[pheno_v5$SNP]
names(geno5) <- pheno_v5$newSNP
write.table(geno5,"geno5.txt", quote = F, row.names = F)

#############
x <- as.numeric(gsub("S6_","",names(geno[["geno6"]])))
x <- data.frame(names(geno[["geno6"]]),as.numeric(gsub("S6_","",names(geno[["geno6"]]))))
names(x) <- c("SNP","ps")
x$chr <- 1



pheno_v2 <- makeGRangesFromDataFrame(x,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chr",
                                     start.field="ps",
                                     end.field="ps")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S6_",pheno_v5$start)

geno6 <- geno[["geno6"]]
geno6 <- geno6[pheno_v5$SNP]
names(geno6) <- pheno_v5$newSNP
write.table(geno6,"geno6.txt", quote = F, row.names = F)
#######
x <- as.numeric(gsub("S8_","",names(geno[["geno8"]])))
x <- data.frame(names(geno[["geno8"]]),as.numeric(gsub("S8_","",names(geno[["geno8"]]))))
names(x) <- c("SNP","ps")
x$chr <- 1



pheno_v2 <- makeGRangesFromDataFrame(x,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chr",
                                     start.field="ps",
                                     end.field="ps")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S8_",pheno_v5$start)

geno8 <- geno[["geno8"]]

geno8 <- geno8[pheno_v5$SNP]
names(geno8) <- pheno_v5$newSNP
geno8[1:4,1:4]
write.table(geno8,"geno8.txt", quote = F, row.names = F)

########3
x <- as.numeric(gsub("S9_","",names(geno[["geno9"]])))
x <- data.frame(names(geno[["geno9"]]),as.numeric(gsub("S9_","",names(geno[["geno9"]]))))
names(x) <- c("SNP","ps")
x$chr <- 1



pheno_v2 <- makeGRangesFromDataFrame(x,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chr",
                                     start.field="ps",
                                     end.field="ps")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S9_",pheno_v5$start)

geno9 <- geno[["geno9"]]

geno9 <- geno9[pheno_v5$SNP]
names(geno9) <- pheno_v5$newSNP
geno9[1:4,1:4]
write.table(geno9,"geno9.txt", quote = F, row.names = F)


#######3x <- as.numeric(gsub("S10_","",names(geno[["geno10"]])))
x <- data.frame(names(geno[["geno10"]]),as.numeric(gsub("S10_","",names(geno[["geno10"]]))))
names(x) <- c("SNP","ps")
x$chr <- 1



pheno_v2 <- makeGRangesFromDataFrame(x,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chr",
                                     start.field="ps",
                                     end.field="ps")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S10_",pheno_v5$start)

geno10 <- geno[["geno10"]]

geno10 <- geno10[pheno_v5$SNP]
names(geno10) <- pheno_v5$newSNP
geno10[1:4,1:4]
write.table(geno10,"geno10.txt", quote = F, row.names = F)
###########

system("ls")
#geno7 <- read.table("chr7_MAF_MAF0.01_pheno.txt", header = T)
geno7 <- vroom("chr7_MAF_MAF0.01_pheno.txt")


x <- as.numeric(gsub("S7_","",names(geno7)))
x <- data.frame(names(geno7),as.numeric(gsub("S7_","",names(geno7))))
names(x) <- c("SNP","ps")
x$chr <- 1



pheno_v2 <- makeGRangesFromDataFrame(x,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqinfo = seqinfo(AGPv2),  # Sanity check
                                     seqnames.field="chr",
                                     start.field="ps",
                                     end.field="ps")


pheno_v4 <- liftOver(pheno_v2,ch) %>% unlist()
pheno_v5 <-  liftOver(pheno_v4,ch2) %>% unlist()
pheno_v5 <- data.frame(pheno_v5)
pheno_v5$newSNP <- paste0("S7_",pheno_v5$start)

#geno7 <- geno[["geno7"]]

geno7 <- geno7[pheno_v5$SNP]
names(geno7) <- pheno_v5$newSNP
geno7[1:4,1:4]
write.table(geno7,"geno7.txt", quote = F, row.names = F)










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
