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
      assign(paste0("geno",i), vroom(paste0("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/geno/RomeroNavarro2017/vcf/chr",i,".hmp.txt")))
      geno <- c(geno, list(get(paste0("geno",i))))
      names(geno)[i] <- paste0("geno",i)
    }
  }
  return(geno)
}
system("ls")
geno <- geno_ld("Zea",10)


###############################################################################
# Overlifting Genome from v2 to v5

geno1 <- as.data.frame(geno[["geno1"]])
geno1 <- as.data.frame(geno1[,c(1,3,4)])


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
