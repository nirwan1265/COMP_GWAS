# Chain files directory
setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data/Overlift")

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


# Loading GWAS files
folder <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data/GWAS/maize"
txt_files <- list.files(folder, pattern = "*.txt")

for (i in 1:length(txt_files)) {
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
