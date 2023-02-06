library(IRanges)
library(GenomicRanges)


setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/ensemblgenomes")
j <- 1
for(i in 1:10 ){
  assign(paste0("db.",sprintf("%02d", j)), read.table(file = paste0("Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.chromosome.",i,".gff3"), sep = "\t", header = TRUE))
  
  j <-  j + 1
}

for(i in paste0("db.", sprintf("%02d", 1:10))){
  d = get(i)
  colnames(d) <- c("Chromosome","Database","Region","Start","End","NA","Strand","NA2","Gene")
  assign(i,d)
}



for(i in paste0("db.", sprintf("%02d", 1:10))){
  d = get(i)
  d$Start <- d$Start-2000
  d$Start[1:4] <- 0
  d$End <- d$End+2000
  assign(i,d)
}




###Annotation the SNPs
##Making GRanges for Database 
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gr.db", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("db.",i))[,"Start"], end = get(paste0("db.",i))[,"End"]), strand = get(paste0("db.",i))[,"Strand"], Region = get(paste0("db.",i))[,"Region"], Gene = get(paste0("db.",i))[,"Gene"]))
}


setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data/GenomicRanges/sorghum")
for(i in sprintf("%02d", 1:10)){
  saveRDS(get(paste0("gr.db", i)), paste0("gr.db", i,".RDS"))
}


