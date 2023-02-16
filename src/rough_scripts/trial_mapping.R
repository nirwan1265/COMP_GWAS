library(GenomicRanges)

setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data/GenomicRanges/sorghum")

for(i in sprintf("%02d",1:10)){
  assign(paste0("gr.db",i), readRDS(paste0("gr.db",i,".RDS")))
}

sol_VL_annot <- sol_VL[,c(2,3,4,6)]

sol_VL_annot$Chr <- sprintf("chr%02s", sol_VL_annot$Chr)

# Create a vector of unique chromosome values
chrs <- unique(sol_VL_annot$Chr)

# Split the data frame based on the chromosome values
split_sol_VL_annots <- split(sol_VL_annot, sol_VL_annot$Chr)

# Assign each split data frame to a separate object
for (i in 1:length(chrs)) {
  assign(paste0("sol_VL_annot_", chrs[i]), split_sol_VL_annots[[i]])
}

# Grange _sol
for(i in sprintf("%02d", 1:10)){
  assign(paste0("gr.q", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = get(paste0("sol_VL_annot_chr",i))[,"Pos"], width = 1,  Marker = get(paste0("sol_VL_annot_chr",i))[,"Marker"],pvalue = get(paste0("sol_VL_annot_chr",i))[,"p"])))
}



# overlap
for(i in sprintf("%02d", 1:10)){
  assign(paste0("common",i), as.data.frame(findOverlapPairs(get(paste0("gr.db",i)), get(paste0("gr.q",i)))))
}


chr3 <- as.data.frame(gr.db03)
