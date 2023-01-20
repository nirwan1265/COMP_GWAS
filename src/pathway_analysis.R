
# Read the pathway database:
# To convert to ENSEMBL or NCBI naming system, change Sobic. to SORBI_3
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum_pathway_database/pythozome/sorghumbicolorcyc/7.0/data")
pathway <- read.csv("pathways_ensembl.csv")
pathway <- as.data.frame(t(pathway[,-1]))
colnames(pathway) <- pathway[1,]
pathway <- pathway[-1,]



#Read the combined pvalue files
#Get the final table for Omni and Magma
sorghum_omni
sorghum_omni$GeneName <- gsub("SORBI","SORBI_", sorghum_omni$GeneName)


#Filtering the pathway genes
x <- as.data.frame(as.matrix(NA))
for(i in 1:ncol(pathway)){
  for(j in 1:nrow(sorghum_omni)){
    if(sorghum_omni[j,1] %in% pathway[,i] == TRUE){
      x[j,i] <- sorghum_omni[j,1]
    }else {
      x[j,i] <- 0
    }  
  }
}

#Sorting the pathway genes and naming the pathways
sorghum.pathway <- as.data.frame(apply(x,2,sort,decreasing = TRUE))
colnames(sorghum.pathway) <- colnames(pathway)

#Removing empty pathways
sorghum.pathway <- sorghum.pathway[,colSums(sorghum.pathway != 0) > 0]

#Sorting by the highest number of genes
sorghum.pathway <- sorghum.pathway[,order(colSums(sorghum.pathway != 0), decreasing = TRUE)]

#SAVING AND READING pathway files
#write.csv(sorghum.pathway,"pathway_OMNI_sorghum.csv")

trial <- sorghum_omni[,c(1,4)]
####################################################################################
####################################################################################

# Adding the pvalues
merged_table <- merge(sorghum.pathway, trial, all.x = TRUE)






