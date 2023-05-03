#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################
# PATHWAY ANALYSIS
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

################################################################################
# Pathway downloaded from
################################################################################
# https://ftp.plantcyc.org/pmn/Pathways/Data_dumps/PMN15.5_January2023/pathways/
#  sorghumbicolorcyc_pathways.20230103
#  corncyc_pathways.20230103



# loading the data
dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Pathways/"

# Preparing the pathway for sorghum
sorghum_pathway <- read.delim(paste0(dir,"sorghumbicolorcyc_pathways.20230103"))
sorghum_pathway <- sorghum_pathway[,c(2,8)]
sorghum_pathway <- as.data.frame(sorghum_pathway[which(sorghum_pathway$Gene.name != "unknown"), ])
setDT(sorghum_pathway)
sorghum_pathway <- as.data.frame(dcast(sorghum_pathway, rowid(Pathway.name) ~ Pathway.name, value.var = "Gene.name", fill = NA_character_))
sorghum_pathway <- sorghum_pathway[,-1]
sorghum_pathway <- sorghum_pathway %>%
  mutate(across(everything(), ~ str_replace(., "^Sobic.", "SORBI_3")))

# Preparing the pathway for maize
maize_pathway <- read.delim(paste0(dir,"corncyc_pathways.20230103"))
maize_pathway <- maize_pathway[,c(2,8)]
maize_pathway <- as.data.frame(maize_pathway[which(maize_pathway$Gene.name != "unknown"), ])
setDT(maize_pathway)
maize_pathway <- as.data.frame(dcast(maize_pathway, rowid(Pathway.name) ~ Pathway.name, value.var = "Gene.name", fill = NA_character_))
maize_pathway <- maize_pathway[,-1]

# NOTRUN
# Read the pathway database:
# To convert to ENSEMBL or NCBI naming system, change Sobic. to SORBI_3
# setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum_pathway_database/pythozome/sorghumbicolorcyc/7.0/data")
# pathway <- read.csv("pathways_ensembl.csv")
# pathway <- as.data.frame(t(pathway[,-1]))
# colnames(pathway) <- pathway[1,]
# pathway <- pathway[-1,]
#

################################################################################
# DATA PREPARATION
################################################################################


#Read the combined magma for maize and sorghum files
dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/Top_MAGMA_hits/"

# Top maize and sorghum hits
# combined_magma_maize <- read.csv(paste0(dir,"maize/Combined_MAGMA_maize_4pheno.csv"))
# combined_magma_sorghum <- read.csv(paste0(dir,"sorghum/Combined_MAGMA_sorghum_4pheno.csv"))

# For MAGMA hits
# combined_magma_maize <- read.csv(paste0(dir,"maize/Top_MAGMA_hits_maize_4pheno.csv"), header= T)
# combined_magma_maize <- combined_magma_maize[which(combined_magma_maize$phenotype == "stp10"), ]
# combined_magma_maize <- combined_magma_maize[which(combined_magma_maize$P_MULTI <= 0.05), ]
# combined_magma_sorghum <- read.csv(paste0(dir,"sorghum/Top_MAGMA_hits_sorghum_4pheno.csv"), header = T)
# combined_magma_sorghum <- combined_magma_sorghum[which(combined_magma_sorghum$phenotype == "stp10"), ]
# combined_magma_sorghum <- combined_magma_sorghum[which(combined_magma_sorghum$P_MULTI <= 0.05), ]
# NOTRUN
#sorghum_omni$GeneName <- gsub("SORBI","SORBI_", sorghum_omni$GeneName)

# Maize
dir_maize <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/MAGMA/maize/MAGMA_maize_2kb_v5_csv/"

# Load MAGMA results
PBR_maize <- read.table(paste0(dir_maize,"PBR1_maize_magma_multi_snpwise.genes.out"), header = T) %>% dplyr::select(c(1,9))
PNZ_maize <- read.table(paste0(dir_maize,"PNZ1_maize_magma_multi_snpwise.genes.out"), header = T) %>% dplyr::select(c(1,9))
POL_maize <- read.table(paste0(dir_maize,"POL1_maize_magma_multi_snpwise.genes.out"), header = T) %>% dplyr::select(c(1,9))
stp_maize <- read.table(paste0(dir_maize,"stp10_maize_magma_multi_snpwise.genes.out"), header = T)%>% dplyr::select(c(1,9))

# File directory
dir_sorghum <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/MAGMA/sorghum/MAGMA_sorghum_2kb_v3_csv/"

# Load MAGMA results
PBR_sorghum <- read.table(paste0(dir_sorghum,"PBR1_magma_multi_snpwise.genes.out"), header = T) %>% dplyr::select(c(1,9))
PNZ_sorghum <- read.table(paste0(dir_sorghum,"PNZ1_magma_multi_snpwise.genes.out"), header = T) %>% dplyr::select(c(1,9))
POL_sorghum <- read.table(paste0(dir_sorghum,"POL1_magma_multi_snpwise.genes.out"), header = T) %>% dplyr::select(c(1,9))
stp_sorghum <- read.table(paste0(dir_sorghum,"stp10_magma_multi_snpwise.genes.out"), header = T) %>% dplyr::select(c(1,9))



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################
# RUNNING PATHWAY ANALYSIS ON SORGHUM
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#Filtering the pathway genes
x <- as.data.frame(as.matrix(NA))
# for(i in 1:ncol(sorghum_pathway)){
#   for(j in 1:nrow(combined_magma_sorghum)){
#     if(combined_magma_sorghum[j,1] %in% sorghum_pathway[,i] == TRUE){
#       x[j,i] <- combined_magma_sorghum[j,1]
#     }else {
#       x[j,i] <- 0
#     }  
#   }
# }
PBR_sorghum <- PBR_sorghum[1:1000,]
for(i in 1:ncol(sorghum_pathway)){
  for(j in 1:nrow(PBR_sorghum)){
    if(PBR_sorghum[j,1] %in% sorghum_pathway[,i] == TRUE){
      x[j,i] <- PBR_sorghum[j,1]
    }else {
      x[j,i] <- 0
    }  
  }
}

#Sorting the pathway genes and naming the pathways
sorghum_pathway_filtered <- as.data.frame(apply(x,2,sort,decreasing = TRUE))
colnames(sorghum_pathway_filtered) <- colnames(sorghum_pathway)

#Removing empty pathways
sorghum_pathway_filtered <- sorghum_pathway_filtered[,colSums(sorghum_pathway_filtered != 0) > 0]

#Sorting by the highest number of genes
sorghum_pathway_filtered <- sorghum_pathway_filtered[,order(colSums(sorghum_pathway_filtered != 0), decreasing = TRUE)]

#SAVING AND READING pathway files
#write.csv(sorghum.pathway,"pathway_OMNI_sorghum.csv")



################################################################################
# GETTING PVALUES
################################################################################
# convert P_MULTI column into long format
combined_magma_sorghum_long <- PBR_sorghum %>%
  tidyr::separate_rows(P_MULTI, sep = ",") %>%
  dplyr::select(GENE, P_MULTI)


# join with sorghum_pathway_filtered based on common genes
joined <- sorghum_pathway_filtered %>%
  tidyr::pivot_longer(cols = everything(), names_to = "pathway", values_to = "gene") %>%
  left_join(combined_magma_sorghum_long, by = c("gene" = "GENE")) %>%
  filter(!is.na(P_MULTI)) # remove rows with NA values in P_MULTI

# pivot wider to get the desired format
result <- joined %>%
  dplyr::select(pathway, P_MULTI) %>%
  dplyr::mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = pathway, values_from = P_MULTI)

result <- result[,-1]
# Find length of longest non-NA vector
max_len <- max(sapply(result, function(x) length(na.omit(x))))

# Create output table
output_table <- data.frame(matrix(NA, nrow = max_len, ncol = ncol(result)))
colnames(output_table) <- colnames(result)

# Fill in output table
for (i in 1:ncol(result)) {
  non_na_vals <- na.omit(result[[i]])
  output_table[1:length(non_na_vals), i] <- non_na_vals
}
output_table <- output_table[,-1]

# Finding highest number of elements in the column
non_na_counts <- apply(output_table, 2, function(x) sum(!is.na(x)))
max_non_na <- max(non_na_counts)


# Replacing NA with 0
output_table[is.na(output_table)] <- 0
output_table <- as.data.frame(apply(output_table,2,sort,decreasing = TRUE))
output_table <- output_table[1:max_non_na,]
output_table[output_table == 0] <- NA


################################################################################
# ACAT
################################################################################

#Doing ACAT:
pathway_sorghum_acat <- do.call(rbind, lapply(output_table, function(col) {
  if (sum(!is.na(col)) > 0) {
    pvalue <- as.numeric(na.omit(col))
    o <- ACAT(Pvals = pvalue)
    return(data.frame(ACAT_result = o))
  } else {
    return(data.frame(ACAT_result = NA))
  }
})) %>% arrange(ACAT_result)

################################################################################
# PLOTTING
################################################################################
str(pathway_sorghum_acat)
x <- pathway_sorghum_acat
x$categories <- rownames(pathway_sorghum_acat)
filtered_data <- x %>% dplyr::filter(ACAT_result <= 0.05) %>% arrange(ACAT_result) #%>% dplyr::slice(1:10)


# wrap long category labels into two lines based on a threshold of 25 characters
filtered_data$categories <- str_wrap(filtered_data$categories, width = 50)

# create bar plot using ggplot2
library(ggplot2)
quartz()

ggplot(filtered_data, aes(x = filtered_data$categories, y = filtered_data$ACAT_result, fill = filtered_data$ACAT_pvalue)) +
  geom_bar(stat = "identity", aes(fill = ACAT_result)) +
  scale_fill_gradient(low = "red", high = "blue", limits = c(0, 0.05), 
                      name = "p-value") +
  labs(title = "Pathway Enrichment through CCT", x = "Biological Processes", y = "p-value") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))


theme(panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12)) 

str(filtered_data)
filtered_data$categories


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################
# RUNNING PATHWAY ANALYSIS ON MAIZE
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#Filtering the pathway genes
x <- as.data.frame(as.matrix(NA))
for(i in 1:ncol(maize_pathway)){
  for(j in 1:nrow(combined_magma_maize)){
    if(combined_magma_maize[j,1] %in% maize_pathway[,i] == TRUE){
      x[j,i] <- combined_magma_maize[j,1]
    }else {
      x[j,i] <- 0
    }  
  }
}

#Sorting the pathway genes and naming the pathways
maize_pathway_filtered <- as.data.frame(apply(x,2,sort,decreasing = TRUE))
colnames(maize_pathway_filtered) <- colnames(maize_pathway)

#Removing empty pathways
maize_pathway_filtered <- maize_pathway_filtered[,colSums(maize_pathway_filtered != 0) > 0]

#Sorting by the highest number of genes
maize_pathway_filtered <- maize_pathway_filtered[,order(colSums(maize_pathway_filtered != 0), decreasing = TRUE)]

#SAVING AND READING pathway files
#write.csv(maize.pathway,"pathway_OMNI_maize.csv")


################################################################################
# GETTING PVALUES
################################################################################
# convert P_MULTI column into long format
combined_magma_maize_long <- combined_magma_maize %>%
  tidyr::separate_rows(P_MULTI, sep = ",") %>%
  dplyr::select(GENE, P_MULTI)


# join with maize_pathway_filtered based on common genes
joined <- maize_pathway_filtered %>%
  tidyr::pivot_longer(cols = everything(), names_to = "pathway", values_to = "gene") %>%
  left_join(combined_magma_maize_long, by = c("gene" = "GENE"))

# pivot wider to get the desired format
result <- joined %>%
  dplyr::select(pathway, P_MULTI) %>%
  dplyr::mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = pathway, values_from = P_MULTI)

result <- result[,-1]
# Find length of longest non-NA vector
max_len <- max(sapply(result, function(x) length(na.omit(x))))

# Create output table
output_table <- data.frame(matrix(NA, nrow = max_len, ncol = ncol(result)))
colnames(output_table) <- colnames(result)

# Fill in output table
for (i in 1:ncol(result)) {
  non_na_vals <- na.omit(result[[i]])
  output_table[1:length(non_na_vals), i] <- non_na_vals
}
output_table <- output_table[,-1]

# Finding highest number of elements in the column
non_na_counts <- apply(output_table, 2, function(x) sum(!is.na(x)))
max_non_na <- max(non_na_counts)


# Replacing NA with 0
output_table[is.na(output_table)] <- 0
output_table <- as.data.frame(apply(output_table,2,sort,decreasing = TRUE))
output_table <- output_table[1:max_non_na,]
output_table[output_table == 0] <- NA

# Saving


################################################################################
################################################################################
# Fisher's Exact Test
################################################################################
################################################################################
library(metap)
# Functions
load_files <- function(dir_path) {
  # Read file names
  data_list <- list.files(dir_path, pattern = "\\.csv")
  
  # Create a named list of file paths
  file_paths <- file.path(dir_path, data_list)
  names(file_paths) <- gsub("\\.csv", "", data_list)
  
  # Load files into a named list of data frames
  data_frames <- lapply(file_paths, read.csv)
  names(data_frames) <- gsub("\\.csv", "", data_list)
  
  return(data_frames)
}
combine_pvalues_fisher <- function(data_frame) {
  combined_pvalues <- do.call(rbind, lapply(data_frame, function(col) {
    non_na_pvalues <- as.numeric(na.omit(col))
    num_pvalues <- length(non_na_pvalues)
    
    if (num_pvalues > 1) {
      o <- sumlog(p = non_na_pvalues)
      return(data.frame(Fishers_result = o$p))
    } else if (num_pvalues == 1) {
      return(data.frame(Fishers_result = non_na_pvalues))
    } else {
      return(data.frame(Fishers_result = NA))
    }
  })) 
  
  return(combined_pvalues)
}


################################################################################
# Maize
################################################################################
##### Loading the files
dir_pathway_maize <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/Pathway/maize/"
maize_data_frames <- load_files(dir_pathway_maize)

# Apply the function to each data frame in the lists
maize_combined_pvalues <- lapply(maize_data_frames, combine_pvalues_fisher)
# Bind the columns from the list of data frames
maize_pvalues_df <- do.call(cbind, maize_combined_pvalues)

# Set the column names to the descriptions of the pathways
colnames(maize_pvalues_df) <- names(maize_combined_pvalues)


################################################################################
# Sorghum
################################################################################
##### Loading the files
dir_pathway_sorghum <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/Pathway/sorghum/"
sorghum_data_frames <- load_files(dir_pathway_sorghum)

# Apply the function to each data frame in the lists
sorghum_combined_pvalues <- lapply(sorghum_data_frames, combine_pvalues_fisher)

equalize_rows <- function(data_frames_list) {
  max_rows <- max(sapply(data_frames_list, nrow))
  
  equalized_data_frames <- lapply(data_frames_list, function(df) {
    diff_rows <- max_rows - nrow(df)
    
    if (diff_rows > 0) {
      new_rows <- data.frame(matrix(NA, nrow = diff_rows, ncol = ncol(df)))
      colnames(new_rows) <- colnames(df)
      return(rbind(df, new_rows))
    } else {
      return(df)
    }
  })
  
  return(equalized_data_frames)
}

# Equalize the number of rows in the sorghum_combined_pvalues list
sorghum_equalized_pvalues <- equalize_rows(sorghum_combined_pvalues)

# Bind the columns from the equalized list of data frames
sorghum_pvalues_df <- do.call(cbind, sorghum_equalized_pvalues)

# Set the column names to the descriptions of the pathways
colnames(sorghum_pvalues_df) <- names(sorghum_equalized_pvalues)


#Save:
write.csv(sorghum_pvalues_df,"sorghum_MLM_pathway_fishers.csv")
write.csv(maize_pvalues_df,"maize_MLM_pathway_fishers.csv")
################################################################################
# ACAT
################################################################################

#Doing ACAT:
pathway_maize_acat <- do.call(rbind, lapply(output_table, function(col) {
  if (sum(!is.na(col)) > 0) {
    pvalue <- as.numeric(na.omit(col))
    o <- ACAT(Pvals = pvalue)
    return(data.frame(ACAT_result = o))
  } else {
    return(data.frame(ACAT_result = NA))
  }
})) %>% arrange(ACAT_result)


# Saving
write.csv(pathway_maize_acat,"pathway_maize_commonPheno_onlysig_acat.csv", row.names = T)
write.csv(pathway_sorghum_acat,"pathway_sorghum_commonPheno_onlysig_acat.csv", row.names = T)

system("pwd")
x <- intersect(rownames(pathway_maize_acat), rownames(pathway_sorghum_acat))
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6237005/
# Interactions of these structures (Phenylalanine) with membrane phospholipids
# have been suggested to cause cell toxicity, in particular in phenylketonuria, 
# by perturbing the phospholipid packing and compromising the membrane structural integrity





###############
#Saving pathways for all
write.csv(PBR_table_pvalue_maize,"PBR_table_pvalue_maize.csv", row.names = F)
write.csv(POL_table_pvalue_maize,"POL_table_pvalue_maize.csv", row.names = F)
write.csv(PNZ_table_pvalue_maize,"PNZ_table_pvalue_maize.csv", row.names = F)
write.csv(stp_table_pvalue_maize,"stp_table_pvalue_maize.csv", row.names = F)

write.csv(PBR_table_pvalue_sorghum,"PBR_table_pvalue_sorghum.csv", row.names = F)
write.csv(POL_table_pvalue_sorghum,"POL_table_pvalue_sorghum.csv", row.names = F)
write.csv(PNZ_table_pvalue_sorghum,"PNZ_table_pvalue_sorghum.csv", row.names = F)
write.csv(stp_table_pvalue_sorghum,"stp_table_pvalue_sorghum.csv", row.names = F)
