#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################
# LMM
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

################################################################################
# SORGHUM
################################################################################


################################################################################
# Phenotype for correlation 
################################################################################
dir <- ("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
pheno <- read.csv(paste0(dir,"/Sorghum_allphospho_africa.csv"))
pheno <- pheno[,-1]

################################################################################
# GWAS for comparison 
################################################################################
# Working directory
dir <- ("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/MAGMA/MAGMA_sorghum_2kb")

# List all the files
file_list <- list.files(dir)

# Loop to get all the gene names with pvalue <= 0.05
for (file_name in file_list){
  file_path <- paste0(dir,"/",file_name)
  file_data <- read.table(file_path, header = TRUE)
  
  # filter based on P_MULTI <= 0.05 and keep only the first column
  filtered_data <- as.data.frame(file_data[file_data$P_MULTI <= 0.05, 1])
  
  # assign the filtered data to a new variable with the same name
  assign(gsub("_magma_multi_snpwise.genes.out","",file_name), filtered_data)
}

rm(file_data)
rm(filtered_data)



# Getting the variables for the Phenotypes
#vars <- c("apa", "ext_P20", "lab", "NPlim",
#          "occ", "org", "PBR1", "PBR2", "PHO1", "PHO2", "PMEH1",
#          "PMEH2", "PNZ1", "PNZ2", "POL1", "POL2", "sec", "sol_Hi",
#          "sol_Lo", "sol", "sol_Mo", "sol_VL", "stp10", "stp100", 
#          "stp20", "stp30", "tot", "TP1", "TP2")
vars <- names(pheno)

result_2 <- data.frame(var1 = character(), length_var1 = numeric(), var2 = character(), length_var2 = numeric(),correlation_var1var2 = numeric(), common_count = numeric(), common_chars = character())

for (i in 1:(length(vars)-1)) {
  for (j in (i+1):length(vars)) {
    v1 <- get(vars[i])
    x <- nrow(v1)
    v2 <- get(vars[j])
    y <- nrow(v2)
    common <- base::intersect(v1$`file_data[file_data$P_MULTI <= 0.05, 1]`, v2$`file_data[file_data$P_MULTI <= 0.05, 1]`)
    common_count <- length(common)
    if (common_count > 0) {
      correlation <- cor(pheno[, vars[i]], pheno[, vars[j]], use = "pairwise.complete.obs")
      result_2 <- rbind(result_2, data.frame(var1 = vars[i],length_var1 = x, var2 = vars[j], length_var2 = y,correlation_var1var2 = correlation,common_count = common_count, common_chars = paste(common, collapse = ", ")))
    }
  }
}

result_2<- result_2[order(-result_2$common_count), ]

system("pwd")
write.csv(result_2, "common_sig_genes_LMM_sorghum_PPhenotype_genes.csv")

summary(result_2)



####
# Intersect of 3 genes


result_3 <- data.frame(var1 = character(), length_var1 = numeric(), var2 = character(), length_var2 = numeric(),var3 = character(), length_var3 = numeric(), common_count = numeric(), common_chars = character())

for (i in 1:(length(vars)-2)) {
  for (j in (i+1):(length(vars)-1)) {
    for(k in (i+2):length(vars)){
      v1 <- get(vars[i])
      x <- nrow(v1)
      v2 <- get(vars[j])
      y <- nrow(v2)
      v3 <- get(vars[k])
      z <- nrow(v3)
      common <- Reduce(base::intersect,list(v1$`file_data[file_data$P_MULTI <= 0.05, 1]`, v2$`file_data[file_data$P_MULTI <= 0.05, 1]`,v3$`file_data[file_data$P_MULTI <= 0.05, 1]`))
      common_count <- length(common)
      if (common_count > 0) {
        result_3 <- rbind(result_3, data.frame(var1 = vars[i],length_var1 = x, var2 = vars[j], length_var2 = y, var3 = vars[k], length_var3 = z, common_count = common_count, common_chars = paste(common, collapse = ", ")))
        
      }
    }
  }
}

result_trial <- result_3[!duplicated(paste(pmin(result_3$var1, result_3$var2), pmax(result_3$var1, result_3$var2))), ]
result_trial <- result_trial[!duplicated(paste(pmin(result_trial$var1, result_trial$var3), pmax(result_trial$var1, result_trial$var3))), ]
result_trial <- result_trial[!duplicated(paste(pmin(result_trial$var2, result_trial$var3), pmax(result_trial$var2, result_trial$var3))), ]


write.csv(result_3, "3pairs_common_sig_genes_LMM_sorghum_PPhenotype_genes.csv",row.names = F)


####
# Intersect of 4 genes


result_4 <- data.frame(var1 = character(), length_var1 = numeric(), var2 = character(), length_var2 = numeric(),var3 = character(), length_var3 = numeric(), var4 = character(), length_var4 = numeric(),common_count = numeric(), common_chars = character())

for (i in 1:(length(vars)-3)) {
  for (j in (i+2):(length(vars)-2)) {
    for(k in (i+3):(length(vars)-1)){
      for(l in (i+4):length(vars)){
        v1 <- get(vars[i])
        x <- nrow(v1)
        v2 <- get(vars[j])
        y <- nrow(v2)
        v3 <- get(vars[k])
        z <- nrow(v3)
        v4 <- get(vars[l])
        a <- nrow(v4)
        common <- Reduce(base::intersect,list(v1$`file_data[file_data$P_MULTI <= 0.05, 1]`, v2$`file_data[file_data$P_MULTI <= 0.05, 1]`,v3$`file_data[file_data$P_MULTI <= 0.05, 1]`,v4$`file_data[file_data$P_MULTI <= 0.05, 1]`))
        common_count <- length(common)
        if (common_count > 0) {
          result_4 <- rbind(result_4, data.frame(var1 = vars[i],length_var1 = x, var2 = vars[j], length_var2 = y, var3 = vars[k], length_var3 = z, var4 = vars[l], length_var4 = a,common_count = common_count, common_chars = paste(common, collapse = ", ")))
        }
      }
    }
  }
}

result_trial <- result_4[!duplicated(paste(pmin(result_4$var1, result_4$var2), pmax(result_4$var1, result_4$var2))), ]
result_trial <- result_trial[!duplicated(paste(pmin(result_trial$var1, result_trial$var3), pmax(result_trial$var1, result_trial$var3))), ]
result_trial <- result_trial[!duplicated(paste(pmin(result_trial$var1, result_trial$var4), pmax(result_trial$var1, result_trial$var4))), ]
result_trial <- result_trial[!duplicated(paste(pmin(result_trial$var2, result_trial$var3), pmax(result_trial$var2, result_trial$var3))), ]
result_trial <- result_trial[!duplicated(paste(pmin(result_trial$var2, result_trial$var4), pmax(result_trial$var2, result_trial$var4))), ]
result_trial <- result_trial[!duplicated(paste(pmin(result_trial$var3, result_trial$var4), pmax(result_trial$var3, result_trial$var4))), ]



write.csv(result_3, "3pairs_common_sig_genes_LMM_sorghum_PPhenotype_genes.csv",row.names = F)






################################################################################
# MAIZE
################################################################################


################################################################################
# Phenotype for correlation 
################################################################################
dir <- ("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
pheno <- read.csv(paste0(dir,"/Maize_allphospho.csv"))
pheno <- pheno[,-1]

################################################################################
# GWAS for comparison 
################################################################################
# Working directory
dir <- ("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GWAS.results/maize_Output/v5/MAGMA_maize_2kb")

# List all the files
file_list <- list.files(dir)

# Loop to get all the gene names with pvalue <= 0.05
for (file_name in file_list){
  file_path <- paste0(dir,"/",file_name)
  file_data <- read.table(file_path, header = TRUE)
  
  # filter based on P_MULTI <= 0.05 and keep only the first column
  filtered_data <- as.data.frame(file_data[file_data$P_MULTI <= 0.05, 1])
  
  # assign the filtered data to a new variable with the same name
  assign(gsub("_magma_multi_snpwise.genes.out","",file_name), filtered_data)
}

rm(file_data)
rm(filtered_data)



# Getting the variables for the Phenotypes
#vars <- c("apa", "ext_P20", "lab", "NPlim",
#          "occ", "org", "PBR1", "PBR2", "PHO1", "PHO2", "PMEH1",
#          "PMEH2", "PNZ1", "PNZ2", "POL1", "POL2", "sec", "sol_Hi",
#          "sol_Lo", "sol", "sol_Mo", "sol_VL", "stp10", "stp100", 
#          "stp20", "stp30", "tot", "TP1", "TP2")
vars <- names(pheno)

result_2 <- data.frame(var1 = character(), length_var1 = numeric(), var2 = character(), length_var2 = numeric(),correlation_var1var2 = numeric(), common_count = numeric(), common_chars = character())

for (i in 1:(length(vars)-1)) {
  for (j in (i+1):length(vars)) {
    v1 <- get(vars[i])
    x <- nrow(v1)
    v2 <- get(vars[j])
    y <- nrow(v2)
    common <- base::intersect(v1$`file_data[file_data$P_MULTI <= 0.05, 1]`, v2$`file_data[file_data$P_MULTI <= 0.05, 1]`)
    common_count <- length(common)
    if (common_count > 0) {
      correlation <- cor(pheno[, vars[i]], pheno[, vars[j]], use = "pairwise.complete.obs")
      result_2 <- rbind(result_2, data.frame(var1 = vars[i],length_var1 = x, var2 = vars[j], length_var2 = y,correlation_var1var2 = correlation,common_count = common_count, common_chars = paste(common, collapse = ", ")))
    }
  }
}

result_2<- result_2[order(-result_2$common_count), ]

system("pwd")
write.csv(result_2, "common_sig_genes_LMM_sorghum_PPhenotype_genes.csv")

summary(result_2)



################################################################################
# MAIZE
################################################################################


################################################################################
# Phenotype for correlation 
################################################################################
dir <- ("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
pheno <- read.csv(paste0(dir,"/Maize_allphospho.csv"))
pheno <- pheno[,-1]

################################################################################
# GWAS for comparison 
################################################################################
# Working directory
dir <- ("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GWAS.results/maize_Output/v5/MAGMA_maize_2kb")

# List all the files
file_list <- list.files(dir)

# Loop to get all the gene names with pvalue <= 0.05
for (file_name in file_list){
  file_path <- paste0(dir,"/",file_name)
  file_data <- read.table(file_path, header = TRUE)
  
  # filter based on P_MULTI <= 0.05 and keep only the first column
  filtered_data <- as.data.frame(file_data[file_data$P_MULTI <= 0.05, 1])
  
  # assign the filtered data to a new variable with the same name
  assign(gsub("_magma_multi_snpwise.genes.out","",file_name), filtered_data)
}

rm(file_data)
rm(filtered_data)



# Getting the variables for the Phenotypes
#vars <- c("apa", "ext_P20", "lab", "NPlim",
#          "occ", "org", "PBR1", "PBR2", "PHO1", "PHO2", "PMEH1",
#          "PMEH2", "PNZ1", "PNZ2", "POL1", "POL2", "sec", "sol_Hi",
#          "sol_Lo", "sol", "sol_Mo", "sol_VL", "stp10", "stp100", 
#          "stp20", "stp30", "tot", "TP1", "TP2")
vars <- names(pheno)

result_2 <- data.frame(var1 = character(), length_var1 = numeric(), var2 = character(), length_var2 = numeric(),correlation_var1var2 = numeric(), common_count = numeric(), common_chars = character())

for (i in 1:(length(vars)-1)) {
  for (j in (i+1):length(vars)) {
    v1 <- get(vars[i])
    x <- nrow(v1)
    v2 <- get(vars[j])
    y <- nrow(v2)
    common <- base::intersect(v1$`file_data[file_data$P_MULTI <= 0.05, 1]`, v2$`file_data[file_data$P_MULTI <= 0.05, 1]`)
    common_count <- length(common)
    if (common_count > 0) {
      correlation <- cor(pheno[, vars[i]], pheno[, vars[j]], use = "pairwise.complete.obs")
      result_2 <- rbind(result_2, data.frame(var1 = vars[i],length_var1 = x, var2 = vars[j], length_var2 = y,correlation_var1var2 = correlation,common_count = common_count, common_chars = paste(common, collapse = ", ")))
    }
  }
}

result_2<- result_2[order(-result_2$common_count), ]

system("pwd")
write.csv(result_2, "common_sig_genes_LMM_maize_PPhenotype_genes.csv")

summary(result_2)






#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################
# GLM
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

################################################################################
# SORGHUM
################################################################################

################################################################################
# Phenotype for correlation 
################################################################################
dir <- ("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
pheno <- read.csv(paste0(dir,"/Sorghum_allphospho_africa.csv"))
pheno <- pheno[,-1]

################################################################################
# GWAS for comparison 
################################################################################
# Working directory
dir <- ("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GBJ/maize/GBJ_maize_2kb_v5_csv")

# List all the files
file_list <- list.files(dir)

# Loop to get all the gene names with pvalue <= 0.05
for (file_name in file_list){
  file_path <- paste0(dir,"/",file_name)
  file_data <- read.csv(file_path, header = TRUE)
  
  # filter based on P_MULTI <= 0.05 and keep only the first column
  filtered_data <- as.data.frame(file_data[file_data$pvalue <= 0.05, 1])
  
  # assign the filtered data to a new variable with the same name
  assign(gsub("_gbj_maize.csv","",file_name), filtered_data)
}

rm(file_data)
rm(filtered_data)



# Getting the variables for the Phenotypes
#vars <- c("apa", "ext_P20", "lab", "NPlim",
#          "occ", "org", "PBR1", "PBR2", "PHO1", "PHO2", "PMEH1",
#          "PMEH2", "PNZ1", "PNZ2", "POL1", "POL2", "sec", "sol_Hi",
#          "sol_Lo", "sol", "sol_Mo", "sol_VL", "stp10", "stp100", 
#          "stp20", "stp30", "tot", "TP1", "TP2")
vars <- c("PNZ1","POL1","stp10","PBR1")

result_2 <- data.frame(var1 = character(), length_var1 = numeric(), var2 = character(), length_var2 = numeric(),correlation_var1var2 = numeric(), common_count = numeric(), common_chars = character())

for (i in 1:(length(vars)-1)) {
  for (j in (i+1):length(vars)) {
    v1 <- get(vars[i])
    x <- nrow(v1)
    v2 <- get(vars[j])
    y <- nrow(v2)
    common <- base::intersect(v1$`file_data[file_data$pvalue <= 0.05, 1]`, v2$`file_data[file_data$pvalue <= 0.05, 1]`)
    common_count <- length(common)
    if (common_count > 0) {
      correlation <- cor(pheno[, vars[i]], pheno[, vars[j]], use = "pairwise.complete.obs")
      result_2 <- rbind(result_2, data.frame(var1 = vars[i],length_var1 = x, var2 = vars[j], length_var2 = y,correlation_var1var2 = correlation,common_count = common_count, common_chars = paste(common, collapse = ", ")))
    }
  }
}

result_2<- result_2[order(-result_2$common_count), ]

system("pwd")
write.csv(result_2, "common_sig_genes_LMM_sorghum_PPhenotype_genes.csv")

summary(result_2)



####
# Intersect of 3 genes


result_3 <- data.frame(var1 = character(), length_var1 = numeric(), var2 = character(), length_var2 = numeric(),var3 = character(), length_var3 = numeric(), common_count = numeric(), common_chars = character())

for (i in 1:(length(vars)-2)) {
  for (j in (i+1):(length(vars)-1)) {
    for(k in (i+2):length(vars)){
      v1 <- get(vars[i])
      x <- nrow(v1)
      v2 <- get(vars[j])
      y <- nrow(v2)
      v3 <- get(vars[k])
      z <- nrow(v3)
      common <- Reduce(base::intersect,list(v1$`file_data[file_data$P_MULTI <= 0.05, 1]`, v2$`file_data[file_data$P_MULTI <= 0.05, 1]`,v3$`file_data[file_data$P_MULTI <= 0.05, 1]`))
      common_count <- length(common)
      if (common_count > 0) {
        result_3 <- rbind(result_3, data.frame(var1 = vars[i],length_var1 = x, var2 = vars[j], length_var2 = y, var3 = vars[k], length_var3 = z, common_count = common_count, common_chars = paste(common, collapse = ", ")))
        
      }
    }
  }
}

result_trial <- result_3[!duplicated(paste(pmin(result_3$var1, result_3$var2), pmax(result_3$var1, result_3$var2))), ]
result_trial <- result_trial[!duplicated(paste(pmin(result_trial$var1, result_trial$var3), pmax(result_trial$var1, result_trial$var3))), ]
result_trial <- result_trial[!duplicated(paste(pmin(result_trial$var2, result_trial$var3), pmax(result_trial$var2, result_trial$var3))), ]


write.csv(result_3, "3pairs_common_sig_genes_LMM_sorghum_PPhenotype_genes.csv",row.names = F)


####
# Intersect of 4 genes


result_4 <- data.frame(var1 = character(), length_var1 = numeric(), var2 = character(), length_var2 = numeric(),var3 = character(), length_var3 = numeric(), var4 = character(), length_var4 = numeric(),common_count = numeric(), common_chars = character())

for (i in 1:(length(vars)-3)) {
  for (j in (i+2):(length(vars)-2)) {
    for(k in (i+3):(length(vars)-1)){
      for(l in (i+4):length(vars)){
        v1 <- get(vars[i])
        x <- nrow(v1)
        v2 <- get(vars[j])
        y <- nrow(v2)
        v3 <- get(vars[k])
        z <- nrow(v3)
        v4 <- get(vars[l])
        a <- nrow(v4)
        common <- Reduce(base::intersect,list(v1$`file_data[file_data$P_MULTI <= 0.05, 1]`, v2$`file_data[file_data$P_MULTI <= 0.05, 1]`,v3$`file_data[file_data$P_MULTI <= 0.05, 1]`,v4$`file_data[file_data$P_MULTI <= 0.05, 1]`))
        common_count <- length(common)
        if (common_count > 0) {
          result_4 <- rbind(result_4, data.frame(var1 = vars[i],length_var1 = x, var2 = vars[j], length_var2 = y, var3 = vars[k], length_var3 = z, var4 = vars[l], length_var4 = a,common_count = common_count, common_chars = paste(common, collapse = ", ")))
        }
      }
    }
  }
}

result_trial <- result_4[!duplicated(paste(pmin(result_4$var1, result_4$var2), pmax(result_4$var1, result_4$var2))), ]
result_trial <- result_trial[!duplicated(paste(pmin(result_trial$var1, result_trial$var3), pmax(result_trial$var1, result_trial$var3))), ]
result_trial <- result_trial[!duplicated(paste(pmin(result_trial$var1, result_trial$var4), pmax(result_trial$var1, result_trial$var4))), ]
result_trial <- result_trial[!duplicated(paste(pmin(result_trial$var2, result_trial$var3), pmax(result_trial$var2, result_trial$var3))), ]
result_trial <- result_trial[!duplicated(paste(pmin(result_trial$var2, result_trial$var4), pmax(result_trial$var2, result_trial$var4))), ]
result_trial <- result_trial[!duplicated(paste(pmin(result_trial$var3, result_trial$var4), pmax(result_trial$var3, result_trial$var4))), ]



write.csv(result_3, "3pairs_common_sig_genes_LMM_sorghum_PPhenotype_genes.csv",row.names = F)






################################################################################
# MAIZE
################################################################################


################################################################################
# Phenotype for correlation 
################################################################################
dir <- ("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
pheno <- read.csv(paste0(dir,"/Maize_allphospho.csv"))
pheno <- pheno[,-1]

################################################################################
# GWAS for comparison 
################################################################################
# Working directory
dir <- ("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GWAS.results/maize_Output/v5/MAGMA_maize_2kb")

# List all the files
file_list <- list.files(dir)

# Loop to get all the gene names with pvalue <= 0.05
for (file_name in file_list){
  file_path <- paste0(dir,"/",file_name)
  file_data <- read.table(file_path, header = TRUE)
  
  # filter based on P_MULTI <= 0.05 and keep only the first column
  filtered_data <- as.data.frame(file_data[file_data$P_MULTI <= 0.05, 1])
  
  # assign the filtered data to a new variable with the same name
  assign(gsub("_magma_multi_snpwise.genes.out","",file_name), filtered_data)
}

rm(file_data)
rm(filtered_data)



# Getting the variables for the Phenotypes
#vars <- c("apa", "ext_P20", "lab", "NPlim",
#          "occ", "org", "PBR1", "PBR2", "PHO1", "PHO2", "PMEH1",
#          "PMEH2", "PNZ1", "PNZ2", "POL1", "POL2", "sec", "sol_Hi",
#          "sol_Lo", "sol", "sol_Mo", "sol_VL", "stp10", "stp100", 
#          "stp20", "stp30", "tot", "TP1", "TP2")
vars <- names(pheno)

result_2 <- data.frame(var1 = character(), length_var1 = numeric(), var2 = character(), length_var2 = numeric(),correlation_var1var2 = numeric(), common_count = numeric(), common_chars = character())

for (i in 1:(length(vars)-1)) {
  for (j in (i+1):length(vars)) {
    v1 <- get(vars[i])
    x <- nrow(v1)
    v2 <- get(vars[j])
    y <- nrow(v2)
    common <- base::intersect(v1$`file_data[file_data$P_MULTI <= 0.05, 1]`, v2$`file_data[file_data$P_MULTI <= 0.05, 1]`)
    common_count <- length(common)
    if (common_count > 0) {
      correlation <- cor(pheno[, vars[i]], pheno[, vars[j]], use = "pairwise.complete.obs")
      result_2 <- rbind(result_2, data.frame(var1 = vars[i],length_var1 = x, var2 = vars[j], length_var2 = y,correlation_var1var2 = correlation,common_count = common_count, common_chars = paste(common, collapse = ", ")))
    }
  }
}

result_2<- result_2[order(-result_2$common_count), ]

system("pwd")
write.csv(result_2, "common_sig_genes_LMM_sorghum_PPhenotype_genes.csv")

summary(result_2)



################################################################################
# MAIZE
################################################################################


################################################################################
# Phenotype for correlation 
################################################################################
dir <- ("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")
pheno <- read.csv(paste0(dir,"/Maize_allphospho.csv"))
pheno <- pheno[,-1]

################################################################################
# GWAS for comparison 
################################################################################
# Working directory
dir <- ("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GWAS.results/maize_Output/v5/MAGMA_maize_2kb")

# List all the files
file_list <- list.files(dir)

# Loop to get all the gene names with pvalue <= 0.05
for (file_name in file_list){
  file_path <- paste0(dir,"/",file_name)
  file_data <- read.table(file_path, header = TRUE)
  
  # filter based on P_MULTI <= 0.05 and keep only the first column
  filtered_data <- as.data.frame(file_data[file_data$P_MULTI <= 0.05, 1])
  
  # assign the filtered data to a new variable with the same name
  assign(gsub("_magma_multi_snpwise.genes.out","",file_name), filtered_data)
}

rm(file_data)
rm(filtered_data)



# Getting the variables for the Phenotypes
#vars <- c("apa", "ext_P20", "lab", "NPlim",
#          "occ", "org", "PBR1", "PBR2", "PHO1", "PHO2", "PMEH1",
#          "PMEH2", "PNZ1", "PNZ2", "POL1", "POL2", "sec", "sol_Hi",
#          "sol_Lo", "sol", "sol_Mo", "sol_VL", "stp10", "stp100", 
#          "stp20", "stp30", "tot", "TP1", "TP2")
vars <- names(pheno)

result_2 <- data.frame(var1 = character(), length_var1 = numeric(), var2 = character(), length_var2 = numeric(),correlation_var1var2 = numeric(), common_count = numeric(), common_chars = character())

for (i in 1:(length(vars)-1)) {
  for (j in (i+1):length(vars)) {
    v1 <- get(vars[i])
    x <- nrow(v1)
    v2 <- get(vars[j])
    y <- nrow(v2)
    common <- base::intersect(v1$`file_data[file_data$P_MULTI <= 0.05, 1]`, v2$`file_data[file_data$P_MULTI <= 0.05, 1]`)
    common_count <- length(common)
    if (common_count > 0) {
      correlation <- cor(pheno[, vars[i]], pheno[, vars[j]], use = "pairwise.complete.obs")
      result_2 <- rbind(result_2, data.frame(var1 = vars[i],length_var1 = x, var2 = vars[j], length_var2 = y,correlation_var1var2 = correlation,common_count = common_count, common_chars = paste(common, collapse = ", ")))
    }
  }
}

result_2<- result_2[order(-result_2$common_count), ]

system("pwd")
write.csv(result_2, "common_sig_genes_LMM_maize_PPhenotype_genes.csv")

summary(result_2)


