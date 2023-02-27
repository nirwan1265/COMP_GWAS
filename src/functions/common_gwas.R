# Working directory
dir <- ("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GWAS.results/sorghum_Output/v3/MAGMA_sorghum_2kb")

# List all the files
file_list <- list.files(dir)

# Loop to get all the gene names with pvalue <= 0.05
for (file_name in file_list){
  file_path <- paste0(dir,"/",file_name)
  file_data <- read.table(file_path, header = TRUE)
  
  # filter based on P_MULTI <= 0.05 and keep only the first column
  filtered_data <- as.data.frame(file_data[file_data$P_MULTI <= 0.05, 1])
  
  # assign the filtered data to a new variable with the same name
  assign(gsub("_multi_snpwise.genes.out","",file_name), filtered_data)
}

rm(file_data)
rm(filtered_data)

ls()

vars <- c("apa_magma", "ext_P20_magma", "lab_magma", "NPlim_magma",
          "occ_magma", "org_magma", "PBR1_magma", "PBR2_magma", "PHO1_magma", "PHO2_magma", "PMEH1_magma",
          "PMEH2_magma", "PNZ1_magma", "PNZ2_magma", "POL1_magma", "POL2_magma", "sec_magma", "sol_Hi_magma",
          "sol_Lo_magma", "sol_magma", "sol_Mo_magma", "sol_VL_magma", "stp10_magma", "stp100_magma", 
          "stp20_magma", "stp30_magma", "tot_magma", "TP1_magma", "TP2_magma")

result <- data.frame(var1 = character(), length_var1 = numeric(), var2 = character(), length_var2 = numeric(), common_count = numeric(), common_chars = character())

for (i in 1:(length(vars)-1)) {
  for (j in (i+1):length(vars)) {
    v1 <- get(vars[i])
    x <- nrow(v1)
    v2 <- get(vars[j])
    y <- nrow(v2)
    common <- base::intersect(v1$`file_data[file_data$P_MULTI <= 0.05, 1]`, v2$`file_data[file_data$P_MULTI <= 0.05, 1]`)
    common_count <- length(common)
    if (common_count > 0) {
      result <- rbind(result, data.frame(var1 = vars[i],length_var1 = x, var2 = vars[j], length_var2 = y,common_count = common_count, common_chars = paste(common, collapse = ", ")))
    }
  }
}

result <- result[order(-result$common_count), ]

result[400:406,1:6]

system("pwd")
write.csv(result, "common_sig_genes_LMM_sorghum_PPhenotype_genes.csv")

summary(result)

intersect(c(2,3,4),c(2,4),2)


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
