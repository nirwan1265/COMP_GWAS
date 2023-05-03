dir <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/GWAS-TASSEL/data/sorghum/"
data <- read.table(paste0(dir,"sorghum_4pheno_30PCA.txt"), header = T) %>% dplyr::select(-1)
str(data)
# Extract the phenotype columns and the PCA columns
phenotypes <- data[,1:4]
PCAs <- data[,5:ncol(data)]

# create a vector with the names of the phenotypes in your data
pheno_names <- c("stp10", "PBR1", "PNZ1", "POL1")  # replace with your actual phenotypes

# create an empty data frame to store the results
results <- data.frame(phenotype = character(),
                      num_PCA = integer())

# loop over the phenotypes
for (pheno in pheno_names) {
  # extract the phenotype vector and the PCA columns from the data
  pheno_vec <- data[, pheno]
  PCA_cols <- colnames(data)[grepl("^EV", colnames(data))]  # assuming your PCA columns start with "EV"
  
  # loop over the PCA columns
  for (PC_col in PCA_cols) {
    # combine the phenotype and the PCA into a data frame
    data_subset <- data.frame(pheno = pheno_vec, PCA = data[, PC_col])
    # run a linear regression
    model <- lm(pheno ~ PCA, data = data_subset)
    # extract the p-value from the model summary
    p_value <- summary(model)$coefficients[2, 4]
    # check if the p-value is significant
    if (p_value > 0.05) {
      # if the p-value is not significant, record the number of PCAs required
      num_PCA <- as.integer(gsub("EV", "", PC_col))
      break  # exit the loop over the PCA columns
    } else {
      num_PCA <- length(PCA_cols)  # if all PCAs are significant, record the maximum number of PCAs
    }
  }
  # add the results for the current phenotype to the results data frame
  results <- rbind(results, data.frame(phenotype = pheno, num_PCA = num_PCA))
}

# print the results
results



# loop through each phenotype
results <- data.frame() # initialize results dataframe
for (pheno in 1:4) { # assuming phenotype columns are 1 to 4
  pheno_name <- colnames(data)[pheno] # get phenotype column name
  sig_pcas <- c() # initialize vector to store significant PCs
  # loop through each PCA
  for (pc in 5:34) { # assuming PCA columns are 5 to 34
    pc_name <- colnames(data)[pc] # get PCA column name
    lm_model <- lm(data[,pheno] ~ data[,pc]) # run linear regression
    p_val <- summary(lm_model)$coefficients[2,4] # get p-value
    if (p_val > 0.05) { # if p-value is not significant
      break # exit loop
    } else {
      sig_pcas <- c(sig_pcas, pc_name) # add significant PCA to vector
    }
  }
  # add results for this phenotype to results dataframe
  results <- rbind(results, data.frame(phenotype = pheno_name, 
                                       significant_pcas = paste(sig_pcas, collapse = ", ")))
}
# print results
results


# loop through each phenotype
results <- data.frame() # initialize results dataframe
for (pheno in 1:4) { # assuming phenotype columns are 1 to 4
  pheno_name <- colnames(data)[pheno] # get phenotype column name
  pheno_data <- data[, c(pheno, 5:34)] # select phenotype and PCA columns
  pheno_data <- na.omit(pheno_data) # remove rows with missing values
  sig_pcas <- c() # initialize vector to store significant PCs
  # loop through each PCA
  for (pc in 2:30) { # assuming PCA columns are 5 to 34
    pc_name <- colnames(pheno_data)[pc] # get PCA column name
    lm_model <- lm(pheno_data[,1] ~ pheno_data[,pc]) # run linear regression
    p_val <- summary(lm_model)$coefficients[2,4] # get p-value
    if (p_val > 0.05) { # if p-value is not significant
      break # exit loop
    } else {
      sig_pcas <- c(sig_pcas, pc_name) # add significant PCA to vector
    }
  }
  # add results for this phenotype to results dataframe
  results <- rbind(results, data.frame(phenotype = pheno_name, 
                                       significant_pcas = paste(sig_pcas, collapse = ", ")))
}
# print results
results

# For all PCAs

# create an empty data frame to store significant PCAs for each phenotype
significant_pcas <- data.frame(phenotype = character(),
                               pcas = character())

# loop through each phenotype
for (i in 1:4) {
  
  # create a subset of data with non-NA values for the current phenotype
  subset_data <- data[!is.na(data[,i]),]
  
  # extract the phenotype vector
  phenotype <- subset_data[,i]
  
  # extract the PCA matrix
  pcas <- subset_data[,5:34]
  
  # loop through each PCA
  for (j in 1:ncol(pcas)) {
    
    # extract the current PCA vector
    pca <- pcas[,j]
    
    # fit a linear regression model
    model <- lm(phenotype ~ pca)
    
    # check if p-value is significant
    if (summary(model)$coef[2,4] < 0.05) {
      
      # add significant PCA to the data frame
      significant_pcas <- rbind(significant_pcas, 
                                data.frame(phenotype = colnames(data)[i], 
                                           pcas = paste0("PCA", j)))
      
    } else {
      
      # if p-value is not significant, break the loop for current phenotype
      break
      
    }
  }
}

# print the significant PCAs for each phenotype
print(significant_pcas)



