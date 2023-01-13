
list.of.packages <- c(
  "foreach",
  "doParallel",
  "ranger",
  "palmerpenguins",
  "tidyverse",
  "kableExtra"
  
)

for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}


# Function for putting everything in a proper list
data_wrangle <- function(path, phenoname, chr, organism){
  preprocess_data <- preprocess(path, phenoname, chr,  organism)
  geno <- geno <- geno_ld(chr)
  pca <- pca_ld()
  return(list(preprocess = preprocess_data,genotype = geno,PCA = pca))
}

gbj_test <- function(path, phenoname, chr, organism){
  #Empty results list to store values for each chromosome
  results <- list()
  
  trial <- data_wrangle(path, phenoname, chr, organism)
  
  for(j in 1:chr){
    
    # Subset the required data
    zstat_df <- as.data.frame(trial$preprocess$Zstat[[j]]) # i
    #zstat_df <- zstat_df[1:100,1:100]
    pvalue_df <- as.data.frame(trial$preprocess$pvalue[[j]])
    #pvalue_df <- pvalue_df[1:100,1:100]
    marker_df <- as.data.frame(trial$preprocess$Marker[[j]])
    #marker_df <- marker_df[1:100,1:100]
    genotype_df <- as.data.frame(trial$genotype[[j]])
    
    ## Function for Subsetting columns with more than one element 
    subset_element <- function (x) length(na.omit(x)) > 1
    
    # Applying the function for each df
    subset_zstat <- sapply(zstat_df, subset_element)
    subset_pvalue <- sapply(pvalue_df, subset_element)
    subset_marker <- sapply(marker_df, subset_element)
    
    # Subsetting the data frame with more than 1 elements
    zstat_df <- zstat_df[, subset_zstat]
    pvalue_df <- pvalue_df[, subset_pvalue]
    marker_df <- marker_df[, subset_marker]
    
    # Subsetting reference genotype
    sub_refgeno <- function(y,z){
      y <- y[!is.na(y)]
      ref_genotype <- as.data.frame(z[y])
      return(ref_genotype)
    }
    
    # Use foreach to apply the function to each column in parallel to subset ref_genotype for each column
    ref_genotype = list()
    ref_genotype <- foreach(i = 1:ncol(marker_df), .combine = c) %dopar% {
      list(sub_refgeno(marker_df[,i], genotype_df))
    }
    
    # PCA 
    pca <- as.data.frame(trial$PCA)
    
    # Use foreach to apply the function to each column in parallel for calculating correlation matrix
    corr_mat = list()
    corr_mat <- foreach(i = 1:20, .combine = c) %dopar% {
      list(GBJ::estimate_ss_cor(ref_pcs = pca, ref_genotypes = as.data.frame(ref_genotype[i]), link_function = 'linear'))
    }
    
    # gbj_analysis = list()
    # gbj_analysis <- foreach(i = 1:20, .combine = c) %dopar% {
    #   list(GBJ::GBJ(test_stats = as.vector(unlist(na.omit(zstat_df[i]))), cor_mat=corr_mat[[i]])$GBJ_pvalue)
    # }
    # results[[j]] <- gbj_analysis
    # names(results)[[j]] <- paste0("chr",j)
    
    omni_analysis = list()
    omni_analysis <- foreach(i = 1:20, .combine = c) %dopar% {
      list(GBJ::OMNI_ss(test_stats = as.vector(unlist(na.omit(zstat_df[i]))), cor_mat=corr_mat[[i]], num_boots = 100)$OMNI_pvalue)
    }
    results[[j]] <- omni_analysis
    names(results)[[j]] <- paste0("chr",j)
    
    for(j in 1:length(zstat_df)){
      names(results[[i]][[j]]) <- colnames(zstat_df)[j]
    }
    
  }
  return(results)
}


#Load genotype
geno_ld <- function(chr){
  geno <- list()
  # it should be in the same folder, but the file is too large to push to git hub 
  #geno <- vroom("Genotype/sorghum/allchrom_africa_filtered.MAF.txt")
  
  #use the above
  for(i in 1:chr){
    assign(paste0("geno",i), vroom(paste0("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/raw/africa.filtered/c",i,"_MAF_sorghum.txt")))
    geno <- c(geno, list(get(paste0("geno",i))))
    names(geno)[i] <- paste0("geno",i)
  }
  return(geno)
}

# Load PCA
pca_ld <- function(pca){
  pca <- vroom("PCA/sorghum/sorghum_PCA.txt") 
  pca <- pca[,-1]
}


#Pre processing step before running GBJ

preprocess <- function(path, phenoname, n, organism){
  
  #Helper lists 
  zstat_list <- list()
  pvalue_list <- list()
  marker_list <- list()
  
  # Set working directory
  setwd(path)
  file_path = paste0(path,"/GWAS/sorghum/")
  
  # Preprocess for Sorghum
  if(organism == "Sorghum bicolor"){
    
    # Load Sorghum Genomic Ranges
    a <- 1
    for(i in sprintf("%02d",1:chr)){
      assign(paste0("gr.db",a), readRDS(paste0("GenomicRanges/sorghum/gr.db",i,".RDS")))
      a = a + 1
    }
    
    # Data processing before GBJ
    a <- 1
    file_list <- list.files(path = file_path, pattern = phenoname)
    for(i in 1:chr){#length(file_list)){
      assign(file_list[i], vroom(paste0(file_path,file_list[i])))
      d = as.data.frame(get(file_list[i]))
      d$zstat = unlist(apply(d,1,zval))
      names(d) <- c("Marker","chr","Start_Position","pvalue","Zvalue")
      d <- d %>% mutate_at(c('chr','Start_Position','pvalue','Zvalue'),as.numeric)
      d <- as.data.frame(d)
      assign(paste0("gr.q", i) , GRanges(seqnames = paste0("chr",sprintf("%02d",i)), ranges = IRanges(start = d[,"Start_Position"], width = 1, zstat = d[,"Zvalue"], Marker = d[,"Marker"],pvalue = d[,"pvalue"])))
      a = a + 1
      assign(paste0("common",i), as.data.frame(findOverlapPairs(get(paste0("gr.db",i)), get(paste0("gr.q",i)))))
      e = get(paste0("common",i))
      e = e[which(e$first.X.Region == "gene"), ]
      e = e[,c(7,16,17,18)]
      colnames(e) = c("Gene","Marker","pvalue","zstat")
      assign(paste0("filter_common", i), e)
      assign(paste0("zstat",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "zstat"), value.var = "zstat"))
      assign(paste0("Marker",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "Marker"), value.var = "Marker"))
      assign(paste0("pvalue",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "pvalue"), value.var = "pvalue"))
      assign(paste0("genename",i),apply(get(paste0("Marker",i)),1,split.names))
      f <- get(paste0("zstat",i))
      f[,1]<- get(paste0("genename",i))
      f <- as.data.frame(t(f))
      colnames(f) <- f[1,]
      f <- f[-1,]
      f <- f %>% mutate_if(is.character,as.numeric, na.rm = T)
      f <- f[mixedsort(row.names(f)), ]
      assign(paste0("zstat",i),f)
      g <- get(paste0("Marker",i))
      g[,1]<- get(paste0("genename",i))
      g <- as.data.frame(t(g))
      colnames(g) <- g[1,]
      g <- g[-1,]
      g <- g[mixedsort(row.names(g)), ]
      assign(paste0("Marker",i),g)
      #return (get(paste0("Marker",i)))
      h <- get(paste0("pvalue",i))
      h[,1]<- get(paste0("genename",i))
      h <- as.data.frame(t(h))
      colnames(h) <- h[1,]
      h <- h[-1,]
      h <- h %>% mutate_if(is.character,as.numeric, na.rm = T)
      h <- h[mixedsort(row.names(h)), ]
      assign(paste0("pvalue",i),h)
    }
  }
  else if(organism == "Zea mays"){
    print("hello")
  }
  zstat_list = list()
  for (i in 1:chr){
    zstat_list <- c(zstat_list, list(get(paste0("zstat",i))))
    names(zstat_list)[i] <- paste0("zstat",i)
  }
  for (i in 1:chr){
    pvalue_list <- c(pvalue_list, list(get(paste0("pvalue",i))))
    names(pvalue_list)[i] <- paste0("pvalue",i)
  }
  for (i in 1:chr){
    marker_list <- c(marker_list, list(get(paste0("Marker",i))))
    names(marker_list)[i] <- paste0("markers",i)
  }
  return (list(Zstat = zstat_list, pvalue = pvalue_list,Marker = marker_list))
}


# Package names
packages <- c("tidyverse","ggplot2", "Rsamtools","GenomicAlignments","rtracklayer","GenomicRanges","AnnotationHub","knitr","gtools","data.table","stringi","GBJ","metap","multtest","Hmisc","devtools","SNPRelate","gdsfmt","dplyr","vcfR","tidyr","AssocTests","SKAT","NCmisc","ACAT","PANTHER.db","UniProt.ws","ape","sp","rgdal","rworldmap","janitor","countrycode","tibble","vroom","gtools","tictoc")

# Install packages not yet installed
# installed_packages <- packages %in% rownames(installed.packages())
# if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages])
# }
# if (any(installed_packages == FALSE)) {
#   BiocManager::install(packages[!installed_packages])
# }
# devtools::install_github("yaowuliu/ACAT")
# remotes::install_github("yaowuliu/ACAT")
# install.packages("/Users/nirwantandukar/Documents/Sorghum root rnaseq data_low phosphorus/org.Sbicolor.eg.db", repos=NULL, type="source")


# Packages loading
suppressMessages(invisible(lapply(packages, library, character.only = TRUE)))


# Subsetting just the genename from the gff3 files

split.names <- function(x,split){
  split.genename <- unlist(strsplit(x, split = ';', fixed = TRUE))[1]
  split.genename2 <- unlist(strsplit(split.genename, split = ":", fixed = TRUE))[2]
  return(split.genename2)
}


# Converting pvalues to Z-scores
zval <- function(x, output){
  pvalue <- unlist(as.numeric(x[4]))
  o <- p.to.Z(pvalue)
  return(o)
}

split.names_m <- function(x,split){
  split.genename <- unlist(strsplit(x, split = ';', fixed = TRUE))[1]
  split.genename2 <- unlist(strsplit(split.genename, split = "=", fixed = TRUE))[2]
  return(split.genename2)
}
