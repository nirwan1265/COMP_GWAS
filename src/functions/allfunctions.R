data_wrangle <- function(path, phenoname, chr, organism){
  preprocess_data <- preprocess(path, phenoname, chr,  organism)
  geno <- geno <- geno_ld(organism,chr)
  pca <- pca_ld(organism)
  return(list(preprocess = preprocess_data,genotype = geno,PCA = pca))
}

geno_ld <- function(organism,chr){
  geno <- list()
  # it should be in the same folder, but the file is too large to push to git hub
  #geno <- vroom("Genotype/sorghum/allchrom_africa_filtered.MAF.txt")
  if(organism == "Sorghum"){
    for(i in 1:chr){
      assign(paste0("geno",i), vroom(paste0("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Lasky.hapmap/raw/africa.filtered/c",i,"_MAF_sorghum.txt")))
      geno <- c(geno, list(get(paste0("geno",i))))
      names(geno)[i] <- paste0("geno",i)
    }
  }
  else {
    for(i in 1:chr){
      assign(paste0("geno",i), vroom(paste0("/Users/nirwantandukar/Documents/RubenLab/Maize/perchr/numerical_MAF_pheno_filtered/geno",i,".txt")))
      geno <- c(geno, list(get(paste0("geno",i))))
      names(geno)[i] <- paste0("geno",i)
    }
  }
  
  return(geno)
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
  #file_path = paste0(path,"/GWAS/maize/")
  # Preprocess for Sorghum
  if(organism == "Sorghum"){
    
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
      assign(paste0("gr.q", i) , GenomicRanges::GRanges(seqnames = paste0("chr",sprintf("%02d",i)), ranges = IRanges(start = d[,"Start_Position"], width = 1, zstat = d[,"Zvalue"], Marker = d[,"Marker"],pvalue = d[,"pvalue"])))
      a = a + 1
      assign(paste0("common",i), as.data.frame(IRanges::findOverlapPairs(get(paste0("gr.db",i)), get(paste0("gr.q",i)))))
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
  else if(organism == "Zea"){
    # Load maize Genomic Ranges
    a <- 1
    for(i in sprintf("%02d",1:chr)){
      assign(paste0("gr.db",a), readRDS(paste0("GenomicRanges/maize/gr.db",i,".RDS")))
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
      assign(paste0("gr.q", i) ,GenomicRanges::GRanges(seqnames = paste0("chr",sprintf("%02d",i)), ranges = IRanges(start = d[,"Start_Position"], width = 1, zstat = d[,"Zvalue"], Marker = d[,"Marker"],pvalue = d[,"pvalue"])))
      a = a + 1
      assign(paste0("common",i), as.data.frame(IRanges::findOverlapPairs(get(paste0("gr.db",i)), get(paste0("gr.q",i)))))
      e = get(paste0("common",i))
      e = e[which(e$first.X.Region == "gene"), ]
      e = e[,c(7,16,17,18)]
      colnames(e) = c("Gene","Marker","pvalue","zstat")
      assign(paste0("filter_common", i), e)
      assign(paste0("zstat",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "zstat"), value.var = "zstat"))
      assign(paste0("Marker",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "Marker"), value.var = "Marker"))
      assign(paste0("pvalue",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "pvalue"), value.var = "pvalue"))
      assign(paste0("genename",i),apply(get(paste0("Marker",i)),1,split.names_m))
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


packages <- c("tidyverse","ggplot2", "Rsamtools","GenomicAlignments","rtracklayer","GenomicRanges","AnnotationHub","knitr","gtools","data.table","stringi","GBJ","metap","multtest","Hmisc","devtools","SNPRelate","gdsfmt","dplyr","vcfR","tidyr","AssocTests","SKAT","NCmisc","ACAT","PANTHER.db","UniProt.ws","ape","sp","rgdal","rworldmap","janitor","countrycode","tibble","vroom","gtools","tictoc")

suppressMessages(invisible(lapply(packages, library, character.only = TRUE)))

list.of.packages <- c(
  "foreach",
  "doParallel",
  "ranger",
  "palmerpenguins",
  "tidyverse",
  "kableExtra",
  "tictoc",
  "parallel"
)

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}


split.names <- function(x,split){
  if(organism == "Sorghum"){
    split.genename <- unlist(strsplit(x, split = ';', fixed = TRUE))[1]
    split.genename2 <- unlist(strsplit(split.genename, split = ":", fixed = TRUE))[2]
  }
  else {
    split.genename <- unlist(strsplit(x, split = ';', fixed = TRUE))[1]
    split.genename2 <- unlist(strsplit(split.genename, split = "=", fixed = TRUE))[2]
  }
  return(split.genename2)
}


zval <- function(x, output){
  pvalue <- unlist(as.numeric(x[4]))
  o <- p.to.Z(pvalue)
  return(o)
}

pca_ld <- function(organism){
  if (organism == "Sorghum"){
    pca <- vroom("PCA/sorghum/sorghum_PCA.txt") 
    pca <- pca[,-1]  
  }
  else{
    pca <- vroom("PCA/maize/maize_PCA.txt") 
    pca <- pca[,-1]  
  }
  return(pca)
}

