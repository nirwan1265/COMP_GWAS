preprocess <- function(path, filename, n, organism){
  # Set working directory
  setwd(path)
  #a <- 1
  file_path = paste0(path,"/GWAS/sorghum/")
  if(organism == "Sorghum bicolor"){
    # Load Sorghum Genomic Ranges
    a <- 1
    for(i in sprintf("%02d",1:10)){
      assign(paste0("gr.db",a), readRDS(paste0("GenomicRanges/sorghum/gr.db",i,".RDS")))
      a = a + 1
    }
    
    # Data processing before GBJ
    a <- 1
    file_list <- list.files(path = file_path, pattern = filename)
    for(i in 1:length(file_list)){
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
       return (get(paste0("Marker",i)))
       h <- get(paste0("pvalue",i))
       h[,1]<- get(paste0("genename",i))
       h <- as.data.frame(t(h))
       colnames(h) <- h[1,]
       h <- h[-1,]
       h <- h %>% mutate_if(is.character,as.numeric, na.rm = T)
       h <- h[mixedsort(row.names(h)), ]
       assign(paste0("pvalue",i),h)
       #print(get(paste0("pvalue",i)))
       #return(get(paste0("pvalue",i)))
    }
  }
}

#system("ls data/GenomicRanges/sorghum")
path = "/Users/nirwan/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/data"
filename <- "tot"
pca <- "pca"
organism <- "Sorghum bicolor"
chr <- 1
preprocess(path, filename, chr,  organism)

x[1:5,1:5]

#Running the Analysis:
#Using the pvalue.combination function:
for(i in sprintf("%02d", 1:10)){
  assign(paste0("pvalue.combine",i), pvalue.combine(get(paste0("gwas",i,".zstat")), get(paste0("gwas",i,".Marker")), get(paste0("gwas",i,".pvalue")), get(paste0("geno",i)), get(paste0("tab",i))))
}
