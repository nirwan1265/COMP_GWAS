preprocess <- function(path, filename, n, organism){
  a <- 1
  for(i in sprintf("%02d",1:10)){
    assign(paste0("gr.db",a), readRDS(paste0("GenomicRanges/sorghum/gr.db",i,".RDS")))
    a = a + 1
    #print(gr.db01)
  }
  print(gr.db2)
  a <- 1
  file_path = paste0(path,"/GWAS/sorghum/")
  if(organism == "Sorghum bicolor"){
    file_list <- list.files(path = file_path, pattern = filename)
    #print(file_list)
     for(i in 1:length(file_list)){
       #print("\n")
       #print(i)
       #setwd(file_path)
       print(file_list[i])
       assign(file_list[i], vroom(paste0(file_path,file_list[i])))
       d = as.data.frame(get(file_list[i]))
       d$zstat = unlist(apply(d,1,zval))
       names(d) <- c("Marker","chr","Start_Position","pvalue","Zvalue")
       d <- d %>% mutate_at(c('chr','Start_Position','pvalue','Zvalue'),as.numeric)
       #print(d)
       assign(paste0("gr.q", i) , GRanges(seqnames = paste0("chr",i), ranges = IRanges(start = d[,"Start_Position"], width = 1, zstat = d[,"Zvalue"], Marker = d[,"Marker"],pvalue = d[,"pvalue"])))
       print(gr.q1)
       assign(paste0("common",i), as.data.frame(findOverlapPairs(get(paste0("gr.db",a)), get(paste0("gr.q",a)))))
       a = a + 1
       print(common01)
       # e = get(paste0("common",i))
       # e = e[which(e$first.X.Region == "gene"), ]
       # e = e[,c(7,16,17,18)]
       # colnames(e) = c("Gene","Marker","pvalue","zstat")
       # print(e[1:4,1:4])
       # assign(paste0("filter_common", i), e)
       # assign(paste0("zstat",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "zstat"), value.var = "zstat"))
       # assign(paste0("Marker",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "Marker"), value.var = "Marker"))
       # assign(paste0("pvalue",i), dcast(setDT(e), Gene~rowid(Gene, prefix = "pvalue"), value.var = "pvalue"))
       # assign(paste0("genename",i),apply(get(paste0("Marker",i)),1,split.names))
       # f <- get(paste0("zstat",i))
       # f[,1]<- get(paste0("genename",i))
       # f <- as.data.frame(t(f))
       # colnames(f) <- f[1,]
       # f <- f[-1,]
       # f <- f %>% mutate_if(is.character,as.numeric, na.rm = T)
       # #print(f[1:3,1:3])
       # assign(paste0("zstat",i),f)
       # #assign(i,f)
       # g <- get(paste0("Marker",i))
       # g[,1]<- get(paste0("genename",i))
       # g <- as.data.frame(t(g))
       # colnames(g) <- g[1,]
       # g <- g[-1,]
       # #print(g[1:3,1:3])
       # assign(paste0("Marker",i),g)
       # #assign(i,g)
       # h <- get(paste0("pvalue",i))
       # h[,1]<- get(paste0("genename",i))
       # h <- as.data.frame(t(h))
       # colnames(h) <- h[1,]
       # h <- h[-1,]
       # h <- h %>% mutate_if(is.character,as.numeric, na.rm = T)
       # print(h[1:3,1:3])
       # assign(paste0("pvalue",i),h)
    }
  }
}

#system("ls data/GenomicRanges/sorghum")
path = "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/scripts/data"
filename <- "tot"
organism <- "Sorghum bicolor"
x <- preprocess(path, filename, 10,  organism)


for(i in sprintf("%02d", as.double(1:length(file_list)))){
 print(typeof(i))
}
