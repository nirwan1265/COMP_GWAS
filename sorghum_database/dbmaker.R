library(knitr)
library(AnnotationForge)


options(width=80)
hook_output = knit_hooks$get('output')
knit_hooks$set(output = function(x, options) {
  # this hook is used only when the linewidth option is not NULL
  if (!is.null(n <- options$linewidth)) {
    x = knitr:::split_lines(x)
    # any lines wider than n should be wrapped
    if (any(nchar(x) > n)) x = strwrap(x, width = n)
    x = paste(x, collapse = '\n')
  }
  hook_output(x, options)
})
#browseVignettes("AnnotationForge")

#setwd("~/Desktop/Practice/Sorghum root rnaseq data_low phosphorus")

setwd("/Users/nirwantandukar/Documents/Sorghum_db")


library(BiocManager)
library(GenomicRanges)
library(AnnotationFilter)
library(AnnotationHub)
library(AnnotationForge)
library(GenomeInfoDb)
library(AnnotationDbi)
library(biomaRt)
library(gower)
library(rlang)
library(tidyr)

makeOrgPackageFromNCBI(version = "0.1",
                       author = "Ntanduk",
                       maintainer = "Ntanduk <ntanduk@ncsu.edu>",
                       outputDir = ".",
                       tax_id = "4558",
                       genus = "Sorghum",
                       species = "bicolor")
?makeOrgPackageFromNCBI()

install.packages("./org.Sbicolor.eg.db", repos=NULL, type="source")
library(org.Sbicolor.eg.db)

R.version


#Annotation hub
hub <- AnnotationHub()
sb <- query(hub, c("sorghum bicolor","orgdb"))
sb <- query(hub, c("sorghum bicolor","NCBI"))
snapshotDate(sb)
organism("Sorghum bicolor", sb)

orgdb_sb <- hub[["AH101176"]]
orgdb_sb$packageName

