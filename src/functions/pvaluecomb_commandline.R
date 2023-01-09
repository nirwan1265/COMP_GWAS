#!/usr/bin/env Rscript 


#Load packages
library(methods)
library(argparser) 
library(optparse)
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
packages <- c("tidyverse","ggplot2", "Rsamtools","GenomicAlignments","rtracklayer","GenomicRanges","AnnotationHub","knitr","gtools","data.table","stringi","GBJ","metap","multtest","Hmisc","devtools","SNPRelate","gdsfmt","dplyr","vcfR","tidyr","AssocTests","SKAT","NCmisc","ACAT","PANTHER.db","UniProt.ws","ape","sp","rgdal","rworldmap","janitor","countrycode","tibble","vroom","gtools","tictoc")

suppressMessages(invisible(lapply(packages, library, character.only = TRUE)))


# OptionParser object
opt_parser <- OptionParser()


# Adding arguments
#p <- arg_parser("GBJ combination with OMNIBUS for GBJ, GHC, and SKAT")

# Add arguments
opt_parser$add_option(c("-p","--path"), action = "store",type = "character", default = NULL, help="Absolute path of the directory")
opt_parser$add_option(c("-f","--filename"), action = "store",type = "character",default = NULL, help="File name of the phenotype without the numbers")
opt_parser$add_option(c( "-n","--chr_no."), action = "store",type = "numeric", default = NULL, help="Number of chromosomes")
opt_parser$add_option(c("-o","--organism"), action = "store",type = "character", default = NULL, help="Scientific name of the organism")

#op <- add_option(op, c("-p","--path"), type = "character", help="Absolute path of the directory")

#op <- add_option(op, c("-f","--filename"), type = "character" help="File name of the phenotype without the numbers")

#op <- add_option(op, c( "-n","--chr_no."),type = "numeric" help="Number of chromosomes")

#op <- add_option(op, c("-o","--organism", type = "character", help="Scientific name of the organism")
                 
                 
# Parsing the command line arguments                 
opt <- parse_args(opt_parser)



# Load the script file with the main function
source("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/src/functions/gbj_test.R")

# Load helper Rscripts with functions
source("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/src/functions/data_wrangle.R")
source("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/src/functions/preprocess.R")
source("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/src/functions/geno_ld.R")
source("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/src/functions/pca_ld.R")
source("/Users/nirwatandukarn/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/src/functions/split_names.R")
source("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/COMP_GWAS/src/functions/zval.R")

# Call the functions 
geno_ld(chr=opt$chr)
preprocess(path=opt$path, phenonameopt$phenoname, n=opt$n, organism=opt$organism)
data_wrangle(path=opt$path, phenoname=opt$phenoname, chr=opt$chr, organism=opt$chromosome)
gbj_test(path=opt$path, phenonameopt$phenoname, chr=opt$chr, organism=pt$organism)


# Extract the values of the arguments
#path <- get_option(op, "path")
#filename <- get_option(op,"filename")
#n <- get_option(op,"number")
#organism <- get_option(op,"organism")


# calling the sccript
#test(path, filename, n, organism)