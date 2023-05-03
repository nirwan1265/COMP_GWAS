
packages <- c("tidyverse","ggplot2","gplots", "Rsamtools","GenomicAlignments","rtracklayer","GenomicRanges","AnnotationHub","knitr","gtools","data.table","stringi","GBJ","metap","multtest","Hmisc","devtools","SNPRelate","gdsfmt","dplyr","vcfR","tidyr","AssocTests","SKAT","NCmisc","ACAT","PANTHER.db","UniProt.ws","ape","sp","rgdal","rworldmap","janitor","countrycode","tibble","vroom","gtools","tictoc","gridExtra","VennDiagram")

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

