# Creating conda environment
# conda create -n orthofinder
# conda activate orthofinder
# conda install -c bioconda orthofinder 
# Download mcscan
# git clone https://github.com/wyp1125/MCScanX.git
# go to folder
# make
# chmod u+x all_files

# to open 
# open -na rstudio


#devtools::install_github("jtlovell/GENESPACE", upgrade = F)
#detach("package:GENESPACE", unload = TRUE)
#devtools::install_github("jtlovell/GENESPACE", upgrade = F)
library(GENESPACE)
library(Biostrings)
library(rtracklayer)


# Assmebly DOwnload
# go to : https://www.ncbi.nlm.nih.gov/assembly/GCF_902167145.1
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000003195.3
# Download genomic gff and translated CDS .faa file for maize and sorghum by clicking the download assembly button 
#setwd("~/Downloads/GENESPACE_data")

# Make bed files:
genomeRepo = "/Users/nirwantandukar/Downloads/GENESPACE/genomeRepo/"
wd = "/Users/nirwantandukar/Downloads/GENESPACE/"
path2mcscanx = "/Users/nirwantandukar/GitHub/MCScanX/"
genespaceWd = "/Users/nirwantandukar/Downloads/GENESPACE/workingDirectory"

# Parse
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = c("maize","sorghum"),
  genomeIDs = c("maize","sorghum"),
  gffString = "gtf",
  faString = "faa",
  presets = "ncbi",
  genespaceWd = genespaceWd)


# -- initalize the run and QC the inputs
gpar <- GENESPACE::init_genespace(
  wd = genespaceWd, 
  path2mcscanx = path2mcscanx)


out <- run_genespace(gpar, overwrite = T)

quartz()
ripd <- plot_riparian(
  gsParam = out,
  refGenome = "maize", 
  useRegions = FALSE)


hits <- read_allBlast(
  filepath = file.path(out$paths$syntenicHits, 
                       "maize_vs_sorghum.allBlast.txt.gz"))
quartz()
ggdotplot(hits = hits, type = "all", verbose = FALSE)

# Region of interest
roi <- data.frame(
  genome = c("maize","maize","sorghum"),
  chr = c("NC_001666.2","NC_001666.2","NC_008602.1"),
  start = c(88, 88,0),
  end = c(Inf, Inf,Inf))

qreturn <- query_hits(gsParam = out, bed = roi, synOnly = TRUE)

genome_blks <- ripd[["blks"]]
roibed <- roi[,c("genome", "chr")]
roibed$color <- c("pink", "cyan","green")
ripd <- plot_riparian(
  gsParam = out, 
  useRegions = FALSE, 
  highlightBed = roibed
  )

ripd <- plot_riparian(
  gsParam = out,
  refGenome = "sorghum", 
  useRegions = FALSE
  )

?plot_riparian()
