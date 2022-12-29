#!/usr/bin/env Rscript --vanilla

library(methods)
library(argparser)


p <- arg_parser("GBJ combination with OMNIBUS for GBJ, GHC, and SKAT")
# Add a positional argument
p <- add_argument(p, "--path", help="Absolute path of the directory", flag = TRUE)
# Add a flag
p <- add_argument(p, "--filename", help="File name of the phenotype without the numbers", flag=TRUE)
# Add another flag
p <- add_argument(p, "--n", help="Number of chromosomes", flag=TRUE)
# Add another flag
p <- add_argument(p, "--organism", help="Scientific name of the organism", flag=TRUE)

argv <- parse_args(p)
