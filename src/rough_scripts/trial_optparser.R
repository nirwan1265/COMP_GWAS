#!/usr/bin/env Rscript --vanilla

# # Load the optparse package
# library(optparse)
# 
# # Create an option object
# option <- make_option(c("-n", "--number"), type="character", help="number")
# 
# # Create an OptionParser object
# op <- OptionParser()
# 
# # Add the option to the OptionParser object
# op <- add_option(op, option)
# 
# # Parse the command line arguments
# op <- parse_args(op)
# 
# # Extract the value of the option
# number <- get_option(op, "number")
# 
# # Print the value of the option
# print(paste("number:", number))



library("optparse")
parser <- OptionParser()
parser <- add_option(parser, c("-v", "--verbose"), action="store_true",
                     default=TRUE, help="Print extra output [default]")
parser <- add_option(parser, c("-q", "--quietly"), action="store_false",
                     dest="verbose", help="Print little output")
parser <- add_option(parser, c("-c", "--count"), type="integer", default=5,
                     help="Number of random normals to generate [default %default]",
                     metavar="number")
parse_args(parser, args = c("--quietly", "--count=15"))
