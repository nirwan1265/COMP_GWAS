
# read in protein.faa file
protein <- readLines("/Users/nirwantandukar/Downloads/GENESPACE/genomeRepo/sorghum/protein.faa")

# create empty data frame
df <- NULL

# loop through lines of protein.faa file and add headers and sequences to data frame
for (i in 1:length(protein)) {
  if (startsWith(protein[i], ">")) {
    header <- protein[i]
    sequence <- ""
  } else {
    sequence <- paste0(sequence, protein[i])
  }
  if (i == length(protein) || startsWith(protein[i+1], ">")) {
    df <- rbind(df, data.frame(headers = header, sequences = sequence, stringsAsFactors = FALSE))
  }
}

protein <- df
protein$id <- gsub("^>\\s*(\\S+).*", "\\1", protein$headers)

# read in cds_from_genomic.fna file
cds <- readLines("/Users/nirwantandukar/Downloads/GENESPACE/genomeRepo/sorghum/cds_from_genomic.fna")

df <- NULL

# loop through lines of cds_from_genomic.fna file and add sequences to data frame based on protein_id match
for (i in 1:length(cds)) {
  if (startsWith(cds[i], ">")) {
    header <- cds[i]
    sequence <- ""
  } else {
    sequence <- paste0(sequence, cds[i])
  }
  if (i == length(cds) || startsWith(cds[i+1], ">")) {
    df <- rbind(df, data.frame(headers = header, sequences = sequence, stringsAsFactors = FALSE))
  }
}

cds = df
library(stringr)
cds$id <- str_extract(cds$headers, "(?<=\\[protein_id=)[^\\]]+")
library(dplyr)
common <- left_join(cds, protein, by = "id")
common <- common[,c(1,5)]
common_t <- t(common)
result <- c(common_t)

writeLines(result, "sorghum.fa")

