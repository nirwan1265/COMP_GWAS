library(dplyr)
library(Biostrings)
library(stringr)

# Maize
maize <- read.table("/Users/nirwantandukar/Downloads/GENESPACE/workingDirectory/maize/maize.bed", sep ="\t")
maize_bed <- maize %>% dplyr::select(V2,V3,V8,V10) %>% filter(V8 == "start_codon")
maize_bed$V10[1]
maize_bed$V10 <- stringr::str_extract(maize_bed$V10, "(?<=protein_id )[^\\s;]+")
maize_bed <- maize_bed[,-3]
maize_bed <- maize_bed[,c(3,1,2)]

maize_bed <- maize_bed %>%
  dplyr::mutate(chr = cumsum(V2 < dplyr::lag(V2, default = 0))) %>%
  dplyr::mutate(chr = paste0("chr",chr)) %>%
  dplyr::select(chr,V2,V3,V10)

# getting fasta files
maize_fasta <- readDNAStringSet("/Users/nirwantandukar/Downloads/GENESPACE/workingDirectory/peptide/maize.fa")
# extract names
maize_names <- names(maize_fasta)

# Create a data frame with the names
names_df <- data.frame(names = maize_names)
# extract common name and version number
names_df <- names_df %>%
  mutate(common_name = str_extract(names, ".*(?=_)"),
         version_num = str_extract(names, "(?<=_)\\w+$"))

# remove rows where version_num is missing
names_df <- names_df %>% filter(!is.na(version_num))

# reorder columns
maize_grouped <- names_df[, c("common_name", "names")]
colnames(maize_grouped) <- c("common_name","multiple_names")



maize_bed2 <- maize_bed %>%
  left_join(maize_grouped, by = c("V9" = "common_name")) %>%
  mutate(V9 = ifelse(!is.na(V4), multiple_names, common_name))
maize_bed2 <- maize_bed2[,-5]


write.table(maize_bed2,"maize.bed",row.names = F,quote = F)
system("pwd")

## SORGHUM
# Getting gff3 file for bed format
sorghum <- read.table("/Users/nirwantandukar/Downloads/GENESPACE/data/sorghum_phytozome/sorghum_v3_gene.gff3")
sorghum_bed <- sorghum %>% dplyr::select(V1,V3,V4,V5,V9) %>% filter(V3 == "gene")
sorghum_bed$V9 <- str_extract(sorghum_bed$V9, "(?<=ID=)[^;]+")
sorghum_bed <- sorghum_bed[,-2]
sorghum_bed$V1 <- str_replace(sorghum_bed$V1, "Chr", "chr")
sorghum_bed$V1 <- str_replace(sorghum_bed$V1, "01", "1")
sorghum_bed$V1 <- str_replace(sorghum_bed$V1, "02", "2")
sorghum_bed$V1 <- str_replace(sorghum_bed$V1, "03", "3")
sorghum_bed$V1 <- str_replace(sorghum_bed$V1, "04", "4")
sorghum_bed$V1 <- str_replace(sorghum_bed$V1, "05", "5")
sorghum_bed$V1 <- str_replace(sorghum_bed$V1, "06", "6")
sorghum_bed$V1 <- str_replace(sorghum_bed$V1, "07", "7")
sorghum_bed$V1 <- str_replace(sorghum_bed$V1, "08", "8")
sorghum_bed$V1 <- str_replace(sorghum_bed$V1, "09", "9")
sorghum_bed$V9 <- str_replace(sorghum_bed$V9, ".v3.2", "")


# getting fasta files
sorghum_fasta <- readDNAStringSet("/Users/nirwantandukar/Downloads/GENESPACE/workingDirectory/peptide/sorghum.fa")
# extract names
sorghum_names <- names(sorghum_fasta)

# Create a data frame with the names
names_df <- data.frame(names = sorghum_names)

# Extract the common name and the version number
names_df <- names_df %>%
  dplyr::mutate(common_name = str_extract(names, "^\\w+\\.\\w+\\d+"),
         version_num = str_extract(names, "\\d+\\.p$"))

# Group by the common name and concatenate the multiple names
sorghum_grouped <- names_df %>%
  dplyr::group_by(common_name) %>%
  dplyr::summarize(multiple_names = paste(names, collapse = ", ")) %>%
  tidyr::separate_rows(multiple_names,sep = ",")

sorghum_bed2 <- sorghum_bed %>%
  left_join(sorghum_grouped, by = c("V9" = "common_name")) %>%
  mutate(V9 = ifelse(!is.na(V4), multiple_names, common_name))
sorghum_bed2 <- sorghum_bed2[,-5]


system("pwd")
write.table(sorghum_bed2,"sorghum.bed",row.names = F,quote = F)
