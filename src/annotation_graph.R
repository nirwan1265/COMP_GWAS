library(dplyr)
library(magrittr)
library(ggplot2)

dir <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/Annotation/"
sorg_PBR_BP <- read.table(paste0(dir,"sorghum_PBR1_BP.txt"), sep ="\t") %>% dplyr::select(c(2,3,5)) %>% magrittr::set_names(c("Biological_Process","No._of_genes","Percentage"))
sorg_PNZ_BP <- read.table(paste0(dir,"sorghum_PNZ1_BP.txt"), sep ="\t") %>% dplyr::select(c(2,3,5)) %>% magrittr::set_names(c("Biological_Process","No._of_genes","Percentage"))
sorg_POL_BP <- read.table(paste0(dir,"sorghum_POL1_BP.txt"), sep ="\t") %>% dplyr::select(c(2,3,5)) %>% magrittr::set_names(c("Biological_Process","No._of_genes","Percentage"))
sorg_stp_BP <- read.table(paste0(dir,"sorghum_stp10_BP.txt"), sep ="\t") %>% dplyr::select(c(2,3,5)) %>% magrittr::set_names(c("Biological_Process","No._of_genes","Percentage"))

# Combine the data frames into a single data frame
df_combined <- rbind(sorg_PBR_BP, sorg_PNZ_BP, sorg_POL_BP, sorg_stp_BP)

#Add a column to indicate the source data frame
# Add a column to indicate the source data frame
df_combined$source <- rep(c("sorg_PBR_BP", "sorg_PNZ_BP", "sorg_POL_BP", "sorg_stp_BP"), times = c(nrow(sorg_PBR_BP), nrow(sorg_PNZ_BP), nrow(sorg_POL_BP), nrow(sorg_stp_BP)))

# Create the bar chart
ggplot(df_combined, aes(x = source, y = No._of_genes, fill = Biological_Process)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Data Frame", y = "Number of Genes", fill = "Biological Process") +
  scale_fill_discrete(name = "Biological Process") +
  theme_minimal()


ggplot(df_combined, aes(x = Biological_Process, y = No._of_genes, fill = source)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ row_group, nrow = 1) +
  labs(x = "Biological Process", y = "Number of Genes", fill = "Data Frame") +
  scale_fill_discrete(name = "Data Frame") +
  theme_minimal()


