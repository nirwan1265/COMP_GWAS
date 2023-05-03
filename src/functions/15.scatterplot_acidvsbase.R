#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################################
# GBJ
################################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

################################################################################################
# MAIZE
################################################################################################

# Set working directory
gbj_result <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GBJ/"

# Maize
maize_gbj_res <- paste0(gbj_result,"/maize/GBJ_maize_2kb_v5_RDS/")
BR_maize <- readRDS(paste0(maize_gbj_res,"PBR1_gbj_maize.RDS")) %>% 
  dplyr::select("GeneName","pvalue") %>% 
  dplyr::rename(PBR1 = pvalue) %>% 
  dplyr::mutate(PBR1=-log10(PBR1))

OL_maize <- readRDS(paste0(maize_gbj_res,"POL1_gbj_maize.RDS")) %>% 
  dplyr::select("GeneName","pvalue") %>% 
  dplyr::rename(POL1 = pvalue) %>%
  dplyr::mutate(POL1=-log10(POL1))

PNZ_maize <- readRDS(paste0(maize_gbj_res,"PNZ1_gbj_maize.RDS")) %>% 
  dplyr::select("GeneName","pvalue") %>% 
  dplyr::rename(POL1 = pvalue) %>%
  dplyr::mutate(POL1=-log10(POL1))

STP_maize <- readRDS(paste0(maize_gbj_res,"stp10_gbj_maize.RDS")) %>% 
  dplyr::select("GeneName","pvalue") %>% 
  dplyr::rename(POL1 = pvalue) %>%
  dplyr::mutate(POL1=-log10(POL1))


# Common gene 
BR_OL_maize <- right_join(BR_maize,OL_maize, by = "GeneName")


# Plot
BR_OL_maize$legend <- ifelse(BR_OL_maize$PBR1 >= 2 & BR_OL_maize$POL1 >= 2, "Both significant",
                             ifelse(BR_OL_maize$PBR1 >= 2 & BR_OL_maize$POL1 < 2, "PBR1",
                                    ifelse(BR_OL_maize$POL1 >= 2 & BR_OL_maize$PBR1 < 2, "POL1", "No significance")))

gbj_maize <- ggplot(na.omit(BR_OL_maize), aes(x = PBR1, y = POL1, color = legend)) +
  geom_point() +
  scale_color_manual(values = c("grey", "purple", "blue", "red")) +
  labs(x = "-log10(PBR1)", y = "-log10(POL1)") + 
  ggtitle("Scatter plot of −log10 pvalues of PBR1 and POL1 for maize using GLM GWAS") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold")) 

################################################################################################
# SORGHUM
################################################################################################

# Set working directory
gbj_result <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GBJ/"

# sorghum
sorghum_gbj_res <- paste0(gbj_result,"/sorghum/GBJ_sorghum_2kb_v3_RDS/")
BR_sorghum <- readRDS(paste0(sorghum_gbj_res,"PBR1_omni_sorghum_2kbext.RDS")) %>% 
  dplyr::select("GeneName","pvalue") %>% 
  dplyr::rename(PBR1 = pvalue) %>% 
  dplyr::mutate(PBR1=-log10(PBR1))

OL_sorghum <- readRDS(paste0(sorghum_gbj_res,"POL1_omni_sorghum_2kbext.RDS")) %>% 
  dplyr::select("GeneName","pvalue") %>% 
  dplyr::rename(POL1 = pvalue) %>%
  dplyr::mutate(POL1=-log10(POL1))

PNZ_sorghum <- readRDS(paste0(sorghum_gbj_res,"PNZ1_omni_sorghum_2kbext.RDS")) %>% 
  dplyr::select("GeneName","pvalue") %>% 
  dplyr::rename(POL1 = pvalue) %>%
  dplyr::mutate(POL1=-log10(POL1))

STP_sorghum <- readRDS(paste0(sorghum_gbj_res,"stp10_omni_sorghum_2kbext.RDS")) %>% 
  dplyr::select("GeneName","pvalue") %>% 
  dplyr::rename(POL1 = pvalue) %>%
  dplyr::mutate(POL1=-log10(POL1))


# Common gene 
BR_OL_PNZ_STP_sorghum <- Reduce(function(x, y) merge(x, y, by = "GeneName"), list(BR_sorghum, OL_sorghum, PNZ_sorghum, STP_sorghum))


# Plot
quartz()
heatmap(as.matrix(BR_OL_PNZ_STP_sorghum[,-1]), Rowv = NA, Colv = NA, col = cm.colors(256))


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################################
# MAGMA
################################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

################################################################################################
# MAIZE
################################################################################################

# Set working directory
magma_result <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/MAGMA/"

# Maize
maize_magma_res <- paste0(magma_result,"maize/MAGMA_maize_2kb_v5_csv/")
BR_maize <- read.table(paste0(maize_magma_res,"PBR1_maize_magma_multi_snpwise.genes.out"), header = T) %>%
  dplyr::select("GENE","P_MULTI") %>% 
  dplyr::rename(PBR1 = P_MULTI) %>% 
  dplyr::mutate(PBR1=-log10(PBR1))

OL_maize <- read.table(paste0(maize_magma_res,"POL1_maize_magma_multi_snpwise.genes.out"), header = T) %>% 
  dplyr::select("GENE","P_MULTI") %>% 
  dplyr::rename(POL1 = P_MULTI) %>%
  dplyr::mutate(POL1=-log10(POL1))

# Common gene 
BR_OL_maize <- right_join(BR_maize,OL_maize, by = "GENE")


# Plot
BR_OL_maize$legend <- ifelse(BR_OL_maize$PBR1 >= 2 & BR_OL_maize$POL1 >= 2, "Both significant",
                             ifelse(BR_OL_maize$PBR1 >= 2 & BR_OL_maize$POL1 < 2, "PBR1",
                                    ifelse(BR_OL_maize$POL1 >= 2 & BR_OL_maize$PBR1 < 2, "POL1", "No significance")))
#quartz()
magma_maize <- ggplot(na.omit(BR_OL_maize), aes(x = PBR1, y = POL1, color = legend)) +
  geom_point() +
  scale_color_manual(values = c("red", "grey", "purple", "blue")) +
  labs(x = "-log10(PBR1)", y = "-log10(POL1)") + 
  ggtitle("Scatter plot of −log10 pvalues of PBR1 and POL1 for maize using LMM GWAS") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold")) 




################################################################################################
# SORGHUM
################################################################################################

# Set working directory
magma_result <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/MAGMA/"

# sorghum
sorghum_magma_res <- paste0(magma_result,"sorghum/MAGMA_sorghum_2kb_v3_csv/")
BR_sorghum <- read.table(paste0(sorghum_magma_res,"PBR1_magma_multi_snpwise.genes.out"), header = T) %>%
  dplyr::select("GENE","P_MULTI") %>% 
  dplyr::rename(PBR1 = P_MULTI) %>% 
  dplyr::mutate(PBR1=-log10(PBR1))

OL_sorghum <- read.table(paste0(sorghum_magma_res,"POL1_magma_multi_snpwise.genes.out"), header = T) %>% 
  dplyr::select("GENE","P_MULTI") %>% 
  dplyr::rename(POL1 = P_MULTI) %>%
  dplyr::mutate(POL1=-log10(POL1))

# Common gene 
BR_OL_sorghum <- right_join(BR_sorghum,OL_sorghum, by = "GENE")


# Plot
BR_OL_sorghum$legend <- ifelse(BR_OL_sorghum$PBR1 >= 2 & BR_OL_sorghum$POL1 >= 2, "Both significant",
                               ifelse(BR_OL_sorghum$PBR1 >= 2 & BR_OL_sorghum$POL1 < 2, "PBR1",
                                      ifelse(BR_OL_sorghum$POL1 >= 2 & BR_OL_sorghum$PBR1 < 2, "POL1", "No significance")))
quartz()
magma_sorghum <- ggplot(na.omit(BR_OL_sorghum), aes(x = PBR1, y = POL1, color = legend)) +
  geom_point() +
  scale_color_manual(values = c("red", "grey", "green", "blue")) +
  labs(x = "-log10(PBR1)", y = "-log10(POL1)") + 
  ggtitle("Acidic VS Basic properties of soil") + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16)) 



# Plotting
quartz()
plot(magma_sorghum)
grid.arrange(gbj_maize,gbj_sorghum,magma_maize,magma_sorghum)





