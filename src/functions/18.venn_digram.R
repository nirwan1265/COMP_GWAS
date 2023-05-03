#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################
# VENN DIAGRAM
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################

################################################################################
################################################################################
# MAGMA
################################################################################
################################################################################

################################################################################
# MAIZE
################################################################################

# File directory
dir_maize <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/MAGMA/maize/MAGMA_maize_2kb_v5_csv/"

# Load MAGMA results
PBR_maize <- read.table(paste0(dir_maize,"PBR1_maize_magma_multi_snpwise.genes.out"), header = T) %>% dplyr::filter(P_MULTI < 0.05) %>% dplyr::select(1) %>% unlist()
PNZ_maize <- read.table(paste0(dir_maize,"PNZ1_maize_magma_multi_snpwise.genes.out"), header = T) %>% dplyr::filter(P_MULTI < 0.05) %>% dplyr::select(1) %>% unlist()
POL_maize <- read.table(paste0(dir_maize,"POL1_maize_magma_multi_snpwise.genes.out"), header = T) %>% dplyr::filter(P_MULTI < 0.05) %>% dplyr::select(1) %>% unlist()
stp_maize <- read.table(paste0(dir_maize,"stp10_maize_magma_multi_snpwise.genes.out"), header = T) %>% dplyr::filter(P_MULTI < 0.05) %>% dplyr::select(1) %>% unlist()

maize_list <- list(PBR_maize = unique(PBR_maize),
               PNZ_maize = unique(PNZ_maize),
               POL_maize = unique(POL_maize),
               stp_maize = unique(stp_maize))

# Plotting
col_list <- c("#249EA0", "#FCBF49", "#E82648", "#09C269")
venn.diagram(maize_list, filename = "Venn_maize_MAGMA.tiff", resolution = 300, imagetype = "tiff", col = "white",lty = "blank", alpha = 0.5, col_list = "col_list",
             fill = col_list, cex = 2, cat.cex = 1.5)



# Common genes in the sets
# Two-way intersections
PBR_PNZ_maize <- intersect(PBR_maize, PNZ_maize)
PBR_POL_maize <- intersect(PBR_maize, POL_maize)
PBR_stp_maize <- intersect(PBR_maize, stp_maize)
PNZ_POL_maize <- intersect(PNZ_maize, POL_maize)
PNZ_stp_maize <- intersect(PNZ_maize, stp_maize)
POL_stp_maize <- intersect(POL_maize, stp_maize)

# Three-way intersections
PBR_PNZ_POL_maize <- intersect(PBR_maize, intersect(PNZ_maize, POL_maize))
PBR_PNZ_stp_maize <- intersect(PBR_maize, intersect(PNZ_maize, stp_maize))
PBR_POL_stp_maize <- intersect(PBR_maize, intersect(POL_maize, stp_maize))
PNZ_POL_stp_maize <- intersect(PNZ_maize, intersect(POL_maize, stp_maize))

# Four-way intersection
PBR_PNZ_POL_stp_maize <- intersect(PBR_maize, intersect(PNZ_maize, intersect(POL_maize, stp_maize)))


combine_maize <- dplyr::bind_rows(data.frame(PBR_PNZ_maize),
                                  data.frame(PBR_POL_maize),
                                  data.frame(PBR_stp_maize),
                                  data.frame(PNZ_POL_maize),
                                  data.frame(PNZ_stp_maize),
                                  data.frame(POL_stp_maize),
                                  data.frame(PBR_PNZ_POL_maize),
                                  data.frame(PBR_PNZ_stp_maize),
                                  data.frame(PBR_POL_stp_maize),
                                  data.frame(PNZ_POL_stp_maize),
                                  data.frame(PBR_PNZ_POL_stp_maize))


# Arrange the columns putting NA/s at the bottom
for (i in seq_along(combine_maize)) {
  # Sort the column, moving NA values to the bottom
  combine_maize[[i]] <- sort(combine_maize[[i]], na.last = TRUE)
}


# Saving
write.csv(combine_maize,"Common_genes_Venn_maize.csv", row.names = F)

################################################################################
# SORGHUM
################################################################################

# File directory
dir_sorghum <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/MAGMA/sorghum/MAGMA_sorghum_2kb_v3_csv/"

# Load MAGMA results
PBR_sorghum <- read.table(paste0(dir_sorghum,"PBR1_magma_multi_snpwise.genes.out"), header = T) %>% dplyr::filter(P_MULTI < 0.05) %>% dplyr::select(1) %>% unlist()
PNZ_sorghum <- read.table(paste0(dir_sorghum,"PNZ1_magma_multi_snpwise.genes.out"), header = T) %>% dplyr::filter(P_MULTI < 0.05) %>% dplyr::select(1) %>% unlist()
POL_sorghum <- read.table(paste0(dir_sorghum,"POL1_magma_multi_snpwise.genes.out"), header = T) %>% dplyr::filter(P_MULTI < 0.05) %>% dplyr::select(1) %>% unlist()
stp_sorghum <- read.table(paste0(dir_sorghum,"stp10_magma_multi_snpwise.genes.out"), header = T) %>% dplyr::filter(P_MULTI < 0.05) %>% dplyr::select(1) %>% unlist()

sorghum_list <- list(PBR_sorghum = unique(PBR_sorghum),
                     PNZ_sorghum = unique(PNZ_sorghum),
                     POL_sorghum = unique(POL_sorghum),
                     stp_sorghum = unique(stp_sorghum))

# Plotting
col_list <- c("#249EA0", "#FCBF49", "#E82648", "#09C269")
venn.diagram(sorghum_list, filename = "Venn_sorghum_MAGMA.tiff", resolution = 300, imagetype = "tiff", col = "white",lty = "blank", alpha = 0.5, col_list = "col_list",
             fill = col_list, cex = 2, cat.cex = 1.5)



# Common genes in the sets
# Two-way intersections
PBR_PNZ_sorghum <- intersect(PBR_sorghum, PNZ_sorghum)
PBR_POL_sorghum <- intersect(PBR_sorghum, POL_sorghum)
PBR_stp_sorghum <- intersect(PBR_sorghum, stp_sorghum)
PNZ_POL_sorghum <- intersect(PNZ_sorghum, POL_sorghum)
PNZ_stp_sorghum <- intersect(PNZ_sorghum, stp_sorghum)
POL_stp_sorghum <- intersect(POL_sorghum, stp_sorghum)

# Three-way intersections
PBR_PNZ_POL_sorghum <- intersect(PBR_sorghum, intersect(PNZ_sorghum, POL_sorghum))
PBR_PNZ_stp_sorghum <- intersect(PBR_sorghum, intersect(PNZ_sorghum, stp_sorghum))
PBR_POL_stp_sorghum <- intersect(PBR_sorghum, intersect(POL_sorghum, stp_sorghum))
PNZ_POL_stp_sorghum <- intersect(PNZ_sorghum, intersect(POL_sorghum, stp_sorghum))

# Four-way intersection
PBR_PNZ_POL_stp_sorghum <- intersect(PBR_sorghum, intersect(PNZ_sorghum, intersect(POL_sorghum, stp_sorghum)))


combine_sorghum <- dplyr::bind_rows(data.frame(PBR_PNZ_sorghum),
                                    data.frame(PBR_POL_sorghum),
                                    data.frame(PBR_stp_sorghum),
                                    data.frame(PNZ_POL_sorghum),
                                    data.frame(PNZ_stp_sorghum),
                                    data.frame(POL_stp_sorghum),
                                    data.frame(PBR_PNZ_POL_sorghum),
                                    data.frame(PBR_PNZ_stp_sorghum),
                                    data.frame(PBR_POL_stp_sorghum),
                                    data.frame(PNZ_POL_stp_sorghum),
                                    data.frame(PBR_PNZ_POL_stp_sorghum))


# Arrange the columns putting NA/s at the bottom
for (i in seq_along(combine_sorghum)) {
  # Sort the column, moving NA values to the bottom
  combine_sorghum[[i]] <- sort(combine_sorghum[[i]], na.last = TRUE)
}

#Saving
write.csv(combine_sorghum,"Common_genes_Venn_sorghum.csv", row.names = F)
