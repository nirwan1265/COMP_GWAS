#Package
library("CMplot")
library(vroom)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##########################################################################################
# LMM
##########################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

##########################################################################################
# MAIZE
##########################################################################################

#Need LMM GWAS file all assoc
# maize
dir <- c("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GWAS.results/maize_Output/v5/LMM_GWAS")
dir <- c("/Users/nirwan/Documents/Alvarez/results/gwas/")
# P Solubility
PBR1 <- vroom(paste0(dir,"PBR1.txt")) #%>% select(1,2,3,5) %>% na.omit()
PBR1 <- PBR1[,c(2,1,3,12)]

# P Retention
PNZ1 <- vroom(paste0(dir,"GLM_sorghum_PNZ1.txt"))# %>% select(1,2,3,5) %>% na.omit()
PNZ1 <- PNZ1[,c(2,1,3,12)]

# Total P
stp10 <- vroom(paste0(dir,"GLM_sorghum_stp10.txt")) #%>% select(1,2,3,5) %>% na.omit()
stp10 <- stp10[,c(2,1,3,12)]

# Alkaline soil P
POL1 <- vroom(paste0(dir,"GLM_sorghum_POL1.txt")) #%>% select(1,2,3,5) %>% na.omit()
POL1 <- POL1[,c(2,1,3,12)]


maize_sol_ret <- cbind(PBR1, PNZ1$p)
names(maize_sol_ret) <- c("SNP","Chromosome","Position","P Solubility","P Retention")

maize_tot_alk <- left_join(stp10, POL1, by = "Marker")
maize_tot_alk <- as.data.frame(maize_tot_alk[,c(1,2,3,4,7)])
names(maize_tot_alk) <- c("SNP","Chromosome","Position","Total P","Alkaline Soil P")

maize_all <- merge(maize_sol_ret,maize_tot_alk)
names(maize_all)

SNPs <- maize_all[
    maize_all$`P Solubility` < 1e-5 |
        maize_all$`P Retention` < 1e-5 |
        maize_all$`Total P` < 1e-5 |
        maize_all$`Alkaline Soil P` < 1e-5 ,1
]

quartz()
setwd("~/Desktop")
# CMplot(maize_all,type="p",plot.type="m",LOG10=TRUE,highlight=SNPs,highlight.type="l",
#        threshold=1e-4,threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),
#        signal.cex=1.2, signal.col="red", highlight.col="grey",highlight.cex=0.7,
#        file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)

CMplot(maize_all, plot.type="m",multracks=TRUE,threshold=c(1e-6,1e-4),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green"),
       signal.cex=0.5, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
       highlight=NULL, highlight.text=NULL, highlight.text.cex=1.4)



##########################################################################################
# SORGHUM
##########################################################################################

#Need LMM GWAS file all assoc
# Sorghum
dir <- c("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GWAS.results/sorghum_Output/v3/LMM_GWAS/GWAS_LMM_results")
system("ls")

# P Solubility
PBR1 <- vroom(paste0(dir,"/PBR1.txt"))
PBR1 <- PBR1[,c(2,1,3,12)]

# P Retention
PNZ1 <- vroom(paste0(dir,"/PNZ1.txt"))
PNZ1 <- PNZ1[,c(2,1,3,12)]

# Total P
stp10 <- vroom(paste0(dir,"/stp10.txt"))
stp10 <- stp10[,c(2,1,3,12)]

# Alkaline soil P
POL1 <- vroom(paste0(dir,"/POL1.txt"))
POL1 <- POL1[,c(2,1,3,12)]




# occ_ret <- vroom("occ_ret.txt")
# occ_ret <- occ_ret[,c(2,1,3,13)]
# tot_tot <- vroom("tot_tot.txt")
# tot_tot <- tot_tot[,c(2,1,3,13)]
# lab_sol <- vroom("lab_sol.txt")
# lab_sol <- lab_sol[,c(2,1,3,13)]
# 

# sorg_all <- cbind(occ_ret, tot_tot$p_wald, lab_sol$p_wald)
# names(sorg_all)[4:6] <- c("Retention","Total","Solubility")
# names(sorg_all)[1:3] <- c("SNP","Chromosome","Position")

sorg_sol_ret <- cbind(PBR1, PNZ1$p_wald)
names(sorg_sol_ret) <- c("SNP","Chromosome","Position","P Solubility","P Retention")

sorg_tot_alk <- left_join(stp10, POL1, by = "rs")
sorg_tot_alk <- as.data.frame(sorg_tot_alk[,c(1,2,3,4,7)])
names(sorg_tot_alk) <- c("SNP","Chromosome","Position","Total P","Alkaline Soil P")

sorg_all <- merge(sorg_sol_ret,sorg_tot_alk)
#sorg_all <- sorg_all[,-c(6,7,8)]
names(sorg_all)


# SNPs <- sorg_all[
#     sorg_all$Retention < 1e-5 |
#     sorg_all$Total < 1e-5 |
#     sorg_all$Solubility < 1e-5, 1
# ]

# SNPs_sol_ret <- sorg_sol_ret[
#     sorg_sol_ret$`P Solubility` < 1e-5 |
#         sorg_sol_ret$`P Retention` < 1e-5, 1
# ]

SNPs <- sorg_all[
    sorg_all$`P Solubility` < 1e-5 |
        sorg_all$`P Retention` < 1e-5 |
        sorg_all$`Total P` < 1e-5 |
        sorg_all$`Alkaline Soil P` < 1e-5 ,1
]

quartz()
setwd("~/Desktop")
# CMplot(sorg_all,type="p",plot.type="m",LOG10=TRUE,highlight=SNPs,highlight.type="l",
#        threshold=1e-4,threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),
#        signal.cex=1.2, signal.col="red", highlight.col="grey",highlight.cex=0.7,
#        file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)

CMplot(sorg_all, plot.type="m",multracks=TRUE,threshold=c(1e-5,1e-4),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green"),
       signal.cex=0.5, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
       highlight=NULL, highlight.text=NULL, highlight.text.cex=1.4)




#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##########################################################################################
# GLM
##########################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

##########################################################################################
# MAIZE
##########################################################################################

maize_dir <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GWAS.results/maize_Output/v5/GLM_GWAS/"
PBR1 <- vroom(paste0(maize_dir,"PBR1.txt"))
PBR1 <- PBR1[,c(2,1,3,5)]
PNZ1 <- vroom(paste0(maize_dir,"PNZ1.txt"))
PNZ1 <- PNZ1[,c(2,1,3,5)]
POL1 <- vroom(paste0(maize_dir,"POL1.txt"))
POL1 <- POL1[,c(2,1,3,5)]
stp10 <- vroom(paste0(maize_dir,"stp10.txt"))
stp10 <- stp10[,c(2,1,3,5)]


maize_sol_ret <- cbind(PBR1, PNZ1$p)
names(maize_sol_ret) <- c("Chromosome","SNP","Position","P Solubility","P Retention")
maize_sol_ret <- maize_sol_ret[,c(2,1,3,4,5)]

maize_tot_alk <- cbind(stp10, POL1$p)
names(maize_tot_alk) <- c("Chromosome","SNP","Position","Total P","Alkaline Soil P")
maize_tot_alk <- maize_tot_alk[,c(2,1,3,4,5)]

maize_all <- merge(maize_sol_ret,maize_tot_alk)
names(maize_all)

SNPs <- maize_all[
    maize_all$`P Solubility` < 1e-5 |
        maize_all$`P Retention` < 1e-5 |
        maize_all$`Total P` < 1e-5 |
        maize_all$`Alkaline Soil P` < 1e-5 ,1
]

quartz()
setwd("~/Desktop")
# CMplot(maize_all,type="p",plot.type="m",LOG10=TRUE,highlight=SNPs,highlight.type="l",
#        threshold=1e-4,threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),
#        signal.cex=1.2, signal.col="red", highlight.col="grey",highlight.cex=0.7,
#        file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)

CMplot(maize_all, plot.type="m",multracks=TRUE,threshold=c(1e-40,1e-20),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green"),
       signal.cex=0.3, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
       highlight=NULL, highlight.text=NULL, highlight.text.cex=1.4)


##########################################################################################
# SORGHUM
##########################################################################################
#30pca

sorghum_dir <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/GWAS-TASSEL/"
sorghum_dir <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GWAS.results/sorghum_Output/v3/GLM_GWAS/GWAS_sorghum_v3_2kb_csv/"
PBR1 <- vroom(paste0(sorghum_dir,"PBR1.txt"))  
PBR1 <- PBR1[,c(2,3,4,6)]
PBR1 <- na.omit(PBR1)
PNZ1 <- vroom(paste0(sorghum_dir,"PNZ1.txt"))
PNZ1 <- PNZ1[,c(2,3,4,6)]
PNZ1 <- na.omit(PNZ1)
POL1 <- vroom(paste0(sorghum_dir,"POL1.txt"))
POL1 <- POL1[,c(2,3,4,6)]
POL1 <- na.omit(POL1)
stp10 <- vroom(paste0(sorghum_dir,"stp10.txt"))
stp10 <- stp10[,c(2,3,4,6)]
stp10 <- na.omit(stp10)


sorghum_sol_ret <- cbind(PBR1, PNZ1$p)
names(sorghum_sol_ret) <- c("SNP","Chromosome","Position","P Solubility","P Retention")

library(dplyr)
sorg_tot_alk <- left_join(stp10, POL1, by = "Marker")
sorg_tot_alk <- as.data.frame(sorg_tot_alk[,c(1,2,3,4,7)])
names(sorg_tot_alk) <- c("SNP","Chromosome","Position","Total P","Alkaline Soil P")


sorghum_all <- base::merge(sorghum_sol_ret,sorg_tot_alk)
names(sorghum_all)

SNPs <- sorghum_all[
    sorghum_all$`P Solubility` < 1e-5 |
        sorghum_all$`P Retention` < 1e-5 |
        sorghum_all$`Total P` < 1e-5 |
        sorghum_all$`Alkaline Soil P` < 1e-5 ,1
]

quartz()
setwd("~/Desktop")
# CMplot(sorghum_all,type="p",plot.type="m",LOG10=TRUE,highlight=SNPs,highlight.type="l",
#        threshold=1e-4,threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),
#        signal.cex=1.2, signal.col="red", highlight.col="grey",highlight.cex=0.7,
#        file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)

CMplot(sorghum_all, plot.type="m",multracks=TRUE,threshold=c(1e-5,1e-4),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green"),
       signal.cex=0.3, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
       highlight=NULL, highlight.text=NULL, highlight.text.cex=1.4)



#################################################################################
# ALL FILES
########### set the directory path where your text files are located
# set the directory path where your text files are located
dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GWAS.results/maize_Output/v5/LMM_GWAS"

# get a list of all the text files in the directory with the .txt extension
file_list <- list.files(dir, pattern = "\\.txt$")

# initialize a list to hold the data frames for each file
data_list <- list()

# loop over the file names and read in the data using vroom
for (file_name in file_list) {
    # read in the file using vroom and add it to the data list
    file_data <- vroom(paste0(dir, "/", file_name))
    data_list[[file_name]] <- file_data[,c("chr", "ps", "p_wald")]
}

# combine the data frames into a single data frame using bind_rows
combined_data <- bind_rows(data_list, .id = "file")

# convert the chromosome column to a factor with levels sorted in ascending order
combined_data$chr <- factor(combined_data$chr, levels = sort(unique(combined_data$chr)))

# calculate the negative log10 of the p-values
combined_data$log_p <- -log10(combined_data$p_wald)

# split the data by file and create a list of ggplot objects
plot_list <- split(combined_data, f = combined_data$file) %>%
    purrr::map(function(df) {
        ggplot2::ggplot(df, ggplot2::aes(x = ps, y = log_p, color = chr)) +
            ggplot2::geom_point(size = 0.5) +
            ggplot2::scale_color_manual(values = rainbow(length(unique(df$chr)))) +
            ggplot2::facet_wrap(. ~ file, nrow = 5, ncol = 6, scales = "free_x") +
            ggplot2::theme_bw() +
            ggplot2::theme(strip.background = ggplot2::element_blank(),
                           strip.text = ggplot2::element_blank())
    })

# plot the list of ggplot objects using gridExtra
library(gridExtra)
quartz()
grid.arrange(grobs = plot_list[1:5], nrow = 5, ncol = 1)




################
#single plot
dir <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Results/GWAS.results/sorghum_Output/v3/LMM_GWAS/GWAS_LMM_results/"
data <- read.table(paste0(dir,"sol_Mo.txt"), header = T) %>% dplyr::select(c(2,1,3,12))
colnames(data) <- c("SNP","Chromosome","Position","sol_Mo")
CMplot(data,type="p",plot.type="m",LOG10=TRUE,threshold=NULL,file="jpg",dpi=300,
       file.output=TRUE,verbose=TRUE,width=14,height=6,chr.labels.angle=45)

CMplot(data, plot.type="m", col=c("grey30","grey60"), LOG10=TRUE, ylim=c(2,12), threshold=c(1e-3,1e-3),
       threshold.lty=c(1,2), threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,
       chr.den.col=NULL, signal.col=c("red","green"), signal.cex=c(1.5,1.5),signal.pch=c(19,19),
       file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
# 'chr.labels.angle': adjust the angle of labels of x-axis (-90 < chr.labels.angle < 90).
# file.name: specify the output file name, the default is corresponding column name when setting ' file.name="" '.