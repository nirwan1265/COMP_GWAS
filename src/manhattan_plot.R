#Package
library("CMplot")
library(vroom)

#Path
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype")

# Data
maize <- vroom("PNZ1_maize.txt")
head(maize)
maize <- maize[,c(2,1,3,11)]
names(maize) <- c("SNP","Chromosome","Position","PNZ1")

sorghum <- vroom("PNZ1_sorghum.txt")
head(sorghum)
names(sorghum) <- c("SNP","Chromosome","Position","PNZ1")

#SNP-density plot
# Need all SNPs
CMplot(maize,type="p",plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300,
       main="Density plot for Maize",file.output=TRUE,verbose=TRUE,width=9,height=6)

CMplot(sorghum,type="p",plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300,
       main="Density plot for Sorghum",file.output=TRUE,verbose=TRUE,width=9,height=6)



# Manhattan GWAS- amplify signals
CMplot(maize, plot.type="m", col=c("grey30","grey60"), LOG10=TRUE, ylim=c(2,22), threshold=c(1e-6,1e-4),
       threshold.lty=c(1,2), threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,
       chr.den.col=NULL, signal.col=c("#F05039","#7CA1CC"), signal.cex=c(1.5,1.5),signal.pch=c(19,19),
       file="tiff",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=25,height=6)


CMplot(sorghum, plot.type="m", col=c("grey30","grey60"), LOG10=TRUE, ylim=c(2,22), threshold=c(1e-6,1e-4),
       threshold.lty=c(1,2), threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,
       chr.den.col=NULL, signal.col=c("#F05039","#7CA1CC"), signal.cex=c(1.5,1.5),signal.pch=c(19,19),
       file="tiff",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=25,height=6)





#Need LML GWAS file all assoc
setwd("~/Desktop")

occ_ret <- vroom("occ_ret.txt")
occ_ret <- occ_ret[,c(2,1,3,13)]
tot_tot <- vroom("tot_tot.txt")
tot_tot <- tot_tot[,c(2,1,3,13)]
lab_sol <- vroom("lab_sol.txt")
lab_sol <- lab_sol[,c(2,1,3,13)]


sorg_all <- cbind(occ_ret, tot_tot$p_wald, lab_sol$p_wald)
names(sorg_all)[4:6] <- c("Retention","Total","Solubility")
names(sorg_all)[1:3] <- c("SNP","Chromosome","Position")

SNPs <- sorg_all[
    sorg_all$Retention < 1e-5 |
    sorg_all$Total < 1e-5 |
    sorg_all$Solubility < 1e-5, 1
]

quartz()
CMplot(sorg_all,type="p",plot.type="m",LOG10=TRUE,highlight=SNPs,highlight.type="l",
       threshold=1e-4,threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),
       signal.cex=1.2, signal.col="red", highlight.col="grey",highlight.cex=0.7,
       file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)

CMplot(sorg_all, plot.type="m",multracks=TRUE,threshold=c(1e-6,1e-4),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green"),
       signal.cex=1, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
       highlight=NULL, highlight.text=NULL, highlight.text.cex=1.4)
