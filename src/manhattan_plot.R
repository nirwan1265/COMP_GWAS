#Package
library("CMplot")

#Path
setwd("~/Desktop")

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



