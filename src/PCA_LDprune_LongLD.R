library(bigsnpr)
library(ggplot2)

NCORES <- nb_cores()

plink.path <- "/Users/nirwantandukar/Documents/plink-1.07-mac-intel/plink2"
prefix.in <- "/Users/nirwantandukar/Documents/plink-1.07-mac-intel/allchrom"


poprescQC.bed <- snp_plinkQC(
  plink.path = plink.path,
  prefix.in = prefix.in,
  file.type = "--vcf",
  geno = 0.04, 
  mind = 0.04, 
  maf = 0.05, 
  #hwe = 1e-10,
  autosome.only = T
)

popresQC2.bed <- snp_plinkIBDQC(poprescQC.bed,
                                plink.path = plink.path,
                                ncores = NCORES)
popresQC.rds <- snp_readBed(popresQC2.bed, "/Users/nirwantandukar/Documents/plink-1.07-mac-intel/ch1")


ch1 <- snp_attach("/Users/nirwantandukar/Documents/plink-1.07-mac-intel/ch1.rds")
G <- ch1$genotypes
CHR <- ch1$map$chromosome
POS <- ch1$map$physical.pos
NCORES <- nb_cores()

# "Verification" there is no missing value
big_counts(G, ind.col = 1:12) # OK





# PCA using LD pruning and removing long distance LD regions
#svd1 <- big_randomSVD(G, snp_scaleBinom(), ncores = NCORES)
ind.keep4 <- snp_clumping(G, CHR, ncores = NCORES,
                          exclude = snp_indLRLDR(CHR, POS))
svd4 <- big_randomSVD(G, snp_scaleBinom(), ncores = NCORES,
                      ind.col = ind.keep4)

plot(svd4, type = "scores", scores = 3:4) +
  aes(color = pop.names[pop]) #+
  labs(color = "Population")
