# GWAS-GBJ
Easy combination of SNP pvalues from GWAS data using the R GBJ package.


## Requirements
1.) GWAS file with the following information
|Markers|Chrom|Position|Pvalue|
|----|----|----|----|
|M1|1|1|0.05|
|M2|2|10|0.04|
|M3|3|25|0.01|

Better if the files are separated into chromosomes. So for 10 chromosomes, you will have 10 files.
Z.State column not required. 


2.) Genotype file
Should be noted in the normal MAF type. Again, faster if it is divided according to the chromosomes.
M1---M2---M3---M4
0---1---2----9

3.) PCA file
THe PCA file that you used to run the model. 10 PCAs are preferred. 





