# GWAS-GBJ
Easy combination of SNP pvalues from GWAS data using the R GBJ package.


## Requirements
1. GWAS text file with the following information
Chromosome 1
|Markers|Chrom|Position|Pvalue|
|----|----|----|----|
|M1|1|1|0.05|
|M2|1|10|0.04|
|M3|1|25|0.01|

Chromosome 2
|Markers|Chrom|Position|Pvalue|
|----|----|----|----|
|M1|2|1|0.05|
|M2|2|10|0.04|
|M3|2|25|0.01|

The files should be separated based on their chromosome. The chromosome should be named as numbers. 

The script to separate chromosome can be found here : 

2. Genotype file
A genotype text file with MAF numbering with markers in the column and subsequent values in the rows. Each row should represent an accesssion/individual.

|M1|M2|M3|M4|
|----|---|---|---|
|0|1|2|-9|
|0|2|0|1|

3. PCA file
THe PCA file that you used to run the model. 10 PCAs are preferred. 





