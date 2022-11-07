# GWAS-GBJ
Easy combination of SNP pvalues from GWAS data using the R GBJ package.


## Requirements
1. GWAS file:

A text file with the following information

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
***
2. Genotype file

A genotype text file with MAF numbering with markers in the column and subsequent values in the rows. Each row should represent an accesssion/individual.

|M1|M2|M3|M4|
|----|---|---|---|
|0|1|2|-9|
|0|2|0|1|
***
3. PCA file

The PCA file that you used to run the linear mixed model. 10 PCAs are preferred. 

|PCA1|PCA2|PCA3|PCA4|PCA5|PCA6|PCA7|PCA8|PCA9|PCA10|
|---|---|---|---|---|---|---|---|---|---|
|0.07|0.78|0.08|0.65|0.48|0.84|0.15|0.26|0.47|0.15|
|0.25|0.19|0.69|0.56|0.47|0.12|0.69|0.47|0.59|0.33|
***







