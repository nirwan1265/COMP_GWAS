#!/bin/sh
for i in *.txt
do awk -F"\t" '{print $2,$1,$3,$11}’ $i > “${i%.*}.filtered.${i##*.}”
done


# $2, $1, $3 and $11 represents markers, chromosome, snp position, pvalues in the gwas file
# run the script in the folder where your gwas results are. 
