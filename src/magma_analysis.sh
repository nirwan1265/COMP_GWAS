#!/bin/zsh

# When everything is in one folder

files=( apa ext_P20 lab NPlim occ org PBR1 PBR2 PHO1 PHO2 PMEH1 PMEH2 PNZ1 PNZ2 POL1 POL2 sec sol_Hi sol_Lo sol_Mo sol sol_VL stp10 stp20 stp30 stp100 tot TP1 TP2 )

for i in $files
do
magma --bfile $i --gene-annot "${i}"_sorghum_LMM.genes.annot --pval "${i}"_sorghum_LMM.pvalue.txt N="$(wc -l < "${i}"_sorghum_LMM.snpfile.txt)" --gene-model multi=snp-wise --out "${i}"_magma_multi_snpwise
done
