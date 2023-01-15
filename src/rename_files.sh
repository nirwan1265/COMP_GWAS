#!/bin/zsh

for f in *.allchrom_maize_pheno.bed; do
    mv $f ${f%.allchrom_maize_pheno.bed}.bed
done

for f in *.allchrom_maize_pheno.bim; do
    mv $f ${f%.allchrom_maize_pheno.bim}.bim
done

for f in *.allchrom_maize_pheno; do
    mv $f ${f%.allchrom_maize_pheno}.fam
done