#!/bin/zsh

for f in *.allchrom.bed; do
    mv $f ${f%.allchrom.bed}.bed
done

for f in *.allchrom.bim; do
    mv $f ${f%.allchrom.bim}.bim
done

for f in *.allchrom.fam; do
    mv $f ${f%.allchrom.fam}.fam
done