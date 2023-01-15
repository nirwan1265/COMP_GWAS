#!/bin/zsh

for n in $(seq 1 $(head -n1 Maize_allphospho.csv | tr "," "\n" | wc -l)); do
    awk -F, -v col="$n" 'NR>1{print $col}' Maize_allphospho.csv > $(head -n1 Maize_allphospho.csv | tr "," "\n" | sed -n ${n}p).txt
done

#Remove "" from file names manually cause I am lazy to write the script for it. 
# Forgive me