#!/bin/zsh


files=( sol.fam sol_Hi.fam sol_Mo.fam sol_VL.fam tot.fam lab.fam org.fam occ.fam sec.fam apa.fam stp10.fam stp100.fam stp20.fam stp30.fam$


for file in "${files[@]}"
do
cp allchr_sorghum_africa.fam $file
done