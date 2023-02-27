#!/bin/zsh
# chr 1 has 1 and 10 so remove 10
# in maize chromosome 10 has a random chr 4's SNPs in the first line so remove that too

head -n 1 tot.txt > header.txt

for i in *.txt
do
head -n 1 $i > "${i%.*}"02.txt
grep "^S2" $i >> "${i%.*}"02.txt
head -n 1 $i > "${i%.*}"03.txt
grep "^S3" $i >> "${i%.*}"03.txt
head -n 1 $i > "${i%.*}"04.txt
grep "^S4" $i >> "${i%.*}"04.txt
head -n 1 $i > "${i%.*}"05.txt
grep "^S5" $i >> "${i%.*}"05.txt
head -n 1 $i > "${i%.*}"06.txt
grep "^S6" $i >> "${i%.*}"06.txt
head -n 1 $i > "${i%.*}"07.txt
grep "^S7" $i >> "${i%.*}"07.txt
head -n 1 $i > "${i%.*}"08.txt
grep "^S8" $i >> "${i%.*}"08.txt
head -n 1 $i > "${i%.*}"09.txt
grep "^S9" $i >> "${i%.*}"09.txt
head -n 1 $i > "${i%.*}"10.txt
grep "^S10" $i >> "${i%.*}"10.txt
head -n 1 $i > "${i%.*}"01.txt
grep "^S1" $i >> "${i%.*}"01.txt
done