#!/bin/sh

for i in *.txt
do
sed '1d' $i > noheader_$(basename $i)
done