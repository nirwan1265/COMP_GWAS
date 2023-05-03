#!/bin/sh

for file in $noheader_*.txt
do
mv $file ${file#noheader_}
done

for file in *_filtered.txt; do mv "$file" "${file/_filtered/}"; done