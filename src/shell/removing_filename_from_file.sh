#!/bin/sh

for file in $noheader_*.txt
do
mv $file ${file#noheader_}
done