#!/bin/bash

for file in *.txt
do
awk '{print $2, $1, $3, $14}' "$file" > "${file%.*}_filtered.txt"
done