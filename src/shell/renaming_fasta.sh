 # renameing fasta file name
 awk '/^>/ {print $1} !/^>/ {print}' sorghum.fa > sorghum_fil.fa