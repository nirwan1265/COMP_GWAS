#!/bin/bash

# Set input and output file names
input_file="input.txt"
output_file="output.txt"

# Loop over each line of the input file
while read col1 col2 col3 col4; do
    # Subtract 2000 from the third column
    new_col3=$((col3 - 2000))
    # Add 2000 to the fourth column
    new_col4=$((col4 + 2000))
    # Output the modified line to the output file
    echo "$col1 $col2 $new_col3 $new_col4" >> $output_file
done < $input_file