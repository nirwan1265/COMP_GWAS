# separate fasta files into headers and sequences:
headers = []
sequences = []
# there should be a space for > genename
with open('protein.faa', 'r') as file:
    lines = file.readlines()
    sequence = ''
    for line in lines:
        if line.startswith('>'):
            headers.append(line.strip())
            if sequence:
                sequences.append(sequence)
                sequence = ''
        else:
            sequence += line.strip()
    sequences.append(sequence)  # add the last sequence to the list

with open('headers.fasta', 'w') as headers_file:
    headers_file.write('\n'.join(headers) + '\n')

with open('sequences.fasta', 'w') as sequences_file:
    sequences_file.write('\n'.join(sequences) + '\n')
    
# extracting the first two elements of the fasta header files:
# in shell
# cut -d ' ' -f 1,2 headers.fasta > gene_names.txt



# Combining the genes name and sequences 
with open('gene_names.txt', 'r') as gene_names_file, open('sequences.fasta', 'r') as sequences_file, open('output.fasta', 'w') as output_file:
    # read the gene names and sequences into two separate lists
    gene_names = gene_names_file.read().splitlines()
    sequences = sequences_file.read().splitlines()

    # write the alternating gene names and sequences to the output file
    for i in range(len(gene_names)):
        output_file.write(gene_names[i] + '\n')
        output_file.write(sequences[i] + '\n')
