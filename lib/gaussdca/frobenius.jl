using GaussDCA
fasta_file=ARGS[1]
output_file=ARGS[2]
FNR = gDCA(fasta_file)
printrank(output_file, FNR)



