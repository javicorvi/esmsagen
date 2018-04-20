using GaussDCA
fasta_file=ARGS[1]
output_file=ARGS[2]
DIR = gDCA(fasta_file, pseudocount = 0.2, score = :DI)
printrank(output_file, DIR)



