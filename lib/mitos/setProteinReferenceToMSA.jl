using MIToS.MSA
fasta_file=ARGS[1]
println(fasta_file)
msa = read(fasta_file, FASTA)
setreference!(msa, "THIO_ECOLI/4-107")
adjustreference!(msa)
write("../data/natural/PF00085_THIO_ECOLI_reference.fasta", msa, FASTA)


