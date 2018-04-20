using MIToS.MSA
using MIToS.Information
fasta_file=ARGS[1]
zmip_result_path=ARGS[2]
println(fasta_file)
msa = read(fasta_file, FASTA)
zmip, mip = buslje09(msa)
writedlm(zmip_result_path,zmip)



