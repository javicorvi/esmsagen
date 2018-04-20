function get_target_nontarget(scores, contacts)
	@assert size(scores) == size(contacts)
	I,J = size(scores)
	@assert I == J
	tar = Float64[]
	non = Float64[]
	for i in 1:(I-1)
		for j in (i+1):I
			if contacts[i,j]
        	        	push!(tar, scores[i,j])
			else
				push!(non, scores[i,j])
			end 
		end
	end
	tar, non
end
using MIToS.MSA
using MIToS.Information
using ROCAnalysis
fasta_file=ARGS[1]
contact_map_file = ARGS[2]
println(fasta_file)
println(contact_map_file)
msa = read(fasta_file, FASTA)
contact_map = readdlm(contact_map_file)
zmip, mip = buslje09(msa,maxgap=1.0)
println(size(zmip))
println(size(contact_map)) 
value_auc = 1.0 - auc(roc(get_target_nontarget(zmip, contact_map)...))
println(value_auc)
