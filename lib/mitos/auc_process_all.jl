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
using MIToS.PDB
using ROCAnalysis
using DataFrames
pdb_name=ARGS[1]
model_name = ARGS[2]
chain_name=ARGS[3]
contact_map_file = ARGS[4]
sequences_paths = ARGS[5] 
result_file = ARGS[6]
result_path_zmip = ARGS[7]
zmip_file_result = "zmip"
println(pdb_name) 
println(model_name) 
println(chain_name) 
println(contact_map_file) 
println(sequences_paths) 
println(result_file) 
contact_map = readdlm(contact_map_file)
contact_map_triu = triu(contact_map,1)	
df = DataFrame(A = AbstractString[], B = Float64[])
for f in filter(x -> (contains(x, pdb_name) && endswith(x, ".cluster")), readdir(sequences_paths))
	println("$sequences_paths$f")           
	beta = f[16:19]
	nsus = f[25:28]
	runs = f[34:search(f,'.',34)-1]
	msa = read("$sequences_paths$f",FASTA)
	zmip, mip = buslje09(msa)
	#zmip_triu = triu(zmip,1)	
	value_auc = 1.0 - auc(roc(get_target_nontarget(zmip, contact_map)...))
	#value_auc_triu = 1.0 - auc(roc(get_target_nontarget(zmip_triu, contact_map_triu)...))
	println(value_auc)
	writedlm("$result_path_zmip$zmip_file_result$f" ,zmip)
	push!(df, [f  value_auc])
end
writetable(result_file, df, separator = ',', header = false)
