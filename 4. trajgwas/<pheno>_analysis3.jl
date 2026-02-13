# Start by loading packages

using DataFrames
using CSV
using TrajGWAS

# SNP index

function block_index(ind, nblocks, size, maxcol)
	if ind > nblocks
		print("Error: index cannot be longer than nblocks")
	else
		first = 1
		last = first + size
	end
	if ind > 1
		first = first + (ind-1) * size
		last = first + size - 1
	end
	if ind == nblocks
		last = first + mod(maxcol, size) - 1
	end
	[first, last]
end

# Create index to match plink files with mobaread4 individuals while fitting a null model then create SNP index

famfileids = CSV.read("/tsd/p805/data/durable/data/genetics/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.fam", DataFrame, header=false)[!,2] # yes this is correct
famfileids = string.(famfileids)

covdf = CSV.read("/tsd/p805/data/durable/projects/ralphp/AP/data/<pheno>.csv", DataFrame)
covdf[!,2] = string.(covdf[!,2])

covrowmask, geneticrowmask = matchindices(@formula(resZ_SCORE ~ 1 + AARGANG + time + sex + genotyping_batch_num + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10), @formula(resZ_SCORE ~ 1), @formula(resZ_SCORE ~ 1 + AARGANG + time + sex + genotyping_batch_num + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10), :IID, covdf, famfileids)

taskid = Base.parse(Int, ENV["SLURM_ARRAY_TASK_ID"]) 
index = block_index(taskid, parse(Int, ARGS[1]), parse(Int, ARGS[2]), parse(Int, ARGS[3])) # SNP index

# Fit model

trajgwas(@formula(resZ_SCORE ~ 1 + AARGANG + time + sex + genotyping_batch_num + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10), @formula(resZ_SCORE ~ 1), @formula(resZ_SCORE ~ 1), :IID, "/tsd/p805/data/durable/projects/ralphp/AP/data/mobaread4.csv", "/cluster/projects/p805/data/genetics/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc", pvalfile = "/tsd/p805/data/durable/projects/ralphp/AP/results/trajread_resid" * string(index[2]), nullfile = "/tsd/p805/data/durable/projects/ralphp/AP/results/trajread_resid_null", covrowinds = covrowmask, geneticrowinds = geneticrowmask, snpinds = index[1]:index[2]) 

exit()
