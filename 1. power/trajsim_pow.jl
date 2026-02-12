using TrajGWAS, CSV, DataFrames, Tables

files = readdir("/Users/ralphp/trajsimdat/power/loglin/", join=true)
pvals_e = zeros(length(files))

for i in 1:length(files)
	dat = CSV.read(files[i], DataFrame)
	pvalpath = "/Users/ralphp/pval.txt"
	nullpath = "/Users/ralphp/null.txt"
	trajgwas(@formula(phe~1+time+snp), @formula(phe~1),  @formula(phe~1+snp), :id, dat, pvalfile=pvalpath, nullfile=nullpath)
	res = readlines("null.txt") # saves as an array
	coefs = res[18:22] # retrieves b and tau coefs incl. pvals
	p_tausnp = coefs[5][50:length(coefs[3])+1]
	if contains(p_tausnp, "<")
		p_tausnp = 0.00000000000000000001
	else p_tausnp = parse(Float64, p_tausnp)
	end
	pvals_e[i] = p_tausnp
end

CSV.write("/Users/ralphp/trajsimdat/power/traj_pvals.csv",  Tables.table(pvals_e), writeheader=false)

files[begin:100:end]
dgm = repeat(["RI","RIRI","RIS"], inner = 100*4) # there are 100 iterations & 4 models
tau = repeat([0, 0.01, 0.05, 0.1], inner = 100, outer = 3)
dgm[begin:100:end]
tau[begin:100:end]

trajres = [dgm tau pvals_e]

CSV.write("/Users/ralphp/trajsimdat/power/trajres.csv",  Tables.table(trajres), writeheader=true)

