Simulation codes and analytic pipeline repository for the manuscript "Longitudinal modelling reveals widespread non-additive genetic effects underlying developmental plasticity". 

Link to preprint: https://www.biorxiv.org/content/10.64898/2025.12.19.695443v1.full.

To re-create the simulations:

1. Download all functions in the 0. dgm, 1. power, and 2. bias folders. Update directories.
2. To check how genetic interactions generate SNP variance effects, refer to Pilot.pdf. Codes are available in Pilot.Rmd.
3. To assess power and bias of our pipeline, run scripts in 1. power and 2. bias folders. You will need to install Julia to compare our pipeline against trajGWAS.

A few notes on the analytic pipeline:

1. Each phenotype per time point was residualized on the GRM using GCTA.
2. The pipeline runs on multiple nodes in parallel using HPC Colossus (Linux).
3. The inputs in the bash are: blocks (nodes), number of SNPs per block, total number of SNPs, analysis dataframe, plink files, and output files.
4. Each output file contains 2000 SNPs. They need to be combined into 1 summary statistics file.
5. Since we relied on LRT to assess SNP mean and variance effect, we repeated the analysis two times for significant SNPs, where we -
   a. added PCs in the residual; and for those that remain significant 
   b. isolated "variance" SNPs by running LRT which compares a mean only against a mean-and-variance model.

