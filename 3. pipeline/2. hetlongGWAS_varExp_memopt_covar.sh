#!/bin/bash
#
#SBATCH --job-name="covar_<pheno>"
#SBATCH --account=p805
#SBATCH --time=00:45:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-61                # THIS IS NOW THE NUMBER OF SNPs IN sigsnps_<pheno>_index.txt

## Set up job enviroment:

source /cluster/bin/jobsetup
module purge
set -o errexit

module load R-bundle-CRAN/2023.12-foss-2023a

Rscript hetlongGWAS_varExp_memopt_covar.R '/tsd/p805/data/durable/projects/ralphp/ANTHRO/data/analysis_<pheno>.Rda' '/tsd/p805/data/durable/projects/ralphp/AP/data/bigsnpr/noimpute/moba.rds' '/tsd/p805/data/durable/projects/ralphp/sigsnps_index2/sigsnps_<pheno>_index.txt' '/tsd/p805/data/durable/projects/ralphp/ANTHRO/results/<outfile>'
