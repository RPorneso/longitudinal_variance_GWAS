#!/bin/bash
#
#SBATCH --job-name="vexp<pheno>"
#SBATCH --account=p805
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1

## Set up job enviroment:

source /cluster/bin/jobsetup
module purge
set -o errexit

module load R-bundle-CRAN/2023.12-foss-2023a

Rscript hetlongGWAS_varExp_memopt_lrtv.R '/tsd/p805/data/durable/projects/ralphp/AP/data/analysis_<pheno>.Rda' '/tsd/p805/data/durable/projects/ralphp/AP/data/bigsnpr/noimpute/moba.rds' '/tsd/p805/data/durable/projects/ralphp/sigsnps_covar_index/exp<pheno>_resid_covarinvar_index.txt' '/tsd/p805/data/durable/projects/ralphp/snps_lrtv_bef_prune/<outfile>'
