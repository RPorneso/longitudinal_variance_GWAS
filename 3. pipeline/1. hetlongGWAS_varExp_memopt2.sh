#!/bin/bash
#
#SBATCH --job-name="<pheno>"
#SBATCH --account=p805
#SBATCH --time=60:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-3491

## Set up job enviroment:

source /cluster/bin/jobsetup
module purge
set -o errexit

module load R-bundle-CRAN/2023.12-foss-2023a

Rscript hetlongGWAS_varExp_memopt2.R 3491 2000 6981748 '/tsd/p805/data/durable/projects/ralphp/ANTHRO/data/analysis_<pheno>.Rda' '/tsd/p805/data/durable/projects/ralphp/AP/data/bigsnpr/noimpute/moba.rds' '/tsd/p805/data/durable/projects/ralphp/ANTHRO/results/<outfile>'
