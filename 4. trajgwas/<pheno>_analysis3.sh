#!/bin/bash
#SBATCH --job-name="traj_<pheno>"
#SBATCH --account=p805_tsd
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-3491

module purge
module load Julia/1.10.0-linux-x86_64
export JULIA_DEPOT_PATH=/tsd/p805/cluster/julia_TrajGWAS

julia <pheno>_analysis3.jl 3491 2000 6981748
