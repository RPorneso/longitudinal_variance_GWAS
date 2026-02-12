#!/usr/bin/bash
#SBATCH --job-name="resid_<grade><pheno>"
#SBATCH --account=p805
#SBATCH --partition=bigmem
#SBATCH --time=80:00:00
#SBATCH --mem-per-cpu=40G
#SBATCH --cpus-per-task=10

/cluster/projects/p805/software/gcta/gcta-1.94.1 \
--reml \
--grm /tsd/p805/data/durable/data/genetics/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc_regenie_500k_snps \
--pheno /tsd/p805/data/durable/projects/ralphp/AP/data/pheno_<grade><pheno>.txt \
--reml-pred-rand \
--out /tsd/p805/data/durable/projects/ralphp/AP/data/resid_<pheno><grade>
