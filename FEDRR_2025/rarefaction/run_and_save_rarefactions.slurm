#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="rarefactions"
#SBATCH --output=rarefactions.out
#SBATCH --partition=shared
#SBATCH --exclude=node117,node118

module purge
module load anaconda/colsa

conda activate r-4.2

Rscript ~/repo/BIL_pathogens/FEDRR_2025/rarefaction/run_and_save_rarefactions.premise.R
