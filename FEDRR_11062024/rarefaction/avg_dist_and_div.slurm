#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="rarefactions"
#SBATCH --output=rarefactions.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate r-4.2

Rscript ~/repo/BIL_pathogens/FEDRR_11062024/rarefactions/avg_dist_and_div.samps.premise.R
