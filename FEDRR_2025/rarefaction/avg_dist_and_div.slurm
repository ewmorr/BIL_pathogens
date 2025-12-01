#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="avg_div"
#SBATCH --output=avg_div_2.out
#SBATCH --partition=shared
#SBATCH --exclude=node117,node118

module purge
module load anaconda/colsa

conda activate r-4.2

Rscript ~/repo/BIL_pathogens/FEDRR_2025/rarefaction/avg_dist_and_div.premise.R
