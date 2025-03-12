#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="avg_div"
#SBATCH --output=avg_div.out
#SBATCH --partition=shared
#SBATCH --exclude=node117,node118

module purge
module load anaconda/colsa

conda activate r-4.2

Rscript ~/repo/BIL_pathogens/FEDRR_all_2024/rarefaction/avg_dist_and_div.funnel_traps_and_spore_traps.premise.R
