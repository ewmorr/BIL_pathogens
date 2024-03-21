#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="itsx"
#SBATCH --output=itsx_unite.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24

module purge
module load linuxbrew/colsa

cd /mnt/home/garnas/ewj4/EDRR_patho/

dir=/mnt/home/garnas/ewj4/blast_dbs/unite_07252023
file=sh_general_release_dynamic_singletons_allEuk_25.07.2023.fasta


ITSx -i $dir/$file -o $dir/sh_general_allEuk_dynamic_singletons_25072023 \
    --preserve T \
    --cpus 24 \
    -p ~/repo/itsxpress_2x_mod/itsxpress/ITSx_db/HMMs
