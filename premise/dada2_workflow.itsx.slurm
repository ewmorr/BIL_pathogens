#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="itsx"
#SBATCH --output=dada2_workflow.itsxpress.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24

module purge
module load linuxbrew/colsa

cd /mnt/home/garnas/ewj4/EDRR_patho/

ITSx -i ASVs.fa -o out_itsx \
    -t F --allow_single_domain 0.00001,10 \
    -E 0.00001 --cpus 24 \
    -p ~/repo/itsxpress_2x_mod/itsxpress/ITSx_db/HMMs

cd /mnt/home/garnas/ewj4/EDRR_patho/itsxpress_mod_dada2_core

ITSx -i ASVs.fa -o out_itsx \
    -t F --allow_single_domain 0.00001,10 \
    -E 0.00001 --cpus 24 \
    -p ~/repo/itsxpress_2x_mod/itsxpress/ITSx_db/HMMs

#ITSx -i itsxpress_out_R2.derep.fasta -o out_itsx_xpress \
#    -t F --allow_single_domain 0.00001,10 \
#    -E 0.00001 --cpus 24 \
#    -p ~/repo/itsxpress_2x_mod/itsxpress/ITSx_db/HMMs

#ITSx -i itsxpress_mod_out_R2.derep.fasta -o out_itsx_xpress_mod \
#    -t F --allow_single_domain 0.00001,10 \
#    -E 0.00001 --cpus 24 \
#    -p ~/repo/itsxpress_2x_mod/itsxpress/ITSx_db/HMMs
