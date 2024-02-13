#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="itsx"
#SBATCH --output=dada2_workflow.itsxpress.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24

module purge
module load linuxbrew/colsa

cd /mnt/home/garnas/ewj4/EDRR_patho/itsx_test

for i in *unique.fasta
do(

    sample=${i%.unique.fasta}

    ITSx -i $i -o $sample \
    -t F --allow_single_domain 1,10 \
    -E 1 --cpus 24 \
    -p ~/repo/itsxpress_2x_mod/itsxpress/ITSx_db/HMMs
)
done

