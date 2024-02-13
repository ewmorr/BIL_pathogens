#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="itsx"
#SBATCH --output=dada2_workflow.itsxpress.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24

module purge
module load linuxbrew/colsa

cd /mnt/home/garnas/ewj4/EDRR_patho/itsx_test

for i in *unique*
do(

    sample=${i%.unique*}

    ITSx -i $i -o $sample \
    -t F --allow_single_domain 1e-5,0 \
    --cpus 24 --graphical F --fasta F \
    -p ~/repo/itsxpress_2x_mod/itsxpress/ITSx_db/HMMs
)
done

