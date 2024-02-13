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
    -p ~/repo/itsxpress_2x_mod/itsxpress/ITSx_db/HMMs \
    --save_regions ITS2 \
    -t F -N 1 \ #N is the minimal number of domains to include a seq. default is 2. need 1
    --cpu 24 --graphical F --fasta F
)
done

