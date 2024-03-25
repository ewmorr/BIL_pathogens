#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="vsearch"
#SBATCH --output=vsearch_UDB.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24

module purge
module load linuxbrew/colsa

cd /mnt/home/garnas/ewj4/blast_dbs/unite_07252023

vsearch --makeudb_usearch sh_general_allEuk_dynamic_singletons_25072023.ITS2.fasta \
    --output sh_general_allEuk_dynamic_singletons_25072023.ITS2.udb \
        --minseqlength 16 --notrunclabels --threads 24

#54 seqs discarded at --minseqlength 16. There were a few more at the default of 32. We'll keep at 16 as these are likely not fungi and our min seq len is 51
