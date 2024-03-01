#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="itsxpress"
#SBATCH --output=dada2_workflow.itsxpress_mod.out
#SBATCH --partition=shared
#SBATCH --cpus-per-task=24

module purge
module load anaconda/colsa
conda activate itsxpress_2x_mod

cd /mnt/home/garnas/ewj4/EDRR_patho

#mkdir itsxpress_mod_out
outdir=itsxpress_mod_out

for i in cutadapt/*R1*
do(

    dir=${i%/*}
    r1File=${i##*/}
    pre=${r1File%R1*}
    post=${r1File##*R1}
    r2File=${pre}R2${post}

    itsxpress --fastq $dir/$r1File --fastq2 $dir/$r2File \
        --reversed_primers --threads 24 --cluster_id 1 --taxa Fungi --region ITS2 \
        --outfile $outdir/$r1File --outfile2 $outdir/$r2File
)
done

conda deactivate itsxpress_2x_mod

Rscript ~/repo/BIL_pathogens/sequence_processing_and_dada2/3_dada2_core.r /mnt/home/garnas/ewj4/EDRR_patho/itsxpress_mod_out /mnt/home/garnas/ewj4/EDRR_patho/itsxpress_mod_dada2_core

module purge
module load linuxbrew/colsa

cd /mnt/home/garnas/ewj4/EDRR_patho/itsxpress_mod_dada2_core

ITSx -i ASVs.fa -o out_itsx_0.01_ \
    -t F --allow_single_domain 0.01,10 \
    -E 0.01 --cpus 24 \
    -p ~/repo/itsxpress_2x_mod/itsxpress/ITSx_db/HMMs

ITSx -i ASVs.fa -o out_itsx_ \
    -t F --allow_single_domain 0.00001,10 \
    -E 0.00001 --cpus 24 \
    -p ~/repo/itsxpress_2x_mod/itsxpress/ITSx_db/HMMs
